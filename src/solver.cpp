#include "solver.hpp"

#include "spectral.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace {
SpectralField makeSpectralField(std::size_t count) { return SpectralField(count); }

void requireFinite(const SpectralField &values, const char *description) {
    if (!std::all_of(values.begin(), values.end(), [](Complex value) {
            return std::isfinite(value.real()) && std::isfinite(value.imag());
        }))
        throw std::runtime_error(std::string(description) + " contains a non-finite coefficient");
}

std::uint64_t resolveRandomSeed(std::uint64_t seed) {
    if (seed != 0)
        return seed;
    return static_cast<std::uint64_t>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());
}

double forcingEnvelope(const Parameters &parameters, std::size_t index) {
    const double k = std::abs(waveNumber(parameters, index));
    if (k == 0.0)
        return 0.0;
    switch (parameters.forcingProfile) {
    case ForcingProfile::annulus:
        return std::abs(k - parameters.forcingWavenumber) < parameters.forcingWidth
                   ? parameters.forcingAmplitude
                   : 0.0;
    case ForcingProfile::gaussian:
        return parameters.forcingAmplitude *
               std::exp(
                   -0.5 *
                   std::pow((k - parameters.forcingWavenumber) / parameters.forcingWidth, 2.0));
    case ForcingProfile::exponential: {
        const double logRatio =
            parameters.forcingShapeOrder * std::log(k / parameters.forcingWavenumber);
        if (!std::isfinite(logRatio) || logRatio > std::log(std::numeric_limits<double>::max()))
            return 0.0;
        const double ratio = std::exp(logRatio);
        return parameters.forcingAmplitude * std::exp(logRatio - ratio);
    }
    case ForcingProfile::logNormal: {
        const double logRatio = std::log(k / parameters.forcingWavenumber);
        return parameters.forcingAmplitude *
               std::exp(-0.5 * std::pow(logRatio / parameters.forcingLogWidth, 2.0));
    }
    case ForcingProfile::singleMode: {
        const long selectedMode =
            std::lround(parameters.forcingWavenumber * parameters.domainLength / (2.0 * sh1dPi));
        return std::abs(signedWave(index, parameters.gridPoints)) == selectedMode
                   ? parameters.forcingAmplitude
                   : 0.0;
    }
    }
    throw std::logic_error("unknown forcing profile");
}
} // namespace

Solver::Solver(Parameters parameters)
    : parameters_(std::move(parameters)), baseTransform_(parameters_.gridPoints),
      nonlinearOperator_(parameters_), linearOperator_(makeSpectralField(parameters_.gridPoints)),
      noise_(parameters_.usesStochasticForcing() ? makeSpectralField(parameters_.gridPoints)
                                                 : SpectralField{}),
      deterministicForcing_(parameters_.usesDeterministicForcing()
                                ? makeSpectralField(parameters_.gridPoints)
                                : SpectralField{}),
      forcingAmplitude_(parameters_.gridPoints),
      stochasticNoiseScale_(parameters_.usesStochasticForcing()
                                ? std::vector<double>(parameters_.gridPoints)
                                : std::vector<double>{}),
      random_(parameters_.randomSeed = resolveRandomSeed(parameters_.randomSeed)) {
    const std::size_t stages = parameters_.integrator == Integrator::etd4 ? 4 : 2;
    for (std::size_t i = 0; i < stages; ++i)
        nonlinearStages_[i] = makeSpectralField(parameters_.gridPoints);
    for (std::size_t i = 1; i < stages; ++i)
        stageStates_[i - 1] = makeSpectralField(parameters_.gridPoints);
    buildLinearOperator();
    buildIntegrationCoefficients();
    if (parameters_.forcingEnabled)
        buildForcing();
}

void Solver::buildLinearOperator() {
    for (std::size_t i = 0; i < parameters_.gridPoints; ++i) {
        const double k = waveNumber(parameters_, i), k2 = k * k;
        const double hamiltonian =
            -parameters_.dispersionCoefficient * k2 + parameters_.chemicalPotential;
        Complex value(0.0, -hamiltonian);
        if (parameters_.hyperviscosity > 0.0)
            value -= parameters_.hyperviscosity * std::pow(k2, parameters_.hyperviscosityOrder);
        if (parameters_.hypoviscosity > 0.0 && (k2 > 0.0 || parameters_.hypoviscosityOrder >= 0.0))
            value -= parameters_.hypoviscosity * std::pow(k2, parameters_.hypoviscosityOrder);
        linearOperator_[i] = value;
    }
    if (parameters_.hypoviscosity > 0.0 && parameters_.hypoviscosityOrder < 0.0)
        linearOperator_[0] = Complex{};
    requireFinite(linearOperator_, "linear operator; reduce dissipation coefficients or orders");
}

void Solver::buildIntegrationCoefficients() {
    const auto phi = [](Complex z, int order) {
        double factorial = 1.0;
        for (int k = 2; k <= order; ++k)
            factorial *= static_cast<double>(k);
        if (std::abs(z) < 2.0) {
            Complex term = 1.0 / factorial;
            Complex sum = term;
            for (int k = 1; k < 96; ++k) {
                term *= z / static_cast<double>(order + k);
                sum += term;
                if (std::abs(term) < 1.e-17 * std::max(1.0, std::abs(sum)))
                    break;
            }
            return sum;
        }
        Complex value = (std::exp(z) - 1.0) / z;
        double previousFactorial = 1.0;
        for (int k = 2; k <= order; ++k) {
            value = (value - 1.0 / previousFactorial) / z;
            previousFactorial *= static_cast<double>(k);
        }
        return value;
    };

    auto &coefficients = coefficients_;
    coefficients.e1 = makeSpectralField(parameters_.gridPoints);
    if (parameters_.integrator != Integrator::integratingFactorRk2) {
        coefficients.q1 = makeSpectralField(parameters_.gridPoints);
        coefficients.f1 = makeSpectralField(parameters_.gridPoints);
    }
    if (parameters_.integrator == Integrator::etd4) {
        coefficients.e2 = makeSpectralField(parameters_.gridPoints);
        coefficients.q2 = makeSpectralField(parameters_.gridPoints);
        coefficients.q3 = makeSpectralField(parameters_.gridPoints);
        coefficients.q4 = makeSpectralField(parameters_.gridPoints);
        coefficients.q5 = makeSpectralField(parameters_.gridPoints);
        coefficients.f2 = makeSpectralField(parameters_.gridPoints);
        coefficients.f3 = makeSpectralField(parameters_.gridPoints);
    }

    const double h = parameters_.timeStep;
    for (std::size_t i = 0; i < parameters_.gridPoints; ++i) {
        const Complex z = h * linearOperator_[i];
        coefficients.e1[i] = std::exp(z);
        if (parameters_.integrator == Integrator::etd2) {
            coefficients.q1[i] = h * phi(z, 1);
            coefficients.f1[i] = h * phi(z, 2);
        } else if (parameters_.integrator == Integrator::etd4) {
            const Complex p1 = phi(z, 1), p2 = phi(z, 2), p3 = phi(z, 3);
            const Complex halfP1 = phi(0.5 * z, 1);
            const Complex halfP2 = phi(0.5 * z, 2);
            coefficients.e2[i] = std::exp(0.5 * z);
            coefficients.q1[i] = 0.5 * h * halfP1;
            coefficients.q2[i] = h * (0.5 * halfP1 - halfP2);
            coefficients.q3[i] = h * halfP2;
            coefficients.q4[i] = h * (p1 - 2.0 * p2);
            coefficients.q5[i] = 2.0 * h * p2;
            coefficients.f1[i] = h * (p1 - 3.0 * p2 + 4.0 * p3);
            coefficients.f2[i] = h * (p2 - 2.0 * p3);
            coefficients.f3[i] = h * (-p2 + 4.0 * p3);
        }
        if (!stochasticNoiseScale_.empty()) {
            const double a = linearOperator_[i].real();
            const double ah = a * h;
            const double variance =
                std::abs(ah) < 1.e-8 ? h * (1.0 + ah) : std::expm1(2.0 * ah) / (2.0 * a);
            if (!(variance >= 0.0) || !std::isfinite(variance))
                throw std::runtime_error("stochastic integration coefficient is non-finite");
            stochasticNoiseScale_[i] = std::sqrt(variance);
        }
    }
    for (const SpectralField *values : coefficients.fields())
        requireFinite(*values, "time-integration coefficients");
}

void Solver::buildForcing() {
    for (std::size_t i = 0; i < parameters_.gridPoints; ++i) {
        const double amplitude = forcingEnvelope(parameters_, i);
        if (!std::isfinite(amplitude))
            throw std::runtime_error("forcing profile contains a non-finite value");
        forcingAmplitude_[i] = amplitude;
        if (amplitude != 0.0) {
            forcedIndices_.push_back(i);
            ++forcedModeCount_;
            waveActionInjectionCoefficient_ += amplitude * amplitude;
            const double signedK = waveNumber(parameters_, i);
            quadraticEnergyInjectionCoefficient_ +=
                (-parameters_.dispersionCoefficient * signedK * signedK +
                 parameters_.chemicalPotential) *
                amplitude * amplitude;
        }
    }
    if (forcedModeCount_ == 0)
        throw std::runtime_error(
            "forcingEnabled is true, but the selected profile contains no modes");
    if (!(waveActionInjectionCoefficient_ > 0.0) ||
        !std::isfinite(waveActionInjectionCoefficient_) ||
        !std::isfinite(quadraticEnergyInjectionCoefficient_))
        throw std::runtime_error("forcing amplitude produces invalid injection coefficients");
    if (parameters_.targetWaveActionInjectionRate > 0.0) {
        const double scale =
            std::sqrt(parameters_.targetWaveActionInjectionRate / waveActionInjectionCoefficient_);
        if (!(scale > 0.0) || !std::isfinite(scale))
            throw std::runtime_error("targetWaveActionInjectionRate cannot be represented");
        for (double &value : forcingAmplitude_)
            value *= scale;
        waveActionInjectionCoefficient_ = parameters_.targetWaveActionInjectionRate;
        quadraticEnergyInjectionCoefficient_ *= scale * scale;
        if (!std::isfinite(quadraticEnergyInjectionCoefficient_))
            throw std::runtime_error("normalized forcing has a non-finite energy coefficient");
    }
    if (parameters_.usesDeterministicForcing())
        for (std::size_t i = 0; i < parameters_.gridPoints; ++i)
            deterministicForcing_[i] = forcingAmplitude_[i];
}

void Solver::generateNoise() {
    std::fill(noise_.begin(), noise_.end(), Complex{});
    constexpr double circularScale = 0.7071067811865475244;
    for (const std::size_t i : forcedIndices_)
        noise_[i] = forcingAmplitude_[i] * circularScale * stochasticNoiseScale_[i] *
                    Complex(normal_(random_), normal_(random_));
}

void Solver::rightHandSide(const SpectralField &input, SpectralField &output) {
    nonlinearOperator_.evaluate(input, output);
    if (parameters_.usesDeterministicForcing())
        for (std::size_t i = 0; i < output.size(); ++i)
            output[i] += deterministicForcing_[i];
}

void Solver::step(SpectralField &wavefunction) {
    const bool stochastic = parameters_.usesStochasticForcing();
    if (stochastic)
        generateNoise();

    SpectralField &stageA = stageStates_[0];
    SpectralField &stageB = stageStates_[1];
    SpectralField &stageC = stageStates_[2];
    SpectralField &nonlinear1 = nonlinearStages_[0];
    SpectralField &nonlinear2 = nonlinearStages_[1];
    SpectralField &nonlinear3 = nonlinearStages_[2];
    SpectralField &nonlinear4 = nonlinearStages_[3];

    rightHandSide(wavefunction, nonlinear1);
    for (std::size_t i = 0; i < wavefunction.size(); ++i)
        stageA[i] = integrationStageA(parameters_.integrator, parameters_.timeStep, i,
                                      coefficients_, wavefunction[i], nonlinear1[i]);
    rightHandSide(stageA, nonlinear2);
    if (parameters_.integrator == Integrator::etd4) {
        for (std::size_t i = 0; i < wavefunction.size(); ++i)
            stageB[i] =
                integrationStageB(i, coefficients_, wavefunction[i], nonlinear1[i], nonlinear2[i]);
        rightHandSide(stageB, nonlinear3);
        for (std::size_t i = 0; i < wavefunction.size(); ++i)
            stageC[i] =
                integrationStageC(i, coefficients_, wavefunction[i], nonlinear1[i], nonlinear3[i]);
        rightHandSide(stageC, nonlinear4);
    }
    for (std::size_t i = 0; i < wavefunction.size(); ++i) {
        const IntegrationFinishValues values{
            .initial = wavefunction[i],
            .stageA = stageA[i],
            .nonlinear1 = nonlinear1[i],
            .nonlinear2 = nonlinear2[i],
            .nonlinear3 = parameters_.integrator == Integrator::etd4 ? nonlinear3[i] : Complex{},
            .nonlinear4 = parameters_.integrator == Integrator::etd4 ? nonlinear4[i] : Complex{}};
        wavefunction[i] = integrationFinish(parameters_.integrator, parameters_.timeStep, i,
                                            coefficients_, values);
        if (stochastic)
            wavefunction[i] += noise_[i];
    }
    enforceStateConstraints(parameters_, wavefunction);
    requireFinite(wavefunction, "wavefunction");
}

void Solver::writeState(const RestartState &state) {
    std::ostringstream randomState, distributionState;
    randomState << random_;
    distributionState << normal_;
    writeWavefunctionAndRestart(parameters_, baseTransform_, state, randomState.str(),
                                distributionState.str());
}

void Solver::validateRunBounds(const RestartState &state) const {
    const long double finalTime =
        static_cast<long double>(state.time) +
        static_cast<long double>(parameters_.timeStep) * parameters_.numberOfSteps;
    if (finalTime > static_cast<long double>(std::numeric_limits<double>::max()))
        throw std::runtime_error("requested run would overflow simulation time");

    const std::uint64_t scheduledOutputs =
        parameters_.numberOfSteps / parameters_.outputIntervalSteps +
        (parameters_.numberOfSteps % parameters_.outputIntervalSteps != 0 ? 1 : 0);
    if (state.frame > std::numeric_limits<std::uint64_t>::max() - scheduledOutputs)
        throw std::runtime_error("requested run would overflow output frame count");
}

void Solver::restoreRandomState(RestartState &state) {
    if (!state.restarting) {
        state.randomSeed = parameters_.randomSeed;
        return;
    }

    std::istringstream engineState(state.randomEngineState);
    if (!(engineState >> random_))
        throw std::runtime_error("cannot restore random-generator state");
    if (!state.normalDistributionState.empty()) {
        std::istringstream distributionState(state.normalDistributionState);
        if (!(distributionState >> normal_))
            throw std::runtime_error("cannot restore normal-distribution state");
    } else {
        normal_.reset();
    }
    parameters_.randomSeed = state.randomSeed;
}

RestartState Solver::prepareRun() {
    std::filesystem::create_directories(parameters_.dataDirectory);
    std::filesystem::create_directories(parameters_.outputDirectory);
    RestartState state = readRestartOrInitial(parameters_, baseTransform_);
    validateRunBounds(state);
    restoreRandomState(state);
    prepareOutput(parameters_, state.restarting, state.frame);
    writeRunRecords(parameters_, state.time, state.frame, forcingAmplitude_, forcedModeCount_,
                    waveActionInjectionCoefficient_, quadraticEnergyInjectionCoefficient_);
    if (!state.restarting)
        writeState(state);

    std::cout << "backend = CPU/FFTW\nmodel = " << modelName(parameters_.model)
              << "\ngridPoints = " << parameters_.gridPoints
              << " timeStep = " << parameters_.timeStep
              << " numberOfSteps = " << parameters_.numberOfSteps
              << " outputIntervalSteps = " << parameters_.outputIntervalSteps
              << "\nrandom seed = " << parameters_.randomSeed << '\n';
    return state;
}

double Solver::writeOutputFrame(const RestartState &state, DiagnosticsAverages &averages) {
    const double energy = appendDiagnostics(parameters_, nonlinearOperator_, state.time,
                                            state.frame, state.wavefunction, averages);
    writeState(state);
    return energy;
}

void Solver::run() {
    RestartState state = prepareRun();
    DiagnosticsAverages averages;
    const auto start = std::chrono::steady_clock::now();
    for (std::uint64_t stepIndex = 0; stepIndex < parameters_.numberOfSteps; ++stepIndex) {
        const std::uint64_t stepNumber = stepIndex + 1;
        step(state.wavefunction);
        state.time += parameters_.timeStep;
        if (stepNumber % parameters_.outputIntervalSteps == 0 ||
            stepNumber == parameters_.numberOfSteps) {
            ++state.frame;
            const double energy = writeOutputFrame(state, averages);
            std::cout << "time = " << state.time << " file = " << state.frame
                      << " Energy = " << energy << '\n';
        }
    }
    const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
    std::cout << "time taken for code is = " << elapsed.count() << '\n';
}
