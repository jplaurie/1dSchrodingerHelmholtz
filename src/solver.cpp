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
SpectralField field(std::size_t count) { return SpectralField(count); }

void requireFinite(const SpectralField &values, const char *description) {
  if (!std::all_of(values.begin(), values.end(), [](Complex value) {
        return std::isfinite(value.real()) && std::isfinite(value.imag());
      }))
    throw std::runtime_error(std::string(description) +
                             " contains a non-finite coefficient");
}

std::uint64_t resolveRandomSeed(std::uint64_t seed) {
  if (seed != 0)
    return seed;
  return static_cast<std::uint64_t>(
      std::chrono::high_resolution_clock::now().time_since_epoch().count());
}

double forcingEnvelope(const Parameters &p, std::size_t index) {
  const double k = std::abs(waveNumber(p, index));
  if (k == 0.0)
    return 0.0;
  switch (p.forcingProfile) {
  case ForcingProfile::annulus:
    return std::abs(k - p.forcingWavenumber) < p.forcingWidth
               ? p.forcingAmplitude
               : 0.0;
  case ForcingProfile::gaussian:
    return p.forcingAmplitude *
           std::exp(-0.5 *
                    std::pow((k - p.forcingWavenumber) / p.forcingWidth, 2.0));
  case ForcingProfile::exponential: {
    const double logRatio =
        p.forcingShapeOrder * std::log(k / p.forcingWavenumber);
    if (!std::isfinite(logRatio) ||
        logRatio > std::log(std::numeric_limits<double>::max()))
      return 0.0;
    const double ratio = std::exp(logRatio);
    return p.forcingAmplitude * std::exp(logRatio - ratio);
  }
  case ForcingProfile::logNormal: {
    const double logRatio = std::log(k / p.forcingWavenumber);
    return p.forcingAmplitude *
           std::exp(-0.5 * std::pow(logRatio / p.forcingLogWidth, 2.0));
  }
  case ForcingProfile::singleMode: {
    const long selectedMode =
        std::lround(p.forcingWavenumber * p.domainLength / (2.0 * sh1dPi));
    return std::abs(signedWave(index, p.gridPoints)) == selectedMode
               ? p.forcingAmplitude
               : 0.0;
  }
  }
  throw std::logic_error("unknown forcing profile");
}
} // namespace

Solver::Solver(Parameters p)
    : p_(std::move(p)), baseTransform_(p_.gridPoints), nonlinear_(p_),
      linear_(field(p_.gridPoints)),
      noise_(p_.usesStochasticForcing() ? field(p_.gridPoints)
                                        : SpectralField{}),
      deterministicForcing_(p_.usesDeterministicForcing() ? field(p_.gridPoints)
                                                          : SpectralField{}),
      forcingAmplitude_(p_.gridPoints),
      noiseScale_(p_.usesStochasticForcing()
                      ? std::vector<double>(p_.gridPoints)
                      : std::vector<double>{}),
      random_(p_.randomSeed = resolveRandomSeed(p_.randomSeed)) {
  const std::size_t stages = p_.integrator == Integrator::etd4 ? 4 : 2;
  for (std::size_t i = 0; i < stages; ++i)
    nonlinearStages_[i] = field(p_.gridPoints);
  for (std::size_t i = 1; i < stages; ++i)
    stageStates_[i - 1] = field(p_.gridPoints);
  buildLinearOperator();
  buildIntegrationCoefficients();
  if (p_.forcingEnabled)
    buildForcing();
}

void Solver::buildLinearOperator() {
  for (std::size_t i = 0; i < p_.gridPoints; ++i) {
    const double k = waveNumber(p_, i), k2 = k * k;
    const double hamiltonian =
        -p_.dispersionCoefficient * k2 + p_.chemicalPotential;
    Complex value(0.0, -hamiltonian);
    if (p_.hyperviscosity > 0.0)
      value -= p_.hyperviscosity * std::pow(k2, p_.hyperviscosityOrder);
    if (p_.hypoviscosity > 0.0 && (k2 > 0.0 || p_.hypoviscosityOrder >= 0.0))
      value -= p_.hypoviscosity * std::pow(k2, p_.hypoviscosityOrder);
    linear_[i] = value;
  }
  if (p_.hypoviscosity > 0.0 && p_.hypoviscosityOrder < 0.0)
    linear_[0] = Complex{};
  requireFinite(linear_,
                "linear operator; reduce dissipation coefficients or orders");
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

  auto &c = coefficients_;
  c.e1 = field(p_.gridPoints);
  if (p_.integrator != Integrator::integratingFactorRk2) {
    c.q1 = field(p_.gridPoints);
    c.f1 = field(p_.gridPoints);
  }
  if (p_.integrator == Integrator::etd4) {
    c.e2 = field(p_.gridPoints);
    c.q2 = field(p_.gridPoints);
    c.q3 = field(p_.gridPoints);
    c.q4 = field(p_.gridPoints);
    c.q5 = field(p_.gridPoints);
    c.f2 = field(p_.gridPoints);
    c.f3 = field(p_.gridPoints);
  }

  const double h = p_.timeStep;
  for (std::size_t i = 0; i < p_.gridPoints; ++i) {
    const Complex z = h * linear_[i];
    c.e1[i] = std::exp(z);
    if (p_.integrator == Integrator::etd2) {
      c.q1[i] = h * phi(z, 1);
      c.f1[i] = h * phi(z, 2);
    } else if (p_.integrator == Integrator::etd4) {
      const Complex p1 = phi(z, 1), p2 = phi(z, 2), p3 = phi(z, 3);
      const Complex halfP1 = phi(0.5 * z, 1);
      const Complex halfP2 = phi(0.5 * z, 2);
      c.e2[i] = std::exp(0.5 * z);
      c.q1[i] = 0.5 * h * halfP1;
      c.q2[i] = h * (0.5 * halfP1 - halfP2);
      c.q3[i] = h * halfP2;
      c.q4[i] = h * (p1 - 2.0 * p2);
      c.q5[i] = 2.0 * h * p2;
      c.f1[i] = h * (p1 - 3.0 * p2 + 4.0 * p3);
      c.f2[i] = h * (p2 - 2.0 * p3);
      c.f3[i] = h * (-p2 + 4.0 * p3);
    }
    if (!noiseScale_.empty()) {
      const double a = linear_[i].real();
      const double ah = a * h;
      const double variance = std::abs(ah) < 1.e-8
                                  ? h * (1.0 + ah)
                                  : std::expm1(2.0 * ah) / (2.0 * a);
      if (!(variance >= 0.0) || !std::isfinite(variance))
        throw std::runtime_error(
            "stochastic integration coefficient is non-finite");
      noiseScale_[i] = std::sqrt(variance);
    }
  }
  for (const SpectralField *values : c.fields())
    requireFinite(*values, "time-integration coefficients");
}

void Solver::buildForcing() {
  for (std::size_t i = 0; i < p_.gridPoints; ++i) {
    const double amplitude = forcingEnvelope(p_, i);
    if (!std::isfinite(amplitude))
      throw std::runtime_error("forcing profile contains a non-finite value");
    forcingAmplitude_[i] = amplitude;
    if (amplitude != 0.0) {
      forcedIndices_.push_back(i);
      ++forcedModeCount_;
      waveActionInjectionCoefficient_ += amplitude * amplitude;
      const double signedK = waveNumber(p_, i);
      quadraticEnergyInjectionCoefficient_ +=
          (-p_.dispersionCoefficient * signedK * signedK +
           p_.chemicalPotential) *
          amplitude * amplitude;
    }
  }
  if (forcedModeCount_ == 0)
    throw std::runtime_error(
        "forcingEnabled is true, but the selected profile contains no modes");
  if (!(waveActionInjectionCoefficient_ > 0.0) ||
      !std::isfinite(waveActionInjectionCoefficient_) ||
      !std::isfinite(quadraticEnergyInjectionCoefficient_))
    throw std::runtime_error(
        "forcing amplitude produces invalid injection coefficients");
  if (p_.targetWaveActionInjectionRate > 0.0) {
    const double scale = std::sqrt(p_.targetWaveActionInjectionRate /
                                   waveActionInjectionCoefficient_);
    if (!(scale > 0.0) || !std::isfinite(scale))
      throw std::runtime_error(
          "targetWaveActionInjectionRate cannot be represented");
    for (double &value : forcingAmplitude_)
      value *= scale;
    waveActionInjectionCoefficient_ = p_.targetWaveActionInjectionRate;
    quadraticEnergyInjectionCoefficient_ *= scale * scale;
    if (!std::isfinite(quadraticEnergyInjectionCoefficient_))
      throw std::runtime_error(
          "normalized forcing has a non-finite energy coefficient");
  }
  if (p_.usesDeterministicForcing())
    for (std::size_t i = 0; i < p_.gridPoints; ++i)
      deterministicForcing_[i] = forcingAmplitude_[i];
}

void Solver::generateNoise() {
  std::fill(noise_.begin(), noise_.end(), Complex{});
  constexpr double circularScale = 0.7071067811865475244;
  for (const std::size_t i : forcedIndices_)
    noise_[i] = forcingAmplitude_[i] * circularScale * noiseScale_[i] *
                Complex(normal_(random_), normal_(random_));
}

void Solver::rightHandSide(const SpectralField &input, SpectralField &output) {
  nonlinear_.evaluate(input, output);
  if (p_.usesDeterministicForcing())
    for (std::size_t i = 0; i < output.size(); ++i)
      output[i] += deterministicForcing_[i];
}

void Solver::step(SpectralField &w) {
  const bool stochastic = p_.usesStochasticForcing();
  if (stochastic)
    generateNoise();
  auto &n = nonlinearStages_;
  auto &stage = stageStates_;
  rightHandSide(w, n[0]);
  for (std::size_t i = 0; i < w.size(); ++i)
    stage[0][i] = integrationStageA(p_.integrator, p_.timeStep, i,
                                    coefficients_, w[i], n[0][i]);
  rightHandSide(stage[0], n[1]);
  if (p_.integrator == Integrator::etd4) {
    for (std::size_t i = 0; i < w.size(); ++i)
      stage[1][i] = integrationStageB(i, coefficients_, w[i], n[0][i], n[1][i]);
    rightHandSide(stage[1], n[2]);
    for (std::size_t i = 0; i < w.size(); ++i)
      stage[2][i] = integrationStageC(i, coefficients_, w[i], n[0][i], n[2][i]);
    rightHandSide(stage[2], n[3]);
  }
  for (std::size_t i = 0; i < w.size(); ++i) {
    const IntegrationFinishValues values{
        .initial = w[i],
        .stageA = stage[0][i],
        .nonlinear1 = n[0][i],
        .nonlinear2 = n[1][i],
        .nonlinear3 = p_.integrator == Integrator::etd4 ? n[2][i] : Complex{},
        .nonlinear4 = p_.integrator == Integrator::etd4 ? n[3][i] : Complex{}};
    w[i] =
        integrationFinish(p_.integrator, p_.timeStep, i, coefficients_, values);
    if (stochastic)
      w[i] += noise_[i];
  }
  enforceStateConstraints(p_, w);
  requireFinite(w, "wavefunction");
}

double Solver::writeFrame(RestartState &state, DiagnosticsAverages *averages) {
  double energy = 0.0;
  if (averages)
    energy = appendDiagnostics(p_, nonlinear_, state.time, state.frame,
                               state.wavefunction, *averages);
  std::ostringstream randomState, distributionState;
  randomState << random_;
  distributionState << normal_;
  writeWavefunctionAndRestart(p_, baseTransform_, state, randomState.str(),
                              distributionState.str());
  return energy;
}

void Solver::run() {
  std::filesystem::create_directories(p_.dataDirectory);
  std::filesystem::create_directories(p_.outputDirectory);
  RestartState state = readRestartOrInitial(p_, baseTransform_);
  if (state.restarting) {
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
    p_.randomSeed = state.randomSeed;
  } else {
    state.randomSeed = p_.randomSeed;
  }
  prepareOutput(p_, state.restarting, state.frame);
  writeRunRecords(p_, state.time, state.frame, forcingAmplitude_,
                  forcedModeCount_, waveActionInjectionCoefficient_,
                  quadraticEnergyInjectionCoefficient_);
  if (!state.restarting)
    writeFrame(state, nullptr);

  std::cout << "backend = CPU/FFTW\nmodel = " << modelName(p_.model)
            << "\ngridPoints = " << p_.gridPoints
            << " timeStep = " << p_.timeStep
            << " numberOfSteps = " << p_.numberOfSteps
            << " outputIntervalSteps = " << p_.outputIntervalSteps
            << "\nrandom seed = " << p_.randomSeed << '\n';

  DiagnosticsAverages averages;
  const auto start = std::chrono::steady_clock::now();
  for (std::uint64_t stepIndex = 0; stepIndex < p_.numberOfSteps; ++stepIndex) {
    const std::uint64_t stepNumber = stepIndex + 1;
    step(state.wavefunction);
    state.time += p_.timeStep;
    if (stepNumber % p_.outputIntervalSteps == 0 ||
        stepNumber == p_.numberOfSteps) {
      ++state.frame;
      const double energy = writeFrame(state, &averages);
      std::cout << "time = " << state.time << " file = " << state.frame
                << " Energy = " << energy << '\n';
    }
  }
  const std::chrono::duration<double> elapsed =
      std::chrono::steady_clock::now() - start;
  std::cout << "time taken for code is = " << elapsed.count() << '\n';
}
