#include "fftw_utils.hpp"
#include "nonlinear.hpp"
#include "parameters.hpp"
#include "spectral.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>

namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}

void testPlaneWave(Model model) {
    Parameters p;
    p.gridPoints = 32;
    p.model = model;
    p.nonlinearityCoefficient = -1.75;
    p.helmholtzParameter = 2.0;
    p.threadCount = 1;
    validateParameters(p);
    NonlinearOperator nonlinear(p);
    SpectralField state(p.gridPoints), result, density;
    const Complex amplitude(0.7, -0.2);
    state[3] = amplitude;
    nonlinear.evaluate(state, result);
    const Complex expected =
        Complex(0.0, -p.nonlinearityCoefficient) * std::norm(amplitude) * amplitude;
    require(std::abs(result[3] - expected) < 2.e-12,
            "plane-wave nonlinear coefficient is incorrect");
    for (std::size_t i = 0; i < result.size(); ++i)
        if (i != 3)
            require(std::abs(result[i]) < 2.e-12, "plane wave generated a spurious nonlinear mode");
    nonlinear.densitySpectrum(state, density);
    require(std::abs(density[0] - Complex(std::norm(amplitude), 0.0)) < 2.e-12,
            "plane-wave density is incorrect");
}

void testSpectralEmbedding() {
    Parameters p;
    p.gridPoints = 16;
    SpectralField base(p.gridPoints), padded(p.paddedGridPoints()), recovered(16);
    for (std::size_t i = 0; i < base.size(); ++i)
        base[i] = Complex(static_cast<double>(i), -static_cast<double>(i));
    base[p.gridPoints / 2] = Complex{};
    embedBaseSpectrum(p, base, padded);
    extractBaseSpectrum(p, padded, recovered);
    require(base == recovered, "spectral embed/extract round trip failed");
}

void testDirectConvolution(Model model) {
    Parameters p;
    p.gridPoints = 32;
    p.model = model;
    p.nonlinearityCoefficient = -0.8;
    p.helmholtzParameter = 0.35;
    validateParameters(p);
    NonlinearOperator nonlinear(p);
    SpectralField state(p.gridPoints), computed;
    const std::map<long, Complex> input{{-2, {0.12, -0.04}},
                                        {-1, {-0.08, 0.07}},
                                        {0, {0.2, 0.03}},
                                        {1, {0.05, -0.11}},
                                        {2, {-0.03, 0.09}}};
    for (const auto &[mode, value] : input) {
        const auto index = mode >= 0
                               ? static_cast<std::size_t>(mode)
                               : static_cast<std::size_t>(static_cast<long>(p.gridPoints) + mode);
        state[index] = value;
    }
    nonlinear.evaluate(state, computed);

    std::map<long, Complex> density;
    for (long q = -4; q <= 4; ++q) {
        Complex value{};
        for (const auto &[mode, amplitude] : input) {
            const auto other = input.find(mode - q);
            if (other != input.end())
                value += amplitude * std::conj(other->second);
        }
        const double k = 2.0 * sh1dPi * static_cast<double>(q) / p.domainLength;
        double multiplier = 1.0;
        if (model == Model::schrodingerHelmholtz)
            multiplier = 1.0 / (1.0 + p.helmholtzParameter * k * k);
        else if (model == Model::longWave)
            multiplier = 1.0 - p.helmholtzParameter * k * k;
        density[q] = multiplier * value;
    }
    for (long mode = -7; mode <= 7; ++mode) {
        Complex product{};
        for (const auto &[waveMode, amplitude] : input) {
            const auto potential = density.find(mode - waveMode);
            if (potential != density.end())
                product += amplitude * potential->second;
        }
        const Complex expected = Complex(0.0, -p.nonlinearityCoefficient) * product;
        const auto index = mode >= 0
                               ? static_cast<std::size_t>(mode)
                               : static_cast<std::size_t>(static_cast<long>(p.gridPoints) + mode);
        require(std::abs(computed[index] - expected) < 3.e-12,
                "dealiased nonlinear term disagrees with direct convolution");
    }
}

void testDensityResponse() {
    Parameters p;
    p.helmholtzParameter = 0.5;
    p.model = Model::schrodingerHelmholtz;
    require(std::abs(densityResponseMultiplier(p, 4.0) - 1.0 / 3.0) < 1.e-15,
            "Schrodinger-Helmholtz density response is incorrect");
    p.model = Model::longWave;
    require(std::abs(densityResponseMultiplier(p, 4.0) + 1.0) < 1.e-15,
            "long-wave density response is incorrect");
    p.model = Model::nonlinearSchrodinger;
    require(densityResponseMultiplier(p, 4.0) == 1.0, "NLS density response is incorrect");
}

void testCutoffConvolution(Model model, bool removeMean, FftwPlanning planning) {
    Parameters p;
    p.gridPoints = 12;
    p.model = model;
    p.fftwPlanning = planning;
    p.nonlinearityCoefficient = 1.3;
    p.helmholtzParameter = 0.35;
    p.hypoviscosity = removeMean ? 0.1 : 0.0;
    configureFftw(p, false);
    NonlinearOperator nonlinear(p);
    SpectralField w(p.gridPoints), actual, density, densityDot;
    const long cutoff = static_cast<long>(p.gridPoints / 2);
    const auto index = [&](long k) {
        return static_cast<std::size_t>(k < 0 ? k + p.gridPoints : k);
    };
    for (std::size_t i = 0; i < w.size(); ++i)
        w[i] = {0.22 * std::cos(0.71 * i), 0.17 * std::sin(0.39 * i)};
    enforceStateConstraints(p, w);
    std::map<long, Complex> rho;
    for (long d = -cutoff + 1; d < cutoff; ++d)
        for (long j = -cutoff + 1; j < cutoff; ++j)
            if (j + d > -cutoff && j + d < cutoff)
                rho[d] += w[index(j + d)] * std::conj(w[index(j)]);
    nonlinear.evaluate(w, actual);
    nonlinear.densitySpectrum(w, density);
    SpectralField rate = actual;
    double massRate = 0.0, momentumRate = 0.0, energyRate = 0.0;
    for (std::size_t i = 0; i < w.size(); ++i) {
        const long k = signedWave(i, w.size());
        Complex expected{};
        if (retainedStateMode(p, i))
            for (const auto &[d, value] : rho)
                if (k - d > -cutoff && k - d < cutoff)
                    expected += Complex(0.0, -p.nonlinearityCoefficient) *
                                densityResponseMultiplier(p, static_cast<double>(d * d)) * value *
                                w[index(k - d)];
        require(std::abs(actual[i] - expected) < 2.e-12,
                "cutoff nonlinear convolution or wavefunction projector is incorrect");
        require(std::abs(density[i] - (k == -cutoff ? Complex{} : rho[k])) < 2.e-12,
                "cutoff density convolution is incorrect");
        const double weight = 0.5 * k * k + 0.2;
        rate[i] += Complex(0.0, -weight) * w[i];
        const double transfer = 2.0 * std::real(std::conj(w[i]) * rate[i]);
        massRate += transfer;
        momentumRate += k * transfer;
        energyRate += weight * transfer;
    }
    nonlinear.densityTimeDerivativeSpectrum(w, rate, densityDot);
    for (std::size_t i = 0; i < w.size(); ++i) {
        const double k = waveNumber(p, i);
        energyRate += p.nonlinearityCoefficient * densityResponseMultiplier(p, k * k) *
                      std::real(std::conj(density[i]) * densityDot[i]);
    }
    require(std::max({std::abs(massRate), std::abs(momentumRate), std::abs(energyRate)}) < 2.e-12,
            "projected mass, momentum, or Hamiltonian transfer does not close");
    // Repeatedly switching diagnostic/evolution calls exercises shared FFT buffers.
    nonlinear.evaluate(w, densityDot);
    for (std::size_t i = 0; i < w.size(); ++i)
        require(std::abs(densityDot[i] - actual[i]) < 2.e-12,
                "diagnostics corrupted nonlinear workspace");
}

void testZeroCoupling() {
    Parameters p;
    p.gridPoints = 8;
    p.nonlinearityCoefficient = 0.0;
    NonlinearOperator nonlinear(p);
    SpectralField w(8, Complex(1.0, 0.2)), result(8, Complex(9.0, 0.0));
    nonlinear.evaluate(w, result);
    require(std::all_of(result.begin(), result.end(), [](Complex v) { return v == Complex{}; }),
            "zero coupling did not clear the nonlinear result");
    w.assign(8, Complex{});
    w[1] = 0.3;
    nonlinear.densitySpectrum(w, result);
    require(std::abs(result[0] - 0.09) < 1.e-14,
            "zero coupling incorrectly disabled density diagnostics");
}

} // namespace

int main(int argc, char **argv) {
    try {
        initializeFftwThreads(argc == 2 ? std::stoi(argv[1]) : 1);
        testPlaneWave(Model::schrodingerHelmholtz);
        testPlaneWave(Model::longWave);
        testPlaneWave(Model::nonlinearSchrodinger);
        testSpectralEmbedding();
        testDirectConvolution(Model::schrodingerHelmholtz);
        testDirectConvolution(Model::longWave);
        testDirectConvolution(Model::nonlinearSchrodinger);
        testDensityResponse();
        for (auto planning : {FftwPlanning::estimate, FftwPlanning::measure})
            for (auto model :
                 {Model::schrodingerHelmholtz, Model::longWave, Model::nonlinearSchrodinger})
                for (bool removeMean : {false, true})
                    testCutoffConvolution(model, removeMean, planning);
        testZeroCoupling();
        finalizeFftwThreads();
        std::cout << "numerical tests passed\n";
        return 0;
    } catch (const std::exception &error) {
        finalizeFftwThreads();
        std::cerr << error.what() << '\n';
        return 1;
    }
}
