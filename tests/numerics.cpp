#include "fftw_utils.hpp"
#include "nonlinear.hpp"
#include "parameters.hpp"
#include "spectral.hpp"

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

} // namespace

int main() {
    try {
        initializeFftwThreads(1);
        testPlaneWave(Model::schrodingerHelmholtz);
        testPlaneWave(Model::longWave);
        testPlaneWave(Model::nonlinearSchrodinger);
        testSpectralEmbedding();
        testDirectConvolution(Model::schrodingerHelmholtz);
        testDirectConvolution(Model::longWave);
        testDirectConvolution(Model::nonlinearSchrodinger);
        testDensityResponse();
        finalizeFftwThreads();
        std::cout << "numerical tests passed\n";
        return 0;
    } catch (const std::exception &error) {
        finalizeFftwThreads();
        std::cerr << error.what() << '\n';
        return 1;
    }
}
