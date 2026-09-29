#include "nonlinear.hpp"

#include "spectral.hpp"

#include <cmath>
#include <stdexcept>

NonlinearOperator::NonlinearOperator(const Parameters &parameters)
    : parameters_(parameters), paddedTransform_(parameters.paddedGridPoints()),
      paddedWavefunction_(parameters.paddedGridPoints()),
      physicalWavefunction_(parameters.paddedGridPoints()), density_(parameters.paddedGridPoints()),
      densitySpectrum_(parameters.paddedGridPoints()), potential_(parameters.paddedGridPoints()),
      nonlinearPhysical_(parameters.paddedGridPoints()),
      nonlinearSpectrum_(parameters.paddedGridPoints()),
      paddedWavefunctionDerivative_(parameters.paddedGridPoints()),
      physicalWavefunctionDerivative_(parameters.paddedGridPoints()),
      densityResponse_(parameters.paddedGridPoints()) {
    const auto cutoff = static_cast<long>(parameters_.gridPoints / 2);
    const auto count = parameters_.paddedGridPoints();
    for (std::size_t i = 0; i < count; ++i) {
        const long mode = signedWave(i, count);
        if (std::abs(mode) >= cutoff)
            continue;
        const double k = 2.0 * sh1dPi * static_cast<double>(mode) / parameters_.domainLength;
        densityResponse_[i] = densityResponseMultiplier(parameters_, k * k);
    }
}

void NonlinearOperator::applyDensityResponse(SpectralField &density) const {
    for (std::size_t i = 0; i < density.size(); ++i)
        density[i] *= densityResponse_[i];
}

void NonlinearOperator::evaluate(const SpectralField &wavefunction, SpectralField &result) {
    if (wavefunction.size() != parameters_.gridPoints)
        throw std::runtime_error("invalid nonlinear input size");
    embedBaseSpectrum(parameters_, wavefunction, paddedWavefunction_);
    paddedTransform_.inverse(paddedWavefunction_, physicalWavefunction_);
    for (std::size_t i = 0; i < density_.size(); ++i)
        density_[i] = std::norm(physicalWavefunction_[i]);
    paddedTransform_.forward(density_, densitySpectrum_);
    applyDensityResponse(densitySpectrum_);
    paddedTransform_.inverse(densitySpectrum_, potential_);
    for (std::size_t i = 0; i < nonlinearPhysical_.size(); ++i)
        nonlinearPhysical_[i] = potential_[i] * physicalWavefunction_[i];
    paddedTransform_.forward(nonlinearPhysical_, nonlinearSpectrum_);
    result.resize(parameters_.gridPoints);
    extractBaseSpectrum(parameters_, nonlinearSpectrum_, result);
    const Complex factor(0.0, -parameters_.nonlinearityCoefficient);
    for (Complex &value : result)
        value *= factor;
}

void NonlinearOperator::densitySpectrum(const SpectralField &wavefunction, SpectralField &result) {
    if (wavefunction.size() != parameters_.gridPoints)
        throw std::runtime_error("invalid density input size");
    embedBaseSpectrum(parameters_, wavefunction, paddedWavefunction_);
    paddedTransform_.inverse(paddedWavefunction_, physicalWavefunction_);
    for (std::size_t i = 0; i < density_.size(); ++i)
        density_[i] = std::norm(physicalWavefunction_[i]);
    paddedTransform_.forward(density_, densitySpectrum_);
    filterPaddedSpectrum(parameters_, densitySpectrum_);
    result.resize(parameters_.gridPoints);
    extractBaseSpectrum(parameters_, densitySpectrum_, result);
}

void NonlinearOperator::densityTimeDerivativeSpectrum(const SpectralField &wavefunction,
                                                      const SpectralField &wavefunctionDot,
                                                      SpectralField &result) {
    if (wavefunction.size() != parameters_.gridPoints ||
        wavefunctionDot.size() != parameters_.gridPoints)
        throw std::runtime_error("invalid density-derivative input size");
    embedBaseSpectrum(parameters_, wavefunction, paddedWavefunction_);
    embedBaseSpectrum(parameters_, wavefunctionDot, paddedWavefunctionDerivative_);
    paddedTransform_.inverse(paddedWavefunction_, physicalWavefunction_);
    paddedTransform_.inverse(paddedWavefunctionDerivative_, physicalWavefunctionDerivative_);
    for (std::size_t i = 0; i < density_.size(); ++i)
        density_[i] = 2.0 * std::real(std::conj(physicalWavefunction_[i]) *
                                      physicalWavefunctionDerivative_[i]);
    paddedTransform_.forward(density_, densitySpectrum_);
    filterPaddedSpectrum(parameters_, densitySpectrum_);
    result.resize(parameters_.gridPoints);
    extractBaseSpectrum(parameters_, densitySpectrum_, result);
}
