#include "nonlinear.hpp"

#include "spectral.hpp"

#include <cmath>
#include <stdexcept>

NonlinearOperator::NonlinearOperator(const Parameters &parameters)
    : parameters_(parameters), paddedWavefunction_(parameters.paddedGridPoints()),
      physicalWavefunction_(parameters.paddedGridPoints()),
      densitySpectrum_(parameters.paddedGridPoints()), potential_(parameters.paddedGridPoints()),
      physicalWavefunctionDerivative_(parameters.paddedGridPoints()),
      densityResponse_(parameters.paddedGridPoints()) {
    inverseWavefunction_ = makeFftwPlan(paddedWavefunction_, physicalWavefunction_, FFTW_BACKWARD);
    forwardDensity_ = makeFftwPlan(potential_, densitySpectrum_, FFTW_FORWARD);
    inverseDensity_ = makeFftwPlan(densitySpectrum_, potential_, FFTW_BACKWARD);
    forwardNonlinear_ = makeFftwPlan(paddedWavefunction_, densitySpectrum_, FFTW_FORWARD);
    inverseDerivative_ =
        makeFftwPlan(paddedWavefunction_, physicalWavefunctionDerivative_, FFTW_BACKWARD);
    const auto cutoff = static_cast<long>(parameters_.gridPoints / 2);
    const auto count = parameters_.paddedGridPoints();
    for (std::size_t i = 0; i < count; ++i) {
        const long mode = signedWave(i, count);
        if (std::abs(mode) >= cutoff)
            continue;
        const double k = 2.0 * sh1dPi * static_cast<double>(mode) / parameters_.domainLength;
        densityResponse_[i] =
            densityResponseMultiplier(parameters_, k * k) / static_cast<double>(count);
    }
}

void NonlinearOperator::transformDensity(const SpectralField &wavefunction) {
    embedBaseSpectrum(parameters_, wavefunction, paddedWavefunction_);
    fftw_execute(inverseWavefunction_.get());
    for (std::size_t i = 0; i < potential_.size(); ++i)
        potential_[i] = std::norm(physicalWavefunction_[i]);
    fftw_execute(forwardDensity_.get());
}

void NonlinearOperator::extractDensity(SpectralField &result) const {
    result.resize(parameters_.gridPoints);
    extractBaseSpectrum(parameters_, densitySpectrum_, result);
    const double scale = 1.0 / static_cast<double>(parameters_.paddedGridPoints());
    for (Complex &value : result)
        value *= scale;
    // Density has its own real, symmetric projector, independent of whether
    // the wavefunction mean is suppressed.
    result[parameters_.gridPoints / 2] = Complex{};
}

void NonlinearOperator::evaluate(const SpectralField &wavefunction, SpectralField &result) {
    if (wavefunction.size() != parameters_.gridPoints)
        throw std::runtime_error("invalid nonlinear input size");
    if (parameters_.nonlinearityCoefficient == 0.0) {
        result.assign(parameters_.gridPoints, Complex{});
        return;
    }
    transformDensity(wavefunction);
    for (std::size_t i = 0; i < densitySpectrum_.size(); ++i)
        densitySpectrum_[i] *= densityResponse_[i];
    fftw_execute(inverseDensity_.get());
    for (std::size_t i = 0; i < paddedWavefunction_.size(); ++i)
        paddedWavefunction_[i] = potential_[i] * physicalWavefunction_[i];
    fftw_execute(forwardNonlinear_.get());
    result.resize(parameters_.gridPoints);
    extractBaseSpectrum(parameters_, densitySpectrum_, result);
    const Complex factor(0.0, -parameters_.nonlinearityCoefficient /
                                  static_cast<double>(parameters_.paddedGridPoints()));
    for (Complex &value : result)
        value *= factor;
    enforceStateConstraints(parameters_, result);
}

void NonlinearOperator::densitySpectrum(const SpectralField &wavefunction, SpectralField &result) {
    if (wavefunction.size() != parameters_.gridPoints)
        throw std::runtime_error("invalid density input size");
    transformDensity(wavefunction);
    extractDensity(result);
}

void NonlinearOperator::densityTimeDerivativeSpectrum(const SpectralField &wavefunction,
                                                      const SpectralField &wavefunctionDot,
                                                      SpectralField &result) {
    if (wavefunction.size() != parameters_.gridPoints ||
        wavefunctionDot.size() != parameters_.gridPoints)
        throw std::runtime_error("invalid density-derivative input size");
    embedBaseSpectrum(parameters_, wavefunction, paddedWavefunction_);
    fftw_execute(inverseWavefunction_.get());
    embedBaseSpectrum(parameters_, wavefunctionDot, paddedWavefunction_);
    fftw_execute(inverseDerivative_.get());
    for (std::size_t i = 0; i < potential_.size(); ++i)
        potential_[i] = 2.0 * std::real(std::conj(physicalWavefunction_[i]) *
                                        physicalWavefunctionDerivative_[i]);
    fftw_execute(forwardDensity_.get());
    extractDensity(result);
}
