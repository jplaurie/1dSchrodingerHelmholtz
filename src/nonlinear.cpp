#include "nonlinear.hpp"

#include "spectral.hpp"

#include <cmath>
#include <stdexcept>

NonlinearOperator::NonlinearOperator(const Parameters &p)
    : p_(p), paddedTransform_(p.paddedGridPoints()),
      paddedWave_(p.paddedGridPoints()), physicalWave_(p.paddedGridPoints()),
      density_(p.paddedGridPoints()), densityHat_(p.paddedGridPoints()),
      potential_(p.paddedGridPoints()),
      nonlinearPhysical_(p.paddedGridPoints()),
      nonlinearHat_(p.paddedGridPoints()),
      paddedDerivative_(p.paddedGridPoints()),
      physicalDerivative_(p.paddedGridPoints()),
      densityResponse_(p.paddedGridPoints()) {
  const auto cutoff = static_cast<long>(p_.gridPoints / 2);
  const auto count = p_.paddedGridPoints();
  for (std::size_t i = 0; i < count; ++i) {
    const long mode = signedWave(i, count);
    if (std::abs(mode) >= cutoff)
      continue;
    const double k = 2.0 * sh1dPi * static_cast<double>(mode) / p_.domainLength;
    densityResponse_[i] = densityResponseMultiplier(p_, k * k);
  }
}

void NonlinearOperator::applyDensityResponse(SpectralField &density) const {
  for (std::size_t i = 0; i < density.size(); ++i)
    density[i] *= densityResponse_[i];
}

void NonlinearOperator::evaluate(const SpectralField &wavefunction,
                                 SpectralField &result) {
  if (wavefunction.size() != p_.gridPoints)
    throw std::runtime_error("invalid nonlinear input size");
  embedBaseSpectrum(p_, wavefunction, paddedWave_);
  paddedTransform_.inverse(paddedWave_, physicalWave_);
  for (std::size_t i = 0; i < density_.size(); ++i)
    density_[i] = std::norm(physicalWave_[i]);
  paddedTransform_.forward(density_, densityHat_);
  applyDensityResponse(densityHat_);
  paddedTransform_.inverse(densityHat_, potential_);
  for (std::size_t i = 0; i < nonlinearPhysical_.size(); ++i)
    nonlinearPhysical_[i] = potential_[i] * physicalWave_[i];
  paddedTransform_.forward(nonlinearPhysical_, nonlinearHat_);
  result.resize(p_.gridPoints);
  extractBaseSpectrum(p_, nonlinearHat_, result);
  const Complex factor(0.0, -p_.nonlinearityCoefficient);
  for (Complex &value : result)
    value *= factor;
}

void NonlinearOperator::densitySpectrum(const SpectralField &wavefunction,
                                        SpectralField &result) {
  if (wavefunction.size() != p_.gridPoints)
    throw std::runtime_error("invalid density input size");
  embedBaseSpectrum(p_, wavefunction, paddedWave_);
  paddedTransform_.inverse(paddedWave_, physicalWave_);
  for (std::size_t i = 0; i < density_.size(); ++i)
    density_[i] = std::norm(physicalWave_[i]);
  paddedTransform_.forward(density_, densityHat_);
  filterPaddedSpectrum(p_, densityHat_);
  result.resize(p_.gridPoints);
  extractBaseSpectrum(p_, densityHat_, result);
}

void NonlinearOperator::densityTimeDerivativeSpectrum(
    const SpectralField &wavefunction, const SpectralField &wavefunctionDot,
    SpectralField &result) {
  if (wavefunction.size() != p_.gridPoints ||
      wavefunctionDot.size() != p_.gridPoints)
    throw std::runtime_error("invalid density-derivative input size");
  embedBaseSpectrum(p_, wavefunction, paddedWave_);
  embedBaseSpectrum(p_, wavefunctionDot, paddedDerivative_);
  paddedTransform_.inverse(paddedWave_, physicalWave_);
  paddedTransform_.inverse(paddedDerivative_, physicalDerivative_);
  for (std::size_t i = 0; i < density_.size(); ++i)
    density_[i] =
        2.0 * std::real(std::conj(physicalWave_[i]) * physicalDerivative_[i]);
  paddedTransform_.forward(density_, densityHat_);
  filterPaddedSpectrum(p_, densityHat_);
  result.resize(p_.gridPoints);
  extractBaseSpectrum(p_, densityHat_, result);
}
