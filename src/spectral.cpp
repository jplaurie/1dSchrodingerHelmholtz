#include "spectral.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

long signedWave(std::size_t index, std::size_t count) {
  return index < (count + 1) / 2
             ? static_cast<long>(index)
             : static_cast<long>(index) - static_cast<long>(count);
}

double waveNumber(const Parameters &p, std::size_t index) {
  return 2.0 * sh1dPi * static_cast<double>(signedWave(index, p.gridPoints)) /
         p.domainLength;
}

double densityResponseMultiplier(const Parameters &p, double k2) {
  switch (p.model) {
  case Model::schrodingerHelmholtz:
    return 1.0 / (1.0 + p.helmholtzParameter * k2);
  case Model::longWave:
    return 1.0 - p.helmholtzParameter * k2;
  case Model::nonlinearSchrodinger:
    return 1.0;
  }
  throw std::logic_error("unknown model in density response");
}

std::size_t paddedIndexForBaseMode(const Parameters &p, std::size_t index) {
  const long mode = signedWave(index, p.gridPoints);
  return mode >= 0 ? static_cast<std::size_t>(mode)
                   : static_cast<std::size_t>(
                         static_cast<long>(p.paddedGridPoints()) + mode);
}

void embedBaseSpectrum(const Parameters &p, const SpectralField &base,
                       SpectralField &padded) {
  if (base.size() != p.gridPoints || padded.size() != p.paddedGridPoints())
    throw std::runtime_error("invalid field size while embedding spectrum");
  std::fill(padded.begin(), padded.end(), Complex{});
  for (std::size_t i = 0; i < base.size(); ++i)
    padded[paddedIndexForBaseMode(p, i)] = base[i];
}

void extractBaseSpectrum(const Parameters &p, const SpectralField &padded,
                         SpectralField &base) {
  if (base.size() != p.gridPoints || padded.size() != p.paddedGridPoints())
    throw std::runtime_error("invalid field size while extracting spectrum");
  for (std::size_t i = 0; i < base.size(); ++i)
    base[i] = padded[paddedIndexForBaseMode(p, i)];
}

void filterPaddedSpectrum(const Parameters &p, SpectralField &field) {
  if (field.size() != p.paddedGridPoints())
    throw std::runtime_error("invalid padded field size");
  const long cutoff = static_cast<long>(p.gridPoints / 2);
  for (std::size_t i = 0; i < field.size(); ++i)
    if (std::abs(signedWave(i, field.size())) >= cutoff)
      field[i] = Complex{};
}

void enforceStateConstraints(const Parameters &p, SpectralField &field) {
  if (field.size() != p.gridPoints)
    throw std::runtime_error("invalid state size");
  field[p.gridPoints / 2] = Complex{};
  if (p.hypoviscosity > 0.0 && p.hypoviscosityOrder < 0.0)
    field[0] = Complex{};
}
