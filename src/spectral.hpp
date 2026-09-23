#pragma once

#include "parameters.hpp"
#include "types.hpp"

#include <cstddef>

constexpr double sh1dPi = 3.141592653589793238462643383279502884;

[[nodiscard]] long signedWave(std::size_t index, std::size_t count);
[[nodiscard]] double waveNumber(const Parameters &p, std::size_t index);
[[nodiscard]] double densityResponseMultiplier(const Parameters &p,
                                               double wavenumberSquared);
[[nodiscard]] std::size_t paddedIndexForBaseMode(const Parameters &p,
                                                 std::size_t index);
void embedBaseSpectrum(const Parameters &p, const SpectralField &base,
                       SpectralField &padded);
void extractBaseSpectrum(const Parameters &p, const SpectralField &padded,
                         SpectralField &base);
void filterPaddedSpectrum(const Parameters &p, SpectralField &field);
void enforceStateConstraints(const Parameters &p, SpectralField &field);
