#pragma once

#include "parameters.hpp"
#include "types.hpp"

#include <cstddef>
#include <span>

constexpr double sh1dPi = 3.141592653589793238462643383279502884;

[[nodiscard]] long signedWave(std::size_t index, std::size_t count);
[[nodiscard]] double waveNumber(const Parameters &parameters, std::size_t index);
[[nodiscard]] bool retainedStateMode(const Parameters &parameters, std::size_t index);
[[nodiscard]] double densityResponseMultiplier(const Parameters &parameters,
                                               double wavenumberSquared);
[[nodiscard]] std::size_t paddedIndexForBaseMode(const Parameters &parameters, std::size_t index);
void embedBaseSpectrum(const Parameters &parameters, std::span<const Complex> base,
                       std::span<Complex> padded);
void extractBaseSpectrum(const Parameters &parameters, std::span<const Complex> padded,
                         std::span<Complex> base);
void filterPaddedSpectrum(const Parameters &parameters, std::span<Complex> field);
void enforceStateConstraints(const Parameters &parameters, SpectralField &field);
