#include "spectral.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

long signedWave(std::size_t index, std::size_t count) {
    return index < (count + 1) / 2 ? static_cast<long>(index)
                                   : static_cast<long>(index) - static_cast<long>(count);
}

double waveNumber(const Parameters &parameters, std::size_t index) {
    return 2.0 * sh1dPi * static_cast<double>(signedWave(index, parameters.gridPoints)) /
           parameters.domainLength;
}

double densityResponseMultiplier(const Parameters &parameters, double k2) {
    switch (parameters.model) {
    case Model::schrodingerHelmholtz:
        return 1.0 / (1.0 + parameters.helmholtzParameter * k2);
    case Model::longWave:
        return 1.0 - parameters.helmholtzParameter * k2;
    case Model::nonlinearSchrodinger:
        return 1.0;
    }
    throw std::logic_error("unknown model in density response");
}

std::size_t paddedIndexForBaseMode(const Parameters &parameters, std::size_t index) {
    const long mode = signedWave(index, parameters.gridPoints);
    return mode >= 0
               ? static_cast<std::size_t>(mode)
               : static_cast<std::size_t>(static_cast<long>(parameters.paddedGridPoints()) + mode);
}

void embedBaseSpectrum(const Parameters &parameters, const SpectralField &base,
                       SpectralField &padded) {
    if (base.size() != parameters.gridPoints || padded.size() != parameters.paddedGridPoints())
        throw std::runtime_error("invalid field size while embedding spectrum");
    std::fill(padded.begin(), padded.end(), Complex{});
    for (std::size_t i = 0; i < base.size(); ++i)
        padded[paddedIndexForBaseMode(parameters, i)] = base[i];
}

void extractBaseSpectrum(const Parameters &parameters, const SpectralField &padded,
                         SpectralField &base) {
    if (base.size() != parameters.gridPoints || padded.size() != parameters.paddedGridPoints())
        throw std::runtime_error("invalid field size while extracting spectrum");
    for (std::size_t i = 0; i < base.size(); ++i)
        base[i] = padded[paddedIndexForBaseMode(parameters, i)];
}

void filterPaddedSpectrum(const Parameters &parameters, SpectralField &field) {
    if (field.size() != parameters.paddedGridPoints())
        throw std::runtime_error("invalid padded field size");
    const long cutoff = static_cast<long>(parameters.gridPoints / 2);
    for (std::size_t i = 0; i < field.size(); ++i)
        if (std::abs(signedWave(i, field.size())) >= cutoff)
            field[i] = Complex{};
}

void enforceStateConstraints(const Parameters &parameters, SpectralField &field) {
    if (field.size() != parameters.gridPoints)
        throw std::runtime_error("invalid state size");
    field[parameters.gridPoints / 2] = Complex{};
    if (parameters.hypoviscosity > 0.0 && parameters.hypoviscosityOrder < 0.0)
        field[0] = Complex{};
}
