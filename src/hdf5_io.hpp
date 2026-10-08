#pragma once

#include "parameters.hpp"
#include "types.hpp"

#include <cstdint>
#include <filesystem>

struct Hdf5Field {
    std::size_t gridPoints = 0;
    double domainLength = 0.0;
    double time = 0.0;
    std::uint64_t frame = 0;
    SpectralField wavefunction;
};

[[nodiscard]] bool hdf5Available();
void writeHdf5Field(const std::filesystem::path &path, const Parameters &parameters, double time,
                    std::uint64_t frame, const SpectralField &physical);
Hdf5Field readHdf5Field(const std::filesystem::path &path);
