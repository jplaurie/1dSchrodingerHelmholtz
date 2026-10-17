#include "hdf5_io.hpp"

#include <chrono>
#include <filesystem>
#include <hdf5.h>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}

void replaceAttribute(const std::filesystem::path &path, const char *name, hid_t type, hid_t space,
                      const void *values) {
    const hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    require(file >= 0, "cannot open HDF5 fixture");
    require(H5Adelete(file, name) >= 0, "cannot delete fixture attribute");
    const hid_t attr = H5Acreate2(file, name, type, space, H5P_DEFAULT, H5P_DEFAULT);
    require(attr >= 0, "cannot create fixture attribute");
    require(H5Awrite(attr, type, values) >= 0, "cannot write fixture attribute");
    require(H5Aclose(attr) >= 0 && H5Fclose(file) >= 0, "cannot close HDF5 fixture");
}

void requireRejected(const std::filesystem::path &path) {
    try {
        (void)readHdf5Field(path);
    } catch (const std::runtime_error &) {
        return;
    }
    throw std::runtime_error("malformed HDF5 attribute was accepted");
}
} // namespace

int main() {
    const auto path =
        std::filesystem::temp_directory_path() /
        ("sh1d_metadata_" +
         std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".h5");
    try {
        Parameters p;
        p.gridPoints = 8;
        const auto reset = [&] { writeHdf5Field(path, p, 0.0, 0, SpectralField(8)); };
        const hsize_t dimensions[] = {32};
        const hid_t array = H5Screate_simple(1, dimensions, nullptr);
        const hid_t scalar = H5Screate(H5S_SCALAR);
        const std::vector<unsigned long long> numbers(32, 8);
        for (const char *name : {"grid_points", "domain_length", "time", "frame", "format"}) {
            reset();
            replaceAttribute(path, name, H5T_NATIVE_ULLONG, array, numbers.data());
            requireRejected(path);
        }
        const hid_t textType = H5Tcopy(H5T_C_S1);
        require(H5Tset_size(textType, 32) >= 0, "cannot size fixture strings");
        std::vector<char> strings(32 * 32, '\0');
        reset();
        replaceAttribute(path, "format", textType, array, strings.data());
        requireRejected(path);
        reset();
        replaceAttribute(path, "format", H5T_NATIVE_ULLONG, scalar, numbers.data());
        requireRejected(path);

        // h5py and other producers may write a scalar variable-length string.
        require(H5Tset_size(textType, H5T_VARIABLE) >= 0, "cannot create variable string");
        const char *format = "sh1d_wavefunction_v1";
        reset();
        replaceAttribute(path, "format", textType, scalar, &format);
        require(readHdf5Field(path).gridPoints == 8, "valid variable string was rejected");
        const std::vector<const char *> variableStrings(32, format);
        replaceAttribute(path, "format", textType, array, variableStrings.data());
        requireRejected(path);
        H5Tclose(textType);
        H5Sclose(scalar);
        H5Sclose(array);
        std::filesystem::remove(path);
        std::cout << "HDF5 metadata validation passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::filesystem::remove(path);
        std::cerr << error.what() << '\n';
        return 1;
    }
}
