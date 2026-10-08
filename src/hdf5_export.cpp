#include "hdf5_io.hpp"
#include "io_utils.hpp"

#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
void writeText(std::ostream &output, const Hdf5Field &field) {
    output << "# x real imaginary\n";
    for (std::size_t i = 0; i < field.gridPoints; ++i) {
        const double x = field.domainLength * static_cast<double>(i) / field.gridPoints;
        output << x << ' ' << field.wavefunction[i].real() << ' ' << field.wavefunction[i].imag()
               << '\n';
    }
}

void writeGnuplot(std::ostream &output, const Hdf5Field &field) {
    output << "# x real imaginary density phase\n";
    for (std::size_t i = 0; i < field.gridPoints; ++i) {
        const double x = field.domainLength * static_cast<double>(i) / field.gridPoints;
        const Complex value = field.wavefunction[i];
        output << x << ' ' << value.real() << ' ' << value.imag() << ' ' << std::norm(value) << ' '
               << std::arg(value) << '\n';
    }
}
} // namespace

int main(int argc, char **argv) {
    try {
        if (argc != 3 && argc != 5)
            throw std::runtime_error(
                "usage: sh1d_hdf5_export INPUT.h5 OUTPUT [--format text|gnuplot]");
        std::string format = "text";
        if (argc == 5) {
            if (std::string(argv[3]) != "--format")
                throw std::runtime_error("expected --format before output format");
            format = argv[4];
        }
        if (format != "text" && format != "gnuplot")
            throw std::runtime_error("export format must be text or gnuplot");

        const Hdf5Field field = readHdf5Field(argv[1]);
        std::ofstream file;
        std::ostream *output = &std::cout;
        if (std::string(argv[2]) != "-") {
            file.open(argv[2]);
            if (!file)
                throw std::runtime_error("cannot create export file: " + std::string(argv[2]));
            output = &file;
        }
        *output << std::scientific << std::setprecision(17);
        if (format == "text")
            writeText(*output, field);
        else
            writeGnuplot(*output, field);
        if (file.is_open())
            closeChecked(file, "failed while writing HDF5 export");
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
