#include "output.hpp"

#include "io_utils.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
class TemporaryDirectory {
  public:
    TemporaryDirectory() {
        const auto suffix = std::chrono::steady_clock::now().time_since_epoch().count();
        path_ = std::filesystem::temp_directory_path() /
                ("sh1d_transaction_tests_" + std::to_string(suffix));
        std::filesystem::create_directory(path_);
    }
    ~TemporaryDirectory() {
        std::error_code ignored;
        std::filesystem::remove_all(path_, ignored);
    }
    [[nodiscard]] const std::filesystem::path &path() const { return path_; }

  private:
    std::filesystem::path path_;
};

void write(const std::filesystem::path &path, const std::string &text,
           std::ios::openmode mode = std::ios::out) {
    std::ofstream output(path, mode);
    output << text;
    output.close();
    if (!output)
        throw std::runtime_error("could not write transaction fixture");
}

std::string read(const std::filesystem::path &path) {
    std::ifstream input(path);
    return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}

void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
} // namespace

int main() {
    try {
        TemporaryDirectory temporary;
        Parameters parameters;
        parameters.dataDirectory = temporary.path() / "data";
        parameters.outputDirectory = temporary.path() / "output";
        std::filesystem::create_directory(parameters.dataDirectory);
        std::filesystem::create_directory(parameters.outputDirectory);
        {
            auto output = numericOutput(temporary.path() / "numeric.txt");
            output << 1.0 << '\n';
            closeChecked(output, "could not close numeric output fixture");
        }
        require(read(temporary.path() / "numeric.txt") == "1.00000000000000000e+00\n",
                "numeric text output is not scientific with round-trip precision");
        for (const char *name : {"diagnostics.csv", "spectra.csv", "fluxes.csv"})
            write(parameters.outputDirectory / name, "header\n");

        beginOutputTransaction(parameters, 0);
        write(parameters.outputDirectory / "diagnostics.csv", "partial\n", std::ios::app);
        write(parameters.dataDirectory / "wavefunction_00000000.dat", "partial\n");
        require(recoverOutputTransaction(parameters), "fresh output recovery was not identified");
        require(read(parameters.outputDirectory / "diagnostics.csv") == "header\n",
                "CSV was not rolled back to its committed size");
        require(!std::filesystem::exists(parameters.dataDirectory / "wavefunction_00000000.dat"),
                "uncommitted frame file was not removed");

        beginOutputTransaction(parameters, 0);
        write(parameters.dataDirectory / "wavefunction_00000000.dat", "complete\n");
        write(parameters.dataDirectory / "restart_state.txt",
              "version 1\ntime 0\nframe 0\ncheckpoint checkpoint_00000000.bin\n");
        require(!recoverOutputTransaction(parameters), "committed output was treated as fresh");
        require(read(parameters.dataDirectory / "wavefunction_00000000.dat") == "complete\n",
                "committed frame file was not retained");
        require(!std::filesystem::exists(parameters.dataDirectory / "output_transaction.txt"),
                "committed transaction journal was not removed");

        bool rejectedCommittedFrame = false;
        try {
            beginOutputTransaction(parameters, 0);
        } catch (const std::runtime_error &) {
            rejectedCommittedFrame = true;
        }
        require(rejectedCommittedFrame, "transaction allowed a committed frame to be rewritten");
        std::cout << "output transaction tests passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
