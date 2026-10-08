#include "parameters.hpp"
#include "hdf5_io.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {
class TemporaryDirectory {
  public:
    TemporaryDirectory() {
        const auto suffix = std::chrono::steady_clock::now().time_since_epoch().count();
        path_ = std::filesystem::temp_directory_path() /
                ("sh1d_parameter_tests_" + std::to_string(suffix));
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

void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}

template <class Function> void requireThrows(Function function, const char *message) {
    try {
        function();
    } catch (const std::runtime_error &) {
        return;
    }
    throw std::runtime_error(message);
}

void writeText(const std::filesystem::path &path, const std::string &text) {
    std::ofstream output(path);
    output << text;
    output.close();
    require(static_cast<bool>(output), "could not create parameter test input");
}

void testParsing() {
    TemporaryDirectory temporary;
    const auto parameterFile = temporary.path() / "test.params";
    writeText(parameterFile, "# Both assignment styles are supported.\n"
                             "gridPoints = 32\n"
                             "domainLength 12.5\n"
                             "numberOfSteps 20\n"
                             "outputIntervalSteps 4\n"
                             "model longWave\n"
                             "integrator rk2\n"
                             "forcingEnabled false\n"
                             "forcingProfile logNormal\n"
                             "fieldOutputFormat text\n"
                             "hdf5CompressionLevel 4\n"
                             "fftwPlanning patient\n"
                             "fftwWisdomFile plans.wisdom\n"
                             "threadCount 2\n");

    const Parameters parameters = readParameters(parameterFile);
    require(parameters.gridPoints == 32 && parameters.domainLength == 12.5,
            "grid settings were not parsed");
    require(parameters.numberOfSteps == 20 && parameters.outputIntervalSteps == 4,
            "run length was not parsed");
    require(parameters.model == Model::longWave, "model was not converted to its enum");
    require(parameters.integrator == Integrator::integratingFactorRk2,
            "integrator was not converted to its enum");
    require(parameters.forcingProfile == ForcingProfile::logNormal,
            "forcing profile was not converted to its enum");
    require(parameters.fieldOutputFormat == FieldOutputFormat::text &&
                parameters.hdf5CompressionLevel == 4 &&
                parameters.fftwPlanning == FftwPlanning::patient &&
                parameters.fftwWisdomFile == "plans.wisdom",
            "output or FFTW settings were not parsed");
    require(!parameters.forcingEnabled && parameters.threadCount == 2,
            "boolean or integer settings were not parsed");

    const auto malformed = temporary.path() / "malformed.params";
    writeText(malformed, "gridPoints 32 extra\n");
    requireThrows([&] { (void)readParameters(malformed); },
                  "malformed parameter line was accepted");

    const auto planning = temporary.path() / "planning.params";
    writeText(planning, "fftwPlanning exhaustive\n");
    requireThrows([&] { (void)readParameters(planning); },
                  "invalid FFTW planning mode was accepted");

    const auto format = temporary.path() / "format.params";
    writeText(format, "fieldOutputFormat netcdf\n");
    requireThrows([&] { (void)readParameters(format); },
                  "invalid field output format was accepted");

    const auto compression = temporary.path() / "compression.params";
    writeText(compression, "hdf5CompressionLevel 10\n");
    requireThrows([&] { (void)readParameters(compression); },
                  "invalid HDF5 compression level was accepted");
}

void testValidation() {
    Parameters parameters;
    parameters.gridPoints = 16'388;
    requireThrows([&] { validateParameters(parameters); },
                  "grid limit above 16384 was not enforced");
    parameters.gridPoints = 16'384;
    parameters.threadCount = 12;
    validateParameters(parameters);
    parameters.threadCount = -1;
    requireThrows([&] { validateParameters(parameters); },
                  "negative thread count was not rejected");

    parameters = Parameters{};
    parameters.forcingEnabled = true;
    parameters.forcingProfile = ForcingProfile::logNormal;
    parameters.forcingLogWidth = 0.0;
    requireThrows([&] { validateParameters(parameters); }, "zero log-normal width was accepted");
    parameters.forcingLogWidth = 0.2;
    parameters.forcingProfile = ForcingProfile::singleMode;
    parameters.forcingWavenumber = 3.25;
    requireThrows([&] { validateParameters(parameters); },
                  "off-grid deterministic forcing was accepted");

    parameters.forcingEnabled = false;
    parameters.dataDirectory = "same-directory";
    parameters.outputDirectory = "./same-directory";
    requireThrows([&] { validateParameters(parameters); },
                  "aliased data and output directories were accepted");

    parameters = Parameters{};
    parameters.fieldOutputFormat = FieldOutputFormat::hdf5;
    if (hdf5Available())
        validateParameters(parameters);
    else
        requireThrows([&] { validateParameters(parameters); },
                      "HDF5 output was accepted by a build without HDF5 support");
}
} // namespace

int main() {
    try {
        testParsing();
        testValidation();
        std::cout << "parameter tests passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
