#include "output.hpp"

#include "io_utils.hpp"
#include "spectral.hpp"

#include <array>
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace {
constexpr std::array csvNames{"diagnostics.csv", "spectra.csv", "fluxes.csv", "modes.csv"};

std::array<std::filesystem::path, 3> frameFiles(const Parameters &parameters, std::uint64_t frame) {
    std::ostringstream suffix;
    suffix << std::setw(8) << std::setfill('0') << frame;
    return {parameters.dataDirectory / ("wavefunction_" + suffix.str() + ".dat"),
            parameters.dataDirectory / ("checkpoint_" + suffix.str() + ".bin"),
            parameters.dataDirectory / ("wavefunction_" + suffix.str() + ".h5")};
}

std::filesystem::path appended(const std::filesystem::path &path, const char *suffix) {
    return path.string() + suffix;
}

std::optional<std::uint64_t> committedFrame(const Parameters &parameters) {
    const auto path = parameters.dataDirectory / "restart_state.txt";
    if (!std::filesystem::exists(path))
        return std::nullopt;
    std::ifstream input(path);
    std::string key, checkpoint;
    std::uint64_t version = 0, frame = 0;
    double time = 0.0;
    if (!(input >> key >> version) || key != "version" || version != 1 || !(input >> key >> time) ||
        key != "time" || !(input >> key >> frame) || key != "frame" ||
        !(input >> key >> checkpoint) || key != "checkpoint")
        throw std::runtime_error("cannot recover output with malformed restart metadata");
    return frame;
}

struct Journal {
    std::uint64_t frame = 0;
    bool previousMetadata = false;
    std::uint64_t previousFrame = 0;
    std::array<bool, csvNames.size()> csvExisted{};
    std::array<std::uintmax_t, csvNames.size()> csvSizes{};
    std::array<bool, 3> frameExisted{};
};

Journal readJournal(const Parameters &parameters) {
    std::ifstream input(parameters.dataDirectory / "output_transaction.txt");
    std::string format, directory;
    Journal journal;
    if (!(input >> format >> std::quoted(directory) >> journal.frame >> journal.previousMetadata >>
          journal.previousFrame) ||
        format != "sh1d_output_transaction_v1")
        throw std::runtime_error("malformed output transaction journal");
    if (std::filesystem::canonical(parameters.outputDirectory) != std::filesystem::path(directory))
        throw std::runtime_error(
            "recover the interrupted run using its original outputDirectory: " + directory);
    for (std::size_t i = 0; i < csvNames.size(); ++i)
        if (!(input >> journal.csvExisted[i] >> journal.csvSizes[i]))
            throw std::runtime_error("malformed CSV offsets in output transaction journal");
    for (bool &existed : journal.frameExisted)
        if (!(input >> existed))
            throw std::runtime_error("malformed frame records in output transaction journal");
    if (!(input >> std::ws).eof())
        throw std::runtime_error("unexpected data in output transaction journal");
    return journal;
}

void removeJournal(const Parameters &parameters) {
    std::filesystem::remove(parameters.dataDirectory / "output_transaction.txt");
    std::filesystem::remove(parameters.dataDirectory / "output_transaction.tmp");
}

std::filesystem::path nextSegmentDirectory(const Parameters &parameters) {
    const auto root = parameters.outputDirectory / "segments";
    std::filesystem::create_directories(root);
    for (std::uint64_t index = 1;; ++index) {
        std::ostringstream name;
        name << "segment_" << std::setw(8) << std::setfill('0') << index;
        const auto candidate = root / name.str();
        if (std::filesystem::create_directory(candidate))
            return candidate;
    }
}

void writeForcingRecords(const Parameters &parameters, const std::filesystem::path &directory,
                         const std::vector<double> &forcingAmplitude, std::size_t forcedModeCount,
                         double waveActionInjectionCoefficient,
                         double quadraticEnergyInjectionCoefficient) {
    auto summary = numericOutput(directory / "forcing_summary.csv");
    summary << "enabled,profile,forced_modes,wave_action_injection_coefficient,"
               "quadratic_energy_injection_coefficient\n"
            << std::boolalpha << parameters.forcingEnabled << ','
            << forcingProfileName(parameters.forcingProfile) << ',' << forcedModeCount << ','
            << waveActionInjectionCoefficient << ',' << quadraticEnergyInjectionCoefficient << '\n';
    closeChecked(summary, "failed while writing forcing_summary.csv");

    auto spectrum = numericOutput(directory / "forcing_spectrum.csv");
    spectrum << "mode,wavenumber,amplitude\n";
    for (std::size_t i = 0; i < forcingAmplitude.size(); ++i)
        spectrum << signedWave(i, parameters.gridPoints) << ',' << waveNumber(parameters, i) << ','
                 << forcingAmplitude[i] << '\n';
    closeChecked(spectrum, "failed while writing forcing_spectrum.csv");
}
} // namespace

bool recoverOutputTransaction(const Parameters &parameters) {
    if (!std::filesystem::exists(parameters.dataDirectory / "output_transaction.txt"))
        return false;
    const Journal journal = readJournal(parameters);
    const auto committed = committedFrame(parameters);
    if (committed && *committed == journal.frame) {
        finishOutputTransaction(parameters);
        return false;
    }
    if (committed.has_value() != journal.previousMetadata ||
        (committed && *committed != journal.previousFrame))
        throw std::runtime_error(
            "restart metadata does not match the interrupted output transaction");

    for (std::size_t i = 0; i < csvNames.size(); ++i) {
        const auto path = parameters.outputDirectory / csvNames[i];
        if (journal.csvExisted[i] && (!std::filesystem::exists(path) ||
                                      std::filesystem::file_size(path) < journal.csvSizes[i]))
            throw std::runtime_error("committed CSV data is missing: " + path.string());
        if (journal.csvExisted[i])
            std::filesystem::resize_file(path, journal.csvSizes[i]);
        else
            std::filesystem::remove(path);
    }
    const auto files = frameFiles(parameters, journal.frame);
    for (std::size_t i = 0; i < files.size(); ++i) {
        const auto backup = appended(files[i], ".previous");
        if (std::filesystem::exists(backup))
            std::filesystem::rename(backup, files[i]);
        else if (!journal.frameExisted[i])
            std::filesystem::remove(files[i]);
        std::filesystem::remove(appended(files[i], ".tmp"));
    }
    std::filesystem::remove(parameters.dataDirectory / "restart_state.tmp");
    removeJournal(parameters);
    return journal.frame == 0 && !committed;
}

void beginOutputTransaction(const Parameters &parameters, std::uint64_t frame) {
    if (std::filesystem::exists(parameters.dataDirectory / "output_transaction.txt"))
        throw std::runtime_error("an output transaction is already active");
    Journal journal;
    journal.frame = frame;
    const auto previous = committedFrame(parameters);
    if (previous && frame <= *previous)
        throw std::runtime_error("output transaction frame must follow the committed frame");
    journal.previousMetadata = previous.has_value();
    journal.previousFrame = previous.value_or(0);
    const auto files = frameFiles(parameters, frame);
    for (std::size_t i = 0; i < files.size(); ++i) {
        journal.frameExisted[i] = std::filesystem::exists(files[i]);
        if (journal.frameExisted[i] &&
            (!parameters.overwriteOutput || !std::filesystem::is_regular_file(files[i])))
            throw std::runtime_error("refusing to overwrite output frame file: " +
                                     files[i].string());
        for (const char *suffix : {".tmp", ".previous"})
            if (std::filesystem::exists(appended(files[i], suffix)))
                throw std::runtime_error("untracked temporary output file: " +
                                         appended(files[i], suffix).string());
    }
    for (std::size_t i = 0; i < csvNames.size(); ++i) {
        const auto path = parameters.outputDirectory / csvNames[i];
        journal.csvExisted[i] = std::filesystem::exists(path);
        if (journal.csvExisted[i])
            journal.csvSizes[i] = std::filesystem::file_size(path);
    }
    const auto temporary = parameters.dataDirectory / "output_transaction.tmp";
    auto output = numericOutput(temporary);
    output << "sh1d_output_transaction_v1\n"
           << std::quoted(std::filesystem::canonical(parameters.outputDirectory).string()) << '\n'
           << journal.frame << ' ' << journal.previousMetadata << ' ' << journal.previousFrame
           << '\n';
    for (std::size_t i = 0; i < csvNames.size(); ++i)
        output << journal.csvExisted[i] << ' ' << journal.csvSizes[i] << '\n';
    for (const bool existed : journal.frameExisted)
        output << existed << '\n';
    closeChecked(output, "cannot write output transaction journal");
    std::filesystem::rename(temporary, parameters.dataDirectory / "output_transaction.txt");
    for (std::size_t i = 0; i < files.size(); ++i)
        if (journal.frameExisted[i])
            std::filesystem::rename(files[i], appended(files[i], ".previous"));
}

void finishOutputTransaction(const Parameters &parameters) {
    const Journal journal = readJournal(parameters);
    if (committedFrame(parameters) != std::optional{journal.frame})
        throw std::runtime_error("cannot finish an uncommitted output transaction");
    for (const auto &path : frameFiles(parameters, journal.frame)) {
        std::filesystem::remove(appended(path, ".previous"));
        std::filesystem::remove(appended(path, ".tmp"));
    }
    removeJournal(parameters);
}

void writeRunRecords(const Parameters &parameters, double startTime, std::uint64_t startFrame,
                     const std::vector<double> &forcingAmplitude, std::size_t forcedModeCount,
                     double waveActionInjectionCoefficient,
                     double quadraticEnergyInjectionCoefficient) {
    const auto segment = nextSegmentDirectory(parameters);
    writeParameterRecord(parameters, segment, "CPU/FFTW");
    writeForcingRecords(parameters, segment, forcingAmplitude, forcedModeCount,
                        waveActionInjectionCoefficient, quadraticEnergyInjectionCoefficient);
    auto invocation = numericOutput(segment / "invocation.txt");
    invocation << "startTime " << startTime << '\n' << "startFrame " << startFrame << '\n';
    closeChecked(invocation, "failed while writing invocation record");

    for (const char *name :
         {"resolved_parameters.txt", "forcing_summary.csv", "forcing_spectrum.csv"}) {
        const auto target = parameters.outputDirectory / name;
        const auto temporary = std::filesystem::path(target.string() + ".tmp");
        std::filesystem::copy_file(segment / name, temporary,
                                   std::filesystem::copy_options::overwrite_existing);
        std::filesystem::rename(temporary, target);
    }
}
