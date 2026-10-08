#include "output.hpp"

#include "hdf5_io.hpp"
#include "io_utils.hpp"
#include "spectral.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace {
constexpr std::array<char, 8> checkpointMagic{'S', 'H', '1', 'D', 'R', 'S', 'T', '1'};
constexpr std::string_view diagnosticsHeader =
    "time,frame,wave_action,linear_energy,chemical_energy,nonlinear_energy,"
    "total_energy,wave_action_dissipation_hypo,wave_action_dissipation_hyper,"
    "quadratic_energy_dissipation_hypo,quadratic_energy_dissipation_hyper";
constexpr std::string_view spectraHeader =
    "time,frame,wavenumber,wave_action_spectrum,quadratic_energy_spectrum,"
    "wave_action_spectrum_average";
constexpr std::string_view fluxesHeader =
    "time,frame,wavenumber,wave_action_flux,energy_flux,wave_action_flux_"
    "average,energy_flux_average";
constexpr std::string_view modesHeader = "time,frame,mode,wavenumber,real,imaginary";

bool finiteField(const SpectralField &field) {
    return std::all_of(field.begin(), field.end(), [](Complex value) {
        return std::isfinite(value.real()) && std::isfinite(value.imag());
    });
}

std::uint64_t parseFrame(std::string_view text, const std::filesystem::path &path) {
    std::uint64_t frame = 0;
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), frame);
    if (error != std::errc{} || end != text.data() + text.size())
        throw std::runtime_error("invalid frame in CSV during recovery: " + path.string());
    return frame;
}

std::string frameName(const std::string &prefix, std::uint64_t frame, const std::string &suffix) {
    std::ostringstream name;
    name << prefix << std::setw(8) << std::setfill('0') << frame << suffix;
    return name.str();
}

template <class T> void writeBinary(std::ostream &out, const T &value) {
    out.write(reinterpret_cast<const char *>(&value), sizeof(value));
}

template <class T> void readBinary(std::istream &in, T &value) {
    in.read(reinterpret_cast<char *>(&value), sizeof(value));
    if (!in)
        throw std::runtime_error("truncated checkpoint");
}

void writeString(std::ostream &out, const std::string &value) {
    const auto size = static_cast<std::uint64_t>(value.size());
    writeBinary(out, size);
    out.write(value.data(), static_cast<std::streamsize>(value.size()));
}

std::string readString(std::istream &in) {
    std::uint64_t size = 0;
    readBinary(in, size);
    if (size > 10'000'000)
        throw std::runtime_error("invalid checkpoint string length");
    std::string value(static_cast<std::size_t>(size), '\0');
    in.read(value.data(), static_cast<std::streamsize>(value.size()));
    if (!in)
        throw std::runtime_error("truncated checkpoint string");
    return value;
}

SpectralField readInitialCondition(const Parameters &parameters, ComplexTransform &transform) {
    SpectralField spectral(parameters.gridPoints), physical(parameters.gridPoints);
    if (parameters.initialConditionFile.empty())
        return spectral;
    const auto extension = parameters.initialConditionFile.extension();
    if (extension == ".h5" || extension == ".hdf5") {
        Hdf5Field field = readHdf5Field(parameters.initialConditionFile);
        const double tolerance =
            1.e-12 * std::max({1.0, std::abs(field.domainLength), parameters.domainLength});
        if (field.gridPoints != parameters.gridPoints ||
            std::abs(field.domainLength - parameters.domainLength) > tolerance)
            throw std::runtime_error(
                "HDF5 initial-condition grid or domain does not match the parameter file");
        physical = std::move(field.wavefunction);
    } else {
        std::ifstream input(parameters.initialConditionFile);
        if (!input)
            throw std::runtime_error("cannot open initial condition: " +
                                     parameters.initialConditionFile.string());
        std::string line;
        std::size_t index = 0, lineNumber = 0;
        while (std::getline(input, line)) {
            ++lineNumber;
            if (const auto comment = line.find('#'); comment != std::string::npos)
                line.erase(comment);
            std::istringstream fields(line);
            std::vector<double> values;
            double value = 0.0;
            while (fields >> value)
                values.push_back(value);
            if (!fields.eof())
                throw std::runtime_error("invalid numeric value on initial-condition line " +
                                         std::to_string(lineNumber));
            if (values.empty())
                continue;
            if ((values.size() != 2 && values.size() != 3) || index >= parameters.gridPoints)
                throw std::runtime_error(
                    "initial-condition rows must contain 'real imag' or 'x real imag'");
            if (!std::all_of(values.begin(), values.end(),
                             [](double entry) { return std::isfinite(entry); }))
                throw std::runtime_error("initial condition contains a non-finite value on line " +
                                         std::to_string(lineNumber));
            const std::size_t offset = values.size() == 3 ? 1 : 0;
            physical[index++] = Complex(values[offset], values[offset + 1]);
        }
        if (index != parameters.gridPoints)
            throw std::runtime_error("initial condition has " + std::to_string(index) +
                                     " rows; expected " + std::to_string(parameters.gridPoints));
    }
    transform.forward(physical, spectral);
    enforceStateConstraints(parameters, spectral);
    return spectral;
}

RestartState readCheckpoint(const Parameters &parameters, const std::filesystem::path &path) {
    std::ifstream in(path, std::ios::binary);
    if (!in)
        throw std::runtime_error("cannot open checkpoint: " + path.string());
    std::array<char, 8> magic{};
    in.read(magic.data(), static_cast<std::streamsize>(magic.size()));
    if (magic != checkpointMagic)
        throw std::runtime_error("invalid checkpoint signature: " + path.string());
    std::uint64_t version = 0, count = 0;
    RestartState state;
    double domainLength = 0.0;
    readBinary(in, version);
    readBinary(in, count);
    readBinary(in, domainLength);
    readBinary(in, state.time);
    readBinary(in, state.frame);
    readBinary(in, state.randomSeed);
    if (version != 1 || count != parameters.gridPoints ||
        std::abs(domainLength - parameters.domainLength) >
            1.e-12 * std::max(1.0, std::abs(parameters.domainLength)))
        throw std::runtime_error("checkpoint grid or domain does not match the parameter file");
    state.randomEngineState = readString(in);
    state.normalDistributionState = readString(in);
    state.wavefunction.resize(parameters.gridPoints);
    for (Complex &coefficient : state.wavefunction) {
        double real = 0.0, imaginary = 0.0;
        readBinary(in, real);
        readBinary(in, imaginary);
        coefficient = Complex(real, imaginary);
    }
    if (!std::isfinite(state.time) || state.time < 0.0 || !finiteField(state.wavefunction) ||
        in.peek() != std::ifstream::traits_type::eof())
        throw std::runtime_error("invalid checkpoint payload: " + path.string());
    state.restarting = true;
    return state;
}

void writeCheckpoint(const Parameters &parameters, const RestartState &state,
                     const std::string &randomEngineState,
                     const std::string &normalDistributionState,
                     const std::filesystem::path &path) {
    writeAtomic(path, [&](const std::filesystem::path &temporary) {
        std::ofstream out(temporary, std::ios::binary);
        if (!out)
            throw std::runtime_error("cannot write checkpoint: " + path.string());
        out.write(checkpointMagic.data(), static_cast<std::streamsize>(checkpointMagic.size()));
        const std::uint64_t version = 1;
        writeBinary(out, version);
        writeBinary(out, static_cast<std::uint64_t>(parameters.gridPoints));
        writeBinary(out, parameters.domainLength);
        writeBinary(out, state.time);
        writeBinary(out, state.frame);
        writeBinary(out, state.randomSeed);
        writeString(out, randomEngineState);
        writeString(out, normalDistributionState);
        for (const Complex coefficient : state.wavefunction) {
            writeBinary(out, coefficient.real());
            writeBinary(out, coefficient.imag());
        }
        out.close();
        if (!out)
            throw std::runtime_error("failed while writing checkpoint: " + path.string());
    });
}

bool managedDataFile(const std::filesystem::path &path) {
    const auto name = path.filename().string();
    return name == "restart_state.txt" || name == "output_transaction.txt" ||
           name.ends_with(".tmp") || name.ends_with(".previous") ||
           name.starts_with("wavefunction_") || name.starts_with("checkpoint_");
}

bool managedOutputFile(const std::filesystem::path &path) {
    const auto name = path.filename().string();
    return name == "diagnostics.csv" || name == "spectra.csv" || name == "fluxes.csv" ||
           name == "modes.csv" || name == "forcing_summary.csv" || name == "forcing_spectrum.csv" ||
           name == "resolved_parameters.txt" || name.ends_with(".tmp");
}

void trimCsv(const std::filesystem::path &path, std::uint64_t committedFrame,
             std::string_view expectedHeader) {
    if (!std::filesystem::exists(path))
        throw std::runtime_error("restart output is missing: " + path.string());
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("cannot read output for recovery: " + path.string());
    std::string line;
    if (!std::getline(input, line) || line != expectedHeader)
        throw std::runtime_error("CSV header does not match this solver version: " + path.string());
    std::vector<std::string> retained{line};
    while (std::getline(input, line)) {
        if (line.empty() || line.front() == '#') {
            retained.push_back(line);
            continue;
        }
        std::istringstream fields(line);
        std::string timeText, frameText;
        if (!std::getline(fields, timeText, ',') || !std::getline(fields, frameText, ','))
            throw std::runtime_error("malformed CSV during recovery: " + path.string());
        const auto frame = parseFrame(frameText, path);
        if (frame <= committedFrame)
            retained.push_back(line);
    }
    input.close();
    writeAtomic(path, [&](const std::filesystem::path &temporary) {
        auto out = numericOutput(temporary);
        for (const auto &entry : retained)
            out << entry << '\n';
        closeChecked(out, "cannot recover output: " + path.string());
    });
}

void createCsv(const std::filesystem::path &path, std::string_view header) {
    auto out = numericOutput(path);
    out << header << '\n';
    closeChecked(out, "cannot create output: " + path.string());
}

} // namespace

RestartState readRestartOrInitial(const Parameters &parameters, ComplexTransform &baseTransform) {
    const auto statePath = parameters.dataDirectory / "restart_state.txt";
    if (!std::filesystem::exists(statePath)) {
        RestartState state;
        state.randomSeed = parameters.randomSeed;
        state.wavefunction = readInitialCondition(parameters, baseTransform);
        return state;
    }
    std::ifstream stateFile(statePath);
    if (!stateFile)
        throw std::runtime_error("cannot open restart state: " + statePath.string());
    std::string key, checkpointName;
    std::uint64_t version = 0;
    std::uint64_t frame = 0;
    double time = 0.0;
    if (!(stateFile >> key >> version) || key != "version" || version != 1 ||
        !(stateFile >> key >> time) || key != "time" || !(stateFile >> key >> frame) ||
        key != "frame" || !(stateFile >> key >> checkpointName) || key != "checkpoint" ||
        !std::isfinite(time) || time < 0.0)
        throw std::runtime_error("invalid restart state: " + statePath.string());
    std::string extra;
    if (stateFile >> extra)
        throw std::runtime_error("unexpected data in restart state: " + statePath.string());
    RestartState state = readCheckpoint(parameters, parameters.dataDirectory / checkpointName);
    if (state.frame != frame || std::abs(state.time - time) > 1.e-12)
        throw std::runtime_error("restart state and checkpoint disagree");
    return state;
}

void prepareOutput(const Parameters &parameters, bool restarting, std::uint64_t committedFrame) {
    std::filesystem::create_directories(parameters.dataDirectory);
    std::filesystem::create_directories(parameters.outputDirectory);
    if (restarting) {
        trimCsv(parameters.outputDirectory / "diagnostics.csv", committedFrame, diagnosticsHeader);
        trimCsv(parameters.outputDirectory / "spectra.csv", committedFrame, spectraHeader);
        trimCsv(parameters.outputDirectory / "fluxes.csv", committedFrame, fluxesHeader);
        const auto modesPath = parameters.outputDirectory / "modes.csv";
        if (std::filesystem::exists(modesPath))
            trimCsv(modesPath, committedFrame, modesHeader);
        else if (parameters.writeModeDiagnostics)
            createCsv(modesPath, modesHeader);
        return;
    }

    bool hasManagedOutput = false;
    for (const auto &entry : std::filesystem::directory_iterator(parameters.dataDirectory))
        hasManagedOutput = hasManagedOutput || managedDataFile(entry.path());
    for (const auto &entry : std::filesystem::directory_iterator(parameters.outputDirectory))
        hasManagedOutput = hasManagedOutput || managedOutputFile(entry.path()) ||
                           entry.path().filename() == "segments";
    if (hasManagedOutput && !parameters.overwriteOutput)
        throw std::runtime_error("output from an existing run is present; choose "
                                 "new directories or set overwriteOutput true");
    if (hasManagedOutput) {
        for (const auto &entry : std::filesystem::directory_iterator(parameters.dataDirectory))
            if (managedDataFile(entry.path()))
                std::filesystem::remove(entry.path());
        for (const auto &entry : std::filesystem::directory_iterator(parameters.outputDirectory)) {
            if (managedOutputFile(entry.path()))
                std::filesystem::remove(entry.path());
            else if (entry.path().filename() == "segments")
                std::filesystem::remove_all(entry.path());
        }
    }
    createCsv(parameters.outputDirectory / "diagnostics.csv", diagnosticsHeader);
    createCsv(parameters.outputDirectory / "spectra.csv", spectraHeader);
    createCsv(parameters.outputDirectory / "fluxes.csv", fluxesHeader);
    if (parameters.writeModeDiagnostics)
        createCsv(parameters.outputDirectory / "modes.csv", modesHeader);
}

void writeWavefunctionAndRestart(const Parameters &parameters, ComplexTransform &baseTransform,
                                 const RestartState &state, const std::string &randomEngineState,
                                 const std::string &normalDistributionState) {
    SpectralField physical;
    baseTransform.inverse(state.wavefunction, physical);
    if (parameters.fieldOutputFormat != FieldOutputFormat::hdf5) {
        const auto wavePath =
            parameters.dataDirectory / frameName("wavefunction_", state.frame, ".dat");
        writeAtomic(wavePath, [&](const std::filesystem::path &temporary) {
            auto out = numericOutput(temporary);
            out << "# x real imaginary\n" << std::setprecision(17);
            const double dx = parameters.domainLength / static_cast<double>(parameters.gridPoints);
            for (std::size_t i = 0; i < physical.size(); ++i)
                out << dx * static_cast<double>(i) << ' ' << physical[i].real() << ' '
                    << physical[i].imag() << '\n';
            closeChecked(out, "failed while writing wavefunction: " + wavePath.string());
        });
    }
    if (parameters.fieldOutputFormat != FieldOutputFormat::text) {
        const auto hdf5Path =
            parameters.dataDirectory / frameName("wavefunction_", state.frame, ".h5");
        writeAtomic(hdf5Path, [&](const std::filesystem::path &temporary) {
            writeHdf5Field(temporary, parameters, state.time, state.frame, physical);
        });
    }

    const auto checkpointName = frameName("checkpoint_", state.frame, ".bin");
    writeCheckpoint(parameters, state, randomEngineState, normalDistributionState,
                    parameters.dataDirectory / checkpointName);
    const auto statePath = parameters.dataDirectory / "restart_state.txt";
    writeAtomic(statePath, [&](const std::filesystem::path &temporary) {
        auto out = numericOutput(temporary);
        out << std::setprecision(17) << "version 1\n"
            << "time " << state.time << '\n'
            << "frame " << state.frame << '\n'
            << "checkpoint " << checkpointName << '\n';
        closeChecked(out, "failed while writing restart state: " + statePath.string());
    });
}
