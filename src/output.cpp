#include "output.hpp"

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
#include <system_error>

namespace {
constexpr std::array<char, 8> checkpointMagic{'S', 'H', '1', 'D',
                                              'R', 'S', 'T', '1'};
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
constexpr std::string_view modesHeader =
    "time,frame,mode,wavenumber,real,imaginary";

bool finiteField(const SpectralField &field) {
  return std::all_of(field.begin(), field.end(), [](Complex value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
  });
}

bool finiteValues(const std::vector<double> &values) {
  return std::all_of(values.begin(), values.end(),
                     [](double value) { return std::isfinite(value); });
}

std::uint64_t parseFrame(std::string_view text,
                         const std::filesystem::path &path) {
  std::uint64_t frame = 0;
  const auto [end, error] =
      std::from_chars(text.data(), text.data() + text.size(), frame);
  if (error != std::errc{} || end != text.data() + text.size())
    throw std::runtime_error("invalid frame in CSV during recovery: " +
                             path.string());
  return frame;
}

void closeChecked(std::ofstream &out, const std::string &message) {
  out.close();
  if (!out)
    throw std::runtime_error(message);
}

std::ofstream numericOutput(const std::filesystem::path &path,
                            std::ios::openmode mode = std::ios::out) {
  std::ofstream out(path, mode);
  if (!out)
    throw std::runtime_error("cannot write output file: " + path.string());
  out << std::setprecision(17);
  return out;
}

std::string frameName(const std::string &prefix, std::uint64_t frame,
                      const std::string &suffix) {
  std::ostringstream name;
  name << prefix << std::setw(8) << std::setfill('0') << frame << suffix;
  return name.str();
}

template <class Writer>
void writeAtomic(const std::filesystem::path &path, Writer writer) {
  const std::filesystem::path temporary = path.string() + ".tmp";
  std::error_code ignored;
  std::filesystem::remove(temporary, ignored);
  try {
    writer(temporary);
    std::filesystem::rename(temporary, path);
  } catch (...) {
    std::filesystem::remove(temporary, ignored);
    throw;
  }
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

SpectralField readInitialCondition(const Parameters &p,
                                   ComplexTransform &transform) {
  SpectralField spectral(p.gridPoints), physical(p.gridPoints);
  if (p.initialConditionFile.empty())
    return spectral;
  std::ifstream input(p.initialConditionFile);
  if (!input)
    throw std::runtime_error("cannot open initial condition: " +
                             p.initialConditionFile.string());
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
      throw std::runtime_error(
          "invalid numeric value on initial-condition line " +
          std::to_string(lineNumber));
    if (values.empty())
      continue;
    if ((values.size() != 2 && values.size() != 3) || index >= p.gridPoints)
      throw std::runtime_error(
          "initial-condition rows must contain 'real imag' or 'x real imag'");
    if (!std::all_of(values.begin(), values.end(),
                     [](double entry) { return std::isfinite(entry); }))
      throw std::runtime_error(
          "initial condition contains a non-finite value on line " +
          std::to_string(lineNumber));
    const std::size_t offset = values.size() == 3 ? 1 : 0;
    physical[index++] = Complex(values[offset], values[offset + 1]);
  }
  if (index != p.gridPoints)
    throw std::runtime_error("initial condition has " + std::to_string(index) +
                             " rows; expected " + std::to_string(p.gridPoints));
  transform.forward(physical, spectral);
  enforceStateConstraints(p, spectral);
  return spectral;
}

RestartState readCheckpoint(const Parameters &p,
                            const std::filesystem::path &path) {
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
  if (version != 1 || count != p.gridPoints ||
      std::abs(domainLength - p.domainLength) >
          1.e-12 * std::max(1.0, std::abs(p.domainLength)))
    throw std::runtime_error(
        "checkpoint grid or domain does not match the parameter file");
  state.randomEngineState = readString(in);
  state.normalDistributionState = readString(in);
  state.wavefunction.resize(p.gridPoints);
  for (Complex &coefficient : state.wavefunction) {
    double real = 0.0, imaginary = 0.0;
    readBinary(in, real);
    readBinary(in, imaginary);
    coefficient = Complex(real, imaginary);
  }
  if (!std::isfinite(state.time) || state.time < 0.0 ||
      !finiteField(state.wavefunction) ||
      in.peek() != std::ifstream::traits_type::eof())
    throw std::runtime_error("invalid checkpoint payload: " + path.string());
  state.restarting = true;
  return state;
}

void writeCheckpoint(const Parameters &p, const RestartState &state,
                     const std::string &randomEngineState,
                     const std::string &normalDistributionState,
                     const std::filesystem::path &path) {
  writeAtomic(path, [&](const std::filesystem::path &temporary) {
    std::ofstream out(temporary, std::ios::binary);
    if (!out)
      throw std::runtime_error("cannot write checkpoint: " + path.string());
    out.write(checkpointMagic.data(),
              static_cast<std::streamsize>(checkpointMagic.size()));
    const std::uint64_t version = 1;
    writeBinary(out, version);
    writeBinary(out, static_cast<std::uint64_t>(p.gridPoints));
    writeBinary(out, p.domainLength);
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
      throw std::runtime_error("failed while writing checkpoint: " +
                               path.string());
  });
}

bool managedDataFile(const std::filesystem::path &path) {
  const auto name = path.filename().string();
  return name == "restart_state.txt" || name.ends_with(".tmp") ||
         name.starts_with("wavefunction_") || name.starts_with("checkpoint_");
}

bool managedOutputFile(const std::filesystem::path &path) {
  const auto name = path.filename().string();
  return name == "diagnostics.csv" || name == "spectra.csv" ||
         name == "fluxes.csv" || name == "modes.csv" ||
         name == "forcing_summary.csv" || name == "forcing_spectrum.csv" ||
         name == "resolved_parameters.txt" || name.ends_with(".tmp");
}

void trimCsv(const std::filesystem::path &path, std::uint64_t committedFrame,
             std::string_view expectedHeader) {
  if (!std::filesystem::exists(path))
    throw std::runtime_error("restart output is missing: " + path.string());
  std::ifstream input(path);
  if (!input)
    throw std::runtime_error("cannot read output for recovery: " +
                             path.string());
  std::string line;
  if (!std::getline(input, line) || line != expectedHeader)
    throw std::runtime_error("CSV header does not match this solver version: " +
                             path.string());
  std::vector<std::string> retained{line};
  while (std::getline(input, line)) {
    if (line.empty() || line.front() == '#') {
      retained.push_back(line);
      continue;
    }
    std::istringstream fields(line);
    std::string timeText, frameText;
    if (!std::getline(fields, timeText, ',') ||
        !std::getline(fields, frameText, ','))
      throw std::runtime_error("malformed CSV during recovery: " +
                               path.string());
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

std::filesystem::path nextSegmentDirectory(const Parameters &p) {
  const auto root = p.outputDirectory / "segments";
  std::filesystem::create_directories(root);
  for (std::uint64_t index = 1;; ++index) {
    std::ostringstream name;
    name << "segment_" << std::setw(8) << std::setfill('0') << index;
    const auto candidate = root / name.str();
    if (!std::filesystem::exists(candidate)) {
      std::filesystem::create_directory(candidate);
      return candidate;
    }
  }
}

} // namespace

RestartState readRestartOrInitial(const Parameters &p,
                                  ComplexTransform &baseTransform) {
  const auto statePath = p.dataDirectory / "restart_state.txt";
  if (!std::filesystem::exists(statePath)) {
    RestartState state;
    state.randomSeed = p.randomSeed;
    state.wavefunction = readInitialCondition(p, baseTransform);
    return state;
  }
  std::ifstream stateFile(statePath);
  if (!stateFile)
    throw std::runtime_error("cannot open restart state: " +
                             statePath.string());
  std::string key, checkpointName;
  std::uint64_t version = 0;
  std::uint64_t frame = 0;
  double time = 0.0;
  if (!(stateFile >> key >> version) || key != "version" || version != 1 ||
      !(stateFile >> key >> time) || key != "time" ||
      !(stateFile >> key >> frame) || key != "frame" ||
      !(stateFile >> key >> checkpointName) || key != "checkpoint" ||
      !std::isfinite(time) || time < 0.0)
    throw std::runtime_error("invalid restart state: " + statePath.string());
  std::string extra;
  if (stateFile >> extra)
    throw std::runtime_error("unexpected data in restart state: " +
                             statePath.string());
  RestartState state = readCheckpoint(p, p.dataDirectory / checkpointName);
  if (state.frame != frame || std::abs(state.time - time) > 1.e-12)
    throw std::runtime_error("restart state and checkpoint disagree");
  return state;
}

void prepareOutput(const Parameters &p, bool restarting,
                   std::uint64_t committedFrame) {
  std::filesystem::create_directories(p.dataDirectory);
  std::filesystem::create_directories(p.outputDirectory);
  if (restarting) {
    trimCsv(p.outputDirectory / "diagnostics.csv", committedFrame,
            diagnosticsHeader);
    trimCsv(p.outputDirectory / "spectra.csv", committedFrame, spectraHeader);
    trimCsv(p.outputDirectory / "fluxes.csv", committedFrame, fluxesHeader);
    const auto modesPath = p.outputDirectory / "modes.csv";
    if (std::filesystem::exists(modesPath))
      trimCsv(modesPath, committedFrame, modesHeader);
    else if (p.writeModeDiagnostics)
      createCsv(modesPath, modesHeader);
    return;
  }

  bool hasManagedOutput = false;
  for (const auto &entry : std::filesystem::directory_iterator(p.dataDirectory))
    hasManagedOutput = hasManagedOutput || managedDataFile(entry.path());
  for (const auto &entry :
       std::filesystem::directory_iterator(p.outputDirectory))
    hasManagedOutput = hasManagedOutput || managedOutputFile(entry.path()) ||
                       entry.path().filename() == "segments";
  if (hasManagedOutput && !p.overwriteOutput)
    throw std::runtime_error("output from an existing run is present; choose "
                             "new directories or set overwriteOutput true");
  if (hasManagedOutput) {
    for (const auto &entry :
         std::filesystem::directory_iterator(p.dataDirectory))
      if (managedDataFile(entry.path()))
        std::filesystem::remove(entry.path());
    for (const auto &entry :
         std::filesystem::directory_iterator(p.outputDirectory)) {
      if (managedOutputFile(entry.path()))
        std::filesystem::remove(entry.path());
      else if (entry.path().filename() == "segments")
        std::filesystem::remove_all(entry.path());
    }
  }
  createCsv(p.outputDirectory / "diagnostics.csv", diagnosticsHeader);
  createCsv(p.outputDirectory / "spectra.csv", spectraHeader);
  createCsv(p.outputDirectory / "fluxes.csv", fluxesHeader);
  if (p.writeModeDiagnostics)
    createCsv(p.outputDirectory / "modes.csv", modesHeader);
}

void writeRunRecords(const Parameters &p, double startTime,
                     std::uint64_t startFrame,
                     const std::vector<double> &forcingAmplitude,
                     std::size_t forcedModeCount,
                     double waveActionInjectionCoefficient,
                     double quadraticEnergyInjectionCoefficient) {
  writeParameterRecord(p, p.outputDirectory, "CPU/FFTW");
  const auto segment = nextSegmentDirectory(p);
  writeParameterRecord(p, segment, "CPU/FFTW");
  {
    auto invocation = numericOutput(segment / "invocation.txt");
    invocation << std::setprecision(17) << "startTime " << startTime << '\n'
               << "startFrame " << startFrame << '\n';
    closeChecked(invocation, "failed while writing invocation record");
  }
  {
    auto summary = numericOutput(p.outputDirectory / "forcing_summary.csv");
    summary << "enabled,profile,forced_modes,wave_action_injection_coefficient,"
               "quadratic_energy_injection_coefficient\n"
            << std::boolalpha << p.forcingEnabled << ','
            << forcingProfileName(p.forcingProfile) << ',' << forcedModeCount
            << ',' << std::setprecision(17) << waveActionInjectionCoefficient
            << ',' << quadraticEnergyInjectionCoefficient << '\n';
    closeChecked(summary, "failed while writing forcing_summary.csv");
  }
  {
    auto spectrum = numericOutput(p.outputDirectory / "forcing_spectrum.csv");
    spectrum << "mode,wavenumber,amplitude\n" << std::setprecision(17);
    for (std::size_t i = 0; i < forcingAmplitude.size(); ++i)
      spectrum << signedWave(i, p.gridPoints) << ',' << waveNumber(p, i) << ','
               << forcingAmplitude[i] << '\n';
    closeChecked(spectrum, "failed while writing forcing_spectrum.csv");
  }
}

void writeWavefunctionAndRestart(const Parameters &p,
                                 ComplexTransform &baseTransform,
                                 const RestartState &state,
                                 const std::string &randomEngineState,
                                 const std::string &normalDistributionState) {
  SpectralField physical;
  baseTransform.inverse(state.wavefunction, physical);
  const auto wavePath =
      p.dataDirectory / frameName("wavefunction_", state.frame, ".dat");
  writeAtomic(wavePath, [&](const std::filesystem::path &temporary) {
    auto out = numericOutput(temporary);
    out << "# x real imaginary\n" << std::setprecision(17);
    const double dx = p.domainLength / static_cast<double>(p.gridPoints);
    for (std::size_t i = 0; i < physical.size(); ++i)
      out << dx * static_cast<double>(i) << ' ' << physical[i].real() << ' '
          << physical[i].imag() << '\n';
    closeChecked(out,
                 "failed while writing wavefunction: " + wavePath.string());
  });

  const auto checkpointName = frameName("checkpoint_", state.frame, ".bin");
  writeCheckpoint(p, state, randomEngineState, normalDistributionState,
                  p.dataDirectory / checkpointName);
  const auto statePath = p.dataDirectory / "restart_state.txt";
  writeAtomic(statePath, [&](const std::filesystem::path &temporary) {
    auto out = numericOutput(temporary);
    out << std::setprecision(17) << "version 1\n"
        << "time " << state.time << '\n'
        << "frame " << state.frame << '\n'
        << "checkpoint " << checkpointName << '\n';
    closeChecked(out,
                 "failed while writing restart state: " + statePath.string());
  });
}

double appendDiagnostics(const Parameters &p, NonlinearOperator &nonlinear,
                         double time, std::uint64_t frame,
                         const SpectralField &wavefunction,
                         DiagnosticsAverages &averages) {
  const std::size_t bins = p.gridPoints / 2;
  std::vector<double> spectrum(bins), quadraticSpectrum(bins);
  std::vector<double> waveTransfer(bins), energyTransfer(bins);
  SpectralField density, nonlinearRhs, wavefunctionDot, densityDot;
  nonlinear.densitySpectrum(wavefunction, density);
  nonlinear.evaluate(wavefunction, nonlinearRhs);
  wavefunctionDot.resize(p.gridPoints);

  double waveAction = 0.0, linearEnergy = 0.0, chemicalEnergy = 0.0;
  double waveHypo = 0.0, waveHyper = 0.0;
  double energyHypo = 0.0, energyHyper = 0.0;
  for (std::size_t i = 0; i < p.gridPoints; ++i) {
    const long mode = signedWave(i, p.gridPoints);
    if (std::abs(mode) >= static_cast<long>(bins))
      continue;
    const std::size_t bin = static_cast<std::size_t>(std::abs(mode));
    const double k = waveNumber(p, i), k2 = k * k;
    const double power = std::norm(wavefunction[i]);
    const double hamiltonian =
        -p.dispersionCoefficient * k2 + p.chemicalPotential;
    waveAction += p.domainLength * power;
    linearEnergy += p.domainLength * (-p.dispersionCoefficient * k2) * power;
    chemicalEnergy += p.domainLength * p.chemicalPotential * power;
    spectrum[bin] += power;
    quadraticSpectrum[bin] += hamiltonian * power;
    wavefunctionDot[i] =
        Complex(0.0, -hamiltonian) * wavefunction[i] + nonlinearRhs[i];
    if (p.hypoviscosity > 0.0 && (k2 > 0.0 || p.hypoviscosityOrder >= 0.0)) {
      const double rate = p.hypoviscosity * std::pow(k2, p.hypoviscosityOrder);
      waveHypo += 2.0 * p.domainLength * rate * power;
      energyHypo += 2.0 * p.domainLength * rate * hamiltonian * power;
    }
    if (p.hyperviscosity > 0.0) {
      const double rate =
          p.hyperviscosity * std::pow(k2, p.hyperviscosityOrder);
      waveHyper += 2.0 * p.domainLength * rate * power;
      energyHyper += 2.0 * p.domainLength * rate * hamiltonian * power;
    }
  }
  nonlinear.densityTimeDerivativeSpectrum(wavefunction, wavefunctionDot,
                                          densityDot);

  double nonlinearEnergy = 0.0;
  for (std::size_t i = 0; i < density.size(); ++i) {
    const double k = waveNumber(p, i);
    nonlinearEnergy += 0.5 * p.nonlinearityCoefficient * p.domainLength *
                       densityResponseMultiplier(p, k * k) *
                       std::norm(density[i]);
  }

  for (std::size_t i = 0; i < p.gridPoints; ++i) {
    const long mode = signedWave(i, p.gridPoints);
    if (std::abs(mode) >= static_cast<long>(bins))
      continue;
    const std::size_t bin = static_cast<std::size_t>(std::abs(mode));
    const double k = waveNumber(p, i);
    const double wave =
        -2.0 * p.domainLength *
        std::real(std::conj(wavefunction[i]) * wavefunctionDot[i]);
    const Complex gradient = Complex(0.0, k) * wavefunction[i];
    const Complex gradientDot = Complex(0.0, k) * wavefunctionDot[i];
    const double energy = 2.0 * p.dispersionCoefficient * p.domainLength *
                              std::real(gradient * std::conj(gradientDot)) +
                          p.chemicalPotential * wave -
                          p.nonlinearityCoefficient * p.domainLength *
                              densityResponseMultiplier(p, k * k) *
                              std::real(density[i] * std::conj(densityDot[i]));
    waveTransfer[bin] += wave;
    energyTransfer[bin] += energy;
  }
  for (std::size_t i = 1; i < bins; ++i) {
    waveTransfer[i] += waveTransfer[i - 1];
    energyTransfer[i] += energyTransfer[i - 1];
  }

  if (averages.spectrum.empty()) {
    averages.spectrum.assign(bins, 0.0);
    averages.waveFlux.assign(bins, 0.0);
    averages.energyFlux.assign(bins, 0.0);
  }
  ++averages.count;
  for (std::size_t i = 0; i < bins; ++i) {
    averages.spectrum[i] += spectrum[i];
    averages.waveFlux[i] += waveTransfer[i];
    averages.energyFlux[i] += energyTransfer[i];
  }

  const double totalEnergy = linearEnergy + chemicalEnergy + nonlinearEnergy;
  if (!std::isfinite(totalEnergy) || !std::isfinite(waveAction) ||
      !std::isfinite(linearEnergy) || !std::isfinite(chemicalEnergy) ||
      !std::isfinite(nonlinearEnergy) || !std::isfinite(waveHypo) ||
      !std::isfinite(waveHyper) || !std::isfinite(energyHypo) ||
      !std::isfinite(energyHyper) || !finiteValues(spectrum) ||
      !finiteValues(quadraticSpectrum) || !finiteValues(waveTransfer) ||
      !finiteValues(energyTransfer) || !finiteValues(averages.spectrum) ||
      !finiteValues(averages.waveFlux) || !finiteValues(averages.energyFlux))
    throw std::runtime_error("diagnostics contain a non-finite value");
  {
    auto out =
        numericOutput(p.outputDirectory / "diagnostics.csv", std::ios::app);
    out << std::setprecision(17) << time << ',' << frame << ',' << waveAction
        << ',' << linearEnergy << ',' << chemicalEnergy << ','
        << nonlinearEnergy << ',' << totalEnergy << ',' << waveHypo << ','
        << waveHyper << ',' << energyHypo << ',' << energyHyper << '\n';
    closeChecked(out, "failed while writing diagnostics.csv");
  }
  {
    auto out = numericOutput(p.outputDirectory / "spectra.csv", std::ios::app);
    const double dk = 2.0 * sh1dPi / p.domainLength;
    for (std::size_t i = 0; i < bins; ++i)
      out << time << ',' << frame << ',' << dk * static_cast<double>(i) << ','
          << spectrum[i] << ',' << quadraticSpectrum[i] << ','
          << averages.spectrum[i] / static_cast<double>(averages.count) << '\n';
    closeChecked(out, "failed while writing spectra.csv");
  }
  {
    auto out = numericOutput(p.outputDirectory / "fluxes.csv", std::ios::app);
    const double dk = 2.0 * sh1dPi / p.domainLength;
    for (std::size_t i = 0; i < bins; ++i)
      out << time << ',' << frame << ',' << dk * static_cast<double>(i) << ','
          << waveTransfer[i] << ',' << energyTransfer[i] << ','
          << averages.waveFlux[i] / static_cast<double>(averages.count) << ','
          << averages.energyFlux[i] / static_cast<double>(averages.count)
          << '\n';
    closeChecked(out, "failed while writing fluxes.csv");
  }
  if (p.writeModeDiagnostics) {
    auto out = numericOutput(p.outputDirectory / "modes.csv", std::ios::app);
    for (std::size_t i = 0; i < p.gridPoints; ++i)
      out << time << ',' << frame << ',' << signedWave(i, p.gridPoints) << ','
          << waveNumber(p, i) << ',' << wavefunction[i].real() << ','
          << wavefunction[i].imag() << '\n';
    closeChecked(out, "failed while writing modes.csv");
  }
  return totalEnergy;
}
