#include "parameters.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <type_traits>

namespace {
bool parseBool(const std::string &text, const std::string &key) {
  if (text == "true" || text == "1")
    return true;
  if (text == "false" || text == "0")
    return false;
  throw std::runtime_error(key + " must be true or false, got: " + text);
}

template <class T>
T parseNumber(const std::string &text, const std::string &key) {
  if constexpr (std::is_unsigned_v<T>)
    if (!text.empty() && text.front() == '-')
      throw std::runtime_error(key + " cannot be negative: " + text);
  std::istringstream input(text);
  T value{};
  input >> value;
  if (!input || !(input >> std::ws).eof())
    throw std::runtime_error("invalid value for " + key + ": " + text);
  return value;
}

std::string trim(const std::string &value) {
  const auto first = value.find_first_not_of(" \t\r\n");
  if (first == std::string::npos)
    return {};
  return value.substr(first, value.find_last_not_of(" \t\r\n") - first + 1);
}

template <class T, std::size_t N>
bool parseMember(
    std::string_view key, const std::string &value, Parameters &parameters,
    const std::pair<std::string_view, T Parameters::*> (&members)[N]) {
  const auto entry = std::find_if(
      std::begin(members), std::end(members),
      [key](const auto &candidate) { return candidate.first == key; });
  if (entry == std::end(members))
    return false;
  if constexpr (std::is_same_v<T, bool>)
    parameters.*entry->second = parseBool(value, std::string(key));
  else if constexpr (std::is_same_v<T, std::filesystem::path>)
    parameters.*entry->second = value;
  else
    parameters.*entry->second = parseNumber<T>(value, std::string(key));
  return true;
}

Model parseModel(const std::string &text) {
  if (text == "schrodingerHelmholtz")
    return Model::schrodingerHelmholtz;
  if (text == "longWave")
    return Model::longWave;
  if (text == "nonlinearSchrodinger")
    return Model::nonlinearSchrodinger;
  throw std::runtime_error("model must be schrodingerHelmholtz, longWave, or "
                           "nonlinearSchrodinger");
}

Integrator parseIntegrator(const std::string &text) {
  if (text == "etd2")
    return Integrator::etd2;
  if (text == "etd4")
    return Integrator::etd4;
  if (text == "rk2")
    return Integrator::integratingFactorRk2;
  throw std::runtime_error("integrator must be etd2, etd4, or rk2");
}

ForcingProfile parseForcingProfile(const std::string &text) {
  if (text == "annulus")
    return ForcingProfile::annulus;
  if (text == "gaussian")
    return ForcingProfile::gaussian;
  if (text == "exponential")
    return ForcingProfile::exponential;
  if (text == "logNormal")
    return ForcingProfile::logNormal;
  if (text == "singleMode")
    return ForcingProfile::singleMode;
  throw std::runtime_error(
      "forcingProfile must be annulus, gaussian, exponential, logNormal, or "
      "singleMode");
}

bool isFinite(double value) { return std::isfinite(value); }
} // namespace

const char *modelName(Model model) {
  switch (model) {
  case Model::schrodingerHelmholtz:
    return "schrodingerHelmholtz";
  case Model::longWave:
    return "longWave";
  case Model::nonlinearSchrodinger:
    return "nonlinearSchrodinger";
  }
  throw std::logic_error("unknown model");
}

const char *integratorName(Integrator integrator) {
  switch (integrator) {
  case Integrator::etd2:
    return "etd2";
  case Integrator::etd4:
    return "etd4";
  case Integrator::integratingFactorRk2:
    return "rk2";
  }
  throw std::logic_error("unknown integrator");
}

const char *forcingProfileName(ForcingProfile profile) {
  switch (profile) {
  case ForcingProfile::annulus:
    return "annulus";
  case ForcingProfile::gaussian:
    return "gaussian";
  case ForcingProfile::exponential:
    return "exponential";
  case ForcingProfile::logNormal:
    return "logNormal";
  case ForcingProfile::singleMode:
    return "singleMode";
  }
  throw std::logic_error("unknown forcing profile");
}

Parameters readParameters(const std::filesystem::path &path) {
  std::ifstream input(path);
  if (!input)
    throw std::runtime_error("cannot open parameter file: " + path.string());
  Parameters p;
  std::string line;
  std::size_t lineNumber = 0;
  while (std::getline(input, line)) {
    ++lineNumber;
    if (const auto comment = line.find('#'); comment != std::string::npos)
      line.erase(comment);
    line = trim(line);
    if (line.empty())
      continue;
    std::replace(line.begin(), line.end(), '=', ' ');
    std::istringstream fields(line);
    std::string key, value, extra;
    fields >> key >> value;
    if (key.empty() || value.empty() || (fields >> extra))
      throw std::runtime_error("invalid parameter line " +
                               std::to_string(lineNumber));

    static constexpr std::pair<std::string_view, std::size_t Parameters::*>
        sizes[]{{"gridPoints", &Parameters::gridPoints}};
    static constexpr std::pair<std::string_view, std::uint64_t Parameters::*>
        counts[]{{"numberOfSteps", &Parameters::numberOfSteps},
                 {"outputIntervalSteps", &Parameters::outputIntervalSteps},
                 {"randomSeed", &Parameters::randomSeed}};
    static constexpr std::pair<std::string_view, double Parameters::*> reals[]{
        {"domainLength", &Parameters::domainLength},
        {"timeStep", &Parameters::timeStep},
        {"dispersionCoefficient", &Parameters::dispersionCoefficient},
        {"nonlinearityCoefficient", &Parameters::nonlinearityCoefficient},
        {"chemicalPotential", &Parameters::chemicalPotential},
        {"helmholtzParameter", &Parameters::helmholtzParameter},
        {"hyperviscosity", &Parameters::hyperviscosity},
        {"hyperviscosityOrder", &Parameters::hyperviscosityOrder},
        {"hypoviscosity", &Parameters::hypoviscosity},
        {"hypoviscosityOrder", &Parameters::hypoviscosityOrder},
        {"forcingWavenumber", &Parameters::forcingWavenumber},
        {"forcingWidth", &Parameters::forcingWidth},
        {"forcingAmplitude", &Parameters::forcingAmplitude},
        {"forcingShapeOrder", &Parameters::forcingShapeOrder},
        {"forcingLogWidth", &Parameters::forcingLogWidth},
        {"targetWaveActionInjectionRate",
         &Parameters::targetWaveActionInjectionRate}};
    static constexpr std::pair<std::string_view, bool Parameters::*> booleans[]{
        {"forcingEnabled", &Parameters::forcingEnabled},
        {"writeModeDiagnostics", &Parameters::writeModeDiagnostics},
        {"overwriteOutput", &Parameters::overwriteOutput}};
    static constexpr std::pair<std::string_view, int Parameters::*> integers[]{
        {"threadCount", &Parameters::threadCount}};
    static constexpr std::pair<std::string_view,
                               std::filesystem::path Parameters::*>
        paths[]{{"initialConditionFile", &Parameters::initialConditionFile},
                {"dataDirectory", &Parameters::dataDirectory},
                {"outputDirectory", &Parameters::outputDirectory}};

    const bool recognized = parseMember(key, value, p, sizes) ||
                            parseMember(key, value, p, counts) ||
                            parseMember(key, value, p, reals) ||
                            parseMember(key, value, p, booleans) ||
                            parseMember(key, value, p, integers) ||
                            parseMember(key, value, p, paths);
    if (key == "model")
      p.model = parseModel(value);
    else if (key == "integrator")
      p.integrator = parseIntegrator(value);
    else if (key == "forcingProfile")
      p.forcingProfile = parseForcingProfile(value);
    else if (!recognized)
      throw std::runtime_error("unknown parameter key on line " +
                               std::to_string(lineNumber) + ": " + key);
  }
  validateParameters(p);
  return p;
}

void validateParameters(const Parameters &p) {
  if (p.gridPoints < 8 || p.gridPoints % 4 != 0 || p.gridPoints > 16'384)
    throw std::runtime_error(
        "gridPoints must be a multiple of four between 8 and 16384");
  if (p.gridPoints >
      2 * (static_cast<std::size_t>(std::numeric_limits<int>::max()) / 3))
    throw std::runtime_error("3/2-rule grid exceeds FFTW integer limits");
  if (!(p.domainLength > 0.0) || !isFinite(p.domainLength))
    throw std::runtime_error("domainLength must be finite and positive");
  if (!(p.timeStep > 0.0) || !isFinite(p.timeStep))
    throw std::runtime_error("timeStep must be finite and positive");
  if (p.numberOfSteps == 0 || p.outputIntervalSteps == 0)
    throw std::runtime_error(
        "numberOfSteps and outputIntervalSteps must be positive");
  if (p.threadCount < 0 || p.threadCount > 2)
    throw std::runtime_error("threadCount must be 0, 1, or 2");
  for (const auto [value, name] :
       {std::pair{p.dispersionCoefficient, "dispersionCoefficient"},
        std::pair{p.nonlinearityCoefficient, "nonlinearityCoefficient"},
        std::pair{p.chemicalPotential, "chemicalPotential"},
        std::pair{p.helmholtzParameter, "helmholtzParameter"},
        std::pair{p.hyperviscosity, "hyperviscosity"},
        std::pair{p.hyperviscosityOrder, "hyperviscosityOrder"},
        std::pair{p.hypoviscosity, "hypoviscosity"},
        std::pair{p.hypoviscosityOrder, "hypoviscosityOrder"}})
    if (!isFinite(value))
      throw std::runtime_error(std::string(name) + " must be finite");
  if (p.helmholtzParameter < 0.0)
    throw std::runtime_error("helmholtzParameter cannot be negative");
  if (p.hyperviscosity < 0.0 || p.hypoviscosity < 0.0 ||
      p.hyperviscosityOrder < 0.0)
    throw std::runtime_error("dissipation coefficients/orders are invalid");
  if (p.forcingEnabled) {
    if (!(p.forcingWavenumber > 0.0) || p.forcingWidth < 0.0 ||
        p.forcingAmplitude < 0.0 || p.forcingShapeOrder <= 0.0 ||
        p.forcingLogWidth <= 0.0 || p.targetWaveActionInjectionRate < 0.0 ||
        !isFinite(p.forcingWavenumber) || !isFinite(p.forcingWidth) ||
        !isFinite(p.forcingAmplitude) || !isFinite(p.forcingShapeOrder) ||
        !isFinite(p.forcingLogWidth) ||
        !isFinite(p.targetWaveActionInjectionRate))
      throw std::runtime_error("forcing parameters are invalid or non-finite");
    if (p.forcingProfile == ForcingProfile::gaussian && p.forcingWidth <= 0.0)
      throw std::runtime_error(
          "forcingWidth must be positive for gaussian forcing");
    if (p.forcingProfile == ForcingProfile::singleMode) {
      const double mode =
          p.forcingWavenumber * p.domainLength / (2.0 * 3.14159265358979323846);
      if (std::abs(mode - std::round(mode)) > 1.e-12 ||
          mode >= static_cast<double>(p.gridPoints) / 2.0)
        throw std::runtime_error("singleMode forcingWavenumber must lie on the "
                                 "Fourier grid below Nyquist");
      if (p.targetWaveActionInjectionRate > 0.0)
        throw std::runtime_error(
            "targetWaveActionInjectionRate applies only to stochastic forcing");
    }
  }
  if (!p.initialConditionFile.empty() &&
      !std::filesystem::exists(p.initialConditionFile))
    throw std::runtime_error("initialConditionFile does not exist: " +
                             p.initialConditionFile.string());
  if (p.dataDirectory.empty() || p.outputDirectory.empty())
    throw std::runtime_error("output directories cannot be empty");
  const auto normalizedData =
      std::filesystem::absolute(p.dataDirectory).lexically_normal();
  const auto normalizedOutput =
      std::filesystem::absolute(p.outputDirectory).lexically_normal();
  if (normalizedData == normalizedOutput ||
      (std::filesystem::exists(normalizedData) &&
       std::filesystem::exists(normalizedOutput) &&
       std::filesystem::equivalent(normalizedData, normalizedOutput)))
    throw std::runtime_error("dataDirectory and outputDirectory must differ");
}

void writeParameterRecord(const Parameters &p,
                          const std::filesystem::path &directory,
                          const std::string &backend) {
  const auto path = directory / "resolved_parameters.txt";
  std::ofstream out(path);
  if (!out)
    throw std::runtime_error("cannot write parameter record: " + path.string());
  out << std::boolalpha << std::setprecision(17) << "backend " << backend
      << '\n'
      << "model " << modelName(p.model) << '\n'
      << "gridPoints " << p.gridPoints << '\n'
      << "domainLength " << p.domainLength << '\n'
      << "timeStep " << p.timeStep << '\n'
      << "numberOfSteps " << p.numberOfSteps << '\n'
      << "outputIntervalSteps " << p.outputIntervalSteps << '\n'
      << "integrator " << integratorName(p.integrator) << '\n'
      << "dispersionCoefficient " << p.dispersionCoefficient << '\n'
      << "nonlinearityCoefficient " << p.nonlinearityCoefficient << '\n'
      << "chemicalPotential " << p.chemicalPotential << '\n'
      << "helmholtzParameter " << p.helmholtzParameter << '\n'
      << "hyperviscosity " << p.hyperviscosity << '\n'
      << "hyperviscosityOrder " << p.hyperviscosityOrder << '\n'
      << "hypoviscosity " << p.hypoviscosity << '\n'
      << "hypoviscosityOrder " << p.hypoviscosityOrder << '\n'
      << "forcingEnabled " << p.forcingEnabled << '\n'
      << "forcingProfile " << forcingProfileName(p.forcingProfile) << '\n'
      << "forcingWavenumber " << p.forcingWavenumber << '\n'
      << "forcingWidth " << p.forcingWidth << '\n'
      << "forcingAmplitude " << p.forcingAmplitude << '\n'
      << "forcingShapeOrder " << p.forcingShapeOrder << '\n'
      << "forcingLogWidth " << p.forcingLogWidth << '\n'
      << "targetWaveActionInjectionRate " << p.targetWaveActionInjectionRate
      << '\n'
      << "randomSeed " << p.randomSeed << '\n'
      << "writeModeDiagnostics " << p.writeModeDiagnostics << '\n'
      << "threadCount " << p.threadCount << '\n'
      << "overwriteOutput " << p.overwriteOutput << '\n'
      << "initialConditionFile " << p.initialConditionFile.string() << '\n'
      << "dataDirectory " << p.dataDirectory.string() << '\n'
      << "outputDirectory " << p.outputDirectory.string() << '\n';
  out.close();
  if (!out)
    throw std::runtime_error("failed while writing parameter record: " +
                             path.string());
}
