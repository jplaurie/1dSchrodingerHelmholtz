#include "parameters.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
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

template <class T> T parseNumber(const std::string &text, const std::string &key) {
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
bool parseMember(std::string_view key, const std::string &value, Parameters &parameters,
                 const std::pair<std::string_view, T Parameters::*> (&members)[N]) {
    const auto entry =
        std::find_if(std::begin(members), std::end(members),
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
    throw std::runtime_error("forcingProfile must be annulus, gaussian, exponential, logNormal, or "
                             "singleMode");
}

struct ParameterSetting {
    std::string key;
    std::string value;
    std::size_t lineNumber;
};

std::optional<ParameterSetting> parseParameterLine(std::string line, std::size_t lineNumber) {
    if (const auto comment = line.find('#'); comment != std::string::npos)
        line.erase(comment);
    line = trim(line);
    if (line.empty())
        return std::nullopt;

    std::replace(line.begin(), line.end(), '=', ' ');
    std::istringstream fields(line);
    ParameterSetting setting{{}, {}, lineNumber};
    std::string extra;
    fields >> setting.key >> setting.value;
    if (setting.key.empty() || setting.value.empty() || (fields >> extra))
        throw std::runtime_error("invalid parameter line " + std::to_string(lineNumber));
    return setting;
}

void applyParameter(const ParameterSetting &setting, Parameters &parameters) {
    static constexpr std::pair<std::string_view, std::size_t Parameters::*> sizes[]{
        {"gridPoints", &Parameters::gridPoints}};
    static constexpr std::pair<std::string_view, std::uint64_t Parameters::*> counts[]{
        {"numberOfSteps", &Parameters::numberOfSteps},
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
        {"targetWaveActionInjectionRate", &Parameters::targetWaveActionInjectionRate}};
    static constexpr std::pair<std::string_view, bool Parameters::*> booleans[]{
        {"forcingEnabled", &Parameters::forcingEnabled},
        {"writeModeDiagnostics", &Parameters::writeModeDiagnostics},
        {"overwriteOutput", &Parameters::overwriteOutput}};
    static constexpr std::pair<std::string_view, int Parameters::*> integers[]{
        {"threadCount", &Parameters::threadCount}};
    static constexpr std::pair<std::string_view, std::filesystem::path Parameters::*> paths[]{
        {"initialConditionFile", &Parameters::initialConditionFile},
        {"dataDirectory", &Parameters::dataDirectory},
        {"outputDirectory", &Parameters::outputDirectory}};

    const std::string &key = setting.key;
    const std::string &value = setting.value;
    const bool recognized =
        parseMember(key, value, parameters, sizes) || parseMember(key, value, parameters, counts) ||
        parseMember(key, value, parameters, reals) ||
        parseMember(key, value, parameters, booleans) ||
        parseMember(key, value, parameters, integers) || parseMember(key, value, parameters, paths);
    if (key == "model")
        parameters.model = parseModel(value);
    else if (key == "integrator")
        parameters.integrator = parseIntegrator(value);
    else if (key == "forcingProfile")
        parameters.forcingProfile = parseForcingProfile(value);
    else if (!recognized)
        throw std::runtime_error("unknown parameter key on line " +
                                 std::to_string(setting.lineNumber) + ": " + key);
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
    Parameters parameters;
    std::string line;
    std::size_t lineNumber = 0;
    while (std::getline(input, line)) {
        ++lineNumber;
        if (const auto setting = parseParameterLine(line, lineNumber))
            applyParameter(*setting, parameters);
    }
    validateParameters(parameters);
    return parameters;
}

void validateParameters(const Parameters &parameters) {
    if (parameters.gridPoints < 8 || parameters.gridPoints % 4 != 0 ||
        parameters.gridPoints > 16'384)
        throw std::runtime_error("gridPoints must be a multiple of four between 8 and 16384");
    if (parameters.gridPoints > 2 * (static_cast<std::size_t>(std::numeric_limits<int>::max()) / 3))
        throw std::runtime_error("3/2-rule grid exceeds FFTW integer limits");
    if (!(parameters.domainLength > 0.0) || !isFinite(parameters.domainLength))
        throw std::runtime_error("domainLength must be finite and positive");
    if (!(parameters.timeStep > 0.0) || !isFinite(parameters.timeStep))
        throw std::runtime_error("timeStep must be finite and positive");
    if (parameters.numberOfSteps == 0 || parameters.outputIntervalSteps == 0)
        throw std::runtime_error("numberOfSteps and outputIntervalSteps must be positive");
    if (parameters.threadCount < 0)
        throw std::runtime_error("threadCount must be nonnegative");
    for (const auto [value, name] :
         {std::pair{parameters.dispersionCoefficient, "dispersionCoefficient"},
          std::pair{parameters.nonlinearityCoefficient, "nonlinearityCoefficient"},
          std::pair{parameters.chemicalPotential, "chemicalPotential"},
          std::pair{parameters.helmholtzParameter, "helmholtzParameter"},
          std::pair{parameters.hyperviscosity, "hyperviscosity"},
          std::pair{parameters.hyperviscosityOrder, "hyperviscosityOrder"},
          std::pair{parameters.hypoviscosity, "hypoviscosity"},
          std::pair{parameters.hypoviscosityOrder, "hypoviscosityOrder"}})
        if (!isFinite(value))
            throw std::runtime_error(std::string(name) + " must be finite");
    if (parameters.helmholtzParameter < 0.0)
        throw std::runtime_error("helmholtzParameter cannot be negative");
    if (parameters.hyperviscosity < 0.0 || parameters.hypoviscosity < 0.0 ||
        parameters.hyperviscosityOrder < 0.0)
        throw std::runtime_error("dissipation coefficients/orders are invalid");
    if (parameters.forcingEnabled) {
        if (!(parameters.forcingWavenumber > 0.0) || parameters.forcingWidth < 0.0 ||
            parameters.forcingAmplitude < 0.0 || parameters.forcingShapeOrder <= 0.0 ||
            parameters.forcingLogWidth <= 0.0 || parameters.targetWaveActionInjectionRate < 0.0 ||
            !isFinite(parameters.forcingWavenumber) || !isFinite(parameters.forcingWidth) ||
            !isFinite(parameters.forcingAmplitude) || !isFinite(parameters.forcingShapeOrder) ||
            !isFinite(parameters.forcingLogWidth) ||
            !isFinite(parameters.targetWaveActionInjectionRate))
            throw std::runtime_error("forcing parameters are invalid or non-finite");
        if (parameters.forcingProfile == ForcingProfile::gaussian && parameters.forcingWidth <= 0.0)
            throw std::runtime_error("forcingWidth must be positive for gaussian forcing");
        if (parameters.forcingProfile == ForcingProfile::singleMode) {
            const double mode = parameters.forcingWavenumber * parameters.domainLength /
                                (2.0 * 3.14159265358979323846);
            if (std::abs(mode - std::round(mode)) > 1.e-12 ||
                mode >= static_cast<double>(parameters.gridPoints) / 2.0)
                throw std::runtime_error("singleMode forcingWavenumber must lie on the "
                                         "Fourier grid below Nyquist");
            if (parameters.targetWaveActionInjectionRate > 0.0)
                throw std::runtime_error(
                    "targetWaveActionInjectionRate applies only to stochastic forcing");
        }
    }
    if (!parameters.initialConditionFile.empty() &&
        !std::filesystem::exists(parameters.initialConditionFile))
        throw std::runtime_error("initialConditionFile does not exist: " +
                                 parameters.initialConditionFile.string());
    if (parameters.dataDirectory.empty() || parameters.outputDirectory.empty())
        throw std::runtime_error("output directories cannot be empty");
    const auto normalizedData =
        std::filesystem::absolute(parameters.dataDirectory).lexically_normal();
    const auto normalizedOutput =
        std::filesystem::absolute(parameters.outputDirectory).lexically_normal();
    if (normalizedData == normalizedOutput ||
        (std::filesystem::exists(normalizedData) && std::filesystem::exists(normalizedOutput) &&
         std::filesystem::equivalent(normalizedData, normalizedOutput)))
        throw std::runtime_error("dataDirectory and outputDirectory must differ");
}

void writeParameterRecord(const Parameters &parameters, const std::filesystem::path &directory,
                          const std::string &backend) {
    const auto path = directory / "resolved_parameters.txt";
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("cannot write parameter record: " + path.string());
    out << std::boolalpha << std::setprecision(17) << "backend " << backend << '\n'
        << "model " << modelName(parameters.model) << '\n'
        << "gridPoints " << parameters.gridPoints << '\n'
        << "domainLength " << parameters.domainLength << '\n'
        << "timeStep " << parameters.timeStep << '\n'
        << "numberOfSteps " << parameters.numberOfSteps << '\n'
        << "outputIntervalSteps " << parameters.outputIntervalSteps << '\n'
        << "integrator " << integratorName(parameters.integrator) << '\n'
        << "dispersionCoefficient " << parameters.dispersionCoefficient << '\n'
        << "nonlinearityCoefficient " << parameters.nonlinearityCoefficient << '\n'
        << "chemicalPotential " << parameters.chemicalPotential << '\n'
        << "helmholtzParameter " << parameters.helmholtzParameter << '\n'
        << "hyperviscosity " << parameters.hyperviscosity << '\n'
        << "hyperviscosityOrder " << parameters.hyperviscosityOrder << '\n'
        << "hypoviscosity " << parameters.hypoviscosity << '\n'
        << "hypoviscosityOrder " << parameters.hypoviscosityOrder << '\n'
        << "forcingEnabled " << parameters.forcingEnabled << '\n'
        << "forcingProfile " << forcingProfileName(parameters.forcingProfile) << '\n'
        << "forcingWavenumber " << parameters.forcingWavenumber << '\n'
        << "forcingWidth " << parameters.forcingWidth << '\n'
        << "forcingAmplitude " << parameters.forcingAmplitude << '\n'
        << "forcingShapeOrder " << parameters.forcingShapeOrder << '\n'
        << "forcingLogWidth " << parameters.forcingLogWidth << '\n'
        << "targetWaveActionInjectionRate " << parameters.targetWaveActionInjectionRate << '\n'
        << "randomSeed " << parameters.randomSeed << '\n'
        << "writeModeDiagnostics " << parameters.writeModeDiagnostics << '\n'
        << "threadCount " << parameters.threadCount << '\n'
        << "overwriteOutput " << parameters.overwriteOutput << '\n'
        << "initialConditionFile " << parameters.initialConditionFile.string() << '\n'
        << "dataDirectory " << parameters.dataDirectory.string() << '\n'
        << "outputDirectory " << parameters.outputDirectory.string() << '\n';
    out.close();
    if (!out)
        throw std::runtime_error("failed while writing parameter record: " + path.string());
}
