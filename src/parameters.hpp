#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>

enum class Model : std::uint8_t {
  schrodingerHelmholtz,
  longWave,
  nonlinearSchrodinger
};
enum class Integrator : std::uint8_t { etd2, etd4, integratingFactorRk2 };
enum class ForcingProfile : std::uint8_t {
  annulus,
  gaussian,
  exponential,
  logNormal,
  singleMode
};

struct Parameters {
  std::size_t gridPoints = 512;
  double domainLength = 6.2831853071795864769;
  double timeStep = 5.0e-5;
  std::uint64_t numberOfSteps = 1000;
  std::uint64_t outputIntervalSteps = 100;

  Model model = Model::schrodingerHelmholtz;
  Integrator integrator = Integrator::etd4;
  double dispersionCoefficient = -0.5;
  double nonlinearityCoefficient = -1.0;
  double chemicalPotential = 0.0;
  double helmholtzParameter = 1.0;

  double hyperviscosity = 0.0;
  double hyperviscosityOrder = 8.0;
  double hypoviscosity = 0.0;
  double hypoviscosityOrder = -2.0;

  bool forcingEnabled = false;
  ForcingProfile forcingProfile = ForcingProfile::annulus;
  double forcingWavenumber = 64.0;
  double forcingWidth = 2.0;
  double forcingAmplitude = 1.0;
  double forcingShapeOrder = 4.0;
  double forcingLogWidth = 0.1;
  double targetWaveActionInjectionRate = 0.0;
  std::uint64_t randomSeed = 1;

  bool writeModeDiagnostics = false;
  int threadCount = 1;
  bool overwriteOutput = false;
  std::filesystem::path initialConditionFile;
  std::filesystem::path dataDirectory = "data";
  std::filesystem::path outputDirectory = "output";

  [[nodiscard]] std::size_t paddedGridPoints() const {
    return 3 * gridPoints / 2;
  }
  [[nodiscard]] bool usesStochasticForcing() const {
    return forcingEnabled && forcingProfile != ForcingProfile::singleMode;
  }
  [[nodiscard]] bool usesDeterministicForcing() const {
    return forcingEnabled && forcingProfile == ForcingProfile::singleMode;
  }
};

[[nodiscard]] const char *modelName(Model model);
[[nodiscard]] const char *integratorName(Integrator integrator);
[[nodiscard]] const char *forcingProfileName(ForcingProfile profile);
Parameters readParameters(const std::filesystem::path &path);
void validateParameters(const Parameters &parameters);
void writeParameterRecord(const Parameters &parameters,
                          const std::filesystem::path &directory,
                          const std::string &backend);
