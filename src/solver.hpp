#pragma once

#include "fftw_utils.hpp"
#include "integrator.hpp"
#include "nonlinear.hpp"
#include "output.hpp"
#include "parameters.hpp"
#include "types.hpp"

#include <array>
#include <random>
#include <vector>

class Solver {
public:
  explicit Solver(Parameters parameters);
  void run();

private:
  void buildLinearOperator();
  void buildIntegrationCoefficients();
  void buildForcing();
  void generateNoise();
  void rightHandSide(const SpectralField &input, SpectralField &output);
  void step(SpectralField &wavefunction);
  double writeFrame(RestartState &state, DiagnosticsAverages *averages);

  Parameters p_;
  ComplexTransform baseTransform_;
  NonlinearOperator nonlinear_;
  SpectralField linear_;
  IntegrationCoefficients coefficients_;
  SpectralField noise_, deterministicForcing_;
  std::array<SpectralField, 4> nonlinearStages_;
  std::array<SpectralField, 3> stageStates_;
  std::vector<double> forcingAmplitude_, noiseScale_;
  std::vector<std::size_t> forcedIndices_;
  std::size_t forcedModeCount_ = 0;
  double waveActionInjectionCoefficient_ = 0.0;
  double quadraticEnergyInjectionCoefficient_ = 0.0;
  std::mt19937_64 random_;
  std::normal_distribution<double> normal_{0.0, 1.0};
};
