#pragma once

#include "fftw_utils.hpp"
#include "nonlinear.hpp"
#include "parameters.hpp"
#include "types.hpp"

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

struct RestartState {
  double time = 0.0;
  std::uint64_t frame = 0;
  std::uint64_t randomSeed = 0;
  SpectralField wavefunction;
  std::string randomEngineState;
  std::string normalDistributionState;
  bool restarting = false;
};

struct DiagnosticsAverages {
  std::uint64_t count = 0;
  std::vector<double> spectrum;
  std::vector<double> waveFlux;
  std::vector<double> energyFlux;
};

RestartState readRestartOrInitial(const Parameters &parameters,
                                  ComplexTransform &baseTransform);
void prepareOutput(const Parameters &parameters, bool restarting,
                   std::uint64_t committedFrame);
void writeRunRecords(const Parameters &parameters, double startTime,
                     std::uint64_t startFrame,
                     const std::vector<double> &forcingAmplitude,
                     std::size_t forcedModeCount,
                     double waveActionInjectionCoefficient,
                     double quadraticEnergyInjectionCoefficient);
void writeWavefunctionAndRestart(const Parameters &parameters,
                                 ComplexTransform &baseTransform,
                                 const RestartState &state,
                                 const std::string &randomEngineState,
                                 const std::string &normalDistributionState);
double appendDiagnostics(const Parameters &parameters,
                         NonlinearOperator &nonlinear, double time,
                         std::uint64_t frame, const SpectralField &wavefunction,
                         DiagnosticsAverages &averages);
