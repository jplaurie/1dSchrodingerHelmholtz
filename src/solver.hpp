#pragma once

#include "fftw_utils.hpp"
#include "host_stepper.hpp"
#include "integrator.hpp"
#include "nonlinear.hpp"
#include "output.hpp"
#include "parameters.hpp"
#include "types.hpp"

#include <cstdint>
#include <random>
#include <vector>

class Solver {
  public:
    explicit Solver(Parameters parameters);
    void run();
    [[nodiscard]] double benchmark(std::uint64_t warmupSteps, std::uint64_t measuredSteps);

  private:
    void buildLinearOperator();
    void buildIntegrationCoefficients();
    void buildForcing();
    void generateNoise();
    void step(SpectralField &wavefunction);
    RestartState prepareRun();
    void validateRunBounds(const RestartState &state) const;
    void restoreRandomState(RestartState &state);
    void writeState(const RestartState &state);
    double writeOutputFrame(const RestartState &state, DiagnosticsAverages &averages);
    [[nodiscard]] SpectralField makeBenchmarkState() const;

    Parameters parameters_;
    ComplexTransform baseTransform_;
    NonlinearOperator nonlinearOperator_;
    SpectralField linearOperator_;
    IntegrationCoefficients coefficients_;
    SpectralField noise_, deterministicForcing_;
    HostIntegrationWorkspace integrationWorkspace_;
    std::vector<double> forcingAmplitude_, stochasticNoiseScale_;
    std::vector<std::size_t> forcedIndices_;
    std::size_t forcedModeCount_ = 0;
    double waveActionInjectionCoefficient_ = 0.0;
    double quadraticEnergyInjectionCoefficient_ = 0.0;
    std::mt19937_64 random_;
    std::normal_distribution<double> normal_{0.0, 1.0};
};
