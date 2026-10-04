#include "fftw_utils.hpp"
#include "parameters.hpp"
#include "solver.hpp"

#include <charconv>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {
template <class Integer> Integer parseInteger(const char *text, const char *name) {
    Integer value{};
    const std::string_view input(text);
    const auto [end, error] = std::from_chars(input.data(), input.data() + input.size(), value);
    if (error != std::errc{} || end != input.data() + input.size())
        throw std::runtime_error(std::string("invalid ") + name + ": " + text);
    return value;
}
} // namespace

int main(int argc, char **argv) {
    try {
        if (argc < 3 || argc > 5)
            throw std::runtime_error(
                "usage: benchmark <resolution> <measured-steps> [warmup-steps] [threads]");

        Parameters parameters;
        parameters.gridPoints = parseInteger<std::size_t>(argv[1], "resolution");
        const std::uint64_t measuredSteps =
            parseInteger<std::uint64_t>(argv[2], "measured step count");
        const std::uint64_t warmupSteps =
            argc >= 4 ? parseInteger<std::uint64_t>(argv[3], "warmup step count") : 3;
        parameters.threadCount = argc >= 5 ? parseInteger<int>(argv[4], "thread count") : 1;
        parameters.model = Model::schrodingerHelmholtz;
        parameters.integrator = Integrator::etd4;
        parameters.timeStep = 1.0e-5;
        parameters.numberOfSteps = measuredSteps;
        parameters.outputIntervalSteps = measuredSteps;
        parameters.dispersionCoefficient = -0.5;
        parameters.nonlinearityCoefficient = -1.0;
        parameters.chemicalPotential = 0.0;
        parameters.helmholtzParameter = 1.0;
        parameters.hyperviscosity = 0.0;
        parameters.hypoviscosity = 0.0;
        parameters.forcingEnabled = false;
        parameters.randomSeed = 1;
        validateParameters(parameters);

        initializeFftwThreads(parameters.threadCount);
        double elapsed = 0.0;
        {
            Solver solver(parameters);
            elapsed = solver.benchmark(warmupSteps, measuredSteps);
        }
        finalizeFftwThreads();
#ifdef SH1D_BENCHMARK_SERIAL
        constexpr const char *backend = "CPU serial";
#else
        constexpr const char *backend = "CPU/OpenMP";
#endif
        std::cout << std::setprecision(17) << "BENCHMARK," << backend << ','
                  << parameters.gridPoints << ',' << parameters.threadCount << ',' << measuredSteps
                  << ',' << elapsed << ',' << elapsed / static_cast<double>(measuredSteps) << '\n';
        return 0;
    } catch (const std::exception &error) {
        finalizeFftwThreads();
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
