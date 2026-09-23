#include "fftw_utils.hpp"
#include "parameters.hpp"
#include "solver.hpp"

#include <exception>
#include <filesystem>
#include <iostream>

int main(int argc, char **argv) {
  try {
    if (argc > 2) {
      std::cerr << "usage: " << argv[0] << " [parameter-file]\n";
      return 2;
    }
    const std::filesystem::path parameterFile =
        argc == 2 ? argv[1] : "params.txt";
    Parameters parameters = readParameters(parameterFile);
    initializeFftwThreads(parameters.threadCount);
    {
      Solver solver(std::move(parameters));
      solver.run();
    }
    finalizeFftwThreads();
    return 0;
  } catch (const std::exception &error) {
    finalizeFftwThreads();
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
