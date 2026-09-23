#pragma once

#include "parameters.hpp"
#include "types.hpp"

#include <fftw3.h>

#include <cstddef>
#include <memory>
#include <vector>

struct FftwDeleter {
  void operator()(fftw_plan_s *plan) const;
};
using FftwPlan = std::unique_ptr<fftw_plan_s, FftwDeleter>;

class ComplexTransform {
public:
  explicit ComplexTransform(std::size_t count);
  ComplexTransform(const ComplexTransform &) = delete;
  ComplexTransform &operator=(const ComplexTransform &) = delete;

  void inverse(const SpectralField &spectrum, SpectralField &physical);
  void forward(const SpectralField &physical, SpectralField &spectrum);

private:
  std::size_t count_;
  SpectralField input_, output_;
  FftwPlan forward_, inverse_;
};

void initializeFftwThreads(int threadCount);
void finalizeFftwThreads();
