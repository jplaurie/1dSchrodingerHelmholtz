#include "fftw_utils.hpp"

#include <algorithm>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {
#ifdef SH1D_HAVE_FFTW_THREADS
bool threadsInitialized = false;
#endif

fftw_complex *data(SpectralField &field) {
  static_assert(sizeof(Complex) == sizeof(fftw_complex));
  return reinterpret_cast<fftw_complex *>(field.data());
}
} // namespace

void FftwDeleter::operator()(fftw_plan_s *plan) const {
  if (plan)
    fftw_destroy_plan(plan);
}

ComplexTransform::ComplexTransform(std::size_t count)
    : count_(count), input_(count), output_(count) {
  forward_.reset(fftw_plan_dft_1d(static_cast<int>(count_), data(input_),
                                  data(output_), FFTW_FORWARD, FFTW_ESTIMATE));
  inverse_.reset(fftw_plan_dft_1d(static_cast<int>(count_), data(input_),
                                  data(output_), FFTW_BACKWARD, FFTW_ESTIMATE));
  if (!forward_ || !inverse_)
    throw std::runtime_error("FFTW could not create transform plans");
}

void ComplexTransform::inverse(const SpectralField &spectrum,
                               SpectralField &physical) {
  if (spectrum.size() != count_)
    throw std::runtime_error("invalid inverse-transform input size");
  input_ = spectrum;
  fftw_execute(inverse_.get());
  physical = output_;
}

void ComplexTransform::forward(const SpectralField &physical,
                               SpectralField &spectrum) {
  if (physical.size() != count_)
    throw std::runtime_error("invalid forward-transform input size");
  input_ = physical;
  fftw_execute(forward_.get());
  spectrum = output_;
  const double scale = 1.0 / static_cast<double>(count_);
  for (Complex &value : spectrum)
    value *= scale;
}

void initializeFftwThreads(int threadCount) {
  int threads = threadCount;
#ifdef _OPENMP
  if (threads == 0)
    threads = std::min(2, omp_get_max_threads());
  omp_set_num_threads(std::max(1, threads));
#else
  if (threads == 0)
    threads = 1;
#endif
  if (threads > 1) {
#ifdef SH1D_HAVE_FFTW_THREADS
    if (!fftw_init_threads())
      throw std::runtime_error("FFTW thread initialization failed");
    threadsInitialized = true;
    fftw_plan_with_nthreads(threads);
#else
    throw std::runtime_error(
        "threadCount > 1 requested, but threaded FFTW is unavailable");
#endif
  }
}

void finalizeFftwThreads() {
#ifdef SH1D_HAVE_FFTW_THREADS
  if (threadsInitialized) {
    fftw_cleanup_threads();
    threadsInitialized = false;
  }
#endif
}
