#include "fftw_utils.hpp"

#include <algorithm>
#include <filesystem>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {
unsigned planningFlags = FFTW_ESTIMATE;
std::filesystem::path wisdomPath;

#ifdef SH1D_HAVE_FFTW_THREADS
bool threadsInitialized = false;
#endif

fftw_complex *fftwData(SpectralField &field) {
    static_assert(sizeof(Complex) == sizeof(fftw_complex));
    return reinterpret_cast<fftw_complex *>(field.data());
}
} // namespace

void configureFftw(const Parameters &parameters, bool importWisdom) {
    switch (parameters.fftwPlanning) {
    case FftwPlanning::estimate:
        planningFlags = FFTW_ESTIMATE;
        break;
    case FftwPlanning::measure:
        planningFlags = FFTW_MEASURE;
        break;
    case FftwPlanning::patient:
        planningFlags = FFTW_PATIENT;
        break;
    }
    wisdomPath = parameters.fftwWisdomFile;
    if (importWisdom && !wisdomPath.empty() && std::filesystem::exists(wisdomPath) &&
        !fftw_import_wisdom_from_filename(wisdomPath.c_str()))
        throw std::runtime_error("cannot import FFTW wisdom: " + wisdomPath.string());
}

void saveFftwWisdom() {
    if (!wisdomPath.empty() && !fftw_export_wisdom_to_filename(wisdomPath.c_str()))
        throw std::runtime_error("cannot export FFTW wisdom: " + wisdomPath.string());
}

unsigned fftwPlanningFlags() { return planningFlags; }

void FftwDeleter::operator()(fftw_plan_s *plan) const {
    if (plan)
        fftw_destroy_plan(plan);
}

ComplexTransform::ComplexTransform(std::size_t count)
    : valueCount_(count), inputBuffer_(count), outputBuffer_(count) {
    forward_.reset(fftw_plan_dft_1d(static_cast<int>(valueCount_), fftwData(inputBuffer_),
                                    fftwData(outputBuffer_), FFTW_FORWARD, fftwPlanningFlags()));
    inverse_.reset(fftw_plan_dft_1d(static_cast<int>(valueCount_), fftwData(inputBuffer_),
                                    fftwData(outputBuffer_), FFTW_BACKWARD, fftwPlanningFlags()));
    if (!forward_ || !inverse_)
        throw std::runtime_error("FFTW could not create transform plans");
}

void ComplexTransform::inverse(const SpectralField &spectrum, SpectralField &physical) {
    if (spectrum.size() != valueCount_)
        throw std::runtime_error("invalid inverse-transform input size");
    inputBuffer_ = spectrum;
    fftw_execute(inverse_.get());
    physical = outputBuffer_;
}

void ComplexTransform::forward(const SpectralField &physical, SpectralField &spectrum) {
    if (physical.size() != valueCount_)
        throw std::runtime_error("invalid forward-transform input size");
    inputBuffer_ = physical;
    fftw_execute(forward_.get());
    spectrum = outputBuffer_;
    const double scale = 1.0 / static_cast<double>(valueCount_);
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
        throw std::runtime_error("threadCount > 1 requested, but threaded FFTW is unavailable");
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
