#pragma once

#include "parameters.hpp"
#include "types.hpp"

#include <fftw3.h>

#include <cstddef>
#include <limits>
#include <memory>
#include <span>
#include <vector>

struct FftwDeleter {
    void operator()(fftw_plan_s *plan) const;
};
using FftwPlan = std::unique_ptr<fftw_plan_s, FftwDeleter>;

template <class T> struct FftwAllocator {
    using value_type = T;
    FftwAllocator() = default;
    template <class U> FftwAllocator(const FftwAllocator<U> &) noexcept {}
    [[nodiscard]] T *allocate(std::size_t count) {
        if (count > std::numeric_limits<std::size_t>::max() / sizeof(T))
            throw std::bad_alloc();
        auto *data = static_cast<T *>(fftw_malloc(count * sizeof(T)));
        if (!data)
            throw std::bad_alloc();
        return data;
    }
    void deallocate(T *data, std::size_t) noexcept { fftw_free(data); }
    template <class U> bool operator==(const FftwAllocator<U> &) const noexcept { return true; }
};
using FftwField = std::vector<Complex, FftwAllocator<Complex>>;

// Buffers must keep their addresses until the returned plan is destroyed.
[[nodiscard]] FftwPlan makeFftwPlan(std::span<Complex> input, std::span<Complex> output, int sign);

class ComplexTransform {
  public:
    explicit ComplexTransform(std::size_t count);
    ComplexTransform(const ComplexTransform &) = delete;
    ComplexTransform &operator=(const ComplexTransform &) = delete;

    void inverse(const SpectralField &spectrum, SpectralField &physical);
    void forward(const SpectralField &physical, SpectralField &spectrum);

  private:
    std::size_t valueCount_;
    SpectralField inputBuffer_, outputBuffer_;
    FftwPlan forward_, inverse_;
};

void configureFftw(const Parameters &parameters, bool importWisdom = true);
void saveFftwWisdom();
[[nodiscard]] unsigned fftwPlanningFlags();
void initializeFftwThreads(int threadCount);
void finalizeFftwThreads();
