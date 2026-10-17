#pragma once

#include "fftw_utils.hpp"
#include "parameters.hpp"
#include "types.hpp"

class NonlinearOperator {
  public:
    explicit NonlinearOperator(const Parameters &parameters);

    // Returns the nonlinear contribution to d(psi_hat)/dt.
    void evaluate(const SpectralField &wavefunction, SpectralField &result);

    // Return the retained Fourier coefficients of |psi|^2 and its time
    // derivative. These raw density spectra are used by the diagnostics.
    void densitySpectrum(const SpectralField &wavefunction, SpectralField &result);
    void densityTimeDerivativeSpectrum(const SpectralField &wavefunction,
                                       const SpectralField &wavefunctionDot, SpectralField &result);

  private:
    void transformDensity(const SpectralField &wavefunction);
    void extractDensity(SpectralField &result) const;

    Parameters parameters_;
    // Reuse the spectral input for the final physical product, and the density
    // spectrum for the final nonlinear spectrum. Plans operate directly on
    // these aligned buffers, avoiding copies through a transform wrapper.
    FftwField paddedWavefunction_, physicalWavefunction_, densitySpectrum_, potential_;
    FftwField physicalWavefunctionDerivative_;
    std::vector<double> densityResponse_;
    FftwPlan inverseWavefunction_, forwardDensity_, inverseDensity_, forwardNonlinear_;
    FftwPlan inverseDerivative_;
};
