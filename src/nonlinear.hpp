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
    void applyDensityResponse(SpectralField &density) const;

    Parameters parameters_;
    ComplexTransform paddedTransform_;
    SpectralField paddedWavefunction_, physicalWavefunction_, density_, densitySpectrum_;
    SpectralField potential_, nonlinearPhysical_, nonlinearSpectrum_;
    SpectralField paddedWavefunctionDerivative_, physicalWavefunctionDerivative_;
    std::vector<double> densityResponse_;
};
