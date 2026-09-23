# 1D Schrödinger–Helmholtz pseudo-spectral solver

A C++20 solver for three related one-dimensional wave equations on a periodic
domain: the nonlocal Schrödinger–Helmholtz equation, its long-wave limit, and
the local cubic nonlinear Schrödinger equation. One `model` flag selects the
equation while the numerical method and output format remain the same.

The solver provides a fully dealiased pseudo-spectral nonlinear evaluation,
ETDRK2, ETDRK4-B, and integrating-factor RK2 time stepping, reproducible
stochastic forcing, atomic checkpoints, automatic restart, CSV diagnostics,
and plotting/movie tools.

## Equations

Let $\rho=|\psi|^2$ on the periodic interval $0\leq x<L$, and let
$k_m=2\pi m/L$. Each model includes the configured forcing and both
dissipation operators. The `model` flag changes only the nonlinear density
response.

### (i) Schrödinger–Helmholtz

Select `model schrodingerHelmholtz`. The full nonlocal response is

```math
i\,\partial_t\psi
=c\,\partial_{xx}\psi
+g\,\psi(1-\beta\partial_{xx})^{-1}|\psi|^2
+\mu\psi
-i\left[\alpha(-\partial_{xx})^q+\nu(-\partial_{xx})^p\right]\psi
+iF(x,t),
```

with

```math
\widehat{(1-\beta\partial_{xx})^{-1}\rho}_k
=\frac{\widehat{\rho}_k}{1+\beta k^2}.
```

The Helmholtz inverse suppresses short-scale density variations. Setting
`helmholtzParameter 0` (that is, $\beta=0$) reduces this response to the local
cubic one.

### (ii) Long-wave limit

Select `model longWave`. Expanding the Helmholtz inverse to first order gives

```math
i\,\partial_t\psi
=c\,\partial_{xx}\psi
+g\,\psi(1+\beta\partial_{xx})|\psi|^2
+\mu\psi
-i\left[\alpha(-\partial_{xx})^q+\nu(-\partial_{xx})^p\right]\psi
+iF(x,t),
```

where

```math
\widehat{(1+\beta\partial_{xx})\rho}_k
=(1-\beta k^2)\widehat{\rho}_k.
```

This approximation is appropriate when $\beta k^2\ll 1$ for the modes that
carry appreciable density. It is the first-order expansion of
$1/(1+\beta k^2)$; it should not be interpreted as an accurate replacement
for the full model at arbitrarily large wavenumber.

### (iii) One-dimensional nonlinear Schrödinger equation

Select `model nonlinearSchrodinger`. The density response is purely local:

```math
i\,\partial_t\psi
=c\,\partial_{xx}\psi
+g|\psi|^2\psi
+\mu\psi
-i\left[\alpha(-\partial_{xx})^q+\nu(-\partial_{xx})^p\right]\psi
+iF(x,t).
```

`helmholtzParameter` is ignored for this model.

### Dissipation and parameter symbols

In Fourier space the damping rate applied to mode $k$ is

```math
D_k=\alpha(k^2)^q+\nu(k^2)^p,
\qquad
\partial_t\widehat\psi_k\big|_{\mathrm{diss}}=-D_k\widehat\psi_k.
```

The symbols in the equations map directly to `params.txt`:

| Symbol | Parameter key | Meaning |
| --- | --- | --- |
| $L$, $N$, $\Delta t$ | `domainLength`, `gridPoints`, `timeStep` | Domain length, grid size, and fixed step |
| $c$ | `dispersionCoefficient` | Linear dispersion coefficient |
| $g$ | `nonlinearityCoefficient` | Nonlinear coupling |
| $\mu$ | `chemicalPotential` | Chemical-potential coefficient |
| $\beta$ | `helmholtzParameter` | Nonlocal length scale squared |
| $\nu,p$ | `hyperviscosity`, `hyperviscosityOrder` | Small-scale damping coefficient and power |
| $\alpha,q$ | `hypoviscosity`, `hypoviscosityOrder` | Large-scale damping coefficient and power |

When $q<0$, the singular zero-mode hypoviscous multiplier is suppressed.
Set both damping coefficients to zero for conservative dynamics.

### Forcing

For the four stochastic profiles, each retained nonzero Fourier mode obeys

```math
d\widehat\psi_k\big|_{\mathrm{force}}=f(k)\,dW_k,
\qquad
\mathbb{E}[dW_k]=0,
\qquad
\mathbb{E}[dW_k\,dW_{k'}^*]=\delta_{kk'}\,dt,
```

where the $W_k$ are independent circular complex Wiener processes. In the
formal physical-space equations above this is denoted by $F(x,t)$. The
available spectral envelopes are

```math
\begin{aligned}
f_{\mathrm{annulus}}(k)
  &=A\,\mathbf{1}_{\{\left|\lvert k\rvert-k_f\right|<\Delta k\}},\\
f_{\mathrm{gaussian}}(k)
  &=A\exp\!\left[-\frac12\left(\frac{|k|-k_f}{\sigma_f}\right)^2\right],\\
f_{\mathrm{exponential}}(k)
  &=A\left(\frac{|k|}{k_f}\right)^s
    \exp\!\left[-\left(\frac{|k|}{k_f}\right)^s\right],\\
f_{\mathrm{logNormal}}(k)
  &=A\exp\!\left[-\frac12
    \left(\frac{\log(|k|/k_f)}{\sigma_{\log}}\right)^2\right].
\end{aligned}
```

Here $A$, $k_f$, $\Delta k=\sigma_f$, $s$, and
$\sigma_{\log}$ are `forcingAmplitude`, `forcingWavenumber`,
`forcingWidth`, `forcingShapeOrder`, and `forcingLogWidth`. The new
`logNormal` option is the same log-space Gaussian envelope used by the 2D
Gross–Pitaevskii solver. A positive `targetWaveActionInjectionRate` rescales
the stochastic envelope so that $\sum_k|f(k)|^2$ equals the requested
spectral injection coefficient.

`singleMode` is deterministic rather than white noise. It sets
$\widehat F_k=A$ at both modes satisfying $|k|=k_f$ and zero elsewhere.
`randomSeed` controls stochastic reproducibility; the complete generator state
is included in checkpoints.

## Conservation laws and energy diagnostics

When forcing and dissipation are disabled (`forcingEnabled false`,
`hyperviscosity 0`, and `hypoviscosity 0`), all three equations conserve wave
action

```math
\mathcal{N} = \int |\psi|^2\,dx
```

and the model-specific Hamiltonian

```math
H = -c\int |\partial_x\psi|^2\,dx + \mu\mathcal{N} + H_{\mathrm{nl}}.
```

The nonlinear energy used by `diagnostics.csv` is, respectively,

```math
H_{\mathrm{nl}}^{\mathrm{SH}}
=\frac{g}{2}\int\rho(1-\beta\partial_{xx})^{-1}\rho\,dx,
```

```math
H_{\mathrm{nl}}^{\mathrm{LW}}
=\frac{g}{2}\int\left[\rho^2-\beta(\partial_x\rho)^2\right]dx,
```

and

```math
H_{\mathrm{nl}}^{\mathrm{NLS}}=\frac{g}{2}\int\rho^2\,dx.
```

The reported `total_energy` is the sum of `linear_energy`, `chemical_energy`,
and this model-specific `nonlinear_energy`. The energy-flux calculation uses
the same response multiplier, so switching `model` changes the evolution,
Hamiltonian, and flux diagnostic together. Automated conservative-run tests
check wave-action and energy drift independently for all three models.

## Spatial and temporal discretization

The spatial discretization uses complex FFTW transforms and a two-pass
3/2-rule evaluation: the wavefunction is embedded on the padded grid,
$|\psi|^2$ is transformed and truncated to the retained band, the Helmholtz
multiplier is applied, and the result is multiplied by $\psi$ before the final
transform. Fourier coefficients use the normalization

```math
\widehat{\psi}_k = \frac{1}{N}\sum_j \psi_j e^{-2\pi i j k/N}.
```

## Requirements

- CMake 3.20 or newer
- A C++20 compiler
- FFTW3 development headers and library
- Python 3 for the end-to-end regression test
- NumPy and Matplotlib for the analysis scripts

OpenMP and FFTW's threads library are optional.

## Build and test

```bash
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/release --parallel
ctest --test-dir build/release --output-on-failure
```

The convenience Makefile provides the same workflow:

```bash
make
make test
```

Use `-DSH1D_OPENMP=OFF` for a strictly serial build.

## Quick start

Run the small forced example from the repository root:

```bash
./build/release/schrodinger_helmholtz examples/quickstart.params
```

Equivalent small examples select the other two equations:

```bash
./build/release/schrodinger_helmholtz examples/long_wave.params
./build/release/schrodinger_helmholtz examples/nonlinear_schrodinger.params
```

It writes physical wavefunctions and checkpoints under `data/quickstart/`,
and diagnostics under `output/quickstart/`. A second invocation with the same
parameter file resumes automatically and performs another `numberOfSteps`.

The executable reads `params.txt` when no parameter path is supplied:

```bash
./build/release/schrodinger_helmholtz
```

## Parameters

Each nonempty line is `key value` or `key = value`; `#` starts a comment.
Unknown keys, invalid values, and extra fields stop the run.

| Key | Purpose |
| --- | --- |
| `gridPoints` | Base-grid size; a multiple of four from 8 through 16384 |
| `domainLength` | Periodic-domain length |
| `timeStep`, `numberOfSteps` | Fixed step and additional steps for this invocation |
| `outputIntervalSteps` | Stored-frame cadence; the final step is always saved |
| `model` | `schrodingerHelmholtz`, `longWave`, or `nonlinearSchrodinger` |
| `integrator` | `etd2`, `etd4`, or `rk2` |
| `dispersionCoefficient` | Coefficient `c` |
| `nonlinearityCoefficient` | Coefficient `g` |
| `chemicalPotential` | Coefficient `mu` |
| `helmholtzParameter` | `beta >= 0`; used by the full and long-wave models |
| `hyperviscosity`, `hyperviscosityOrder` | Small-scale damping coefficient and order |
| `hypoviscosity`, `hypoviscosityOrder` | Large-scale damping coefficient and order |
| `forcingEnabled` | Enable deterministic or stochastic forcing |
| `forcingProfile` | `annulus`, `gaussian`, `exponential`, `logNormal`, or `singleMode` |
| `forcingWavenumber`, `forcingWidth` | Forcing location and width in physical wavenumber |
| `forcingAmplitude`, `forcingShapeOrder` | Forcing magnitude and exponential-profile order |
| `forcingLogWidth` | Log-space standard deviation for `logNormal` forcing |
| `targetWaveActionInjectionRate` | Rescale stochastic forcing when positive |
| `randomSeed` | Reproducible 64-bit seed; zero selects and records a time-based seed |
| `initialConditionFile` | Optional `real imag` or `x real imag` physical field |
| `threadCount` | `1` or `2`; zero chooses automatically but is capped at two |
| `writeModeDiagnostics` | Write every complex Fourier mode to `modes.csv` |
| `overwriteOutput` | Replace managed output from an existing non-restart run |
| `dataDirectory`, `outputDirectory` | State and diagnostic locations |

Stochastic profiles use independent circular complex Gaussian increments.
Their covariance is integrated exactly through the linear damping over one
time step, and both random-generator states are stored in every checkpoint.
`singleMode` is deterministic and acts at both signs of the selected physical
wavenumber.

## Output and restart

Every fresh run saves frame zero, followed by the requested cadence and final
step. Output frames are numbered with eight digits.

| Location | Contents |
| --- | --- |
| `dataDirectory/wavefunction_NNNNNNNN.dat` | `x real imaginary` physical field |
| `dataDirectory/checkpoint_NNNNNNNN.bin` | Spectral state, time, frame, and random-generator state |
| `dataDirectory/restart_state.txt` | Latest committed checkpoint |
| `outputDirectory/diagnostics.csv` | Wave action, energy components, and dissipation rates |
| `outputDirectory/spectra.csv` | Instantaneous and running-average wave-action spectra |
| `outputDirectory/fluxes.csv` | Instantaneous and running-average wave-action/energy fluxes |
| `outputDirectory/modes.csv` | Optional complex Fourier modes |
| `outputDirectory/forcing_*.csv` | Resolved forcing profile and injection coefficients |
| `outputDirectory/segments/` | Per-invocation resolved parameters and start state |

Snapshots, checkpoints, and restart-state files are written through temporary
files and renamed atomically. On restart, CSV rows beyond the last committed
frame are removed before output resumes. A detected restart always takes
precedence over `initialConditionFile`; use new output directories for an
independent run.

## Plotting and movies

The [`scripts/`](scripts/) directory contains notebooks for the wavefunction,
spectrum, and time diagnostics, plus a command-line movie renderer. See
[`scripts/README.md`](scripts/README.md) for examples.

## Parallel execution

The solver is intentionally CPU-only and uses at most two FFTW/OpenMP threads.
`gridPoints` is limited to 16384. A one-dimensional run at the usual
$N=512$ performs small FFTs, so thread-management overhead outweighs the
saved arithmetic; one thread is the recommended default. Two threads can help
near the resolution limit.

For parameter scans or ensemble statistics, run independent simulations in
parallel with distinct `dataDirectory` and `outputDirectory` values. This
coarse-grained approach should scale much better than distributing a single
small 1D FFT. MPI and CUDA backends are intentionally outside the scope of
this solver.

## Repository layout

```text
src/          C++20 solver, parameter handling, FFTs, output, and restart
examples/     Small reproducible parameter file
tests/        Numerical and end-to-end restart regression tests
scripts/      Plotting notebooks, shared readers, and movie renderer
params.txt    Representative forced/dissipated run
```

