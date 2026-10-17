# 1D Schrödinger–Helmholtz audit

Date: 2026-10-17. Reviewed commit: `62d60c4`.

Scope: numerical projection, integrators, forcing, diagnostics, restart/output,
HDF5, plotting, and performance-sensitive paths. The findings below describe
the original code; their locations refer to the reviewed commit. P1 means high priority; P2 means
normal priority. Performance opportunities are distinguished from reproduced
correctness defects below.

## Correction status

All six findings and both performance opportunities have been addressed:

- Every integration-stage RHS is projected; nonlinear results and forcing use
  the same admissible wavefunction band.
- HDF5 readers validate scalar dataspaces and numeric/string types before
  reading and handle variable-length strings safely. Plotting decodes bytes.
- Version-2 checkpoints persist full-precision diagnostic sums. Legacy
  checkpoints reconstruct sums from instantaneous CSV history once.
- CSV recovery streams with bounded memory and avoids rewriting clean files.
- Zero coupling skips nonlinear FFTs. Nonlinear FFT plans now operate on
  aligned reusable buffers directly, with normalization fused into response
  and final extraction scaling.
- The README documents the projected Hamiltonian, mean removal, restart
  averages, checkpoint compatibility, and trajectory changes.

Regression coverage includes dense cutoff convolutions and invariant transfer
closure for all three models, FFTW estimate/measure planning, independent
nonuniform temporal convergence for all integrators with/without mean removal,
discarded forcing support and normalization, v1/v2 restart averages, streamed
CSV trimming and clean-file preservation, malformed HDF5 scalar/string
attributes, and Python plotting of C++ HDF5-only snapshots.

Post-fix verification: all **9 Release tests**, **8 ASan/UBSan tests** (including
LeakSanitizer outside the sandbox), and **5 tests with HDF5 and OpenMP disabled**
passed. Nonuniform ETDRK4 orders range from 3.998 to 4.004 across all models
and mean-removal settings; ETD2/RK2 orders range from 1.996 to 2.003.

A local one-thread Release comparison used the existing complete-timestep
benchmark, three alternating before/after trials per resolution, and median
times. Planning, startup, diagnostics, and output were excluded:

| Grid points | Before, microseconds/step | After, microseconds/step | Throughput speedup |
| --- | ---: | ---: | ---: |
| 512 | 28.64 | 25.67 | 1.12x |
| 4096 | 394.18 | 354.27 | 1.11x |
| 16384 | 2720.86 | 2080.29 | 1.31x |

Zero-coupling nonlinear evaluation now costs approximately 0.027, 0.29, and
1.05 microseconds at these resolutions, versus the original approximately
6, 94, and 640 microseconds. These are local measurements, not portable
speedup guarantees. The original published benchmark figure remains a record
of its historical build.

## Verification

- A fresh Release build in `/tmp/sh1d-audit-build` completed successfully.
- All six existing CTest tests passed.
- Independent, non-wrapping convolutions tested densely populated spectra up to
  the retained-band boundary for all three models.
- End-to-end time convergence was compared with a separately implemented
  projected Fourier ODE integrated with SciPy DOP853, using tight tolerances.
- A separate AddressSanitizer/UndefinedBehaviorSanitizer build reproduced an
  HDF5 attribute buffer overflow.
- Small isolated solver runs reproduced the forcing, plotting, and average
  issues below. Existing simulation output was not used or modified.

## Projection: the intended density treatment is sound

The local `../2dGrossPitaevskii` implementation and its documentation explicitly
use the density-projected formulation. The analogous 1D nonlinear RHS is

\[
  F(\psi)=-ig\,P_\psi\left[\psi\,R\,P_\rho(|\psi|^2)\right],
\]

where the response multiplier is respectively
`1/(1 + beta*k^2)`, `1 - beta*k^2`, or `1`. Both density cutoff endpoints are
excluded, preserving a real density response. The density filter itself should
be retained. This formulation is described for truncated GP in
[Krstulovic and Brachet, Eq. 18](https://arxiv.org/pdf/1010.0116).

For valid 1D states, the computed nonlinear coefficients in retained modes
agreed with independent density-projected convolution to within `3.2e-17` in
the tested dense eight-point case. With the wavefunction projector applied to
the RHS, wave-action, momentum, and projected-Hamiltonian rates closed to
within `2.5e-16` before multiplying by the domain length, across all three models.

Unlike the current 2DGP implementation, 1D explicitly excludes its wavefunction
Nyquist mode; it also excludes the mean when negative-order hypoviscosity is
active. These restrictions must hold during intermediate stages as well.

## Findings

### F1 — P1: omitted stage projection reduces temporal order to one

Locations: `src/nonlinear.cpp:48`, `src/host_stepper.hpp:44`,
`src/host_stepper.hpp:71`, `src/spectral.cpp:65`.

`NonlinearOperator::evaluate` returns a generally nonzero Nyquist coefficient.
The stepper uses that coefficient to construct intermediate states and only
calls `enforceConstraints` after completing the step. The forbidden mode then
participates in later nonlinear evaluations and changes retained modes before
being discarded. The same mechanism affects the mean when hypoviscosity
requires it to vanish.

End-to-end ETDRK4 test: `N=8`, `L=2*pi`, `g=1.3`, `beta=0.35`, `c=-0.5`,
`mu=0`, no forcing or damping, final time `0.4`. Initial Fourier coefficients
are `0.22*cos(0.71*i) + 0.17j*sin(0.39*i)` for indices `i=0,...,7`, with
Nyquist set to zero. The maximum spectral errors were:

| Model | 16 steps | 32 steps | 64 steps | 128 steps | Final observed order |
| --- | ---: | ---: | ---: | ---: | ---: |
| Schrödinger–Helmholtz | 3.269e-5 | 1.636e-5 | 8.187e-6 | 4.094e-6 | 1.000 |
| Long-wave | 1.106e-5 | 5.509e-6 | 2.749e-6 | 1.373e-6 | 1.001 |
| NLS | 7.197e-5 | 3.607e-5 | 1.805e-5 | 9.031e-6 | 0.999 |

A separate `N=32` SH test populated only modes `-2,-1,1,2` and used
`hypoviscosity=0.1`, `hypoviscosityOrder=-2`. ETDRK4 again showed order `0.999`.
This second trigger is relevant to the default parameter file, which removes
the mean. ETD2 and integrating-factor RK2 also exhibited order loss.

In an isolated zero-linear-operator test using the actual host stepper,
projecting each RHS restored ETDRK4: errors fell from `2.291e-10` at 16 steps
to `1.433e-11` at 32 steps and `8.956e-13` at 64 steps. The unchanged stepper
gave `4.857e-5`, `2.428e-5`, and `1.214e-5` on the same test.

Correction: use the same wavefunction projector consistently for nonlinear
RHS values, forcing, and stage states. Enforcing constraints on each stage
before its RHS evaluation also prevents the feedback. Preserve the intended
1D band convention rather than copying the different 2D wavefunction band
without updating the rest of the solver. Add cutoff-populated convolution and
nonuniform temporal-convergence tests, including mean removal.

### F2 — P1: HDF5 attribute shapes can cause memory corruption

Locations: `src/hdf5_io.cpp:82` and `src/hdf5_io.cpp:91`.

`readScalarAttribute` passes the address of one scalar to `H5Aread` without
checking the attribute dataspace. HDF5 writes every attribute element into
that buffer. Replacing `grid_points` with a 32-element unsigned-integer
attribute produced an AddressSanitizer-confirmed **256-byte write into an
8-byte stack variable** at line 86. This affects the exporter and HDF5 initial
condition loading.

The string reader has the analogous element-count issue: it allocates space
for one datatype element without checking that the attribute is scalar. It
also assumes a fixed-length string datatype.

Correction: inspect and validate the attribute dataspace before reading, and
validate string type/class and variable-length handling. Reject malformed
attributes cleanly. Add malformed scalar-array and string-array fixtures.

### F3 — P2: forcing normalization includes a discarded Nyquist mode

Location: `src/solver.cpp:184`, interacting with `src/spectral.cpp:65`.

`buildForcing` evaluates every base-grid index, including the excluded
wavefunction Nyquist mode. Its amplitude contributes to `forcedModeCount`
and target-injection normalization, but the resulting noise is deleted by
the state constraint.

Reproduction: `N=8`, `L=2*pi`, zero initial state, `g=c=mu=0`, no damping,
annulus forcing with `forcingWavenumber=4`, `forcingWidth=0.1`, amplitude 1,
and target injection 1. The run succeeds and records one forced mode and
injection coefficient 1, while the field and wave action remain exactly zero.
Broad Gaussian and other profiles can also allocate some forcing to this
discarded mode.

Correction: apply the admissible-wavefunction mask before counting modes,
normalizing amplitudes, and constructing noise. Reject profiles whose
admissible support is empty.

### F4 — P2: plotting rejects the solver's HDF5 snapshots

Location: `scripts/sh1d_plotting.py:45`.

The C++ writer stores `format` as a fixed-length HDF5 string. In the installed
h5py, it loads as `np.bytes_(b'sh1d_wavefunction_v1')`. Comparing this directly
with the Python string always fails. `load_wavefunction` raised
`ValueError: unsupported HDF5 wavefunction format` on a freshly generated
HDF5-only snapshot. Notebooks and movies using the helper inherit the failure;
`both` output can hide it by loading the text snapshot first.

Correction: decode byte-valued attributes before comparison and test the
Python helper on a file generated by the C++ solver.

### F5 — P2: running averages silently restart with each invocation

Locations: `src/solver.cpp:351`, `src/diagnostics.cpp:102`.

Every invocation constructs empty `DiagnosticsAverages`, although CSV output
is appended to the previous run. Consequently the advertised running-average
columns change their averaging window at restart. The checkpoint does not
store the sums or sample count.

Reproduction: zero state, `g=c=mu=0`, no damping, deterministic forcing of
amplitude 1 at modes `+/-1`, `dt=0.1`, one output per step. Two one-step
invocations give instantaneous shell spectra `0.02` then `0.08`; frame two
reports average `0.08`, whereas an uninterrupted two-step run's average is
`0.05`.

Correction: restore sums/counts from checkpointed statistics or committed
output. If segment-local averages are intended, explicitly label and document
that scope in the data. Add restart-equivalence checks for averages as well
as the wavefunction.

### F6 — P2: restart recovery scales memory with entire CSV size

Location: `src/output.cpp:224`.

Every restart calls `trimCsv`, which loads all retained rows into
`vector<string>` and rewrites the complete file, even when there are no rows
to remove. Long spectra/mode histories can require many gigabytes of memory
just to resume a small simulation. This is established by the implementation;
an out-of-memory run was not attempted.

Correction: stream rows into the atomic temporary output using bounded
memory, and avoid rewriting unchanged files where safe. Keep header and frame
validation and transaction recovery guarantees.

## Performance opportunities

1. **Skip nonlinear evaluation when `g=0`.** `src/nonlinear.cpp:34` performs
   all four padded transforms before multiplying the result by zero. A local
   one-thread Release microbenchmark measured approximately 6, 94, and 640
   microseconds per unnecessary evaluation at N=512, 4096, and 16384.
   ETDRK4 performs four such evaluations per step. Return a zero RHS directly;
   diagnostic density transforms remain necessary when requested.

2. **Reduce transform buffer copies.** `src/fftw_utils.cpp:65` and `:73` copy
   input into an internal buffer and output back to the caller for every
   transform. The four-transform nonlinearity therefore adds eight full
   padded-field copies per evaluation, or 32 per ETDRK4 step. The local 2DGP
   backend demonstrates plans bound to reusable working buffers. That
   approach, with aligned FFTW allocations and fused normalization/filter
   passes, is worth benchmarking here. If using new-array execution instead,
   preserve the FFTW plan's layout and alignment requirements described in the
   [FFTW documentation](https://www.fftw.org/fftw3_doc/New_002darray-Execute-Functions.html).
   No speedup estimate is claimed without an implementation benchmark.

## Scientific documentation and test coverage

The finite-resolution nonlinear Hamiltonian actually evaluated is

\[
 H_{\rm nl,h}=\frac{gL}{2}\sum_{|m|<N/2}
   R(k_m)\,\left|\widehat{P_\rho|\psi|^2}_m\right|^2.
\]

The continuum energy expressions in the README should be accompanied by this
projected form and an explicit distinction between density projection and
projecting the full cubic only once, as in the 2DGP documentation. The two
methods differ at finite resolution, even when both are correctly dealiased.
Also reconcile the README's statement that negative-order hypoviscosity
suppresses the singular zero-mode multiplier with `params.txt` and the code,
which remove the wavefunction mean itself.

Existing convolution tests populate only modes -2 through +2 at N=32; their
density support never reaches the cutoff. The convergence test uses a uniform
solution. These explain why the current suite misses F1. Retain these tests
and add dense spectra, explicit forbidden-mode checks, nonuniform convergence,
restart-statistics checks, malformed HDF5 metadata, and HDF5 plotting coverage.

## Local audit artifacts

The following temporary artifacts were created during this review:

- `/tmp/sh1d_projection_probe.cpp`: exact projected convolution, invariant
  checks, and isolated stage-projection experiment.
- `/tmp/sh1d_solver_order_audit.py`: independent Fourier ODE reference and
  end-to-end solver convergence runs.
- `/tmp/sh1d_audit_cases.py`: forcing, HDF5 plotting, metadata, mean-removal,
  and restart-average fixtures.
- `/tmp/sh1d_nonlinear_benchmark.cpp`: zero-coupling evaluation timing.
- `/tmp/sh1d-audit-cases/nonscalar-attribute.h5`: malformed attribute fixture.
- `/tmp/sh1d-audit-build` and `/tmp/sh1d-audit-asan`: isolated builds.

The Python fixture scripts used fixed temporary directories on their initial
execution; use fresh directories when repeating them so automatic restart
does not alter the intended experiments.

Recommended correction order: F1 and F2 first, then F3–F5, followed by restart
memory use and measured performance improvements.
