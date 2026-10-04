# CPU thread-scaling benchmark

This benchmark times complete ETD4-B timesteps for the Schrödinger--Helmholtz model. Each step
performs four nonlinear evaluations, and each evaluation contains four complex FFTs on the
3/2-padded grid.

Construction, FFT planning, allocation, coefficient creation, state initialization, and warm-up
happen before the timed interval. Diagnostics and file output are excluded. The serial reference
is compiled without OpenMP or threaded FFTW; the other configurations use the production
OpenMP/threaded-FFTW build with the requested thread count.

```bash
cmake -S . -B build/benchmarks -DCMAKE_BUILD_TYPE=Release
cmake --build build/benchmarks --parallel \
  --target sh1d_benchmark_cpu_serial sh1d_benchmark_cpu
python3 benchmarks/run_benchmarks.py --overwrite
python3 benchmarks/plot_benchmarks.py
```

Raw trials are in `results.csv`, with machine/build metadata and exact timing exclusions in
`system.json`. The runner uses the Python standard library; plot generation also requires
Matplotlib.
