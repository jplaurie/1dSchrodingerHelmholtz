#!/usr/bin/env python3
"""Check temporal order for uniform and independently projected Fourier ODEs."""

from __future__ import annotations

import argparse
import cmath
import math
import subprocess
import tempfile
from pathlib import Path


def run_case(
    executable: Path,
    root: Path,
    method: str,
    steps: int,
    final_time: float,
    amplitude: complex,
) -> float:
    case = root / f"{method}-{steps}"
    case.mkdir()
    initial = case / "initial.dat"
    initial.write_text(
        "".join(f"{amplitude.real:.17g} {amplitude.imag:.17g}\n" for _ in range(8)),
        encoding="utf-8",
    )
    parameters = case / "run.params"
    parameters.write_text(
        f"""gridPoints 8
domainLength {2.0 * math.pi:.17g}
timeStep {final_time / steps:.17g}
numberOfSteps {steps}
outputIntervalSteps {steps}
model nonlinearSchrodinger
integrator {method}
dispersionCoefficient -0.5
nonlinearityCoefficient 1.3
chemicalPotential 0.2
helmholtzParameter 0
hyperviscosity 0
hypoviscosity 0
forcingEnabled false
randomSeed 1
threadCount 1
initialConditionFile {initial}
dataDirectory {case / 'data'}
outputDirectory {case / 'output'}
""",
        encoding="utf-8",
    )
    result = subprocess.run(
        [str(executable), str(parameters)], text=True, capture_output=True, check=False
    )
    if result.returncode:
        raise RuntimeError(
            f"{method}/{steps} failed:\n{result.stdout}\n{result.stderr}"
        )
    values = (case / "data/wavefunction_00000001.dat").read_text(
        encoding="utf-8"
    ).splitlines()
    _, real, imaginary = values[1].split()
    actual = complex(float(real), float(imaginary))
    frequency = 0.2 + 1.3 * abs(amplitude) ** 2
    exact = amplitude * cmath.exp(-1j * frequency * final_time)
    return abs(actual - exact)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path, required=True)
    args = parser.parse_args()
    expected_orders = {"rk2": 2, "etd2": 2, "etd4": 4}
    amplitude = complex(0.7, 0.2)
    steps = (4, 8, 16, 32)
    with tempfile.TemporaryDirectory(prefix="sh1d-convergence-") as temporary:
        root = Path(temporary)
        for method, expected_order in expected_orders.items():
            errors = [
                run_case(args.executable, root, method, count, 0.8, amplitude)
                for count in steps
            ]
            observed = math.log(errors[-2] / errors[-1], 2)
            minimum = expected_order - 0.35
            if not math.isfinite(observed) or observed < minimum:
                raise AssertionError(
                    f"{method} observed order {observed:.3f}, expected at least "
                    f"{minimum:.2f}; errors={errors}"
                )
            print(f"{method}: order={observed:.3f}, errors={errors}")
        for model in ("schrodingerHelmholtz", "longWave", "nonlinearSchrodinger"):
            for remove_mean in (False, True):
                initial, reference = projected_reference(model, remove_mean)
                for method, expected_order in expected_orders.items():
                    errors = [
                        run_nonuniform(args.executable, root, method, model, remove_mean,
                                       count, initial, reference)
                        for count in (8, 16, 32)
                    ]
                    observed = math.log2(errors[-2] / errors[-1])
                    if not math.isfinite(observed) or observed < expected_order - 0.35:
                        raise AssertionError(f"{model}/{method}/remove_mean={remove_mean}: "
                                             f"order={observed:.3f}, errors={errors}")
                    print(f"{model}/{method}/remove_mean={remove_mean}: order={observed:.3f}")


def projected_reference(model: str, remove_mean: bool) -> tuple[list[complex], list[complex]]:
    """Direct non-wrapping convolutions plus a fine, independent RK4 solve."""
    modes = (0, 1, 2, 3, -4, -3, -2, -1)
    indices = {k: i for i, k in enumerate(modes) if k != -4 and not (remove_mean and k == 0)}
    initial = [complex(.22 * math.cos(.71 * i), .17 * math.sin(.39 * i))
               if k in indices else 0j for i, k in enumerate(modes)]
    pairs = {d: [(i, indices[k - d]) for k, i in indices.items() if k - d in indices]
             for d in range(-3, 4)}
    response = {d: (1 / (1 + .35 * d * d) if model == "schrodingerHelmholtz"
                    else 1 - .35 * d * d if model == "longWave" else 1)
                for d in pairs}
    products = {i: [(d, indices[k - d]) for d in pairs if k - d in indices]
                for k, i in indices.items()}
    linear = [-.5j * k * k - (.1 / k ** 4 if remove_mean and k else 0) for k in modes]

    def rhs(w):
        rho = {d: response[d] * sum(w[i] * w[j].conjugate() for i, j in entries)
               for d, entries in pairs.items()}
        return [linear[i] * w[i] - 1.3j * sum(rho[d] * w[j] for d, j in products[i])
                if i in products else 0j for i in range(8)]

    w = initial.copy()
    h = .4 / 4096
    for _ in range(4096):
        k1 = rhs(w)
        k2 = rhs([v + h / 2 * a for v, a in zip(w, k1)])
        k3 = rhs([v + h / 2 * a for v, a in zip(w, k2)])
        k4 = rhs([v + h * a for v, a in zip(w, k3)])
        w = [v + h / 6 * (a + 2 * b + 2 * c + d)
             for v, a, b, c, d in zip(w, k1, k2, k3, k4)]
    return initial, w


def run_nonuniform(executable, root, method, model, remove_mean, steps, initial, reference):
    case = root / f"{model}-{method}-{remove_mean}-{steps}"
    case.mkdir()
    physical = [sum(a * cmath.exp(2j * math.pi * i * j / 8) for i, a in enumerate(initial))
                for j in range(8)]
    (case / "initial.dat").write_text("".join(f"{v.real:.17g} {v.imag:.17g}\n" for v in physical))
    settings = dict(gridPoints=8, domainLength=2 * math.pi, timeStep=.4 / steps,
                    numberOfSteps=steps, outputIntervalSteps=steps, model=model, integrator=method,
                    dispersionCoefficient=-.5, nonlinearityCoefficient=1.3, helmholtzParameter=.35,
                    hypoviscosity=.1 if remove_mean else 0, hypoviscosityOrder=-2,
                    forcingEnabled="false", threadCount=1, initialConditionFile=case / "initial.dat",
                    dataDirectory=case / "data", outputDirectory=case / "output")
    params = case / "run.params"
    params.write_text("".join(f"{k} {v}\n" for k, v in settings.items()))
    completed = subprocess.run([str(executable), str(params)], capture_output=True, text=True)
    if completed.returncode:
        raise RuntimeError(completed.stderr)
    physical = [complex(*map(float, line.split()[1:]))
                for line in (case / "data/wavefunction_00000001.dat").read_text().splitlines()
                if not line.startswith("#")]
    actual = [sum(v * cmath.exp(-2j * math.pi * i * j / 8) for j, v in enumerate(physical)) / 8
              for i in range(8)]
    return max(abs(v - exact) for v, exact in zip(actual, reference))


if __name__ == "__main__":
    main()
