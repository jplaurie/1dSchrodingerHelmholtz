#!/usr/bin/env python3
"""Check temporal order against the exact spatially uniform solution."""

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


if __name__ == "__main__":
    main()
