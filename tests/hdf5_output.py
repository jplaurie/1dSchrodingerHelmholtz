#!/usr/bin/env python3
"""Exercise HDF5 snapshots, export, FFTW wisdom, and HDF5 input."""

from __future__ import annotations

import argparse
import math
import subprocess
import tempfile
from pathlib import Path


def run(command: list[Path | str]) -> None:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode:
        raise RuntimeError(
            f"command failed: {' '.join(map(str, command))}\n"
            f"{result.stdout}\n{result.stderr}"
        )


def rows(path: Path) -> list[list[float]]:
    return [
        [float(value) for value in line.split()]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--exporter", type=Path, required=True)
    args = parser.parse_args()

    with tempfile.TemporaryDirectory(prefix="sh1d-hdf5-") as temporary:
        root = Path(temporary)
        initial = root / "initial.dat"
        with initial.open("w", encoding="utf-8") as stream:
            for index in range(8):
                x = 2.0 * math.pi * index / 8
                value = complex(0.4 + 0.03 * math.cos(x), 0.2 * math.sin(2 * x))
                stream.write(f"{x:.17g} {value.real:.17g} {value.imag:.17g}\n")

        wisdom = root / "fftw.wisdom"
        parameters = root / "run.params"
        parameters.write_text(
            f"""gridPoints 8
domainLength {2.0 * math.pi:.17g}
timeStep 0.001
numberOfSteps 1
outputIntervalSteps 1
model nonlinearSchrodinger
integrator etd4
dispersionCoefficient -0.5
nonlinearityCoefficient 1
chemicalPotential 0
helmholtzParameter 0
hyperviscosity 0
hypoviscosity 0
forcingEnabled false
randomSeed 1
threadCount 1
fieldOutputFormat both
hdf5CompressionLevel 1
fftwPlanning measure
fftwWisdomFile {wisdom}
initialConditionFile {initial}
dataDirectory {root / 'data'}
outputDirectory {root / 'output'}
""",
            encoding="utf-8",
        )
        run([args.executable, parameters])
        if not wisdom.is_file():
            raise AssertionError("FFTW wisdom file was not written")
        for frame in (0, 1):
            stem = f"wavefunction_{frame:08d}"
            text_field = root / "data" / f"{stem}.dat"
            hdf5_field = root / "data" / f"{stem}.h5"
            exported = root / f"{stem}_exported.dat"
            if not text_field.is_file() or not hdf5_field.is_file():
                raise AssertionError(f"missing text or HDF5 field for frame {frame}")
            run([args.exporter, hdf5_field, exported])
            expected, actual = rows(text_field), rows(exported)
            if len(expected) != 8 or len(actual) != 8:
                raise AssertionError("unexpected exported field length")
            error = max(
                abs(left - right)
                for expected_row, actual_row in zip(expected, actual, strict=True)
                for left, right in zip(expected_row, actual_row, strict=True)
            )
            if error > 1.0e-12:
                raise AssertionError(f"HDF5 export differs from text output: {error}")

        gnuplot = root / "gnuplot.dat"
        run(
            [
                args.exporter,
                root / "data/wavefunction_00000001.h5",
                gnuplot,
                "--format",
                "gnuplot",
            ]
        )
        if len(rows(gnuplot)) != 8 or any(len(row) != 5 for row in rows(gnuplot)):
            raise AssertionError("unexpected gnuplot export format")

        run([args.executable, parameters])
        if not (root / "data/wavefunction_00000002.h5").is_file():
            raise AssertionError("HDF5 output was not preserved across restart")

        imported_root = root / "from-hdf5"
        imported_parameters = root / "from-hdf5.params"
        imported_parameters.write_text(
            f"""gridPoints 8
domainLength {2.0 * math.pi:.17g}
timeStep 0.001
numberOfSteps 1
outputIntervalSteps 1
model nonlinearSchrodinger
integrator etd2
dispersionCoefficient -0.5
nonlinearityCoefficient 0
chemicalPotential 0
helmholtzParameter 0
hyperviscosity 0
hypoviscosity 0
forcingEnabled false
randomSeed 1
threadCount 1
fieldOutputFormat hdf5
initialConditionFile {root / 'data/wavefunction_00000001.h5'}
dataDirectory {imported_root / 'data'}
outputDirectory {imported_root / 'output'}
""",
            encoding="utf-8",
        )
        run([args.executable, imported_parameters])
        if not (imported_root / "data/wavefunction_00000000.h5").is_file():
            raise AssertionError("HDF5 initial condition run did not write HDF5 output")
        if (imported_root / "data/wavefunction_00000000.dat").exists():
            raise AssertionError("HDF5-only mode unexpectedly wrote text output")
        run([args.executable, imported_parameters])
        if not (imported_root / "data/wavefunction_00000002.h5").is_file():
            raise AssertionError("HDF5-only output did not restart cleanly")

        resolved = (root / "output/resolved_parameters.txt").read_text(encoding="utf-8")
        for key in ("solverVersion 0.3.0", "gitCommit ", "gitDirty "):
            if key not in resolved:
                raise AssertionError(f"resolved parameters omit provenance field: {key}")
        print("HDF5 snapshots, export, import, and FFTW wisdom passed")


if __name__ == "__main__":
    main()
