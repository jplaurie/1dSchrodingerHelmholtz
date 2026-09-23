#!/usr/bin/env python3
"""End-to-end output and restart regression test for the CPU solver."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import tempfile
from pathlib import Path


def write_initial(path: Path, count: int) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for index in range(count):
            x = 2.0 * math.pi * index / count
            phase_1 = complex(math.cos(x), math.sin(x))
            phase_2 = complex(math.cos(-2.0 * x), math.sin(-2.0 * x))
            phase_3 = complex(math.cos(3.0 * x), math.sin(3.0 * x))
            value = 0.08 * (phase_1 + 0.7 * phase_2 + 0.45 * phase_3)
            handle.write(f"{x:.17g} {value.real:.17g} {value.imag:.17g}\n")


def write_parameters(
    path: Path,
    initial: Path,
    data: Path,
    output: Path,
    steps: int,
    integrator: str = "etd4",
    model: str = "schrodingerHelmholtz",
    forcing_enabled: bool = True,
    forcing_profile: str = "annulus",
    helmholtz_parameter: float = 1.0,
    output_interval: int = 10,
) -> None:
    path.write_text(
        f"""gridPoints 32
domainLength {2.0 * math.pi:.17g}
timeStep 0.0001
numberOfSteps {steps}
outputIntervalSteps {output_interval}
model {model}
integrator {integrator}
dispersionCoefficient -0.5
nonlinearityCoefficient -1
chemicalPotential 0
helmholtzParameter {helmholtz_parameter:.17g}
hyperviscosity 0
hyperviscosityOrder 2
hypoviscosity 0
hypoviscosityOrder -2
forcingEnabled {str(forcing_enabled).lower()}
forcingProfile {forcing_profile}
forcingWavenumber 4
forcingWidth 1
forcingAmplitude 0.01
forcingShapeOrder 4
forcingLogWidth 0.2
targetWaveActionInjectionRate 0
randomSeed 1234
writeModeDiagnostics true
threadCount 1
overwriteOutput false
initialConditionFile {initial}
dataDirectory {data}
outputDirectory {output}
""",
        encoding="utf-8",
    )


def run(executable: Path, parameters: Path) -> None:
    subprocess.run(
        [str(executable), str(parameters)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def require_run_failure(executable: Path, parameters: Path) -> None:
    result = subprocess.run(
        [str(executable), str(parameters)],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if result.returncode == 0:
        raise AssertionError(f"invalid run unexpectedly succeeded: {parameters}")


def read_wavefunction(path: Path) -> list[complex]:
    values: list[complex] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            _, real, imaginary = line.split()
            values.append(complex(float(real), float(imaginary)))
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path, required=True)
    args = parser.parse_args()

    conservation_results: list[tuple[str, str, float]] = []
    with tempfile.TemporaryDirectory(prefix="sh1d-regression-") as temporary:
        root = Path(temporary)
        initial = root / "initial.dat"
        write_initial(initial, 32)

        malformed_initial = root / "malformed-initial.dat"
        malformed_initial.write_text("0 1 0 trailing-text\n", encoding="utf-8")
        malformed_parameters = root / "malformed-initial.params"
        write_parameters(
            malformed_parameters,
            malformed_initial,
            root / "malformed-data",
            root / "malformed-output",
            1,
            forcing_enabled=False,
            output_interval=1,
        )
        require_run_failure(args.executable, malformed_parameters)

        split_parameters = root / "split.params"
        split_data, split_output = root / "split-data", root / "split-output"
        write_parameters(split_parameters, initial, split_data, split_output, 20)
        run(args.executable, split_parameters)
        run(args.executable, split_parameters)

        full_parameters = root / "full.params"
        full_data, full_output = root / "full-data", root / "full-output"
        write_parameters(full_parameters, initial, full_data, full_output, 40)
        run(args.executable, full_parameters)

        split = read_wavefunction(split_data / "wavefunction_00000004.dat")
        full = read_wavefunction(full_data / "wavefunction_00000004.dat")
        if len(split) != 32 or len(full) != 32:
            raise AssertionError("unexpected wavefunction length")
        error = max(abs(a - b) for a, b in zip(split, full, strict=True))
        if error > 2.0e-12:
            raise AssertionError(f"split restart differs from full run: {error}")

        with (split_output / "diagnostics.csv").open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        if [int(row["frame"]) for row in rows] != [1, 2, 3, 4]:
            raise AssertionError("diagnostics did not append cleanly across restart")
        for name in (
            "spectra.csv",
            "fluxes.csv",
            "modes.csv",
            "forcing_summary.csv",
            "forcing_spectrum.csv",
            "resolved_parameters.txt",
        ):
            if not (split_output / name).is_file():
                raise AssertionError(f"missing output file: {name}")
        if len(list((split_output / "segments").iterdir())) != 2:
            raise AssertionError("restart did not create a second segment record")

        corrupt_parameters = root / "corrupt-checkpoint.params"
        corrupt_data = root / "corrupt-checkpoint-data"
        corrupt_output = root / "corrupt-checkpoint-output"
        write_parameters(
            corrupt_parameters,
            initial,
            corrupt_data,
            corrupt_output,
            10,
        )
        run(args.executable, corrupt_parameters)
        with (corrupt_data / "checkpoint_00000001.bin").open("ab") as handle:
            handle.write(b"unexpected trailing data")
        require_run_failure(args.executable, corrupt_parameters)

        header_parameters = root / "bad-header.params"
        header_data = root / "bad-header-data"
        header_output = root / "bad-header-output"
        write_parameters(
            header_parameters,
            initial,
            header_data,
            header_output,
            10,
        )
        run(args.executable, header_parameters)
        diagnostics_path = header_output / "diagnostics.csv"
        diagnostics_text = diagnostics_path.read_text(encoding="utf-8")
        diagnostics_path.write_text(
            diagnostics_text.replace("time,frame,", "bad,header,", 1),
            encoding="utf-8",
        )
        require_run_failure(args.executable, header_parameters)

        for method, model in (
            ("etd2", "nonlinearSchrodinger"),
            ("rk2", "longWave"),
        ):
            parameters = root / f"{method}.params"
            data, output = root / f"{method}-data", root / f"{method}-output"
            write_parameters(
                parameters,
                initial,
                data,
                output,
                10,
                integrator=method,
                model=model,
            )
            run(args.executable, parameters)
            if not (data / "wavefunction_00000001.dat").is_file():
                raise AssertionError(f"{method} did not write its final frame")

        lognormal_parameters = root / "lognormal.params"
        lognormal_data = root / "lognormal-data"
        lognormal_output = root / "lognormal-output"
        write_parameters(
            lognormal_parameters,
            initial,
            lognormal_data,
            lognormal_output,
            10,
            forcing_profile="logNormal",
        )
        run(args.executable, lognormal_parameters)
        with (lognormal_output / "forcing_spectrum.csv").open(
            encoding="utf-8"
        ) as handle:
            forcing_rows = list(csv.DictReader(handle))
        forcing_by_mode = {
            int(row["mode"]): float(row["amplitude"]) for row in forcing_rows
        }
        expected_mode_2 = 0.01 * math.exp(
            -0.5 * (math.log(2.0 / 4.0) / 0.2) ** 2
        )
        if not math.isclose(forcing_by_mode[4], 0.01, rel_tol=1.0e-13):
            raise AssertionError("log-normal forcing does not peak at k_f")
        if not math.isclose(
            forcing_by_mode[2], expected_mode_2, rel_tol=1.0e-13
        ):
            raise AssertionError("log-normal forcing envelope is incorrect")
        resolved = (lognormal_output / "resolved_parameters.txt").read_text(
            encoding="utf-8"
        )
        resolved_parameters = dict(
            line.split(maxsplit=1) for line in resolved.splitlines() if line
        )
        profile_recorded = resolved_parameters.get("forcingProfile") == "logNormal"
        width_recorded = math.isclose(
            float(resolved_parameters.get("forcingLogWidth", "nan")), 0.2
        )
        if not profile_recorded or not width_recorded:
            raise AssertionError("resolved log-normal parameters were not recorded")

        for profile in ("gaussian", "exponential", "singleMode"):
            parameters = root / f"forcing-{profile}.params"
            data = root / f"forcing-{profile}-data"
            output = root / f"forcing-{profile}-output"
            write_parameters(
                parameters,
                initial,
                data,
                output,
                10,
                forcing_profile=profile,
            )
            run(args.executable, parameters)
            if not (data / "wavefunction_00000001.dat").is_file():
                raise AssertionError(f"{profile} forcing did not finish")

        for model in (
            "schrodingerHelmholtz",
            "longWave",
            "nonlinearSchrodinger",
        ):
            parameters = root / f"conservation-{model}.params"
            data = root / f"conservation-{model}-data"
            output = root / f"conservation-{model}-output"
            write_parameters(
                parameters,
                initial,
                data,
                output,
                2000,
                model=model,
                forcing_enabled=False,
                helmholtz_parameter=0.02,
                output_interval=100,
            )
            run(args.executable, parameters)
            with (output / "diagnostics.csv").open(encoding="utf-8") as handle:
                diagnostics = list(csv.DictReader(handle))
            if len(diagnostics) != 20:
                raise AssertionError(f"{model} conservation run lost output frames")
            with (output / "fluxes.csv").open(encoding="utf-8") as handle:
                fluxes = list(csv.DictReader(handle))
            terminal_fluxes = fluxes[15::16]
            if len(terminal_fluxes) != len(diagnostics):
                raise AssertionError(f"{model} flux output has incomplete frames")
            for quantity in ("wave_action_flux", "energy_flux"):
                closure_error = max(
                    abs(float(row[quantity])) for row in terminal_fluxes
                )
                if closure_error > 1.0e-10:
                    raise AssertionError(
                        f"{model} {quantity} does not close: {closure_error:.3e}"
                    )
            for quantity in ("wave_action", "total_energy"):
                values = [float(row[quantity]) for row in diagnostics]
                scale = max(abs(values[0]), 1.0e-14)
                relative_drift = max(abs(value - values[0]) for value in values) / scale
                conservation_results.append((model, quantity, relative_drift))
                if relative_drift > 2.0e-8:
                    raise AssertionError(
                        f"{model} {quantity} relative drift is {relative_drift:.3e}"
                    )

    for model, quantity, drift in conservation_results:
        print(f"{model} {quantity} maximum relative drift: {drift:.3e}")
    print("regression tests passed")


if __name__ == "__main__":
    main()
