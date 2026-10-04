#!/usr/bin/env python3
"""Run calibrated serial and threaded 1D Schrodinger--Helmholtz benchmarks."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import platform
import re
import subprocess
import sys


def output(command: list[str]) -> str:
    try:
        return subprocess.check_output(command, text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def physical_cores() -> int:
    cores: set[tuple[str, str]] = set()
    for cpu in Path("/sys/devices/system/cpu").glob("cpu[0-9]*"):
        try:
            cores.add(((cpu / "topology/physical_package_id").read_text().strip(),
                       (cpu / "topology/core_id").read_text().strip()))
        except OSError:
            pass
    return len(cores) or (os.cpu_count() or 1)


def parse_result(stdout: str) -> dict[str, str]:
    line = next((item for item in reversed(stdout.splitlines())
                 if item.startswith("BENCHMARK,")), "")
    fields = next(csv.reader([line]), [])
    if len(fields) != 7:
        raise RuntimeError(f"could not parse benchmark output:\n{stdout}")
    return dict(zip(("record", "backend", "resolution", "threads", "steps", "seconds",
                     "seconds_per_step"), fields))


def run_once(executable: Path, resolution: int, steps: int, warmup: int,
             threads: int) -> dict[str, str]:
    environment = os.environ.copy()
    environment.update({"OMP_NUM_THREADS": str(threads), "OMP_PROC_BIND": "close",
                        "OMP_PLACES": "cores"})
    completed = subprocess.run(
        [str(executable), str(resolution), str(steps), str(warmup), str(threads)],
        text=True, capture_output=True, env=environment
    )
    if completed.returncode:
        raise RuntimeError(f"benchmark failed ({executable.name}, {threads} threads):\n"
                           f"{completed.stdout}{completed.stderr}")
    return parse_result(completed.stdout)


def calibrate(executable: Path, resolution: int, warmup: int, threads: int,
              target_seconds: float, minimum_steps: int) -> int:
    steps = 1
    while True:
        result = run_once(executable, resolution, steps, warmup, threads)
        seconds = float(result["seconds"])
        if seconds >= .1 or steps >= 1_000_000:
            break
        steps *= 4
    return max(minimum_steps, min(1_000_000, math.ceil(target_seconds * steps / seconds)))


def metadata(repo: Path, args: argparse.Namespace) -> dict[str, object]:
    cpu_model = "unknown"
    try:
        match = re.search(r"^model name\s*:\s*(.+)$", Path("/proc/cpuinfo").read_text(),
                          re.MULTILINE)
        if match:
            cpu_model = match.group(1).strip()
    except OSError:
        pass
    return {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": platform.node(), "platform": platform.platform(),
        "cpu_model": cpu_model, "logical_cpus": os.cpu_count(),
        "physical_cores_detected": physical_cores(), "thread_counts": args.threads,
        "compiler": output(["c++", "--version"]).splitlines()[0],
        "fftw": output(["pkg-config", "--modversion", "fftw3"]),
        "git_commit": output(["git", "-C", str(repo), "rev-parse", "HEAD"]),
        "git_dirty": output(["git", "-C", str(repo), "status", "--porcelain"])
                     not in ("", "unknown"),
        "resolutions": args.resolutions, "trials": args.trials,
        "target_timed_seconds_per_trial": args.target_seconds,
        "minimum_timed_steps": args.minimum_timed_steps,
        "warmup_steps": args.warmup_steps,
        "model": "Schrodinger--Helmholtz", "integrator": "ETD4-B",
        "dealiasing": "two-pass 3/2 padding",
        "timing_excludes": ["process startup", "FFTW thread initialization", "FFT planning",
                            "allocation", "coefficient construction", "state initialization",
                            "warm-up steps", "diagnostics", "file output"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=Path("build/benchmarks"))
    parser.add_argument("--output", type=Path, default=Path("benchmarks/results.csv"))
    parser.add_argument("--metadata", type=Path, default=Path("benchmarks/system.json"))
    parser.add_argument("--resolutions", nargs="+", type=int,
                        default=[512, 1024, 2048, 4096, 8192, 16384])
    parser.add_argument("--threads", nargs="+", type=int, default=[1, 2, 4, 8, 12])
    parser.add_argument("--trials", type=int, default=3)
    parser.add_argument("--target-seconds", type=float, default=.75)
    parser.add_argument("--warmup-steps", type=int, default=3)
    parser.add_argument("--minimum-timed-steps", type=int, default=3)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if (any(value <= 0 for value in [*args.resolutions, *args.threads, args.trials,
                                     args.warmup_steps, args.minimum_timed_steps]) or
            args.target_seconds <= 0):
        parser.error("resolutions, threads, counts, and durations must be positive")
    if len(set(args.threads)) != len(args.threads):
        parser.error("thread counts must not contain duplicates")

    repo = Path(__file__).resolve().parents[1]
    build = (repo / args.build_dir).resolve() if not args.build_dir.is_absolute() else args.build_dir
    serial = build / "sh1d_benchmark_cpu_serial"
    threaded = build / "sh1d_benchmark_cpu"
    for executable in (serial, threaded):
        if not executable.is_file():
            parser.error(f"missing {executable}; build the benchmark targets first")
    result_path = (repo / args.output).resolve() if not args.output.is_absolute() else args.output
    metadata_path = ((repo / args.metadata).resolve()
                     if not args.metadata.is_absolute() else args.metadata)
    if result_path.exists() and not args.overwrite:
        parser.error(f"{result_path} exists; pass --overwrite to replace it")
    result_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(metadata(repo, args), indent=2) + "\n")

    columns = ["backend", "resolution", "trial", "threads", "steps", "warmup_steps",
               "seconds", "seconds_per_step", "grid_points_per_second"]
    configurations = [("CPU serial", serial, 1)] + [
        ("CPU/OpenMP", threaded, threads) for threads in args.threads
    ]
    with result_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for resolution in args.resolutions:
            for backend, executable, threads in configurations:
                steps = calibrate(executable, resolution, args.warmup_steps, threads,
                                  args.target_seconds, args.minimum_timed_steps)
                suffix = "thread" if threads == 1 else "threads"
                label = backend if backend == "CPU serial" else f"{backend} ({threads} {suffix})"
                print(f"{label:>24} N={resolution:<5} timed_steps={steps} "
                      f"warmup_steps={args.warmup_steps}", flush=True)
                for trial in range(1, args.trials + 1):
                    result = run_once(executable, resolution, steps, args.warmup_steps, threads)
                    step_time = float(result["seconds_per_step"])
                    writer.writerow({"backend": backend, "resolution": resolution,
                                     "trial": trial, "threads": threads,
                                     "steps": result["steps"],
                                     "warmup_steps": args.warmup_steps,
                                     "seconds": result["seconds"],
                                     "seconds_per_step": step_time,
                                     "grid_points_per_second": resolution / step_time})
                    stream.flush()
                    print(f"  trial {trial}: {step_time:.6g} s/step", flush=True)
    print(f"wrote {result_path}\nwrote {metadata_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
