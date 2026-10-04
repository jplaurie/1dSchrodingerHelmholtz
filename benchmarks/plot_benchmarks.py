#!/usr/bin/env python3
"""Plot complete-timestep timing and serial-normalized thread scaling."""

import argparse
import csv
import json
from pathlib import Path
from statistics import median

import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator


COLORS = {"serial": "#577590", 1: "#277da1", 2: "#4d908e", 4: "#90be6d",
          8: "#f8961e", 12: "#d1495b"}
MARKERS = {"serial": "D", 1: "o", 2: "s", 4: "^", 8: "v", 12: "P"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("benchmarks/results.csv"))
    parser.add_argument("--metadata", type=Path, default=Path("benchmarks/system.json"))
    parser.add_argument("--output-prefix", type=Path,
                        default=Path("benchmarks/thread_scaling"))
    args = parser.parse_args()
    with args.input.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        parser.error("the input file has no benchmark rows")
    metadata = json.loads(args.metadata.read_text())

    samples: dict[object, dict[int, list[float]]] = {}
    for row in rows:
        key: object = "serial" if row["backend"] == "CPU serial" else int(row["threads"])
        samples.setdefault(key, {}).setdefault(int(row["resolution"]), []).append(
            float(row["seconds_per_step"]))
    order: list[object] = ["serial", *metadata["thread_counts"]]
    labels = {"serial": "CPU serial (no OpenMP)"}
    labels.update({threads: f"OpenMP/FFTW ({threads} thread{'s' if threads != 1 else ''})"
                   for threads in metadata["thread_counts"]})
    figure, (timing, speedup) = plt.subplots(1, 2, figsize=(10.8, 4.3))
    medians: dict[object, dict[int, float]] = {}
    for key in order:
        if key not in samples:
            continue
        resolutions = sorted(samples[key])
        values = [median(samples[key][n]) for n in resolutions]
        medians[key] = dict(zip(resolutions, values))
        lower = [value - min(samples[key][n]) for n, value in zip(resolutions, values)]
        upper = [max(samples[key][n]) - value for n, value in zip(resolutions, values)]
        timing.errorbar(resolutions, values, yerr=[lower, upper], label=labels[key],
                        color=COLORS[key], marker=MARKERS[key], linewidth=2, capsize=3)
    timing.set(xscale="log", yscale="log", xlabel="Spectral resolution, N",
               ylabel="Time per ETD4 step (s)", title="Complete timestep time")
    timing.grid(True, which="both", alpha=.25)
    timing.legend(frameon=False, fontsize=8.5)

    baseline = medians["serial"]
    for key in order:
        if key not in medians:
            continue
        resolutions = sorted(set(baseline) & set(medians[key]))
        speedup.plot(resolutions, [baseline[n] / medians[key][n] for n in resolutions],
                     label=labels[key], color=COLORS[key], marker=MARKERS[key], linewidth=2)
    speedup.axhline(1, color=".45", linewidth=1, linestyle="--")
    speedup.set(xscale="log", yscale="log", xlabel="Spectral resolution, N",
                ylabel="Speedup over serial CPU", title="Thread-count speedup")
    speedup.grid(True, which="both", alpha=.25)
    resolutions = sorted({n for values in samples.values() for n in values})
    for axis in (timing, speedup):
        axis.set_xscale("log", base=2)
        axis.set_xticks(resolutions, labels=[f"{n:,}" for n in resolutions])
        axis.xaxis.set_minor_locator(NullLocator())
    figure.suptitle(f"1D Schrodinger--Helmholtz ETD4 benchmark\n{metadata['cpu_model']}",
                    fontsize=11)
    figure.tight_layout()
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    for extension in ("svg", "png"):
        path = args.output_prefix.with_suffix(f".{extension}")
        figure.savefig(path, dpi=180, bbox_inches="tight")
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
