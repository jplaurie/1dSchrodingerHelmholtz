"""Readers and plotting helpers for 1D Schrödinger-Helmholtz output."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

FRAME_PATTERN = re.compile(r"wavefunction_(\d{8})\.dat$")


def available_frames(data_directory: str | Path) -> list[int]:
    """Return sorted frame numbers found in a solver data directory."""
    frames: list[int] = []
    for path in Path(data_directory).glob("wavefunction_*.dat"):
        match = FRAME_PATTERN.match(path.name)
        if match:
            frames.append(int(match.group(1)))
    return sorted(frames)


def load_wavefunction(
    data_directory: str | Path, frame: int
) -> tuple[np.ndarray, np.ndarray]:
    """Load x and complex psi arrays for one frame."""
    path = Path(data_directory) / f"wavefunction_{frame:08d}.dat"
    values = np.loadtxt(path)
    if values.ndim != 2 or values.shape[1] != 3:
        raise ValueError(f"unexpected wavefunction format: {path}")
    return values[:, 0], values[:, 1] + 1j * values[:, 2]


def load_csv(path: str | Path) -> dict[str, np.ndarray]:
    """Load a solver CSV as a dictionary of float arrays."""
    with Path(path).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return {}
    return {
        key: np.asarray([float(row[key]) for row in rows], dtype=float)
        for key in rows[0]
    }


def frame_rows(table: dict[str, np.ndarray], frame: int) -> dict[str, np.ndarray]:
    """Select the rows belonging to one output frame."""
    if not table:
        raise ValueError("the table is empty")
    mask = table["frame"].astype(int) == frame
    if not np.any(mask):
        raise ValueError(f"frame {frame} is absent from the table")
    return {key: values[mask] for key, values in table.items()}


def plot_wavefunction(
    data_directory: str | Path, frame: int, axes=None
):
    """Plot real(psi), imaginary(psi), and intensity for one frame."""
    x, psi = load_wavefunction(data_directory, frame)
    if axes is None:
        _, axes = plt.subplots(2, 1, sharex=True, figsize=(8, 6))
    axes[0].plot(x, psi.real, label=r"$\mathrm{Re}\,\psi$")
    axes[0].plot(x, psi.imag, label=r"$\mathrm{Im}\,\psi$")
    axes[0].set_ylabel("wavefunction")
    axes[0].legend()
    axes[1].plot(x, np.abs(psi) ** 2, color="black")
    axes[1].set_xlabel(r"$x$")
    axes[1].set_ylabel(r"$|\psi|^2$")
    return axes


def plot_spectrum(output_directory: str | Path, frame: int, axis=None):
    """Plot the instantaneous and running-average wave-action spectrum."""
    table = frame_rows(load_csv(Path(output_directory) / "spectra.csv"), frame)
    if axis is None:
        _, axis = plt.subplots(figsize=(7, 4))
    positive = table["wavenumber"] > 0.0
    axis.loglog(
        table["wavenumber"][positive],
        table["wave_action_spectrum"][positive],
        label="instantaneous",
    )
    axis.loglog(
        table["wavenumber"][positive],
        table["wave_action_spectrum_average"][positive],
        label="running average",
    )
    axis.set_xlabel(r"$|k|$")
    axis.set_ylabel(r"$n_k$")
    axis.legend()
    return axis


def plot_diagnostics(output_directory: str | Path, axes=None):
    """Plot total energy and wave action versus time."""
    table = load_csv(Path(output_directory) / "diagnostics.csv")
    if not table:
        raise ValueError("diagnostics.csv contains no frames")
    if axes is None:
        _, axes = plt.subplots(2, 1, sharex=True, figsize=(7, 6))
    axes[0].plot(table["time"], table["total_energy"])
    axes[0].set_ylabel("total energy")
    axes[1].plot(table["time"], table["wave_action"])
    axes[1].set_ylabel("wave action")
    axes[1].set_xlabel("time")
    return axes

