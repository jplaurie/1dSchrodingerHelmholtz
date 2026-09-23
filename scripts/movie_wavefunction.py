#!/usr/bin/env python3
"""Render wavefunction and intensity frames as an MP4 or GIF."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

from sh1d_plotting import available_frames, load_wavefunction


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data_directory", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--fps", type=int, default=20)
    parser.add_argument("--stride", type=int, default=1)
    args = parser.parse_args()
    if args.stride < 1:
        parser.error("--stride must be positive")

    frames = available_frames(args.data_directory)[:: args.stride]
    if not frames:
        parser.error("no wavefunction frames were found")
    samples = [load_wavefunction(args.data_directory, frame) for frame in frames]
    x = samples[0][0]
    amplitude = max(float(np.max(np.abs(psi))) for _, psi in samples)
    intensity = max(float(np.max(np.abs(psi) ** 2)) for _, psi in samples)

    figure, axes = plt.subplots(2, 1, sharex=True, figsize=(9, 6))
    (real_line,) = axes[0].plot([], [], label=r"$\mathrm{Re}\,\psi$")
    (imaginary_line,) = axes[0].plot([], [], label=r"$\mathrm{Im}\,\psi$")
    (intensity_line,) = axes[1].plot([], [], color="black")
    axes[0].legend(loc="upper right")
    axes[0].set_ylabel("wavefunction")
    axes[1].set_ylabel(r"$|\psi|^2$")
    axes[1].set_xlabel(r"$x$")
    axes[0].set_xlim(float(x[0]), float(x[-1]))
    axes[0].set_ylim(-1.05 * amplitude, 1.05 * amplitude)
    axes[1].set_ylim(0.0, 1.05 * intensity)
    title = axes[0].set_title("")

    def update(position: int):
        _, psi = samples[position]
        real_line.set_data(x, psi.real)
        imaginary_line.set_data(x, psi.imag)
        intensity_line.set_data(x, np.abs(psi) ** 2)
        title.set_text(f"frame {frames[position]}")
        return real_line, imaginary_line, intensity_line, title

    movie = animation.FuncAnimation(
        figure, update, frames=len(frames), interval=1000 / args.fps, blit=False
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    writer = "pillow" if args.output.suffix.lower() == ".gif" else "ffmpeg"
    movie.save(args.output, writer=writer, fps=args.fps)


if __name__ == "__main__":
    main()

