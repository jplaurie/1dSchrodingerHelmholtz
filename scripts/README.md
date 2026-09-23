# Analysis scripts

The plotting helpers read the modern solver's whitespace-delimited
wavefunction files and CSV diagnostics. They require Python 3, NumPy, and
Matplotlib.

From the repository root, start Jupyter and open one of the notebooks:

```bash
jupyter lab scripts/
```

The notebooks expect the quick-start output in `data/quickstart/` and
`output/quickstart/`. Change `DATA_DIRECTORY`, `OUTPUT_DIRECTORY`, or `FRAME`
in the first configuration cell for another run.

To render every stored wavefunction frame:

```bash
python scripts/movie_wavefunction.py data/quickstart output/quickstart/wavefunction.mp4
python scripts/movie_wavefunction.py data/quickstart output/quickstart/wavefunction.gif --fps 12
```

MP4 output requires an FFmpeg installation; GIF output uses Pillow. Use
`--stride N` to render every Nth stored frame.

