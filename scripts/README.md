# Analysis scripts

The plotting helpers read both whitespace-delimited and HDF5 wavefunction
files, together with the CSV diagnostics. They require Python 3, NumPy, and
Matplotlib; reading HDF5 snapshots additionally requires `h5py`. When a frame
exists in both formats, the text file is preferred.

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
