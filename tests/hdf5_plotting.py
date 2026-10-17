"""Load C++ HDF5-only output through the plotting helper."""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from sh1d_plotting import load_wavefunction

parser = argparse.ArgumentParser()
parser.add_argument("--executable", type=Path, required=True)
args = parser.parse_args()
with tempfile.TemporaryDirectory(prefix="sh1d-plotting-") as temporary:
    root = Path(temporary)
    initial = root / "initial.dat"
    initial.write_text("0.7 0.2\n" * 8)
    params = root / "run.params"
    params.write_text(f"""gridPoints 8
numberOfSteps 1
outputIntervalSteps 1
nonlinearityCoefficient 0
dispersionCoefficient 0
fieldOutputFormat hdf5
initialConditionFile {initial}
dataDirectory {root / 'data'}
outputDirectory {root / 'output'}
""")
    subprocess.run([str(args.executable), str(params)], check=True, capture_output=True)
    x, psi = load_wavefunction(root / "data", 1)
    np.testing.assert_allclose(psi, .7 + .2j, atol=1.e-14)
    np.testing.assert_allclose(x, 2 * np.pi * np.arange(8) / 8, atol=1.e-14)
    with h5py.File(root / "data/wavefunction_00000001.h5", "r+") as handle:
        del handle.attrs["format"]
        handle.attrs["format"] = "sh1d_wavefunction_v1"
    np.testing.assert_array_equal(load_wavefunction(root / "data", 1)[1], psi)
print("HDF5 plotting accepts fixed and variable strings")
