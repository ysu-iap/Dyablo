import h5py
import matplotlib.pyplot as plt
import argparse
from pathlib import Path
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("filenames", nargs="+")
parser.add_argument("--output-folder", default=".")
parser.add_argument("--grid-size", default=64, type=int)

args = parser.parse_args()

out_folder = Path(args.output_folder)
out_folder.mkdir(parents=True, exist_ok=True)  # create output folder if it doesn't exist

def plot_file(fname: Path):
    with h5py.File(fname, "r") as f:
        x, y, z = f["coordinates"][:].T
        m = f["mass"][:]

    bins = np.linspace(0, 1, args.grid_size + 1)

    fig, axes = plt.subplots(1, 1, figsize=(5, 5))
    norm = plt.matplotlib.colors.Normalize()
    axes.hist2d(x, y, bins=bins, norm=norm, weights=m)

    fig.tight_layout()
    fig.savefig(out_folder / f"{fname.stem}.png")
    plt.close(fig)

for fname in args.filenames:
    plot_file(Path(fname))
