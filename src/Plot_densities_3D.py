#!/usr/bin/env python3
"""
Slicer/scatter for DF3 volumes.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt


def load_df3(path):
    """Load a DF3 file and return a float64 cube normalized to [0,1].

    DF3 format: 6-byte header (nx,ny,nz as big-endian 16-bit words), then
    nx*ny*nz uint8 voxels in column-major (Fortran) order.
    """
    with open(path, "rb") as f:
        h = f.read(6)
        if len(h) != 6:
            raise ValueError("DF3 header too short")
        nx = (h[0] << 8) | h[1]
        ny = (h[2] << 8) | h[3]
        nz = (h[4] << 8) | h[5]
        raw = f.read()
    arr = np.frombuffer(raw, dtype=np.uint8)
    if arr.size != nx * ny * nz:
        raise ValueError(f"Expected {nx*ny*nz} voxels, got {arr.size}")
    cube = arr.reshape((nx, ny, nz), order="F").astype(np.float64) / 255.0
    return cube


def parse_args(argv):
    out = {"file": None}
    if len(argv) == 0:
        out["file"] = "density_version5.df3"
        return out
    out["file"] = argv[0]
    return out


def main():
    parsed = parse_args(sys.argv[1:])
    fname = parsed["file"]
    cube = load_df3(fname)

    # Choose center slices for XY, XZ, YZ views
    nx, ny, nz = cube.shape
    ix = nx // 2
    iy = ny // 2
    iz = nz // 2

    # Prepare figure with three side-by-side images
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    axes[0].imshow(cube[:, :, iz].T, origin="lower", cmap="magma")
    axes[0].set_title(f"XY (Z index {iz})")
    axes[0].set_xlabel("X")
    axes[0].set_ylabel("Y")

    axes[1].imshow(cube[:, iy, :].T, origin="lower", cmap="magma")
    axes[1].set_title(f"XZ (Y index {iy})")
    axes[1].set_xlabel("X")
    axes[1].set_ylabel("Z")

    axes[2].imshow(cube[ix, :, :].T, origin="lower", cmap="magma")
    axes[2].set_title(f"YZ (X index {ix})")
    axes[2].set_xlabel("Y")
    axes[2].set_ylabel("Z")

    # Auto-save PNG next to the DF3 file (same base name)
    base = fname.rsplit(".", 1)[0]
    out_png = base + ".png"
    fig.savefig(out_png, dpi=200)
    print("Saved:", out_png)


if __name__ == "__main__":
    main()
