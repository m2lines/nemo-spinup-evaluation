"""Prune restart.nc and mesh_mask.nc into the A/B fixture dirs.

Time-evolving grid files (grid_T_2D, grid_T_3D, grid_U_3D, grid_V_3D) are
handled by CDO in build_fixtures.sh — this script only handles the two files
that have no meaningful time axis to slice/shift.

restart.nc gets a deterministic Gaussian-noise variant in the B dir so that
diff-mode comparisons have something to do.

Usage:
    python tests/utils/create_test_data.py SRC_DIR DST_A DST_B

SRC_DIR must contain restart.nc and mesh_mask.nc.
"""

import argparse
from pathlib import Path

import numpy as np
import xarray as xr

restart_vars_to_keep = [
    "nav_lon",
    "nav_lat",
    "nav_lev",
    "time_counter",
    "un",
    "vn",
    "tn",
    "sn",
    "rhop",
    "sshn",
]

mesh_mask_vars_to_keep = [
    "e1t",
    "e2t",
    "e3t_0",
    "tmask",
    "umask",
    "e3u_0",
    "e2u",
    "e3v_0",
    "e3w_0",
    "e1v",
    "vmask",
    "nav_lat",
    "nav_lon",
    "nav_lev",
]


def add_noise(ds, exclude_vars=None, noise_level=1e-3, seed=42):
    """Add small deterministic Gaussian noise to numeric data variables."""
    rng = np.random.default_rng(seed)
    exclude_vars = set(exclude_vars or [])
    ds_out = ds.copy()
    for v in ds.data_vars:
        if v in exclude_vars:
            continue
        arr = ds[v].values
        if not np.issubdtype(arr.dtype, np.number):
            continue
        scale = np.nanstd(arr)
        if not np.isfinite(scale) or scale == 0:
            scale = 1.0
        noise = rng.normal(0.0, noise_level * scale, size=arr.shape)
        ds_out[v] = xr.DataArray(
            arr + noise, dims=ds[v].dims, attrs=ds[v].attrs, name=v
        )
    return ds_out


def _subset(src, vars_to_keep):
    ds = xr.open_dataset(src)
    present = [v for v in vars_to_keep if v in ds.variables]
    return ds[present]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("src_dir", type=Path)
    p.add_argument("dst_a", type=Path)
    p.add_argument("dst_b", type=Path)
    args = p.parse_args()

    args.dst_a.mkdir(parents=True, exist_ok=True)
    args.dst_b.mkdir(parents=True, exist_ok=True)

    # mesh_mask: identical in both dirs
    mesh = _subset(args.src_dir / "mesh_mask.nc", mesh_mask_vars_to_keep)
    mesh.to_netcdf(args.dst_a / "mesh_mask.nc")
    mesh.to_netcdf(args.dst_b / "mesh_mask.nc")
    print(f"wrote mesh_mask.nc to {args.dst_a} and {args.dst_b}")

    # restart: clean -> A, noisy -> B
    restart = _subset(args.src_dir / "restart.nc", restart_vars_to_keep)
    restart.to_netcdf(args.dst_a / "restart.nc")
    add_noise(
        restart,
        exclude_vars={"nav_lon", "nav_lat", "nav_lev", "time_counter"},
    ).to_netcdf(args.dst_b / "restart.nc")
    print(f"wrote restart.nc (clean) to {args.dst_a} and (noisy) to {args.dst_b}")


if __name__ == "__main__":
    main()
