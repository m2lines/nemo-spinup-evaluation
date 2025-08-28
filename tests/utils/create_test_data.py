from pathlib import Path

import xarray as xr

# Source and destination directories
src_dir = Path("tests/data/DINO")
dst_dir = Path("tests/data/DINO_subsampled_3")
dst_dir.mkdir(exist_ok=True)

files = [
    "grid_T_3D.nc",
    "grid_T_2D.nc",
    "grid_U_3D.nc",
    "grid_V_3D.nc",
    "mesh_mask.nc",
    "restart.nc",
]

# Variables to keep in restart.nc
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


for fname in files:
    src_file = src_dir / fname
    dst_file = dst_dir / fname
    if not src_file.exists():
        print(f"File {src_file} not found, skipping.")
        continue

    ds = xr.open_dataset(src_file)
    if fname == "restart.nc":
        # Only keep selected variables
        vars_present = [v for v in restart_vars_to_keep if v in ds.variables]
        ds = ds[vars_present]
        # No time subsampling needed if only one time point
    elif fname == "mesh_mask.nc":
        vars_present = [v for v in mesh_mask_vars_to_keep if v in ds.variables]
        ds = ds[vars_present]
    # Subsample first 2 time steps if time or time_counter exists
    elif "time" in ds.dims:
        ds = ds.isel(time=slice(0, 2))
    elif "time_counter" in ds.dims:
        ds = ds.isel(time_counter=slice(0, 2))
        # else: no time dimension, copy as is

    ds.to_netcdf(dst_file)
    print(f"Processed {fname} written to {dst_file}")

print("Subsampling and variable selection complete.")
