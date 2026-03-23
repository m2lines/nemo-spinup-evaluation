"""Restructure Samudra output to NEMO-like 4D format."""

import re

import xarray as xr

# --- Input and output paths ---
input_path = "/home/sg2147/Samudra/rollout/2025-10-28-rollout-om4_samudra/"
"predictions.nc"
output_path = "/home/sg2147/Samudra/rollout/2025-10-28-rollout-om4_samudra/"
"predictions_nemo_format.nc"

# --- Open dataset ---
ds = xr.open_dataset(input_path)


def combine_levels(ds, base_name, depth_name="deptht"):
    """Combine split variables (e.g., so_0 ... so_18) into one 4D variable."""
    vars_found = sorted(
        [v for v in ds.data_vars if re.match(rf"{base_name}_\d+", v)],
        key=lambda x: int(x.split("_")[1]),
    )
    if not vars_found:
        msg = f"No variables found for base name '{base_name}'"
        raise ValueError(msg)

    # Concatenate along new depth dimension
    da = xr.concat([ds[v] for v in vars_found], dim=depth_name)
    da = da.assign_coords({depth_name: range(da.sizes[depth_name])})
    da.name = base_name
    return da


# --- Combine so and thetao into 4D arrays ---
so_4d = combine_levels(ds, "so")
thetao_4d = combine_levels(ds, "thetao")

# --- Add zos if it exists (2D surface variable) ---
data_vars = {"so": so_4d, "thetao": thetao_4d}

data_vars["zos"] = ds["zos"]

# --- Create NEMO-like dataset ---
ds_combined = xr.Dataset(data_vars)

# --- Copy coordinate variables ---
for coord in ["time", "lat", "lon"]:
    if coord in ds.coords:
        ds_combined = ds_combined.assign_coords({coord: ds[coord]})

# --- Save final structured dataset ---
ds_combined.to_netcdf(output_path)
