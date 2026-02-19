"""Reformat OM4 Zarr dataset to combine level-wise variables into 4D variables."""

# Standard library imports
import re

import xarray as xr

# Open the Zarr dataset
input_path = "/home/sg2147/Samudra_data/data.zarr"
output_path = "/home/sg2147/Samudra_data/OM4-reformatted/om4_combined.nc"


ds = xr.open_zarr(input_path, consolidated=True)


# Function to combine level-wise variables
def combine_levelwise_variables(ds, base_name):
    """
    Combine OM4-style levelwise variables into a single 4D variable.

    For example: so_lev_10_0, so_lev_50_0, ... → so(time, lev, y, x)
    """
    pattern = re.compile(f"^{base_name}_lev_([0-9]+)_0$")
    level_vars, levels = [], []

    for var in ds.data_vars:
        match = pattern.match(var)
        if match:
            level = float(match.group(1))
            levels.append(level)
            level_vars.append(ds[var])

    if not level_vars:
        print(f"No variables found for '{base_name}'")
        return None

    # Sort by depth level
    sorted_vars = [v for _, v in sorted(zip(levels, level_vars), key=lambda x: x[0])]
    sorted_levels = sorted(levels)

    # Stack along a new 'lev' dimension
    combined = xr.concat(sorted_vars, dim="lev")
    combined = combined.assign_coords(lev=("lev", sorted_levels))

    # Reorder dimensions: (time, lev, y, x)
    combined = combined.transpose("time", "lev", "y", "x")

    return combined


# Combine salinity and temperature level-wise variables
so_4d = combine_levelwise_variables(ds, "so")  # salinity
thetao_4d = combine_levelwise_variables(ds, "thetao")  # temperature

# Get sea surface height variable
zos = ds.get("zos", None)

# create new dataset with combined variables
output_vars = {}
if so_4d is not None:
    output_vars["so"] = so_4d
if thetao_4d is not None:
    output_vars["thetao"] = thetao_4d
if zos is not None:
    output_vars["zos"] = zos

ds_combined = xr.Dataset(output_vars)

# Add global attributes
ds_combined.attrs["title"] = "OM4 data reformatted to NEMO-like structure"
ds_combined.attrs["source"] = "Converted from OM4 Zarr"
ds_combined.attrs["created_by"] = "convert_om4_zarr_to_nemo_like_nc.py"

# Save to NetCDF
print(f"Saving combined dataset to {output_path}")
ds_combined.to_netcdf(output_path)
print("File saved successfully.")
