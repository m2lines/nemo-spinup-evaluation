"""Reformat OM4 Zarr dataset to combine level-wise variables into 4D variables."""

# Standard library imports
import re

import xarray as xr

# Open the Zarr dataset
input_path = "/home/sg2147/Samudra_data/data.zarr"
output_path = "/home/sg2147/Samudra_data/OM4_reformatted/om4_combined.nc"


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
    levels, level_vars = zip(*sorted(zip(levels, level_vars), key=lambda x: x[0]))

    # Concatenate along new dimension
    combined = xr.concat(level_vars, dim="lev")

    # Assign coordinate for lev
    combined = combined.assign_coords(lev=("lev", list(levels)))

    print(f"Combined {len(levels)} levels for '{base_name}'")
    return combined


# Combine salinity and temperature level-wise variables
so_4d = combine_levelwise_variables(ds, "so")  # salinity
thetao_4d = combine_levelwise_variables(ds, "thetao")  # temperature

# Get sea surface height variable
zos = ds.get("zos", None)
