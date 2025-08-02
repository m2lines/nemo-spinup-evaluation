import xarray as xr
import os
import glob
from src.variable_aliases import VARIABLE_ALIASES, standardize_variables
from xarray import open_dataset

import os
import glob
import xarray as xr
from typing import Dict


def load_output_variables(
    var_to_file_map: dict[str, str], path: str
) -> dict[str, xr.DataArray]:
    """
    Load output variables as a flat dictionary: variable -> DataArray.
    """
    file_cache: dict[str, xr.Dataset] = {}
    variable_data: dict[str, xr.DataArray] = {}

    for var, relpath in var_to_file_map.items():
        fpath = os.path.join(path, relpath)
        if not os.path.exists(fpath):
            raise FileNotFoundError(f"File {relpath} not found at {fpath}")

        if relpath not in file_cache:
            file_cache[relpath] = xr.open_dataset(fpath)

        # Patch time_counter into ssh if missing
        if var == "ssh" and "time_counter" not in file_cache[relpath]:
            t_file = next(
                (p for v, p in var_to_file_map.items() if v == "temperature"), None
            )
            if t_file and t_file in file_cache:
                file_cache[relpath]["time_counter"] = file_cache[t_file].get(
                    "time_counter", None
                )

        variable_data[var] = file_cache[relpath][var]

    return variable_data


def load_mesh_mask(path: str) -> xr.Dataset:
    """Load the NEMO mesh mask file."""
    return xr.open_dataset(path)


def standardize_variables(ds_or_da, aliases: dict):
    """Standardize variable names using provided aliases."""
    if isinstance(ds_or_da, xr.Dataset):
        rename_dict = {k: v for k, v in aliases.items() if k in ds_or_da.data_vars}
        return ds_or_da.rename(rename_dict)
    elif isinstance(ds_or_da, xr.DataArray):
        name = ds_or_da.name
        if name in aliases:
            return ds_or_da.rename(aliases[name])
        return ds_or_da
    else:
        raise TypeError("Input must be an xarray.Dataset or xarray.DataArray")


def load_dino_outputs(
    mode: str, path: str, setup: dict, VARIABLE_ALIASES: dict
) -> dict[str, xr.Dataset]:
    """
    Load DINO model outputs based on YAML configuration.

    Parameters
    ----------
    mode : str
        One of "restart", "output", or "both".
    path : str
        Base directory where model files are stored.
    setup : dict
        Parsed YAML configuration for the DINO model.
    VARIABLE_ALIASES : dict
        Mapping from raw variable names to standardized names.

    Returns
    -------
    dict[str, xr.Dataset or dict[str, xr.DataArray]]
        A dictionary with keys like 'mesh_mask', 'restart', 'output'.
    """
    data = {}

    # Mesh mask is always required
    if "mesh_mask" not in setup:
        raise ValueError("Mesh mask file is required in the setup configuration.")
    mesh_mask_path = os.path.join(path, setup["mesh_mask"])
    data["mesh_mask"] = load_mesh_mask(mesh_mask_path)

    # Optional restart file
    if mode in ["restart", "both"]:
        restart_file = setup.get("restart_file")
        if restart_file:
            restart_path = os.path.join(path, restart_file)
        else:
            restart_candidates = glob.glob(os.path.join(path, "*restart.nc"))
            if not restart_candidates:
                raise FileNotFoundError("No restart file found matching '*restart.nc'")
            restart_path = restart_candidates[0]
        data["restart"] = xr.open_dataset(restart_path)

    # Output files
    if mode in ["output", "both"] and "output_variables" in setup:
        raw_output = load_output_variables(setup["output_variables"], path)
        standardized_output = {
            std_name: standardize_variables(da, VARIABLE_ALIASES)
            for std_name, da in raw_output.items()
        }
        data["output"] = standardized_output

    # Standardize mesh_mask and restart datasets
    for key in ["mesh_mask", "restart"]:
        if key in data:
            data[key] = standardize_variables(data[key], VARIABLE_ALIASES)

    return data
