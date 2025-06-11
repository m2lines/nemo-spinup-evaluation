"""Functions to load NEMO model output and restart files based on YAML configuration."""

import glob
import os

import xarray as xr

from src.standardise_inputs import VARIABLE_ALIASES, standardise


def load_output_files(setup: dict, path: str) -> dict[str, xr.Dataset]:
    """Load output files based on the setup configuration.

    This function is used to load output files for the NEMO model based
    on the provided setup. It does not load restart files.

    Parameters
    ----------
    setup : dict
        The setup configuration containing output file mappings.
    path : str
        The base path to the output files.

    Returns
    -------
    dict[str, xr.Dataset]
        A dictionary mapping output variable names to their xarray datasets.
    """
    output_data = {}
    for field, filename in setup.items():
        # TODO: check if the file exists
        if not os.path.exists(os.path.join(path, filename)):
            msg = f"File {filename} not found in path {path}."
            raise FileNotFoundError(msg)
        if field == "T_sampled":
            # Special case for T_sampled, which is sampled from T grid
            data_grid_T = xr.open_dataset(os.path.join(path, setup["T"]))
            data_grid_T_sampled = xr.open_dataset(os.path.join(path, filename))
            data_grid_T_sampled["time_counter"] = data_grid_T["time_counter"]
            output_data[field] = data_grid_T_sampled
        else:
            output_data[field] = xr.open_dataset(os.path.join(path, filename))
        # make a correction here: if field is 'T_sampled'

    return output_data


def load_mesh_mask(path: str) -> xr.Dataset:
    """Load the NEMO mesh mask file."""
    ds = xr.open_dataset(path)
    required_vars = ["tmask", "e1t", "e2t", "e3t_0"]
    missing = [var for var in required_vars if var not in ds.variables]
    if missing:
        msg = f"Mesh mask file {path} is missing required variables: {missing}"
        raise ValueError(msg)
    return ds


def load_dino_outputs(mode, path, setup: dict) -> dict[str, xr.Dataset]:
    """Load DINO model inputs based on YAML configuration.

    Parameters
    ----------
    mode : str
        Specifies which data to load: 'restart', 'output', or 'both'.
    path : str
        The base directory path where files are located.
    setup : dict
        Parsed YAML config containing mesh_mask, restart_files, output_files, etc.

    Returns
    -------
    dict[str, xr.Dataset]
        Dictionary mapping source (e.g., 'mesh_mask', 'restart',
    'output') to xarray datasets or variables.
    """
    data = {}

    # Always required
    # TODO: check if the mesh_mask file exists
    if "mesh_mask" not in setup:
        msg = "Mesh mask file is required in the setup configuration."
        raise ValueError(msg)
    mesh_mask_path = setup["mesh_mask"]
    data["mesh_mask"] = load_mesh_mask(os.path.join(path, mesh_mask_path))

    # Optional: restart file
    if mode in ["restart", "both"]:
        # get the restart path using glob in path - it ends with restart.nc
        # open a file that ends with restart.nc
        restart_path = glob.glob(os.path.join(path, "*restart.nc"))[0]
        data["restart"] = xr.open_dataset(restart_path)

    # Optional: grid T, U, V, and sampled T
    if mode in ["output", "both"] and "output_variables" in setup:
        # check if output_variables is
        setup_output = setup["output_variables"]
        data["output"] = load_output_files(setup_output, path)

    # now standardise the data
    for key, dataset in data.items():
        if isinstance(dataset, xr.Dataset):
            data[key] = standardise(dataset, VARIABLE_ALIASES)
        elif isinstance(dataset, dict):  # e.g., output files
            data[key] = {
                k: standardise(v, VARIABLE_ALIASES) for k, v in dataset.items()
            }

    # print(data)
    return data
