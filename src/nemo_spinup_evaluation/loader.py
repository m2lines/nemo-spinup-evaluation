"""Functions to load NEMO model output and restart files based on YAML configuration."""

import os
from pathlib import Path
from typing import Dict, Mapping, cast

import xarray as xr

VarSpec = Mapping[str, Mapping[str, str]]
VarMap = Mapping[str, list[str]]


def _open_cached(cache: Dict[str, xr.Dataset], base: str, relpath: str) -> xr.Dataset:
    """Open a dataset once, caching by relative path."""
    if relpath not in cache:
        full = os.path.join(base, relpath)
        if not os.path.exists(full):
            msg = f"File not found: {full}"
            raise FileNotFoundError(msg)
        cache[relpath] = xr.open_dataset(full)
    return cache[relpath]


MAX_DISPLAYED_VARIABLES = 20


def _check_required_coords(
    data: xr.DataArray | xr.Dataset, required: tuple[str, ...], name: str
):
    """
    Check for expected coordinates in a dataset.

    Parameters
    ----------
    data : xr.DataArray | xr.Dataset
       The dataset to check for required coordinates.
    required : tuple[str, ...]
        List of required coordinate names.
    name : str
        Dataset name to use in error messages.

    Raises
    ------
    KeyError
        When a specified coordinate is missing.
    """
    missing = [coord for coord in required if coord not in data.coords]
    if missing:
        msg = f"Required coordinate(s) missing in {name}: {', '.join(missing)}"
        raise KeyError(msg)


def _check_grid_time_alignment(grid_data: Mapping[str, xr.DataArray | xr.Dataset]):
    """
    Check grid variables have aligned time coordinates to ensure temporal consistency.

    Parameters
    ----------
    grid_data : Mapping[str, xr.DataArray | xr.Dataset]
        A dictionary of all grid variables and corresponding DataArray.

    Raises
    ------
    ValueError
        When the time_counter of any variable differs from the others.
    """
    # Use the first variable as a reference
    first_var = next(iter(grid_data))
    ref_time = grid_data[first_var]["time_counter"]

    # Compare all variables to the reference
    for name, da in grid_data.items():
        # .equals is not used here because grid variables can optionally include
        # time_centered as an auxiliary coordinate which would cause the comparison to
        # fail regardless of time_counter equality
        if (
            len(da.time_counter) != len(ref_time.time_counter)
            or not (ref_time.values == da.time_counter.values).all()
        ):
            msg = (
                "Time coordinates for grid variable "
                f"{name} differ from {first_var}. "
                "All grid files must have identical time steps "
                "with the same start time and frequency."
            )
            raise ValueError(msg)


def resolve_file_path(file_name: str, sim_path: str) -> Path:
    """Resolve a file path, handling absolute and relative paths."""
    p = Path(file_name)
    candidate = p if p.is_absolute() else Path(sim_path) / file_name
    if not candidate.is_file():
        hint = (
            "Ensure `mesh_mask` and `restart` exist under --sim-path, "
            "or use absolute paths in your YAML."
        )
        msg = f"File not found: {candidate}. {hint}"
        raise FileNotFoundError(msg)

    return candidate


def load_mesh_mask(path: Path) -> xr.Dataset:
    """Load the NEMO mesh mask file and validate required fields."""
    if not path.exists():
        msg = f"Mesh mask file not found: {path}"
        raise FileNotFoundError(msg)
    ds = xr.open_dataset(path)
    required_vars = ["tmask", "e1t", "e2t", "e3t_0"]
    missing = [v for v in required_vars if v not in ds.variables]
    if missing:
        msg = f"Mesh mask file {path} is missing required variables: {missing}"
        raise ValueError(msg)
    return ds


def load_grid_variables(
    base: str, output_specs: VarSpec, files_cache: Dict[str, xr.Dataset]
) -> Dict[str, xr.DataArray]:
    """
    Build a dict of {canonical_name: DataArray} with a single open per file.

    Parameters
    ----------
    base : str
        The base directory for loading data files.
    output_specs : VarSpec
        The variable specifications for the output data.
    files_cache : Dict[str, xr.Dataset]
        A cache of opened xarray datasets, keyed by file path.

    Returns
    -------
    Dict[str, xr.DataArray]
        A dictionary mapping canonical variable names to their DataArray objects.
    """
    # Pull the arrays
    out: Dict[str, xr.DataArray] = {}
    for canon, spec in output_specs.items():
        ds = _open_cached(files_cache, base, spec["file"])
        # Select specified variable
        out[canon] = ds[spec["var"]]

    return out


def standardise_vars(
    data: xr.DataArray | xr.Dataset, variable_map: VarMap
) -> xr.DataArray | xr.Dataset:
    """
    Rename variables/coords/dims to the canonical field names.

    The single variable in a DataArray is not renamed.
    All variables in Datasets are renamed.
    Coords are renamed in both DataArray and Dataset.
    nav_lat and nav_lon are promoted to coordinates.

    Parameters
    ----------
    data : xr.DataArray | xr.Dataset
        Input data to be standardised.
    variable_map : VarMap
        Mapping of canonical field names and their variations for conversion.

    Returns
    -------
    xr.DataArray | xr.Dataset
        Standardised data, maintaining the same type as the input.
    """
    rename_map = {}
    for std, aliases in variable_map.items():
        for alias in aliases:
            if hasattr(data, "variables") and alias in data.variables:
                rename_map[alias] = std
                break
            if alias in data.coords:
                rename_map[alias] = std
                break

    data = data.rename(rename_map)

    # Promote nav_lat and nav_lon to coordinates to enable .sel() indexing.
    # Previously, they were not inherited by dataarray subsets, causing thresholding
    # to use integer indices instead of degrees.
    # See PR [#76](https://github.com/m2lines/nemo-spinup-evaluation/pull/76
    # for more details.
    for name in ("nav_lat", "nav_lon"):
        if name in data and name not in data.coords:
            data = data.set_coords(name)  # zero-copy promotion to coordinate

    return data


def load_dino_data(
    mode: str,
    base: str,
    setup: Mapping[str, object],
) -> Dict[str, object]:
    """
    Load DINO data according to YAML setup.

    Parameters
    ----------
    mode : str
       The mode of operation (e.g., "restart", "output", "both").
    base : str
       The base directory for loading data files.
    setup : Mapping[str, object]
       A mapping containing the YAML configuration for data loading.

    Returns
    -------
    dict
        A dictionary with the following keys:

        - ``mesh_mask``: xr.Dataset
        - ``restart``: xr.Dataset or None
        - ``grid``: dict mapping canonical name to xr.DataArray
        - ``files``: dict mapping relative path to xr.Dataset
        - ``paths``: dict with keys ``base``, ``mesh_mask``, ``restart``,
          ``output_files``

    """
    if mode not in ["output", "restart", "both"]:
        msg = "Mode must be one of 'output', 'restart', 'both'"
        raise ValueError(msg)

    base = os.path.abspath(base)

    data: Dict[str, object] = {}
    files_cache: Dict[str, xr.Dataset] = {}

    # Initialize paths dictionary
    paths: Dict[str, object] = {
        "base": base,
        "mesh_mask": None,
        "restart": None,
        "output_files": [],
    }

    # mesh mask (required)
    if "mesh_mask" not in setup:
        msg = "setup must specify 'mesh_mask'."
        raise ValueError(msg)

    # Resolve mesh mask path
    mesh_mask_path = resolve_file_path(str(setup["mesh_mask"]), base)
    paths["mesh_mask"] = str(mesh_mask_path)
    data["mesh_mask"] = load_mesh_mask(mesh_mask_path)

    # restart (optional / controlled by mode)
    data["restart"] = None
    restart_file = str(setup.get("restart") or "")

    if mode in ("restart", "both"):
        restart_path = resolve_file_path(restart_file, base)
        paths["restart"] = restart_path
        data["restart"] = xr.open_dataset(restart_path)

    # outputs (optional / controlled by mode)
    data["grid"] = {}
    if mode in ("output", "both"):
        if "output_variables" not in setup:
            msg = "Setup file missing output_variables section."
            raise KeyError(msg)

        # Load grid variables and store paths of each file loaded
        var_specs = cast(VarSpec, setup["output_variables"])
        data["grid"] = load_grid_variables(base, var_specs, files_cache)
        paths["output_files"] = [os.path.join(base, relpath) for relpath in files_cache]

        # Check grid variables for temporal alignment
        _check_grid_time_alignment(data["grid"])

        # Restart file is optional in output mode
        if mode == "output" and restart_file:
            restart_path = resolve_file_path(restart_file, base)
            paths["restart"] = restart_path
            data["restart"] = xr.open_dataset(restart_path)

    # expose the file cache
    data["files"] = files_cache
    data["paths"] = paths

    # standardise names after loading has taken place
    if "variable_map" in setup:
        var_map = cast(VarMap, setup["variable_map"])

        # 1) mesh_mask
        data["mesh_mask"] = standardise_vars(data["mesh_mask"], var_map)

        # 2) restart if present
        if data["restart"] is not None:
            data["restart"] = standardise_vars(data["restart"], var_map)

        # 3) each requested variable DataArray
        data["grid"] = {
            name: standardise_vars(da, var_map) for name, da in data["grid"].items()
        }

    # Check for expected coordinates after all other processing
    required_coords = ("time_counter", "depth", "nav_lat", "nav_lon")
    _check_required_coords(data["mesh_mask"], required_coords, "mesh mask")

    if data["restart"] is not None:
        _check_required_coords(data["restart"], required_coords, "restart file")

    for name, da in data["grid"].items():
        _check_required_coords(da, ("time_counter", "nav_lat", "nav_lon"), name)

    return data
