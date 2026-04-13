"""Functions to load NEMO model output and restart files based on YAML configuration."""

import glob
import logging
import os
from pathlib import Path
from typing import Dict, Mapping, Optional, Union

import xarray as xr

from nemo_spinup_evaluation.standardise_inputs import VARIABLE_ALIASES, standardise

logger = logging.getLogger(__name__)

VarSpec = Mapping[str, Union[str, Mapping[str, str]]]


def _open_cached(
    cache: Dict[str, xr.Dataset], base: str, relpath: str, lazy: bool = True
) -> xr.Dataset:
    """Open a dataset once, caching by relative path."""
    if relpath not in cache:
        full = os.path.join(base, relpath)
        if not os.path.exists(full):
            msg = f"File not found: {full}"
            raise FileNotFoundError(msg)
        open_kw = {"chunks": {}} if lazy else {}
        cache[relpath] = xr.open_dataset(full, **open_kw)
    return cache[relpath]


MAX_DISPLAYED_VARIABLES = 20


def _infer_var_name(ds: xr.Dataset, canon: str) -> str:
    """
    Choose the actual variable name inside a dataset given a canonical name.

    Uses VARIABLE_ALIASES (e.g., {"temperature": ["toce","thetao","temp", ...], ...})
    and falls back to the canonical name if present.
    """
    # exact match first
    if canon in ds.variables:
        return canon

    # alias match next
    aliases = VARIABLE_ALIASES.get(canon, [])

    for name in aliases:
        if name in ds.variables:
            return name

    # fail loudly if we can't find a matching variable
    available = list(ds.variables)
    msg = (
        f"Could not find a variable for '{canon}'. "
        f"Checked aliases {aliases!r} and '{canon}'. "
        f"Available in file: {available[:MAX_DISPLAYED_VARIABLES]}"
        f"{'...' if len(available) > MAX_DISPLAYED_VARIABLES else ''}"
    )
    raise KeyError(msg)


def _normalise_var_specs(
    var_specs: VarSpec, file_cache: Dict[str, xr.Dataset], base: str, lazy: bool = True
) -> Dict[str, Dict[str, str]]:
    """
    Normalise variable specifications in yaml file.

    Accepts either:
      simple form:  {"temperature": "grid_T_3D.nc"}
      rich form:    {"temperature": {"file": "grid_T_3D.nc", "var": "toce"}}

    Returns canonical mapping:
      {canon: {"file": str, "var": str}}

    The file_cache is used to avoid reopening the same dataset multiple times.
    It is built up as the function processes each variable specification.

    Parameters
    ----------
    var_specs : VarSpec
        The variable specifications to normalise.
    file_cache : Dict[str, xr.Dataset]
        A cache of opened xarray datasets, keyed by file path.
    base : str
        The base directory for resolving relative file paths.
    lazy : bool
        Whether to use dask-backed lazy loading.

    Returns
    -------
    Dict[str, Dict[str, str]]
        A mapping of canonical variable names to their specifications.
    """
    normalised: Dict[str, Dict[str, str]] = {}
    for canon, spec in var_specs.items():
        if isinstance(spec, str):
            # simple form → infer variable name
            ds = _open_cached(file_cache, base, spec, lazy=lazy)
            var_name = _infer_var_name(ds, canon)
            normalised[canon] = {"file": spec, "var": var_name}

        elif isinstance(spec, Mapping):
            # rich form should include 'file' and 'var'
            # the 'var' corresponds to the variable name in the dataset
            # therefore we don't need to infer it
            if "file" not in spec:
                msg = f"Missing 'file' for variable '{canon}'."
                raise ValueError(msg)
            fname = str(spec["file"])
            ds = _open_cached(file_cache, base, fname, lazy=lazy)

            if "var" in spec:
                var_name = str(spec["var"])  # user-specified → trust it
            else:
                var_name = _infer_var_name(ds, canon)

            entry = {"file": fname, "var": var_name}
            normalised[canon] = entry

        else:
            msg = f"Bad spec for '{canon}': {spec!r}"
            raise TypeError(msg)

    return normalised


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


def resolve_mesh_mask(mesh_mask: str, sim_path: str) -> Path:
    """Resolve the mesh mask path, handling absolute and relative paths."""
    p = Path(mesh_mask)
    candidate = p if p.is_absolute() else Path(sim_path) / mesh_mask
    if not candidate.exists():
        hint = (
            "Set `mesh_mask` to an absolute path in your YAML, "
            "or ensure it exists under --sim-path."
        )
        msg = f"Mesh mask file not found: {candidate}. {hint}"
        raise FileNotFoundError(msg)

    return candidate


def load_mesh_mask(path: Path, lazy: bool = True) -> xr.Dataset:
    """Load the NEMO mesh mask file and validate required fields."""
    if not path.exists():
        msg = f"Mesh mask file not found: {path}"
        raise FileNotFoundError(msg)
    open_kw = {"chunks": {}} if lazy else {}
    ds = xr.open_dataset(path, **open_kw)
    required_vars = ["tmask", "e1t", "e2t", "e3t_0"]
    missing = [v for v in required_vars if v not in ds.variables]
    if missing:
        msg = f"Mesh mask file {path} is missing required variables: {missing}"
        raise ValueError(msg)
    return ds


def get_restart_file_path(base: str, restart_hint: Optional[str]) -> Optional[str]:
    """
    Get the restart file path based on base directory and hint.

    If restart_hint ends with '.nc', it is treated as a direct file name.
    Otherwise, we search for files matching '*{hint}*.nc' in the base directory.

    Parameters
    ----------
    base : str
        The base directory to search for restart files.
    restart_hint : Optional[str]
        An optional hint for the restart file name (e.g., 'restart').

    Returns
    -------
    Optional[str]
        The path to the restart file if found, otherwise None.
    """
    if restart_hint and restart_hint.endswith(".nc"):
        candidate = os.path.join(base, restart_hint)
        return candidate if os.path.exists(candidate) else None

    pattern_core = restart_hint if restart_hint else "restart"
    matches = sorted(glob.glob(os.path.join(base, f"*{pattern_core}*.nc")))
    return matches[0] if matches else None


def load_grid_variables(
    base: str,
    output_specs: VarSpec,
    files_cache: Dict[str, xr.Dataset],
    lazy: bool = True,
) -> Dict[str, xr.DataArray]:
    """
    Build a dict of {canonical_name: DataArray} with a single open per file.

    Supports simple and rich variable specs (see _normalise_var_specs).

    Parameters
    ----------
    base : str
        The base directory for loading data files.
    output_specs : VarSpec
        The variable specifications for the output data.
    files_cache : Dict[str, xr.Dataset]
        A cache of opened xarray datasets, keyed by file path.
    lazy : bool
        Whether to use dask-backed lazy loading.

    Returns
    -------
    Dict[str, xr.DataArray]
        A dictionary mapping canonical variable names to their DataArray objects.
    """
    norm = _normalise_var_specs(output_specs, files_cache, base, lazy=lazy)

    # Pull the arrays
    out: Dict[str, xr.DataArray] = {}
    for canon, spec in norm.items():
        ds = _open_cached(files_cache, base, spec["file"], lazy=lazy)
        out[canon] = ds[spec["var"]]

    return out


def load_dino_data(
    mode: str,
    base: str,
    setup: Mapping[str, object],
    do_standardise: bool = True,
    lazy: bool = True,
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
    do_standardise : bool
       Whether to standardise variable names using aliases.
    lazy : bool
       Whether to use dask-backed lazy loading.

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

    if "mesh_mask" not in setup:
        msg = "setup must specify 'mesh_mask'."
        raise ValueError(msg)

    open_kw = {"chunks": {}} if lazy else {}
    base = os.path.abspath(base)

    # Resolve paths upfront
    mesh_mask_path = resolve_mesh_mask(str(setup["mesh_mask"]), base)
    restart_hint = str(setup.get("restart_files") or "")
    restart_path = get_restart_file_path(base, restart_hint)

    if restart_path is None and mode in ("restart", "both"):
        msg = "No restart file found matching pattern."
        raise FileNotFoundError(msg)
    if restart_path is None and mode == "output":
        logger.warning(
            "No restart file found. The check_density_computed metric will be skipped."
        )

    # Load data
    data: Dict[str, object] = {}
    files_cache: Dict[str, xr.Dataset] = {}

    data["mesh_mask"] = load_mesh_mask(mesh_mask_path, lazy=lazy)

    if restart_path is not None:
        data["restart"] = xr.open_dataset(restart_path, **open_kw)
    else:
        data["restart"] = None

    data["grid"] = {}
    if mode in ("output", "both") and "output_variables" in setup:
        var_specs: VarSpec = setup["output_variables"]
        vars_map = load_grid_variables(base, var_specs, files_cache, lazy=lazy)
        data["grid"].update(vars_map)

    # Build paths metadata
    paths: Dict[str, object] = {
        "base": base,
        "mesh_mask": str(mesh_mask_path),
        "restart": restart_path,
        "output_files": [os.path.join(base, relpath) for relpath in files_cache],
    }

    data["files"] = files_cache
    data["paths"] = paths

    # standardise names after loading has taken place
    if do_standardise:
        # 1) mesh_mask
        data["mesh_mask"] = standardise(data["mesh_mask"], VARIABLE_ALIASES)

        # 2) restart if present
        if data["restart"] is not None:
            data["restart"] = standardise(data["restart"], VARIABLE_ALIASES)

        # 3) each requested variable DataArray
        data["grid"] = {
            name: standardise(da, VARIABLE_ALIASES) for name, da in data["grid"].items()
        }

    # Check for expected coordinates after all other processing
    required_coords = ("time_counter", "depth", "nav_lat", "nav_lon")
    _check_required_coords(data["mesh_mask"], required_coords, "mesh mask")

    if data["restart"] is not None:
        _check_required_coords(data["restart"], required_coords, "restart file")

    for name, da in data["grid"].items():
        _check_required_coords(da, ("time_counter", "nav_lat", "nav_lon"), name)

    # Check that grid variables are temporally aligned
    if data["grid"]:
        _check_grid_time_alignment(data["grid"])

    return data
