"""Functions to load NEMO model output and restart files based on YAML configuration."""

import glob
import os
from typing import Dict, Mapping, Optional, Union

import xarray as xr

from spinup_evaluation.standardise_inputs import VARIABLE_ALIASES, standardise

VarSpec = Mapping[str, Union[str, Mapping[str, str]]]


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
    var_specs: VarSpec, file_cache: Dict[str, xr.Dataset], base: str
) -> Dict[str, Dict[str, str]]:
    """
    Normalise variable specifications in yaml file.

    Accepts either:
      simple form:  {"temperature": "grid_T_3D.nc"}
      rich form:    {"temperature": {"file": "grid_T_3D.nc", "var": "teoce",
      "time_from" : "density"}

    Returns canonical mapping:
      {canon: {"file": str, "var": str, "time_from": Optional[str]}}

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

    Returns
    -------
    Dict[str, Dict[str, str]]
        A mapping of canonical variable names to their specifications.
    """
    normalised: Dict[str, Dict[str, str]] = {}
    for canon, spec in var_specs.items():
        if isinstance(spec, str):
            # simple form → infer variable name
            ds = _open_cached(file_cache, base, spec)
            var_name = _infer_var_name(ds, canon)
            normalised[canon] = {"file": spec, "var": var_name}

        elif isinstance(spec, Mapping):
            # rich form should include 'file', 'var', and optionally 'time_from'
            # the 'var' corresponds to the variable name in the dataset
            # therefore we don't need to infer it
            if "file" not in spec:
                msg = f"Missing 'file' for variable '{canon}'."
                raise ValueError(msg)
            fname = str(spec["file"])
            ds = _open_cached(file_cache, base, fname)

            if "var" in spec:
                var_name = str(spec["var"])  # user-specified → trust it
            else:
                var_name = _infer_var_name(ds, canon)

            entry = {"file": fname, "var": var_name}
            if "time_from" in spec:
                entry["time_from"] = str(spec["time_from"])
            normalised[canon] = entry

        else:
            msg = f"Bad spec for '{canon}': {spec!r}"
            raise TypeError(msg)

    return normalised


def load_mesh_mask(path: str) -> xr.Dataset:
    """Load the NEMO mesh mask file and validate required fields."""
    ds = xr.open_dataset(path)
    required_vars = ["tmask", "e1t", "e2t", "e3t_0"]
    missing = [v for v in required_vars if v not in ds.variables]
    if missing:
        msg = f"Mesh mask file {path} is missing required variables: {missing}"
        raise ValueError(msg)
    return ds


def _find_restart_ds(base: str, restart_hint: Optional[str]) -> Optional[xr.Dataset]:
    """
    Find and open a restart dataset.

    If restart_hint is provided (e.g., 'restart'), we look for *{hint}*.nc
    Otherwise we try '*restart*.nc'

    Parameters
    ----------
    base : str
        The base directory to search for restart files.
    restart_hint : Optional[str]
        An optional hint for the restart file name (e.g., 'restart').

    Returns None if not found.
    """
    pattern_core = restart_hint if restart_hint else "restart"
    matches = sorted(glob.glob(os.path.join(base, f"*{pattern_core}*.nc")))
    return xr.open_dataset(matches[0]) if matches else None


def load_grid_variables(
    base: str, output_specs: VarSpec, files_cache: Dict[str, xr.Dataset]
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

    Returns
    -------
    Dict[str, xr.DataArray]
        A dictionary mapping canonical variable names to their DataArray objects.
    """
    norm = _normalise_var_specs(output_specs, files_cache, base)

    # Pull the arrays
    out: Dict[str, xr.DataArray] = {}
    for canon, spec in norm.items():
        ds = _open_cached(files_cache, base, spec["file"])
        out[canon] = ds[spec["var"]]

    # NOTE: temporary hard-coded time donor to avoid label-mismatch  empties; remove
    # when YAML 'time_from' returns

    recipient = "ssh"
    donor = "temperature"
    if donor:
        if donor not in out:
            msg = f"time_from donor '{donor}' not found among requested variables."
            raise KeyError(msg)
        donor_time = out[donor].coords.get("time_counter", None)
        if donor_time is None:
            msg = f"Donor '{donor}' has no 'time_counter' coordinate to copy."
            raise ValueError(msg)

        out[recipient] = out[recipient].assign_coords(time_counter=donor_time)

    return out


def load_dino_data(
    mode: str,
    base: str,
    setup: Mapping[str, object],
    do_standardise: bool = True,
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

    Returns
    -------
      {
        "mesh_mask": xr.Dataset,
        "restart": xr.Dataset | None,
        "grid": {canon: xr.DataArray, ...},
        "files": {relative_path: xr.Dataset, ...}  # the shared cache
      }
    """
    data: Dict[str, object] = {}
    files_cache: Dict[str, xr.Dataset] = {}

    # mesh mask (required)
    if "mesh_mask" not in setup:
        msg = "setup must specify 'mesh_mask'."
        raise ValueError(msg)
    mesh_mask_path = os.path.join(base, str(setup["mesh_mask"]))
    data["mesh_mask"] = load_mesh_mask(mesh_mask_path)

    # restart (optional / controlled by mode)
    restart_hint = str(setup.get("restart_files") or "")
    if mode in ("restart", "both"):
        data["restart"] = _find_restart_ds(base, restart_hint)
        if data["restart"] is None:
            msg = "No restart file found matching pattern."
            raise FileNotFoundError(msg)

    # outputs (optional / controlled by mode)
    data["grid"] = {}
    if mode in ("output", "both") and "output_variables" in setup:
        var_specs: VarSpec = setup["output_variables"]  # simple or rich form accepted
        # Load variables; this will populate a cache of files so
        # we do not have to keep reopening files we already opened
        data["restart"] = _find_restart_ds(base, restart_hint)
        vars_map = load_grid_variables(base, var_specs, files_cache)
        data["grid"].update(vars_map)

    # expose the file cache
    data["files"] = files_cache
    # dict_keys(['grid_T_3D.nc', 'grid_T_2D.nc', 'grid_U_3D.nc', 'grid_V_3D.nc'])

    # standardise names after loading has taken place
    if do_standardise:
        # 1) mesh_mask
        data["mesh_mask"] = standardise(data["mesh_mask"], VARIABLE_ALIASES)

        # 2) restart if present
        if data.get("restart") is not None:
            data["restart"] = standardise(data["restart"], VARIABLE_ALIASES)

        # 3) each requested variable DataArray
        data["grid"] = {
            name: standardise(da, VARIABLE_ALIASES) for name, da in data["grid"].items()
        }

    return data
