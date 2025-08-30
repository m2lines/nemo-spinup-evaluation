"""Compute module for spinup evaluation.

This module provides functions for retrieving, resampling, aligning, and caching
xarray-based variable data at various temporal granularities for metric analysis.
"""

import xarray as xr

from spinup_evaluation.standardise_inputs import VARIABLE_ALIASES, standardise

from .cache import (
    load_aligned_from_disk,
    load_resampled_from_disk,
    save_aligned_to_disk,
    save_resampled_to_disk,
)
from .fix_time import ensure_time, find_time_dimension, resample_like, same_time_axis

GRANULARITY_ORDER = ["10d", "1m", "3m", "1y"]


# USED (1) and (2) public interface
def get_data(
    var,
    granularity,
    variable_file_map,
    cache,
    analysis,
    allow_resampling=True,
    disk_cache_dir="./resampled_cache",
    save_to_cache=True,
):
    """
    Retrieve data for a variable at a given granularity.

    We are using memory and disk cache, and resampling from the closest
    available finer granularity if needed.
    """
    key = (var, granularity)
    if key in cache:
        print("MEMORY CACHE HIT", key)
        return cache[key]

    # aligned artifact takes precedence (skips re-binning later)
    aligned = load_aligned_from_disk(var, granularity, disk_cache_dir)
    if aligned is not None:
        cache[key] = aligned
        return aligned

    # 1. Try direct file first
    if granularity in analysis["direct_files"].get(var, []):
        file_entry = next(
            e for e in variable_file_map[var] if e["granularity"] == granularity
        )
        print(f"Loading {var}@{granularity} directly from {file_entry['file']}")
        ds = xr.open_dataset(file_entry["file"])
        ds = standardise(ds, VARIABLE_ALIASES)
        if var in ds:
            data = ds[var]
            cache[key] = data
            return data
        else:
            msg = f"Variable {var} not found in {file_entry['file']}"
            raise ValueError(msg)

    # 2. Otherwise, resample from closest available finer granularity
    if allow_resampling:
        gran_idx = GRANULARITY_ORDER.index(granularity)
        finer_grans = [
            g
            for g in analysis["direct_files"].get(var, [])
            if GRANULARITY_ORDER.index(g) < gran_idx
        ]
        if not finer_grans:
            msg = (
                f"No finer granularity available for {var} to resample to {granularity}"
            )
            raise ValueError(msg)

        # NOTE: Pick the closest (highest index < gran_idx).
        # Prefer to get finest for accuracy
        # source_gran = max(finer_grans, key=lambda g: GRANULARITY_ORDER.index(g))

        # Pick the finest granularity
        source_gran = min(finer_grans, key=lambda g: GRANULARITY_ORDER.index(g))
        source_file = next(
            e["file"] for e in variable_file_map[var] if e["granularity"] == source_gran
        )

        # Try disk cache for resampled data
        cached_resampled = load_resampled_from_disk(
            var, source_gran, granularity, source_file, disk_cache_dir
        )
        if cached_resampled is not None:
            cache[key] = cached_resampled
            return cached_resampled

        print(f"Resampling {var}: {source_gran} → {granularity}")

        # Get source data (recursive call)
        finer_data = get_data(
            var,
            source_gran,
            variable_file_map,
            cache,
            analysis,
            allow_resampling,
            disk_cache_dir,
            save_to_cache,
        )

        # Resample
        time_dim = find_time_dimension(finer_data)
        if time_dim is None:
            msg = f"No time dimension found in {var}. Dims: {finer_data.dims}"
            raise ValueError(msg)

        # Use fallback_freq for resampling frequency string
        freq_map = {"10d": "10D", "1m": "1ME", "3m": "3ME", "1y": "1YE"}
        resample_kwargs = {time_dim: freq_map[granularity]}
        resampled = finer_data.resample(
            **resample_kwargs, label="right", closed="right"
        ).mean(skipna=True)

        # Optionally save to disk cache
        if save_to_cache:
            save_resampled_to_disk(
                resampled.compute(),
                var,
                source_gran,
                granularity,
                source_file,
                disk_cache_dir,
            )
            # Reload as lazy for memory efficiency
            cached_data = load_resampled_from_disk(
                var, source_gran, granularity, source_file, disk_cache_dir
            )
            cache[key] = cached_data
            return cached_data
        else:
            print("  Not saving to cache (save_to_cache=False)")
            cache[key] = resampled  # Keep as lazy
            return resampled

    msg = f"{var} not available at {granularity} and resampling not allowed"
    raise ValueError(msg)


# Maybe USED (1) - called from run_all_metrics_with_cache and
# run_metrics_intelligently_with_cache
# NOTE: This doesn't do alignment - so it will fail.
def run_metric_with_cache(
    metric_name,
    metric_function,
    required_vars,
    granularity,
    variable_file_map,
    analysis,
    cache=None,
    down_sample=True,
    disk_cache_dir="./resampled_cache",
    save_to_cache=True,
):
    """Run_metric with disk caching."""
    try:
        print(f"Running {metric_name} at {granularity}")
        inputs = [
            get_data(
                var,
                granularity,
                variable_file_map,
                cache,
                analysis,
                down_sample,
                disk_cache_dir,
                save_to_cache,
            )
            for var in required_vars
        ]
        result = metric_function(*inputs)
        print(f"✓ Success: {metric_name}")
        return result
    except (KeyError, ValueError, OSError) as e:
        print(f"✗ Failed: {metric_name} - {e}")
        return None


# USED (1) - called from preload_and_align_all_variables
def load_all_variables_at_granularity(
    gran, variable_file_map, analysis, cache, disk_cache_dir, save_to_cache
):
    """
    Load all variables available at a given granularity (direct or resampled).

    Returns a dict: {var: xarray.DataArray}
    """
    available_vars = [
        var for var, grans in analysis["variable_availability"].items() if gran in grans
    ]
    loaded_vars = {}
    for var in available_vars:
        try:
            data = get_data(
                var,
                gran,
                variable_file_map,
                cache,
                analysis=analysis,
                allow_resampling=True,
                disk_cache_dir=disk_cache_dir,
                save_to_cache=save_to_cache,
            )
            loaded_vars[var] = data
        except (KeyError, ValueError, OSError) as e:
            print(f"✗ Failed to load {var}@{gran}: {e}")
    return loaded_vars


# USED (1) - called from preload_and_align_all_variables
def align_all_to_reference(
    loaded_vars: dict[str, xr.DataArray],
    ref_var: str,
    fallback_freq: str | None = None,
) -> dict[str, xr.DataArray]:
    """
    For each variable, return a version whose time axis exactly matches reference.

    Only resamples when needed; otherwise returns the data unchanged.
    """
    ref = ensure_time(loaded_vars[ref_var])
    aligned = {}
    for name, da in loaded_vars.items():
        da_ensured = ensure_time(da)
        if same_time_axis(da_ensured, ref):
            aligned[name] = da_ensured
        else:
            aligned[name] = resample_like(da_ensured, ref, fallback_freq=fallback_freq)
    return aligned


# USED (1) - public interface
def preload_and_align_all_variables(
    variable_file_map,
    granularities,
    analysis,
    disk_cache_dir="./resampled_cache",
    save_to_cache=False,
    persist_aligned=False,
):
    """Preload and align all variables to a common reference."""
    cache = {}
    aligned_vars_by_gran = {}

    for gran in granularities:
        print(f"\n=== GRANULARITY: {gran} ===")
        loaded_vars = load_all_variables_at_granularity(
            gran, variable_file_map, analysis, cache, disk_cache_dir, save_to_cache
        )

        direct_vars = analysis["direct_from_file_by_granularity"].get(gran, [])
        ref_var = next((var for var in direct_vars if var in loaded_vars), None)
        if ref_var is None:
            print(f"No reference variable for {gran}")
            aligned_vars_by_gran[gran] = loaded_vars
            continue
        print(f"Using '{ref_var}' as reference for alignment at {gran}")

        # Get calendar and units from the reference variable's time coordinate
        # ref_data = loaded_vars[ref_var]
        # time_dim = find_time_dimension(ref_data)
        # calendar = ref_data[time_dim].attrs.get("calendar", "360_day")
        # units = ref_data[time_dim].attrs.get("units", "days since 0001-01-01")

        # Align all variables to the reference bins
        # aligned_vars = align_all_variables_to_reference_bins(
        #     loaded_vars, ref_var, calendar, units, time_dim=time_dim,
        #     fallback_freq=gran
        # )

        # NOTE: There's enough information in loaded vars
        # or in aligned_vars or even cache at this given granularity
        # to get the granularity of the data he target was sampled from.
        aligned_vars = align_all_to_reference(loaded_vars, ref_var, fallback_freq=gran)

        # Optional saving of aligned files to disk
        # Note that these are different to
        if persist_aligned:
            for var, da in aligned_vars.items():
                save_aligned_to_disk(da, var, gran, ref_var, cache_dir=disk_cache_dir)

        for var, data in aligned_vars.items():
            cache[(var, gran)] = data

    return cache
