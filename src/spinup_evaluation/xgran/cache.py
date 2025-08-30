"""Cache utilities for resampled and aligned xarray data.

This module provides functions to save, load, and manage cached
resampled and time-aligned xarray DataArrays, including metadata
handling and cache statistics.
"""

import glob
import hashlib
import json
import os
from pathlib import Path

import pandas as pd
import xarray as xr

from .fix_time import ensure_time


# USED (1) - called from get_data
def save_resampled_to_disk(
    data, var, source_gran, target_gran, source_file, cache_dir="./resampled_cache"
):
    """Save resampled data to cache."""
    os.makedirs(cache_dir, exist_ok=True)

    cache_file = get_cache_filename(var, source_gran, target_gran, cache_dir)
    metadata_file = cache_file.replace(".nc", "_metadata.json")

    try:
        # Save data
        print(f"Caching {var} {source_gran}→{target_gran} to {cache_file}")
        data.to_netcdf(cache_file)

        # Save metadata
        metadata = {
            "variable": var,
            "source_granularity": source_gran,
            "target_granularity": target_gran,
            "source_file": source_file,
            "source_hash": get_file_hash(source_file),
            "created": pd.Timestamp.now().isoformat(),
            "shape": list(data.shape),
            "chunks": str(data.chunks) if hasattr(data, "chunks") else None,
        }

        with open(metadata_file, "w") as f:
            json.dump(metadata, f, indent=2)

        print(f"✓ Cached {var} {source_gran}→{target_gran}")

    except (OSError, TypeError) as e:
        print(f"Cache save failed: {e}")


# USED (1) - called from get_data
def load_resampled_from_disk(
    var, source_gran, target_gran, source_file, cache_dir="./resampled_cache"
):
    """Load resampled data from cache if available and valid."""
    cache_file = get_cache_filename(var, source_gran, target_gran, cache_dir)
    metadata_file = cache_file.replace(".nc", "_metadata.json")

    if not (os.path.exists(cache_file) and os.path.exists(metadata_file)):
        return None

    # Check if source file has changed
    try:
        with open(metadata_file, "r") as f:
            metadata = json.load(f)

        current_hash = get_file_hash(source_file)
        if metadata.get("source_hash") != current_hash:
            print(
                (
                    f"Source file changed, cache invalid for {var} "
                    f"{source_gran}→{target_gran}"
                )
            )
            return None

        print(f"Loading {var} {source_gran}→{target_gran} from cache")
        return xr.open_dataarray(cache_file, chunks={"time": 100})

    except (OSError, json.JSONDecodeError) as e:
        print(f"Cache load failed: {e}")
        return None


# USED (1)
def save_aligned_to_disk(
    data: xr.DataArray,
    var: str,
    target_gran: str,
    ref_var: str,
    cache_dir: str = "./resampled_cache",
) -> str:
    """Save time-aligned data to cache."""
    Path(cache_dir).mkdir(parents=True, exist_ok=True)
    th = _time_hash(data)
    fname = f"{var}_ALIGNED_to_{target_gran}_{th}.nc"
    fpath = os.path.join(cache_dir, fname)
    ensure_time(data).to_netcdf(fpath)

    meta = {
        "variable": var,
        "target_granularity": target_gran,
        "aligned": True,
        "reference_variable": ref_var,
        "time_hash": th,
        "calendar": ensure_time(data)["time"].attrs.get("calendar"),
        "units": ensure_time(data)["time"].attrs.get("units"),
    }
    with open(fpath.replace(".nc", "_metadata.json"), "w") as f:
        json.dump(meta, f, indent=2)
    return fpath


# USED (1) (mainly) and (2) - need to check how this affects (2)
def load_aligned_from_disk(var, target_gran, cache_dir="./resampled_cache"):
    """Load aligned data from cache if available."""
    pattern = os.path.join(cache_dir, f"{var}_ALIGNED_to_{target_gran}_*.nc")
    hits = sorted(glob.glob(pattern))
    if hits:
        print(f"Loading {var} ALIGNED→{target_gran} from cache")
        return xr.open_dataarray(hits[-1], chunks={"time": 100})
    return None


# USED (1)
def _time_hash(da: xr.DataArray) -> str:
    t = ensure_time(da)["time"].values
    return hashlib.sha1("".join(map(str, t)).encode()).hexdigest()[:16]  # noqa 5324


def get_cache_filename(var, source_gran, target_gran, cache_dir="./resampled_cache"):
    """Generate consistent cache filename."""
    cache_key = f"{var}_{source_gran}_to_{target_gran}"
    return os.path.join(cache_dir, f"{cache_key}.nc")


def get_file_hash(filepath):
    """Get hash of source file to detect changes."""
    with open(filepath, "rb") as f:
        # Hash first 1MB for speed
        return hashlib.md5(f.read(1024 * 1024)).hexdigest()  # noqa 5324


def write_metadata(cache_file, meta: dict):
    """Write metadata for cached file."""
    meta_file = cache_file.replace(".nc", "_metadata.json")
    with open(meta_file, "w") as f:
        json.dump(meta, f, indent=2)


def show_cache_stats(cache_dir="./resampled_cache"):
    """Show cache statistics."""
    if not os.path.exists(cache_dir):
        print("No disk cache found")
        return

    nc_files = [f for f in os.listdir(cache_dir) if f.endswith(".nc")]
    if not nc_files:
        print("Cache directory exists but is empty")
        return

    total_size = 0
    cached_items = []

    for nc_file in nc_files:
        file_path = os.path.join(cache_dir, nc_file)
        size = os.path.getsize(file_path)
        total_size += size

        # Parse filename: "temperature_10d_to_1m.nc"
        name_parts = nc_file[:-3].split("_to_")
        PARTS = 2
        if len(name_parts) == PARTS:
            source_parts = name_parts[0].split("_")
            if len(source_parts) >= PARTS:
                var = "_".join(source_parts[:-1])
                source_gran = source_parts[-1]
                target_gran = name_parts[1]
                cached_items.append((var, source_gran, target_gran, size))

    print("\n=== DISK CACHE STATISTICS ===")
    print(f"Cache directory: {cache_dir}")
    print(f"Total files: {len(nc_files)}")
    print(f"Total size: {total_size / 1024**3:.2f} GB")

    if cached_items:
        print("\nCached resampled data:")
        for var, source_gran, target_gran, size in sorted(cached_items):
            print(f"  {var}: {source_gran}→{target_gran} ({size / 1024**2:.1f} MB)")
    else:
        print("No valid cached items found")
