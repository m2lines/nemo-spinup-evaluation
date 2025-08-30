"""Utility functions for handling and aligning time dimensions in xarray DataArrays."""

import cftime
import numpy as np
import xarray as xr


def ensure_time(da: xr.DataArray) -> xr.DataArray:
    """Ensure the DataArray has a 'time' dimension, renaming if needed."""
    for cand in ("time", "time_counter", "t", "T", "Time"):
        if cand in da.dims:
            return da if cand == "time" else da.rename({cand: "time"})
    for d in da.dims:
        if "time" in d.lower():
            return da.rename({d: "time"})
    msg = "No time-like dimension found"
    raise ValueError(msg)


def same_time_axis(a: xr.DataArray, b: xr.DataArray) -> bool:
    """Check if two DataArrays have the same time axis."""
    a, b = ensure_time(a), ensure_time(b)
    if a.sizes["time"] != b.sizes["time"]:
        return False
    return np.array_equal(a["time"].values, b["time"].values)


def resample_like(
    data: xr.DataArray, ref: xr.DataArray, fallback_freq=None
) -> xr.DataArray:
    """Bin `data` to the bins implied by `ref.time` and label with `ref.time`."""
    data, ref = ensure_time(data), ensure_time(ref)
    ref_time = ref["time"].values
    cal = ref["time"].attrs.get("calendar", "360_day")
    units = ref["time"].attrs.get("units", "days since 0001-01-01")

    # build edges from ref labels (centre-stamped compatible)
    t_num = cftime.date2num(ref_time, units, cal)

    MIN_LENGTH = 2
    if len(t_num) < MIN_LENGTH:
        msg = "Reference time axis too short to define bins"
        raise ValueError(msg)
    mids = 0.5 * (t_num[:-1] + t_num[1:])
    first_edge = t_num[0] - (mids[0] - t_num[0])
    last_edge = t_num[-1] + (t_num[-1] - mids[-1])
    edges = cftime.num2date(
        np.concatenate([[first_edge], mids, [last_edge]]), units, cal
    )

    try:
        out = (
            data.groupby_bins("time", edges, labels=ref_time, right=False)
            .mean(skipna=True)
            .rename({"time_bins": "time"})
            .assign_coords(time=("time", ref_time))
        )
    except Exception:
        if fallback_freq is None:
            raise
        out = data.resample(time=fallback_freq).mean(skipna=True)
    return out


def find_time_dimension(data):
    """
    Find the time dimension in the dataset.

    Different climate models use different names for time dimensions
    """
    possible_time_dims = ["time", "time_counter", "t", "T", "Time"]

    for dim in data.dims:
        if dim in possible_time_dims:
            return dim

    # If no exact match, look for anything with 'time' in the name
    for dim in data.dims:
        if "time" in dim.lower():
            return dim

    return None


# def overwrite_time_coordinate(data, ref_data):
#     """
#     Overwrite the time coordinate of `data` with that of `ref_data`.
#     """
#     time_dim = find_time_dimension(data)
#     ref_time_dim = find_time_dimension(ref_data)
#     if time_dim and ref_time_dim:
#         # Only overwrite if lengths match
#         if data.sizes[time_dim] == ref_data.sizes[ref_time_dim]:
#             return data.assign_coords({time_dim: ref_data[ref_time_dim].values})
#         else:
#             print(f"✗ Cannot overwrite time: length mismatch ({data.sizes[time_dim]}
# vs {ref_data.sizes[ref_time_dim]})")
#     return data

# def align_variables_by_overwriting_time(loaded_vars, ref_var):
#     """
#     Overwrites the time coordinate of all variables in loaded_vars with ref_var.
#     """
#     ref_data = loaded_vars[ref_var]
#     aligned_vars = {}
#     for var, data in loaded_vars.items():
#         print(data.time.values)
#         print(ref_data.time.values)
#         aligned_vars[var] = overwrite_time_coordinate(data, ref_data)
#     return aligned_vars

# def resample_to_reference_bins(
#     data, ref_time, calendar, units, time_dim="time_counter", fallback_freq=None
# ):
#     """
#     Resample data to bins defined by ref_time (using cftime logic).
#     If this fails, fall back to xarray's .resample() with fallback_freq.
#     """
#     import cftime
#     import numpy as np
#     import xarray as xr

#     try:
#         numeric = cftime.date2num(ref_time, units, calendar)
#         if len(numeric) < 2:
#             raise ValueError("Reference time axis too short for binning.")

#         midpoints = (numeric[:-1] + numeric[1:]) / 2
#         first_edge = numeric[0] - (midpoints[0] - numeric[0])
#         last_edge = numeric[-1] + (numeric[-1] - midpoints[-1])
#         bin_edges_numeric = np.concatenate([[first_edge], midpoints, [last_edge]])
#         bin_edges = cftime.num2date(bin_edges_numeric, units, calendar)

#         resampled = data.groupby_bins(
#             time_dim, bin_edges, labels=ref_time, right=False
#         ).mean()
#         resampled = resampled.rename({f"{time_dim}_bins": time_dim})
#         resampled = resampled.assign_coords({time_dim: ref_time})
#         return resampled

#     except Exception as e:
#         print(f"[resample_to_reference_bins] Binning failed: {e}")
#         if fallback_freq is not None:
#             print(f"[resample_to_reference_bins] Falling back to \
# .resample({fallback_freq})")
#             resampled = data.resample({time_dim: fallback_freq}).mean()
#             return resampled
#         else:
#             raise RuntimeError(
#                 "Resampling to reference bins failed and no fallback_freq provided."
#             )
