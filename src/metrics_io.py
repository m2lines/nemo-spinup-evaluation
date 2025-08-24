"""
Provide functions for writing metric results to a CSV file.

Includes handling both time series and static metrics using xarray DataArrays.
"""

import csv
from collections import defaultdict

import xarray as xr


def write_metric_results(results, output_filepath):
    """
    Output the results of the metrics to a CSV file with time_counter as index.

    Parameters
    ----------
    results : dict
        The results of the metrics (some can be time series).
    output_filepath : str
        The path to the output CSV file.
    """
    indexed_data = defaultdict(dict)
    time_labels = []
    static_metrics = {}

    # Find any one DataArray with time_counter to use its time index
    time_array = None
    for result in results.values():
        if isinstance(result, xr.DataArray) and "time_counter" in result.dims:
            time_array = result["time_counter"]
            break

    if time_array is not None:
        time_labels = [t.values for t in time_array]

    for name, result in results.items():
        if isinstance(result, xr.DataArray) and "time_counter" in result.dims:
            for i, val in enumerate(result.values):
                indexed_data[i][name] = f"{val:.6f}"
        else:
            val = result.item() if hasattr(result, "item") else result
            static_metrics[name] = f"{val:.6f}" if isinstance(val, float) else str(val)

    # Define header
    header = [
        "timestamp",
        "check_density_from_file",
        "check_density_computed",
        "temperature_500m_30NS_metric",
        "temperature_BWbox_metric",
        "temperature_DWbox_metric",
        "ACC_Drake_metric",
        "ACC_Drake_metric_2",
        "NASTG_BSF_max",
    ]

    with open(output_filepath, mode="w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)

        for i in sorted(indexed_data.keys()):
            row = [time_labels[i] if i < len(time_labels) else ""]
            for metric in header[1:]:
                row.append(indexed_data[i].get(metric, static_metrics.get(metric, "")))
            writer.writerow(row)

    print(f"\n Successfully wrote metrics to '{output_filepath}'")
