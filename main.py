"""Evaluate grid files using physical metrics."""

import argparse
import csv
import os
import sys
from collections import defaultdict

import xarray as xr

from src.metrics import (
    ACC_Drake_metric,
    NASTG_BSF_max,
    check_density,
    temperature_500m_30NS_metric,
    temperature_BWbox_metric,
    temperature_DWbox_metric,
)
from src.utils import get_density, get_depth
from variable_aliases import VARIABLE_ALIASES, standardize_variables


def read_data(datafilepath, maskfilepath, variable_dict):
    """Read the data and mask files and standardize variable names."""
    filename = os.path.basename(datafilepath)

    print(f"The provided file is: {filename}\n")

    # # Determine whether to decode CF conventions
    # decode_cf = not any(
    #     key in filename for key in ["grid_T", "grid_T_sampled", "grid_U", "grid_V"]
    # )

    # Open the dataset
    data = xr.open_dataset(datafilepath)

    data = standardize_variables(data, variable_dict)

    # Open and standardize the mesh mask
    mesh_mask = xr.open_dataset(maskfilepath)
    mesh_mask = standardize_variables(mesh_mask, variable_dict)

    return data, mesh_mask


def apply_metrics_restart(data, mask):
    """
    Apply metrics to the dataset and mask.

    Parameters
    ----------
    data : xarray.Dataset
        The standardized dataset containing ocean model variables.
    mask : xarray.Dataset
        The standardized mask dataset.

    Returns
    -------
    dict
        A dictionary containing the results of the metrics.
    """
    results = {}
    metric_functions = {
        "check_density_from_file": lambda d: check_density(d["density"][0]),
        "check_density_computed": lambda d: check_density(
            get_density(
                d["temperature"], d["salinity"], get_depth(restart, mask), mask["tmask"]
            )[0]
        ),
        "temperature_500m_30NS_metric": lambda d: temperature_500m_30NS_metric(
            d["temperature"][0], mask
        ),
        "temperature_BWbox_metric": lambda d: temperature_BWbox_metric(
            d["temperature"][0], mask
        ),
        "temperature_DWbox_metric": lambda d: temperature_DWbox_metric(
            d["velocity_zonal"][0], mask
        ),
        "ACC_Drake_metric": lambda d: ACC_Drake_metric(d["velocity_zonal"][0], mask),
        "NASTG_BSF_max": lambda d: NASTG_BSF_max(
            d["velocity_meridional"][0], d["ssh"], mask
        ),
    }

    for name, func in metric_functions.items():
        try:
            results[name] = func(data)
        except (ValueError, TypeError, KeyError, AttributeError, ImportError) as e:
            results[name] = f"Error: {e}"
            print(f"Error in metric {name}: {e}")

    return results


def apply_metrics_grid(
    data_grid_T, data_grid_T_sampled, data_grid_U, data_grid_V, restart, mask
):
    """
    Apply metrics to ocean model grid datasets.

    Parameters
    ----------
    data_grid_T : xarray.Dataset
        The standardized T-grid dataset containing temperature, salinity,
        and density variables.
    data_grid_T_sampled : xarray.Dataset
        The sampled T-grid dataset containing SSH (sea surface height) data.
    data_grid_U : xarray.Dataset
        The standardized U-grid dataset containing zonal velocity data.
    data_grid_V : xarray.Dataset
        The standardized V-grid dataset containing meridional velocity data.
    restart : xarray.Dataset
        The restart dataset containing model state information.
    mask : xarray.Dataset
        The standardized mask dataset containing grid masks.

    Returns
    -------
    dict
        A dictionary containing the results of the applied metrics, with
        metric names as keys and computed values or error messages as values.
    """
    results = {}
    metric_functions = {
        "check_density_from_file": lambda: check_density(data_grid_T["density"]),
        "check_density_computed": lambda: check_density(
            get_density(
                data_grid_T["temperature"],
                data_grid_T["salinity"],
                get_depth(restart, mask),
                mask["tmask"],
            )[0]
        ),
        "temperature_500m_30NS_metric": lambda: temperature_500m_30NS_metric(
            data_grid_T["temperature"], mask
        ),
        "temperature_BWbox_metric": lambda: temperature_BWbox_metric(
            data_grid_T["temperature"], mask
        ),
        "temperature_DWbox_metric": lambda: temperature_DWbox_metric(
            data_grid_U["velocity_zonal"], mask
        ),
        "ACC_Drake_metric": lambda: ACC_Drake_metric(
            data_grid_U["velocity_zonal"], mask
        ),
        "NASTG_BSF_max": lambda: NASTG_BSF_max(
            data_grid_V["velocity_meridional"], data_grid_T_sampled["ssh"], mask
        ),
    }

    for name, func in metric_functions.items():
        try:
            result = func()
            results[name] = result
            # if hasattr(result, "plot"):
            #     result.plot()
        except (ValueError, TypeError, KeyError, AttributeError, ImportError) as e:
            results[name] = f"Error: {e}"
            print(f"Error in metric {name}: {e}")

    return results


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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute climate model diagnostics from a restart file and \
        mesh_mask."
    )

    parser.add_argument(
        "--restart",
        type=str,
        help="Path to the model restart file (e.g., restart.nc)",
    )
    parser.add_argument(
        "--grid_T",
        type=str,
        help="Path to the NEMO grid_T file (e.g., _grid_T.nc)",
    )
    parser.add_argument(
        "--grid_T_sampled",
        type=str,
        help="Path to save output metric values (e.g., _grid_T_sampled.nc)",
    )
    parser.add_argument(
        "--grid_U",
        type=str,
        help="Path to the NEMO grid_U file (e.g., _grid_U.nc)",
    )
    parser.add_argument(
        "--grid_V",
        type=str,
        help="Path to the NEMO grid_V file (e.g., _grid_V.nc)",
    )
    parser.add_argument(
        "--mesh-mask",
        type=str,
        required=True,
        help="Path to the NEMO mesh mask file (e.g., mesh_mask.nc)",
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Path to save output metric values (default: metrics_results.csv)",
    )

    args = parser.parse_args()

    # Ensure at least one of the inputs is provided
    if not args.restart and not (
        args.grid_T and args.grid_T_sampled and args.grid_U and args.grid_V
    ):
        print("Error: You must give path to either restart file or grid files.")
        parser.print_help()
        sys.exit(1)

    # Collect files paths based on arguments provided
    data_files = []
    if args.restart:
        data_files.append(args.restart)
    if args.grid_T:
        data_files.append(args.grid_T)
    if args.grid_T_sampled:
        data_files.append(args.grid_T_sampled)
    if args.grid_U:
        data_files.append(args.grid_U)
    if args.grid_V:
        data_files.append(args.grid_V)

    if args.restart and not args.grid_T:
        restart, mesh_mask = read_data(args.restart, args.mesh_mask, VARIABLE_ALIASES)
        results = apply_metrics_restart(restart, mesh_mask)
        write_metric_results(results, "results/metrics_results_restart.csv")
    else:
        restart, mesh_mask = read_data(args.restart, args.mesh_mask, VARIABLE_ALIASES)
        data_grid_T, mesh_mask = read_data(
            args.grid_T, args.mesh_mask, VARIABLE_ALIASES
        )
        data_grid_T_sampled, mesh_mask = read_data(
            args.grid_T_sampled, args.mesh_mask, VARIABLE_ALIASES
        )
        data_grid_T_sampled["time_counter"] = data_grid_T["time_counter"]
        data_grid_U, mesh_mask = read_data(
            args.grid_U, args.mesh_mask, VARIABLE_ALIASES
        )
        data_grid_V, mesh_mask = read_data(
            args.grid_V, args.mesh_mask, VARIABLE_ALIASES
        )
        results = apply_metrics_grid(
            data_grid_T,
            data_grid_T_sampled,
            data_grid_U,
            data_grid_V,
            restart,
            mesh_mask,
        )
        write_metric_results(results, "results/metrics_results_grid.csv")


##################################
### Command to run grid files ####
##################################
# python3 main.py --restart ../nc_files/nc_files/DINO_00576000_restart.nc
#                 --grid_T ../nc_files/nc_files/DINO_1y_grid_T.nc
#                 --grid_T_sampled ../nc_files/nc_files/DINO_1m_To_1y_grid_T.nc
#                 --grid_U ../nc_files/nc_files/DINO_1y_grid_U.nc
#                 --grid_V ../nc_files/nc_files/DINO_1y_grid_V.nc
#                 --mesh-mask ../nc_files/nc_files/mesh_mask.nc

##################################
### Command to run restart files ####
##################################
# python3 main.py --restart ../nc_files/nc_files/DINO_00576000_restart.nc
#                 --mesh-mask ../nc_files/nc_files/mesh_mask.nc
