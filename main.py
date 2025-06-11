"""Evaluate grid files using physical metrics."""

import argparse
import csv
import os
import sys
import glob
from collections import defaultdict
import yaml
import pandas as pd

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


def load_output_files(setup: dict, path: str) -> dict[str, xr.Dataset]:
    """Load output files based on the setup configuration."""
    output_data = {}
    for field, filename in setup.items():
        # TODO: check if the file exists
        if not os.path.exists(os.path.join(path, filename)):
            raise FileNotFoundError(f"File {filename} not found in path {path}.")
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


# these could be general functions for loading any mesh_mask
# could extract the requried vars from the setup.yaml or VARIABLE_ALIASES file.
def load_mesh_mask(path: str) -> xr.Dataset:
    """Load the NEMO mesh mask file."""
    ds = xr.open_dataset(path)
    # required_vars = ["tmask", "e1t", "e2t", "e3t"]  # TODO: update with what you need
    # missing = [var for var in required_vars if var not in ds.variables]
    # if missing:
    #     raise ValueError(f"Mesh mask file {path} is missing required variables: {missing}")
    return ds


def load_dino_outputs(
    mode, path, setup: dict, VARIABLE_ALIASES: dict
) -> dict[str, xr.Dataset]:
    """Load DINO model inputs based on YAML configuration.

    Parameters
    ----------
    config (dict)
        Parsed YAML config containing mesh_mask, restart_files,
    output_files, etc.

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
        raise ValueError("Mesh mask file is required in the setup configuration.")
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
            data[key] = standardize_variables(dataset, VARIABLE_ALIASES)
        elif isinstance(dataset, dict):  # e.g., output files
            data[key] = {
                k: standardize_variables(v, VARIABLE_ALIASES)
                for k, v in dataset.items()
            }

    # print(data)
    return data


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
                d["temperature"], d["salinity"], get_depth(data, mask), mask["tmask"]
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


def apply_metrics_output(
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
        description="Compute climate model diagnostics from restart or output files."
    )
    parser.add_argument(
        "--sim-path",
        type=str,
        help="Path to the NEMO simulation directory",
    )
    parser.add_argument(
        "--ref-path",
        type=str,
        help="Path to the reference NEMO simulation directory",
    )
    parser.add_argument(
        "--mode",
        type=str,
        help="Process mode: 'output' for output files, 'restart' for restart files,\
         'both' for both",
        choices=["output", "restart", "both"],
        default="both",
    )
    parser.add_argument(
        "--results",
        type=str,
        help="Path to save output metric values (default: metrics_results.csv)",
    )

    args = parser.parse_args()

    # Ensure at least one of the inputs is provided
    if not args.sim_path and not args.ref_path:
        print("Error: You must give a path to a directory.")
        parser.print_help()
        sys.exit(1)

    # Collect files paths based on arguments provided

    data_files = []
    # Logic to collect data files based on the provided paths
    # read in contents of yaml file DINO-sim.yaml which maps
    # variables to grid files.

    if args.sim_path:
        # Load DINO-sim.yaml to get expected DINO files
        yaml_file = "DINO-setup.yaml"
        if not os.path.exists(yaml_file):
            print(f"Error: The file {yaml_file} does not exist.")
            sys.exit(1)
        with open(yaml_file, "r") as file:
            dino_setup = yaml.safe_load(file)
        # Collect files based on the yaml mapping

    data = load_dino_outputs(args.mode, args.sim_path, dino_setup, VARIABLE_ALIASES)

    # If ref_sim_path is provided, load reference outputs
    comparison = args.ref_path is not None

    if comparison:
        data_ref = load_dino_outputs(
            args.mode, args.ref_path, dino_setup, VARIABLE_ALIASES
        )

    if args.mode in ["restart", "both"]:
        results = apply_metrics_restart(data["restart"], data["mesh_mask"])
        if comparison:
            ref_results = apply_metrics_restart(
                data_ref["restart"], data_ref["mesh_mask"]
            )

    if args.mode in ["output", "both"]:
        # For output mode, we expect grid_T, grid_T_sampled, grid_U, grid_V
        grid_output = data["output"]
        # use magic operator to pass to function.
        results = apply_metrics_output(
            grid_output["T"],
            grid_output["T_sampled"],
            grid_output["u"],
            grid_output["v"],
            data["restart"],
            data["mesh_mask"],
        )

        if comparison:
            grid_output_ref = data_ref["output"]
            results_ref = apply_metrics_output(
                grid_output_ref["T"],
                grid_output_ref["T_sampled"],
                grid_output_ref["u"],
                grid_output_ref["v"],
                data_ref["restart"],
                data_ref["mesh_mask"],
            )
        # write_metric_results(results, "results/metrics_results_grid.csv")

        results_ref = {f"ref_{key}": value for key, value in results_ref.items()}

        # compute differences between results and reference results
        results_diff = {}
        results_stats = {}
        for key in results:
            # compute MAE and RMSE using xarrays
            results_diff[f"diff_{key}_ae"] = results[key] - results_ref[f"ref_{key}"]
            results_stats[f"diff_{key}_mae"] = (
                results[key].mean() - results_ref.get(f"ref_{key}", 0).mean()
            )
            results_stats[f"diff_{key}_rmse"] = (
                (results[key] - results_ref.get(f"ref_{key}", 0)) ** 2
            ).mean() ** 0.5

        # now append all results together
        results.update(results_ref)
        results.update(results_diff)

        time_array = grid_output["T"].get("time_counter", None)
        # convert results to dataframe
        results_df = pd.DataFrame(results)
        # add time_counter as index if available
        if time_array is not None:
            results_df["timestamp"] = time_array.values
            results_df.set_index("timestamp", inplace=True)

        results_df.to_csv("results/metrics_results_grid_2.csv", index=True)


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
