"""Evaluate grid files using physical metrics."""

import argparse
import csv
import os
import sys
import glob
from collections import defaultdict
import yaml
import pandas as pd


from src.metrics import (
    ACC_Drake_metric,
    NASTG_BSF_max,
    check_density,
    temperature_500m_30NS_metric,
    temperature_BWbox_metric,
    temperature_DWbox_metric,
)
from src.utils import get_density, get_depth, sanitize_metrics_dict
from src.metrics_io import write_metrics_to_csv
from src.variable_aliases import VARIABLE_ALIASES

from src.loader import load_dino_outputs


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

    if restart is not None:
        metric_functions["check_density_computed"] = lambda: check_density(
            get_density(
                data_grid_T["temperature"],
                data_grid_T["salinity"],
                get_depth(restart, mask),
                mask["tmask"],
            )[0]
        )

    for name in sorted(metric_functions):
        try:
            result = metric_functions[name]()
            results[name] = result
        except (ValueError, TypeError, KeyError, AttributeError, ImportError) as e:
            results[name] = f"Error: {e}"
            print(f"Error in metric {name}: {e}")

    return results


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

    dino_setup = {}

    if args.sim_path:
        # Load DINO-setup.yaml to get expected DINO files
        # point to the config file in the configs directory
        config_file = os.path.join("configs", "DINO-setup.yaml")
        if not os.path.exists(config_file):
            print(f"Error: The file {config_file} does not exist.")
            sys.exit(1)
        with open(config_file, "r") as file:
            dino_setup = yaml.safe_load(file)
        # Collect files based on the yaml mapping

    data = load_dino_outputs(args.mode, args.sim_path, dino_setup, VARIABLE_ALIASES)

    # If ref_sim_path is provided, load reference outputs
    comparison = args.ref_path is not None

    data_ref = {}
    if comparison:
        data_ref = load_dino_outputs(
            args.mode, args.ref_path, dino_setup, VARIABLE_ALIASES
        )

    results_ref = {}
    if args.mode in ["restart", "both"]:
        print("Applying metrics to restart files...")
        results = sanitize_metrics_dict(
            apply_metrics_restart(data["restart"], data["mesh_mask"])
        )
        if comparison:
            results_ref = sanitize_metrics_dict(
                apply_metrics_restart(data_ref["restart"], data_ref["mesh_mask"])
            )

            results_ref = {f"ref_{key}": value for key, value in results_ref.items()}

        output_filepath = args.results or "metrics_results_restart.csv"
        write_metrics_to_csv(
            results,
            output_filepath,
            results_ref=results_ref if comparison else None,
            time_array=data["restart"].get("time_counter", None),
        )

    if args.mode in ["output", "both"]:
        print("Applying metrics to output files...")
        # For output mode, we expect grid_T, grid_T_sampled, grid_U, grid_V
        grid_output = data["output"]
        # use magic operator to pass to function.
        restart = None
        # data["restart"] might not exist therefore check if it is None
        if "restart" in data and data["restart"] is not None:
            restart = data["restart"]

        results = apply_metrics_output(
            grid_output["T"],
            grid_output["T_sampled"],
            grid_output["u"],
            grid_output["v"],
            restart,
            data["mesh_mask"],
        )
        results_ref = {}
        if comparison:
            grid_output_ref = data_ref["output"]

            if "restart" in data_ref and data_ref["restart"] is not None:
                restart = data_ref["restart"]

            results_ref = apply_metrics_output(
                grid_output_ref["T"],
                grid_output_ref["T_sampled"],
                grid_output_ref["u"],
                grid_output_ref["v"],
                restart,
                data_ref["mesh_mask"],
            )
            # write_metric_results(results, "results/metrics_results_grid.csv")
            # rename the results_ref keys to start with ref_ - this avoids collision
            results_ref = {f"ref_{key}": value for key, value in results_ref.items()}

        from src.metrics_io import write_metrics_to_csv

        output_filepath = args.results or "metrics_results.csv"
        write_metrics_to_csv(
            results,
            output_filepath,
            results_ref=results_ref if comparison else None,
            time_array=grid_output["T"].get("time_counter", None),
        )


##################################
### Command to run ####
##################################
# python3 main.py --sim data/DINO_EXP00/
#                 --ref data/DINO_EXP00_ref/
#                 --mode {output, restart, both}
#                 --results metrics_results.csv
