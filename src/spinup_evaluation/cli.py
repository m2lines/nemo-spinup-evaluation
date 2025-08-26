"""Evaluate grid files using physical metrics."""

import argparse
import os
import sys

import yaml

from spinup_evaluation.loader import load_dino_data
from spinup_evaluation.metrics_io import write_metric_results
from spinup_evaluation.utils import get_density, get_depth

try:
    # Try relative import first (when running as part of package)
    from .metrics import (
        ACC_Drake_metric,
        ACC_Drake_metric_2,
        NASTG_BSF_max,
        check_density,
        temperature_500m_30NS_metric,
        temperature_BWbox_metric,
        temperature_DWbox_metric,
    )
except ImportError:
    # Fallback to absolute import (for development/testing)
    from spinup_evaluation.metrics import (
        ACC_Drake_metric,
        ACC_Drake_metric_2,
        NASTG_BSF_max,
        check_density,
        temperature_500m_30NS_metric,
        temperature_BWbox_metric,
        temperature_DWbox_metric,
    )


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
            d["velocity_u"][0], mask
        ),
        "ACC_Drake_metric": lambda d: ACC_Drake_metric(d["velocity_u"][0], mask),
        "ACC_Drake_metric_2": lambda d: ACC_Drake_metric_2(
            d["velocity_u"][0], d["ssh"][0], mask
        ),
        "NASTG_BSF_max": lambda d: NASTG_BSF_max(d["velocity_v"][0], d["ssh"], mask),
    }

    for name, func in metric_functions.items():
        try:
            results[name] = func(data)
        except (ValueError, TypeError, KeyError, AttributeError, ImportError) as e:
            results[name] = f"Error: {e}"
            print(f"Error in metric {name}: {e}")

    return results


def apply_metrics_output(grid_output, restart, mask):
    """
    Apply metrics to ocean model grid datasets.

    Parameters
    ----------
    grid_output : xarray.Dataset
        The standardized grid output dataset containing model variables.
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
        "check_density_from_file": lambda: check_density(
            grid_output["density"]
        ),  # Updated
        "temperature_500m_30NS_metric": lambda: temperature_500m_30NS_metric(
            grid_output["temperature"], mask
        ),
        "temperature_BWbox_metric": lambda: temperature_BWbox_metric(
            grid_output["temperature"], mask
        ),
        "temperature_DWbox_metric": lambda: temperature_DWbox_metric(
            grid_output["velocity_u"], mask
        ),
        "ACC_Drake_metric": lambda: ACC_Drake_metric(grid_output["velocity_u"], mask),
        "ACC_Drake_metric_2": lambda: ACC_Drake_metric_2(
            grid_output["velocity_u"], grid_output["ssh"], mask
        ),
        "NASTG_BSF_max": lambda: NASTG_BSF_max(
            grid_output["velocity_v"], grid_output["ssh"], mask
        ),
    }

    if restart is not None:
        metric_functions["check_density_computed"] = lambda: check_density(
            get_density(
                grid_output["temperature"],
                grid_output["salinity"],
                get_depth(restart, mask),
                mask["tmask"],
            )[0]
        )

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


def main():
    """
    Command-line interface for computing climate model diagnostics.

    Parses arguments, loads configuration and data, applies metrics, and writes results.
    """
    parser = argparse.ArgumentParser(
        description="Compute climate model diagnostics from restart or output files."
    )
    parser.add_argument(
        "--sim-path",
        type=str,
        help="Path to the NEMO simulation directory",
    )
    parser.add_argument(
        "--ref-sim-path",
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
        "--results-dir",
        type=str,
        help="Path to save output metric values",
        default="results/",
    )
    parser.add_argument(
        "--result-file-prefix",
        type=str,
        help="Prefix to save the final result file",
        default="metrics_results",
    )
    args = parser.parse_args()

    # Ensure at least one of the inputs is provided
    if not args.sim_path and not args.ref_sim_path:
        print("Error: You must give a path to a directory.")
        parser.print_help()
        sys.exit(1)

    # Logic to collect data files based on the provided paths
    # read in contents of yaml file DINO-setup.yaml which maps
    # variables to grid files.

    if args.sim_path:
        # Load DINO-setup.yaml to get expected DINO files
        # point to the config file in the configs directory
        config_file = os.path.join("configs", "DINO-setup.yaml")
        if not os.path.exists(config_file):
            print(f"Error: The file {config_file} does not exist.")
            sys.exit(1)
        with open(config_file, "r") as file:
            dino_setup = yaml.safe_load(file)

    # Ensure the results directory exists
    results_dir = os.path.dirname(args.results_dir)
    if results_dir and not os.path.exists(results_dir):
        os.makedirs(results_dir, exist_ok=True)
    prefix = args.result_file_prefix

    data = load_dino_data(args.mode, args.sim_path, dino_setup, do_standardise=True)

    if args.mode in ["restart", "both"]:
        results = apply_metrics_restart(data["restart"], data["mesh_mask"])
        write_metric_results(results, f"{results_dir}/{prefix}_restart.csv")
    if args.mode in ["output", "both"]:
        grid_output = data["grid"]

        restart = None
        # data["restart"] might not exist therefore check if it is None
        if "restart" in data and data["restart"] is not None:
            restart = data["restart"]

        # use magic operator to pass to function.
        results = apply_metrics_output(
            grid_output,
            restart,
            data["mesh_mask"],
        )
        write_metric_results(results, f"{results_dir}/{prefix}_grid.csv")


##################################
### Running spinup-evaluation ####
##################################
# python3 main.py --sim-path /path/to/simulation_directory
#                 --mode output
#
# The code now uses a configuration YAML (e.g., configs/DINO-setup.yaml)
# to determine which files to load and which variables to map.

if __name__ == "__main__":
    main()
