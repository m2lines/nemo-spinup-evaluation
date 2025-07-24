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
from src.utils import get_density, get_depth, sanitize_metrics_dict
from src.metrics_io import write_metrics_to_csv
from src.variable_aliases import VARIABLE_ALIASES, standardize_variables
from src.loader import (
    load_dino_outputs,
    load_mesh_mask,
    load_output_files,
)  # noqa: E501


metric_functions = {
    "check_density_from_file": lambda ctx: check_density(ctx["density"]),
    "check_density_computed": lambda ctx: (
        check_density(
            get_density(
                ctx["temperature"],
                ctx["salinity"],
                get_depth(ctx["restart"], ctx["mask"]),
                ctx["mask"]["tmask"],
            )[0]
        )
        if ctx.get("restart") is not None
        else "Skipped: restart missing"
    ),
    "temperature_500m_30NS_metric": lambda ctx: temperature_500m_30NS_metric(
        ctx["temperature"], ctx["mask"]
    ),
    "temperature_BWbox_metric": lambda ctx: temperature_BWbox_metric(
        ctx["temperature"], ctx["mask"]
    ),
    "temperature_DWbox_metric": lambda ctx: temperature_DWbox_metric(
        ctx["velocity_zonal"], ctx["mask"]
    ),
    "ACC_Drake_metric": lambda ctx: ACC_Drake_metric(
        ctx["velocity_zonal"], ctx["mask"]
    ),
    "NASTG_BSF_max": lambda ctx: NASTG_BSF_max(
        ctx["velocity_meridional"], ctx["ssh"], ctx["mask"]
    ),
}


def apply_metrics(context):
    results = {}

    for name in sorted(metric_functions):
        try:
            results[name] = metric_functions[name](context)
        except Exception as e:
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
    mode = args.mode

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

    # RESTART MODE
    if mode in ["restart", "both"]:
        print("Applying metrics to restart files...")
        restart_context = build_restart_context(data)

        results = sanitize_metrics_dict(apply_metrics(restart_context))
        if comparison:
            restart_context_ref = build_restart_context(data_ref)
            results_ref = sanitize_metrics_dict(apply_metrics(restart_context))
            results_ref = {f"ref_{key}": value for key, value in results_ref.items()}

        output_filepath = args.results or "metrics_results_restart.csv"
        write_metrics_to_csv(
            results,
            output_filepath,
            results_ref=results_ref if comparison else None,
            time_array=data["restart"].get("time_counter", None),
        )

    # OUTPUT MODE
    if mode in ["output", "both"]:
        print("Applying metrics to output files...")
        # For output mode, we expect grid_T, grid_T_sampled, grid_U, grid_V
        # use magic operator to pass to function.

        output_context = build_output_context(data)
        results = apply_metrics(output_context)
        results_ref = {}
        if comparison:
            output_context_ref = build_output_context(data_ref)
            results_ref = apply_metrics(output_context_ref)
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
