"""Run metrics on data file."""

import argparse
import os

import xarray as xr

from src.metrics import (
    ACC_Drake_metric,
    NASTG_BSF_max,
    check_density,
    temperature_500m_30NS_metric,
    temperature_BWbox_metric,
    temperature_DWbox_metric,
)


def read_data(datafilepath, maskfilepath):
   """Read the data and mask files."""
   filename = os.path.basename(datafilepath)
   if "restart" in filename:
       print("The provided file is a 'restart' file\n")
       data = xr.open_dataset(datafilepath).rename(
           {"nav_lev":"depth","y":"nav_lat","x":"nav_lon"})
   elif "grid_T" in filename:
       print("The provided file is a 'grid_T' file\n")
       data = xr.open_dataset(datafilepath, decode_cf=False).rename(
       {"deptht": "depth", "y": "nav_lat", "x": "nav_lon"})
   elif "grid_U" in filename:
       print("The provided file is a 'grid_U' file\n")
       data=xr.open_dataset(datafilepath,decode_cf=False).rename(
           {"depthu":"depth","y":"nav_lat","x":"nav_lon"})
   elif "grid_V" in filename:
       print("The provided file is a 'grid_V' file\n")
       data=xr.open_dataset(datafilepath,decode_cf=False).rename(
           {"depthv":"depth","y":"nav_lat","x":"nav_lon"})


   print(f"Successfully loaded dataset from {filename}")
   mask = xr.open_dataset(maskfilepath).rename(
       {"nav_lev": "depth", "y": "nav_lat", "x": "nav_lon"}
   )
   print(f"Successfully loaded mesh mask for {filename}")


   return data, mask


def apply_metrics_restart(restart, mask):
    """
    Apply metrics to restart netcdf and mask.

    Parameters
    ----------
    restart  : xarray.Dataset
        The restart dataset containing ocean model variables.
    mask  : xarray.Dataset
        The dataset containing mask variables.

    Returns
    -------
        dict: A dictionary containing the results of the metrics.
    """
    metrics_fns = {
        "check_density": lambda restart: check_density(restart.rhop),
        "temperature_500m_30NS_metric": lambda restart: temperature_500m_30NS_metric(
            restart.tn, mask
        ),
        "temperature_BWbox_metric": lambda restart: temperature_BWbox_metric(
            restart.tn, mask
        ),
        "temperature_DWbox_metric": lambda restart: temperature_DWbox_metric(
            restart.un, mask
        ),
        "ACC_Drake_metric": lambda restart: ACC_Drake_metric(restart.un, mask),
        "NASTG_BSF_max": lambda restart: NASTG_BSF_max(restart.vn, restart.sshn, mask),
    }

    results = {}

    for name, func in metrics_fns.items():
        try:
            results[name] = func(restart)
            print(f"Successfully ran the metric `{name}` on the given file")
        except (ValueError, KeyError, TypeError) as e:
            results[name] = f"Error: {e!s}"
            print(f"Error running metric `{name}` on the given file")

    return results


def output_metrics(results, output_filepath):
    """
    Output the results of the metrics to a file.

    Parameters
    ----------
    results : dict
        The results of the metrics.
    output_filepath : str
        The path to the output file.
    """
    with open(output_filepath, "w") as f:
        for name, result in results.items():
            if isinstance(result.item(), (int, float)):
                f.write(f"{name}: {result.item():.6f}\n")
                print(f"{name}: {result.item():.6f}\n")
            else:
                f.write(f"{name}: {result}\n")
                print(f"{name}: {result}\n")
        print("Successfully wrote metrics to file")


if __name__ == "__main__":
    # Example use:
    # python main.py --restart ../../spinup-data/nemo-raw/DINO_00576000_restart.nc \
    # --mesh-mask ../../spinup-data/nemo-raw/mesh_mask.nc

    parser = argparse.ArgumentParser(
        description="Compute climate model diagnostics from a restart file and \
        mesh_mask."
    )

    parser.add_argument(
        "--restart",
        type=str,
        required=True,
        help="Path to the model restart file (e.g., restart.nc)",
    )
    parser.add_argument(
        "--mesh-mask",
        type=str,
        required=True,
        help="Path to the NEMO mesh mask file (e.g., mesh_mask.nc)",
    )
    parser.add_argument(
        "--grid_T",
        type=str,
        help="Path to the NEMO grid_T file (e.g., _grid_T.nc)",
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
        "--output",
        type=str,
        default="metrics_results.txt",
        help="Path to save output metric values (default: metrics_results.txt)",
    )

    args = parser.parse_args()
    
    # Ensure at least one of the inputs is provided
    if not args.restart and not args.grid:
        print("Error: You must specify at least one of --restart or --grid.")
        parser.print_help()
        sys.exit(1)

    restart_data, mesh_mask = read_data(args.restart, args.mesh_mask)
    results = apply_metrics_restart(restart_data, mesh_mask)
    output_filepath = args.output
    output_metrics(results, output_filepath)
