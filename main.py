"""Run metrics on data file."""

import os

import xarray as xr

# from utils import get_density, get_depth
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
    elif "grid" in filename:
        print("The provided file is a 'grid' file\n")

    data = xr.open_dataset(datafilepath, decode_cf=False).rename(
        {"deptht": "depth", "y": "nav_lat", "x": "nav_lon"}
    )
    print(f"Successfully loaded dataset from {filename}")
    mask = xr.open_dataset(maskfilepath).rename(
        {"nav_lev": "depth", "y": "nav_lat", "x": "nav_lon"}
    )
    print(f"Successfully loaded mesh mask for {filename}")

    return data, mask


def apply_metrics(data, mask):
    """
    Apply metrics to the data and mask.

    Parameters
    ----------
    data  : xarray.Dataset
        The dataset containing ocean model variables.
    mask  : (xarray.Dataset)
        The dataset containing mask variables.

    Returns
    -------
        dict: A dictionary containing the results of the metrics.
    """
    metrics_dict = {
        "check_density": check_density,
        "temperature_500m_30NS_metric": temperature_500m_30NS_metric,
        "temperature_BWbox_metric": temperature_BWbox_metric,
        "temperature_DWbox_metric": temperature_DWbox_metric,
        "ACC_Drake_metric": ACC_Drake_metric,
        "NASTG_BSF_max": NASTG_BSF_max,
    }

    results = {}

    for name, func in metrics_dict.items():
        try:
            results[name] = func(data, mask)
        except Exception as e:
            results[name] = f"Error: {e!s}"

    print(f"Successfully ran the metric `{name}` on the given file")


## Example
# filepath = "/home/sg2147/nc_files/nc_files/DINO_00576000_restart.nc"
filepath = "data/nemo-raw/DINO_1y_grid_T.nc"
maskfile = "data/nemo-raw/mesh_mask.nc"
data, mask = read_data(filepath, maskfile)
results = apply_metrics(data, mask)
