"""Run metrics on data file."""

import os

import xarray as xr
import numpy as np

# from utils import get_density, get_depth
from src.metrics import (
    ACC_Drake_metric,
    NASTG_BSF_max,
    check_density,
    temperature_500m_30NS_metric,
    temperature_BWbox_metric,
    temperature_DWbox_metric,
    wrap_temperature_numpy_as_xarray,
)


def read_data(datafilepath, maskfilepath):
    """Read the data and mask files."""
    filename = os.path.basename(datafilepath)
    if "restart" in filename:
        print("The provided file is a 'restart' file\n")
    elif "grid" in filename:
        print("The provided file is a 'grid' file\n")

    # this fails on restart files because deptht is not a variable

    if "grid" in filename:
        data = xr.open_dataset(datafilepath, decode_cf=False).rename(
            {"deptht": "depth", "y": "nav_lat", "x": "nav_lon"}
        )
    elif "restart" in filename:
        data = xr.open_dataset(datafilepath).rename(
            {"nav_lev": "depth", "y": "nav_lat", "x": "nav_lon"}
        )

    print(f"Successfully loaded dataset from {filename}")

    mask = xr.open_dataset(maskfilepath).rename(
        {"nav_lev": "depth", "y": "nav_lat", "x": "nav_lon"}
    )
    print(f"Successfully loaded mesh mask for {filename}")

    return data, mask


def read_prediction_data(datafilepath):

    # Read prediction data in dir : pred_so.npy, pred_thetao.npy, pred_zos.npy
    pred_so = np.load(os.path.join(datafilepath, "pred_so.npy"))
    pred_thetao = np.load(os.path.join(datafilepath, "pred_thetao.npy"))
    pred_zos = np.load(os.path.join(datafilepath, "pred_zos.npy"))

    return pred_so, pred_thetao, pred_zos


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



    # for name, func in metrics_dict.items():
    #     try:
    #         results[name] = func(data, mask)
    #     except Exception as e:
    #         results[name] = f"Error: {e!s}"

    # print(f"Successfully ran the metric `{name}` on the given file")


## Example
# filepath = "/home/sg2147/nc_files/nc_files/DINO_00576000_restart.nc"
# raw_path = "data/nemo-raw/DINO_1y_grid_T.nc"


# check that the restart density is the same as the density computed
# from the projected data (pred_so, pred_thetao, pred_zos)
# def check_density(restart, mask):
#     """
#     Check the density of the restart file and the projected data.
#
#     Parameters
#     ----------
#     restart  : xarray.Dataset
#         The dataset containing ocean model variables.
#     mask     : xarray.Dataset
#         The dataset containing mask variables.
#
#     Returns


    # print(pred_so.rhop)j
from src.utils import get_deptht, get_density


def apply_metrics_on_restart():
    # Read the restart file
    restart_path = "../../spinup-data/nemo-raw/DINO_00576000_restart.nc"
    maskfile = "../../spinup-data/nemo-raw/mesh_mask.nc"
    data, mask = read_data(restart_path, maskfile)
    print(check_density(data.rhop))
    
    return 

def check_density_equivalence_between_restart_files(restart_old, restart_new):
    """
    Check the density of the restart file and the projected data.

    Parameters
    ----------
    restart  : xarray.Dataset
        The dataset containing ocean model variables.
    mask     : xarray.Dataset
        The dataset containing mask variables.

    Returns
    -------
    bool: True if the densities are equal, False otherwise.
    """
    # Check if the two arrays are equal
    return np.allclose(restart_old.rhop, restart_new.rhop)


def load_all_data():
    prediction_path = "../../spinup-data/prediction-out/"
    pred_so, pred_thetao, pred_zos = read_prediction_data(prediction_path)

    maskfilepath = "../../spinup-data/nemo-raw/mesh_mask.nc"
    mask = xr.open_dataset(maskfilepath).rename(
    {"nav_lev": "depth", "y": "nav_lat", "x": "nav_lon"})

    restart_path = "../../spinup-data/restart/NEW_DINO_00576000_restart.nc"
    restart_path_old = "../../spinup-data/nemo-raw/DINO_00576000_restart.nc"

    maskfile = "../../spinup-data/nemo-raw/mesh_mask.nc"
    restart_data, mask = read_data(restart_path, maskfile)
    restart_data_old, mask = read_data(restart_path_old, maskfile)


    return pred_so, pred_thetao, pred_zos, mask, restart_data, restart_data_old

def apply_metrics_on_preds():
    """
    Apply metrics to the prediction data.

    It checks that the check_density function is giving the same results
    on the new restart file and the projected data (.npz).

    This means that we don't have to regenerate the restart file.

    Parameters
    ----------

    """

    pred_so, pred_thetao, pred_zos, mask, restart_data, restart_data_old = load_all_data()
    # rhop = get_density(pred_thetao[-1], pred_so[-1], pred_zos[-1], mask.tmask)
    # for curiosity bring in the restart file's deptht
    print(check_density_equivalence_between_restart_files(restart_data, restart_data_old))
    rhop_1 = get_density(pred_thetao[-1], pred_so[-1], get_deptht(restart_data, mask), mask.tmask)[0].fillna(0)
    
    rhop_2 = get_density(pred_thetao[-1], pred_so[-1], pred_zos[-1], mask.tmask)[0].fillna(0)
    print("density ", check_density(restart_data.rhop), check_density(restart_data_old.rhop), check_density(rhop_1), check_density(rhop_2))

    # I'm happy with the metrics for density
    # now temp- this is hard, we don't have functions like .sel as it it's an xarray, it's a numpy array. 
    # tpred = temperature_BWbox_metric_numpy(pred_thetao[-1], mask.depth, mask.tmask)

    thetao = wrap_temperature_numpy_as_xarray(pred_thetao[-1], mask)
    tpred = temperature_500m_30NS_metric(thetao, mask)
    tpred_r = temperature_500m_30NS_metric(restart_data.tn, mask)
    print(tpred, tpred_r)

    tpBW = temperature_BWbox_metric(thetao, mask)
    tpBW_r = temperature_BWbox_metric(restart_data.tn, mask)
    print(tpBW, tpBW_r)

    tpDW = temperature_DWbox_metric(thetao, mask)
    tpDW_r = temperature_DWbox_metric(restart_data.tn, mask)
    print(tpDW, tpDW_r)
    
    # ACC_Drake_metric
    uo = restart_data.uo.squeeze()
    acc = ACC_Drake_metric(uo, mask)

    return rhop_1

def apply_metrics_on_grid():
    # Read the grid file
    grid_path = "../../spinup-data/nemo-raw/DINO_1y_grid_T.nc"
    maskfile = "../../spinup-data/nemo-raw/mesh_mask.nc"
    data, mask = read_data(grid_path, maskfile)
    # check_density(data.rhop)
    return data.rhop

def check_equivalence(rhop, restart_rhop):
    """
    Check if the density from the restart file and the projected data are equal.

    Parameters
    ----------
    rhop : numpy.array
        The density of the projected data.
    restart_rhop : numpy.array
        The density of the restart file.

    Returns
    -------
    bool: True if the densities are equal, False otherwise.
    """
    # Check if the two arrays are equal
    return np.allclose(rhop, restart_rhop)


# check if ssh from the updated restart file is the same as the ssh from the pre_ssh.so

if __name__ == "__main__":
    # data, mask = read_data(raw_path, maskfile)
    # results = apply_metrics(data, mask)
    rhop_restart = apply_metrics_on_restart()

    rhop_pred = apply_metrics_on_preds()
    # equivalence = check_equivalence(rhop_pred, rhop_restart)
    # print(f"Are the densities equal? {equivalence}")
    # print(results)
