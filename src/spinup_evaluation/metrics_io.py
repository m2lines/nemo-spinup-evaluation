"""
Provide functions for writing metric results to a CSV file.

Includes handling both time series and static metrics using xarray DataArrays.
"""
# spinup_evaluation/metrics_io.py

from __future__ import annotations

import os
from typing import Any, Mapping

import pandas as pd
import xarray as xr  # type: ignore


def write_metric_results(results: Mapping[str, Any], out_path: str) -> None:
    """
    Write metric results to CSV files.

    Writes two CSVs:
      - {out_path}_timeseries.csv : timestamp + sim/ref/diff/diff_ae series
      - {out_path}_summary.csv    : metric, mae, rmse (from diff_*_mae / diff_*_rmse)


    Parameters
    ----------
    results : Mapping[str, Any]
        The metric results to write.
    out_path : str
        The base output path (without file extension).
    """
    # compute diff (on timeseries if present)
    ts_cols = {}
    timestamp = None
    for k, v in results.items():
        if isinstance(v, xr.DataArray) and "time_counter" in v.dims:
            if timestamp is None:  # pick canonical time axis
                timestamp = [str(t) for t in v["time_counter"].values]
            ts_cols[k] = v.values
        elif isinstance(v, xr.DataArray) and v.ndim == 0:
            ts_cols[k] = [v.item()]
        elif not k.endswith(("_mae", "_rmse")):
            ts_cols[k] = [v]

    if timestamp is None:  # restart-only case
        ts_df = pd.DataFrame(ts_cols)
        ts_df.insert(0, "timestamp", [None])
    else:
        ts_df = pd.DataFrame(ts_cols, index=timestamp).reset_index()
        ts_df = ts_df.rename(columns={"index": "timestamp"})

    # compute summary statistics
    rows = []
    for k, v in results.items():
        if k.startswith("diff_") and k.endswith(("_mae", "_rmse")):
            if k.endswith("_mae"):
                metric = k[len("diff_") : -4]  # strip "_mae"
                stat = "mae"
            elif k.endswith("_rmse"):
                metric = k[len("diff_") : -5]  # strip "_rmse"
                stat = "rmse"
            else:
                continue
            rows.append((metric, stat, v))
    summary_df = (
        pd.DataFrame(rows, columns=["metric", "stat", "value"])
        .pivot(index="metric", columns="stat", values="value")
        .reset_index()
    )

    ts_path = f"{out_path}_timeseries.csv"
    sm_path = f"{out_path}_summary.csv"
    os.makedirs(os.path.dirname(os.path.abspath(ts_path)), exist_ok=True)
    print(f"Writing timeseries metrics to '{ts_path}'")
    ts_df.to_csv(ts_path, index=False)
    print(f"Writing summary metrics to '{sm_path}'")
    summary_df.to_csv(sm_path, index=False)
