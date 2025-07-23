import pandas as pd
import xarray as xr
from collections import defaultdict
import numpy as np

import pandas as pd
import xarray as xr
import numpy as np
from collections import defaultdict


def write_metrics_to_csv(results, output_filepath, results_ref=None, time_array=None):
    """
    Write evaluation metrics to two CSV files:
      1. Time-varying metrics indexed by timestamp
      2. Scalar summary metrics (MAE and RMSE only)

    Parameters
    ----------
    results : dict
        Dictionary of predicted metrics (xarray.DataArray or scalar).
    output_filepath : str
        Path to time-varying CSV file. Summary will be saved as *_summary.csv.
    results_ref : dict, optional
        Dictionary of reference metrics to compare against.
    time_array : xarray.DataArray or None, optional
        Optional time index to use; inferred from results if not provided.
    """
    results_all = results.copy()
    static_metrics = {}

    # Compute AE and summary stats (MAE, RMSE)
    if results_ref:
        for key in results:
            ref_key = f"ref_{key}"
            if ref_key in results_ref:
                pred = results[key]
                ref = results_ref[ref_key]
                try:
                    results_all[f"diff_{key}_ae"] = pred - ref
                    static_metrics[f"diff_{key}_mae"] = abs(pred - ref).mean().item()
                    static_metrics[f"diff_{key}_rmse"] = (
                        ((pred - ref) ** 2).mean() ** 0.5
                    ).item()
                except Exception as e:
                    print(f"Could not compute diffs for '{key}': {e}")

        results_all.update(results_ref)
        # sort so that diff_* metrics come last
        results_all = dict(
            sorted(results_all.items(), key=lambda x: x[0].startswith("diff_"))
        )

    # Try to get time array if not provided
    if time_array is None:
        for val in results_all.values():
            if isinstance(val, xr.DataArray) and "time_counter" in val.dims:
                time_array = val["time_counter"]
                break

    # Separate time-varying metrics only
    df_data = defaultdict(dict)
    time_vals = time_array.values if time_array is not None else []

    for key, val in results_all.items():
        if isinstance(val, xr.DataArray) and "time_counter" in val.dims:
            # breakpoint()
            for i, v in enumerate(val.values):
                df_data[i][key] = v

    # Write time-varying metrics
    df = pd.DataFrame.from_dict(df_data, orient="index")
    if len(time_vals) and len(df) > 0:
        df["timestamp"] = time_vals[: len(df)]
        df.set_index("timestamp", inplace=True)

    df.to_csv(output_filepath, index=True)
    print(f"Metrics written to '{output_filepath}'")

    # Write only MAE/RMSE summary
    if static_metrics:
        summary_df = pd.DataFrame([static_metrics])
        summary_path = output_filepath.replace(".csv", "_summary.csv")
        summary_df.to_csv(summary_path, index=False)
        print(f"Scalar summary metrics (MAE/RMSE) written to '{summary_path}'")
    else:
        print("No summary metrics written (MAE/RMSE not found)")
