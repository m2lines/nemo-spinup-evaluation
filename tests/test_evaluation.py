from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from spinup_evaluation.cli import main


@pytest.fixture
def dummy_input_file(tmp_path):
    """
    Fixture to create a dummy NetCDF input file.
    """
    rng = np.random.default_rng()
    dummy_data = xr.Dataset(
        {"temperature": (("lat", "lon"), rng.random((10, 10)))},
        coords={
            "lat": np.linspace(-90, 90, 10),
            "lon": np.linspace(0, 360, 10, endpoint=False),
        },
    )
    input_file = tmp_path / "input.nc"
    dummy_data.to_netcdf(input_file)
    return input_file


@pytest.fixture
def real_input_dir():
    """
    Fixture to provide the path to the directory with real NetCDF input files.
    """
    return Path("tests/data/DINO_subsampled_A/")


@pytest.fixture
def output_file_path(tmp_path):
    """
    Fixture to provide an output file path.
    """
    return Path(tmp_path / "results")


def assert_csv_equal(path_a: str, path_b: str, float_tol: float = 1e-12) -> None:
    df_a = pd.read_csv(path_a)
    df_b = pd.read_csv(path_b)

    # same column sets
    assert set(df_a.columns) == set(df_b.columns), (
        f"Different columns: only in A={set(df_a.columns) - set(df_b.columns)}, "
        f"only in B={set(df_b.columns) - set(df_a.columns)}"
    )

    # reorder B's columns to match A
    df_b = df_b[df_a.columns]

    # same row count
    assert len(df_a) == len(df_b), f"Different row counts: {len(df_a)} vs {len(df_b)}"

    # numeric vs non-numeric comparison
    for col in df_a.columns:
        a, b = df_a[col], df_b[col]
        if col.lower() in {"timestamp", "time_counter", "time", "time_centered"}:
            # Require identical values (string compare is fine)
            assert a.equals(b), f"Time column {col} does not match."
            continue
        if not np.allclose(
            a.to_numpy(), b.to_numpy(), atol=float_tol, rtol=0, equal_nan=True
        ):
            diffs = np.abs(a.to_numpy() - b.to_numpy())
            msg = (
                f"Column {col} differs "
                f"(max abs diff {np.nanmax(diffs)} > tol {float_tol})"
            )
            raise AssertionError(msg)
        elif not a.equals(b):
            msg = f"Column {col} has non-numeric mismatches."
            raise AssertionError(msg)


@pytest.mark.parametrize("run", ["no_ref", "with_ref", "self_diff"])
def test_integration_with_real_data(real_input_dir, output_file_path, run):
    args = [
        "--sim-path",
        str(real_input_dir),
        "--mode",
        "both",
        "--results-dir",
        str(output_file_path),
    ]

    if run == "with_ref":
        args += ["--ref-sim-path", "tests/data/DINO_subsampled_B"]
    elif run == "self_diff":
        args += ["--ref-sim-path", str(real_input_dir)]

    exit_code = main(args)
    assert exit_code == 0

    out_grid = output_file_path / "metrics_results_grid.csv"
    out_restart = output_file_path / "metrics_results_restart.csv"
    assert out_grid.exists()
    assert out_restart.exists()

    if run in {"no_ref", "with_ref"}:
        base = Path(
            "tests/data/results/diff"
            if run == "with_ref"
            else "tests/data/results/original"
        )
        ref_grid = base / "metrics_results_grid.csv"
        ref_restart = base / "metrics_results_restart.csv"
        print(out_grid)
        assert_csv_equal(str(out_grid), str(ref_grid))
        assert_csv_equal(str(out_restart), str(ref_restart))

    elif run == "self_diff":
        # Load outputs and check all diffs are zero
        grid_df = pd.read_csv(out_grid)
        restart_df = pd.read_csv(out_restart)

        # Any column starting with "diff_" should be zero
        for df in (grid_df, restart_df):
            diff_cols = [c for c in df.columns if c.startswith("diff_")]
            for col in diff_cols:
                assert np.allclose(df[col].fillna(0).values, 0.0, atol=1e-12), (
                    f"Nonzero in {col}"
                )

        # For summary file: all MAE/RMSE should be zero
        summary_grid = pd.read_csv(
            out_grid.with_name("metrics_results_grid_summary.csv")
        )
        summary_restart = pd.read_csv(
            out_restart.with_name("metrics_results_restart_summary.csv")
        )
        for df in (summary_grid, summary_restart):
            assert np.allclose(df[["mae", "rmse"]].fillna(0).values, 0.0, atol=1e-12)
