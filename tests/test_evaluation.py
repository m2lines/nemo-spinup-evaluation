import subprocess
from pathlib import Path

import numpy as np
import pytest
import xarray as xr


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
    return Path("tests/data/subsampled/")


@pytest.fixture
def output_file_path(tmp_path):
    """
    Fixture to provide an output file path.
    """
    return tmp_path / "metrics_results_grid.csv"


def test_integration_with_real_data(real_input_dir, tmp_path):
    """
    Integration test using actual data from tests/data.
    Compares generated results to reference files.
    """
    # Output file paths
    output_grid = tmp_path / "metrics_results_grid.csv"
    output_restart = tmp_path / "metrics_results_restart.csv"

    # Run CLI for 'both' mode to generate both files
    result = subprocess.run(  # noqa: S603
        [
            str(Path(__import__("sys").executable)),
            "-m",
            "spinup_evaluation.cli",
            "--sim-path",
            str(real_input_dir),
            "--mode",
            "both",
            "--results",
            str(tmp_path / "metrics_results"),
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, f"stderr: {result.stderr}"
    assert output_grid.exists()

    # Compare output files to reference files
    ref_grid = Path("tests/data/results_resampled/metrics_results_grid.csv")
    ref_restart = Path("tests/data/results_resampled/metrics_results_restart.csv")

    with open(output_grid) as out, open(ref_grid) as ref:
        assert out.read() == ref.read(), "Grid results do not match reference."

    with open(output_restart) as out, open(ref_restart) as ref:
        assert out.read() == ref.read(), "Restart results do not match reference."
