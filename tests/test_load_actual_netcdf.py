import os

import xarray as xr


def test_actual_load():
    """
    Test loading an actual NetCDF file to ensure compatibility.
    """
    test_file = os.path.join(
        os.path.dirname(__file__), "data/subsampled", "grid_T_2D.nc"
    )
    data = xr.open_dataset(test_file)
    assert data is not None
    assert "ssh" in data.variables
    assert "sst" in data.variables
