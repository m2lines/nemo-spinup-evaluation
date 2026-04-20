from pathlib import Path

import cftime
import pytest
import xarray as xr
import yaml

from nemo_spinup_evaluation.loader import (
    _check_grid_time_alignment,
    _check_required_coords,
    load_dino_data,
    load_grid_variables,
    standardise_vars,
)

"""
Parameters found in these tests mimic loader.load_dino_data:

Parameters
----------
test_data_path : fixture
    The base directory for loading data files.
dino_setup : fixture
    A mapping containing the YAML configuration for data loading.
mode : parametrize
    The mode of operation, one of "restart", "output" or "both"
"""


@pytest.fixture
def test_data_path() -> Path:
    """Fixture path to the test data directory."""
    return Path("tests/data/DINO_subsampled_A")


@pytest.fixture
def dino_setup():
    """Fixture to load the DINO setup configuration."""
    config_path = Path("tests/DINO-test-setup.yaml")
    with config_path.open("r") as f:
        return yaml.safe_load(f)


def test_load_dino_data_invalid_mode(test_data_path, dino_setup):
    """Test error handling when an invalid mode option is provided."""

    with pytest.raises(ValueError, match=r"'output', 'restart', 'both'"):
        load_dino_data("invalid_mode", test_data_path, dino_setup)


@pytest.mark.parametrize("mode", ["restart", "output", "both"])
def test_load_dino_data_missing_mesh_mask(test_data_path, mode):
    """Test error handling when mesh_mask is missing from setup."""

    # Setup with missing mesh_mask file
    bad_setup = {
        "mesh_mask": "missing_mesh_mask.nc",
        "restart": "restart.nc",
        "output_variables": {"temperature": "grid_T_3D.nc"},
    }

    with pytest.raises(FileNotFoundError, match=r"mesh_mask"):
        load_dino_data(mode, test_data_path, bad_setup)

    # Setup without mesh_mask
    bad_setup = {
        "restart": "restart.nc",
        "output_variables": {"temperature": "grid_T_3D.nc"},
    }

    with pytest.raises(ValueError, match=r"mesh_mask"):
        load_dino_data(mode, test_data_path, bad_setup)


@pytest.mark.parametrize("mode", ["restart", "output", "both"])
def test_load_dino_data_missing_restart_file(test_data_path, mode):
    """Test error handling when restart file is missing."""

    # Setup with non-existent restart file
    bad_setup = {
        "mesh_mask": "mesh_mask.nc",
        "restart": "missing_restart.nc",
        "output_variables": {"temperature": {"file": "grid_T_3D.nc", "var": "toce"}},
    }

    with pytest.raises(FileNotFoundError, match=r"restart"):
        load_dino_data(mode, test_data_path, bad_setup)


def test_load_dino_data_no_restart_output_mode(test_data_path):
    """Test that restart file is optional in output mode."""

    # Setup with no restart entry
    setup = {
        "mesh_mask": "mesh_mask.nc",
        "output_variables": {"temperature": {"file": "grid_T_3D.nc", "var": "toce"}},
        "variable_map": {"depth": ["nav_lev"]},
    }

    data = load_dino_data("output", test_data_path, setup)
    assert data["restart"] is None


@pytest.mark.parametrize("mode", ["output", "both"])
def test_load_dino_data_missing_output_file(test_data_path, mode):
    """Test error handling when output file is missing."""

    # Setup with non-existent output file
    bad_setup = {
        "mesh_mask": "mesh_mask.nc",
        "restart": "restart.nc",
        "output_variables": {"temperature": {"file": "missing.nc", "var": "toce"}},
    }

    with pytest.raises(FileNotFoundError, match=r"missing.nc"):
        load_dino_data(mode, test_data_path, bad_setup)


@pytest.mark.parametrize("mode", ["restart", "output", "both"])
def test_load_dino_data_modes(test_data_path, dino_setup, mode):
    """Test full load of DINO data in each mode."""

    # Load data
    data = load_dino_data(mode, test_data_path, dino_setup)

    # Verify structure
    assert "mesh_mask" in data
    assert "restart" in data
    assert "grid" in data
    assert "files" in data
    assert "paths" in data

    # Verify mesh_mask loaded
    assert isinstance(data["mesh_mask"], xr.Dataset)
    assert "tmask" in data["mesh_mask"]

    # Verify restart loaded
    assert isinstance(data["restart"], xr.Dataset)

    # Verify grid variables
    if mode == "restart":
        assert not data["grid"]
    else:
        assert isinstance(data["grid"], dict)
        assert isinstance(data["grid"]["temperature"], xr.DataArray)
        assert isinstance(data["grid"]["salinity"], xr.DataArray)
        assert isinstance(data["grid"]["density"], xr.DataArray)
        assert isinstance(data["grid"]["ssh"], xr.DataArray)
        assert isinstance(data["grid"]["velocity_u"], xr.DataArray)
        assert isinstance(data["grid"]["velocity_v"], xr.DataArray)

    # Verify paths
    assert isinstance(data["paths"], dict)
    assert "base" in data["paths"]
    assert "mesh_mask" in data["paths"]
    assert "restart" in data["paths"]

    if mode in ["output", "both"]:
        assert "output_files" in data["paths"]
        assert len(data["paths"]["output_files"]) > 0


@pytest.mark.parametrize("mode", ["restart", "output", "both"])
def test_load_dino_data_nav_latlon_promotion(test_data_path, dino_setup, mode):
    """Test that nav_lat and nav_lon are promoted to coordinates."""

    # Load data in output mode
    data = load_dino_data(mode, test_data_path, dino_setup)

    assert isinstance(data["restart"], xr.Dataset)
    assert isinstance(data["restart"]["temperature"], xr.DataArray)

    # Verify nav_lat and nav_lon are coordinates
    assert "nav_lat" in data["restart"]["temperature"].coords
    assert "nav_lon" in data["restart"]["temperature"].coords

    # Check a selection to ensure promoted coords are accessible
    assert len(data["restart"]["temperature"].nav_lat) == len(data["restart"].y)

    DEGREES = -65.0
    NUM_BELOW = 806
    sel = data["restart"]["temperature"].nav_lat < DEGREES
    num_true = sel.sum()
    assert num_true == NUM_BELOW


def test_load_grid_variables(test_data_path):
    """Test loading grid variables."""

    grid_specs = {
        "temperature": {"file": "grid_T_3D.nc", "var": "toce"},
        "salinity": {"file": "grid_T_3D.nc", "var": "soce"},
        "density": {"file": "grid_T_3D.nc", "var": "rhop"},
        "ssh": {"file": "grid_T_2D.nc", "var": "ssh"},
    }

    files_cache = {}
    grid_vars = load_grid_variables(str(test_data_path), grid_specs, files_cache)

    # Verify grid variable names and datasets
    assert isinstance(grid_vars, dict)

    assert grid_vars["temperature"].name == "toce"
    assert grid_vars["salinity"].name == "soce"
    assert grid_vars["density"].name == "rhop"
    assert grid_vars["ssh"].name == "ssh"

    assert isinstance(grid_vars["temperature"], xr.DataArray)
    assert isinstance(grid_vars["salinity"], xr.DataArray)
    assert isinstance(grid_vars["density"], xr.DataArray)
    assert isinstance(grid_vars["ssh"], xr.DataArray)


def test_load_grid_variables_missing_var(test_data_path):
    """
    Test error handling for non-existent variables when loading grid files.
    """

    grid_specs = {
        "temperature": {"file": "grid_T_3D.nc", "var": "toce"},
        "salinity": {"file": "grid_T_3D.nc", "var": "no_var"},
        "density": {"file": "grid_T_3D.nc", "var": "rhop"},
        "ssh": {"file": "grid_T_2D.nc", "var": "no_var"},
    }

    files_cache = {}
    with pytest.raises(KeyError, match=r"No variable named"):
        load_grid_variables(str(test_data_path), grid_specs, files_cache)


def test_check_required_coords():
    """Test that _check_required_coords validates coordinates."""

    ds = xr.Dataset(
        {"temperature": (("time_counter", "y", "x"), [[[0.0]]])},
        coords={
            "time_counter": [0],
            "nav_lat": (("y", "x"), [[0.0]]),
            "nav_lon": (("y", "x"), [[0.0]]),
        },
    )

    # Check that coords present are recognised
    _check_required_coords(ds, ("time_counter", "nav_lat", "nav_lon"), "test_dataset")
    _check_required_coords(ds, ("nav_lat", "nav_lon"), "test_dataset")

    # Check that missing coords cause a KeyError: for one missing
    with pytest.raises(KeyError, match=r"missing in test_dataset: depth"):
        _check_required_coords(ds, ("time_counter", "depth"), "test_dataset")

    # ...and multiple missing
    with pytest.raises(KeyError, match=r"missing in test_dataset: missing, depth"):
        _check_required_coords(ds, ("missing", "depth", "time_counter"), "test_dataset")


def test_check_grid_time_alignment(test_data_path):
    """Test that _check_grid_time_alignment validates time_counters."""

    grid_specs = {
        "temperature": {"file": "grid_T_3D.nc", "var": "toce"},
        "salinity": {"file": "grid_T_3D.nc", "var": "soce"},
        "density": {"file": "grid_T_3D.nc", "var": "rhop"},
        "ssh": {"file": "grid_T_2D.nc", "var": "ssh"},
    }

    files_cache = {}
    grid_vars = load_grid_variables(test_data_path, grid_specs, files_cache)

    # Test data with correct temporal alignment
    _check_grid_time_alignment(grid_vars)

    # Test different times by changing the first entry
    orig_coords = grid_vars["ssh"].coords["time_counter"].values.copy()
    new_coords = orig_coords.copy()
    new_coords[0] = cftime.Datetime360Day(1, 1, 1, has_year_zero=True)
    grid_vars["ssh"] = grid_vars["ssh"].assign_coords(time_counter=new_coords)

    with pytest.raises(ValueError, match=r"ssh differ from temperature"):
        _check_grid_time_alignment(grid_vars)

    ## Reset coords
    grid_vars["ssh"] = grid_vars["ssh"].assign_coords(time_counter=orig_coords)

    # Test different length time_counters by adding a new time step
    extra = grid_vars["ssh"].isel(time_counter=-1)
    extra = extra.assign_coords(
        time_counter=cftime.Datetime360Day(3, 7, 1, has_year_zero=True)
    )
    grid_vars["ssh"] = xr.concat([grid_vars["ssh"], extra], dim="time_counter")

    with pytest.raises(ValueError, match=r"ssh differ from temperature"):
        _check_grid_time_alignment(grid_vars)


def test_standardise_vars(test_data_path, dino_setup):
    """Test renaming of loaded variables."""

    # Load all data (mode, both)
    data = load_dino_data("both", test_data_path, dino_setup)

    # Check that mesh mask variables were renamed
    assert isinstance(data["mesh_mask"], xr.Dataset)
    assert all(v in data["mesh_mask"] for v in ("time_counter", "depth"))

    # Check that restart variables were renamed
    assert isinstance(data["restart"], xr.Dataset)
    assert all(
        v in data["restart"]
        for v in (
            "velocity_u",
            "velocity_v",
            "temperature",
            "salinity",
            "density",
            "ssh",
        )
    )

    # Grid variables are not renamed
    assert isinstance(data["grid"], dict)
    assert data["grid"]["temperature"].name == "toce"


def test_standardise_vars_nav_latlon():
    """Test that nav_lat and nav_lon variables are renamed and promoted."""

    # Simplified dataset with lat and lon as unexpected variables (not coords)
    ds = xr.Dataset(
        {
            "temperature": (("time", "y", "x"), [[[0.0]]]),
            "latitude": (("y", "x"), [[0.0]]),
            "longitude": (("y", "x"), [[0.0]]),
        },
        coords={
            "time": [0],
        },
    )

    # Confirm pre-standardisation layout
    assert "latitude" in ds.variables
    assert "longitude" in ds.variables
    assert "time" in ds.coords

    var_map = {
        "time_counter": ["time"],
        "nav_lat": ["latitude"],
        "nav_lon": ["longitude"],
    }

    ds_std = standardise_vars(ds, var_map)

    assert "time_counter" in ds_std.coords
    assert "nav_lat" in ds_std.coords
    assert "nav_lon" in ds_std.coords
