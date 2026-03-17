from pathlib import Path

import pytest
import xarray as xr
import yaml

from spinup_evaluation.loader import load_dino_data, load_grid_variables

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
    config_path = Path("configs/DINO-setup.yaml")
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
        "restart_files": "restart.nc",
        "output_variables": {"temperature": "grid_T_3D.nc"},
    }

    with pytest.raises(FileNotFoundError, match=r"Mesh mask"):
        load_dino_data(mode, test_data_path, bad_setup)

    # Setup without mesh_mask
    bad_setup = {
        "restart_files": "restart",
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
        "restart_files": "missing_restart",
        "output_variables": {"temperature": "grid_T_3D.nc"},
    }

    with pytest.raises(FileNotFoundError, match=r"restart file"):
        load_dino_data(mode, test_data_path, bad_setup)


@pytest.mark.parametrize("mode", ["output", "both"])
def test_load_dino_data_missing_output_file(test_data_path, mode):
    """Test error handling when output file is missing."""

    # Setup with non-existent output file
    bad_setup = {
        "mesh_mask": "mesh_mask.nc",
        "restart_files": "restart",
        "output_variables": {"temperature": "missing.nc"},
    }

    with pytest.raises(FileNotFoundError, match=r"missing.nc"):
        load_dino_data(mode, test_data_path, bad_setup)


@pytest.mark.parametrize("mode", ["restart", "output", "both"])
def test_load_dino_data_modes(test_data_path, dino_setup, mode):
    """Test full load of DINO data in each mode."""

    # Load data in restart mode
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
    data = load_dino_data(mode, test_data_path, dino_setup, do_standardise=True)

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


def test_load_grid_variables_simple_format(test_data_path):
    """Test loading grid variables using simple format (variable: filename)."""

    # Simple format - just filename
    simple_specs = {
        "temperature": "grid_T_3D.nc",
        "salinity": "grid_T_3D.nc",
        "density": "grid_T_3D.nc",
        "ssh": "grid_T_2D.nc",
    }

    files_cache = {}
    grid_vars = load_grid_variables(str(test_data_path), simple_specs, files_cache)

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


def test_load_grid_variables_rich_format(test_data_path):
    """Test loading grid variables using rich format (variable: {file: .., var: ..})."""

    # Rich format - explicit file and variable names
    # Custom names are used to test that the explicit variable names
    # are used instead of substitutions based on the key name
    rich_specs = {
        "temperature": {"file": "grid_T_3D.nc", "var": "toce"},
        "custom_salinity": {"file": "grid_T_3D.nc", "var": "soce"},
        "custom_density": {"file": "grid_T_3D.nc", "var": "rhop"},
        "ssh": {"file": "grid_T_2D.nc", "var": "ssh"},
    }

    files_cache = {}
    grid_vars = load_grid_variables(str(test_data_path), rich_specs, files_cache)

    # Verify grid variable names and datasets
    assert isinstance(grid_vars, dict)

    assert grid_vars["temperature"].name == "toce"
    assert grid_vars["custom_salinity"].name == "soce"
    assert grid_vars["custom_density"].name == "rhop"
    assert grid_vars["ssh"].name == "ssh"

    assert isinstance(grid_vars["temperature"], xr.DataArray)
    assert isinstance(grid_vars["custom_salinity"], xr.DataArray)
    assert isinstance(grid_vars["custom_density"], xr.DataArray)
    assert isinstance(grid_vars["ssh"], xr.DataArray)


def test_load_grid_variables_rich_format_missing_var(test_data_path):
    """
    Test error handling for non-existent variables when loading grid in rich format.
    """

    # Rich format - explicit file and variable names, some invalid
    rich_specs = {
        "temperature": {"file": "grid_T_3D.nc", "var": "toce"},
        "custom_salinity": {"file": "grid_T_3D.nc", "var": "no_var"},
        "custom_density": {"file": "grid_T_3D.nc", "var": "rhop"},
        "ssh": {"file": "grid_T_2D.nc", "var": "no_var"},
    }

    files_cache = {}
    with pytest.raises(KeyError, match=r"No variable named"):
        load_grid_variables(str(test_data_path), rich_specs, files_cache)
