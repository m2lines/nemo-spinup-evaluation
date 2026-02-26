from pathlib import Path

import pytest
import xarray as xr
import yaml

from spinup_evaluation.loader import load_dino_data, load_grid_variables


@pytest.fixture
def dino_setup():
    """Load the DINO setup configuration."""
    config_path = Path("configs/DINO-setup.yaml")
    with config_path.open("r") as f:
        return yaml.safe_load(f)


@pytest.fixture
def test_data_path() -> Path:
    """Path to the test data directory."""
    return Path("tests/data/DINO_subsampled_A")


def test_load_dino_data_invalid_mode(test_data_path, dino_setup):
    """Test error handling when an invalid mode is provided."""

    with pytest.raises(ValueError, match=r"'output', 'restart', 'both'"):
        load_dino_data("invalid_mode", test_data_path, dino_setup)


def test_load_dino_data_missing_mesh_mask(test_data_path):
    """Test error handling when mesh_mask is missing from setup."""

    # Setup without mesh_mask
    bad_setup = {
        "restart_files": "restart",
        "output_variables": {"temperature": "grid_T_3D.nc"},
    }

    with pytest.raises(ValueError, match=r"mesh_mask"):
        load_dino_data("restart", test_data_path, bad_setup)


def test_load_dino_data_missing_restart_file(test_data_path):
    """Test error handling when restart file is missing."""

    # Setup with non-existent restart file
    bad_setup = {
        "mesh_mask": "mesh_mask.nc",
        "restart_files": "missing_restart",
        "output_variables": {"temperature": "grid_T_3D.nc"},
    }

    with pytest.raises(FileNotFoundError, match=r"restart file"):
        load_dino_data("restart", test_data_path, bad_setup)

    with pytest.raises(FileNotFoundError, match=r"restart file"):
        load_dino_data("output", test_data_path, bad_setup)

    with pytest.raises(FileNotFoundError, match=r"restart file"):
        load_dino_data("both", test_data_path, bad_setup)


def test_load_dino_data_missing_output_file(test_data_path):
    """Test error handling when output file is missing."""

    # Setup with non-existent output file
    bad_setup = {
        "mesh_mask": "mesh_mask.nc",
        "restart_files": "restart",
        "output_variables": {"temperature": "missing.nc"},
    }

    with pytest.raises(FileNotFoundError, match=r"missing.nc"):
        load_dino_data("output", test_data_path, bad_setup)


def test_load_dino_data_restart_mode(test_data_path, dino_setup):
    """Test loading data in restart mode."""

    # Load data in restart mode
    data = load_dino_data("restart", test_data_path, dino_setup)

    # Verify the structure
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

    # Verify grid is empty in restart mode
    assert not data["grid"]

    # Verify paths
    assert isinstance(data["paths"], dict)
    assert "base" in data["paths"]
    assert "mesh_mask" in data["paths"]
    assert "restart" in data["paths"]


def test_load_dino_data_output_mode(test_data_path, dino_setup):
    """Test loading data in output mode."""

    # Load data in output mode
    data = load_dino_data("output", test_data_path, dino_setup)

    # Verify the structure
    assert "mesh_mask" in data
    assert "restart" in data
    assert "grid" in data
    assert "files" in data
    assert "paths" in data

    # Verify mesh_mask loaded
    assert isinstance(data["mesh_mask"], xr.Dataset)

    # Verify restart loaded
    assert isinstance(data["restart"], xr.Dataset)

    # Verify grid variables
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
    assert "output_files" in data["paths"]
    assert len(data["paths"]["output_files"]) > 0


def test_load_dino_data_both_mode(test_data_path, dino_setup):
    """Test loading data in both mode."""

    # Load data in both mode
    data = load_dino_data("both", test_data_path, dino_setup)

    # Verify the structure
    assert "mesh_mask" in data
    assert "restart" in data
    assert "grid" in data
    assert "files" in data
    assert "paths" in data

    # Verify mesh_mask loaded
    assert isinstance(data["mesh_mask"], xr.Dataset)

    # Verify restart loaded
    assert isinstance(data["restart"], xr.Dataset)

    # Verify grid variables
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
    assert "output_files" in data["paths"]
    assert len(data["paths"]["output_files"]) > 0


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
    # Note that custom names for temperature and ssh
    # will currently fail time donor allocation
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
