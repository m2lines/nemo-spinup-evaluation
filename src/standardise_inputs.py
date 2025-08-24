"""Lookup table for mapping standard variables with different variants."""

import xarray as xr

VARIABLE_ALIASES = {
    "temperature": [
        "toce",
        "tn",
        "temperature",
    ],  # Example: can be called 'toce' or 'tn' in various datasets
    "velocity_u": [
        "un",
        "u",
        "uoce",
        "zonal_velocity",
        "velocity_u",
    ],  # Zonal velocity can be 'un', 'u', or 'zonal_velocity'
    "velocity_v": [
        "vn",
        "v",
        "voce",
        "meridional_velocity",
        "velocity_v",
    ],  # Meridional velocity
    "depth": [
        "depth",
        "nav_lev",
        "deptht",
        "depthu",
        "depthv",
    ],  # Depth can be 'depth', 'nav_lev', or 'deptht'
    "latitude": ["nav_lat", "y"],  # Latitude can be 'nav_lat' or 'y'
    "longitude": ["nav_lon", "x"],  # Longitude can be 'nav_lon' or 'x'
    "ssh": ["sshn", "ssh"],  # Sea surface height could be 'sshn' or 'ssh'
    "time_counter": [
        "time_counter",
        "time",
    ],  # Time variable may be 'time_counter' or 'time'
    "salinity": ["so", "sn", "soce", "salinity"],  # Example for salinity variable
    "density": ["rhop", "density"],  # Density might be named 'rhop' or 'density'
    # Add more mappings if necessary depending on your datasets
}


def standardise(dataset, variable_dict):
    """
    Rename variables/coords/dims using aliases.

    Returns the same type as input (DataArray or Dataset).
    """
    is_da = isinstance(dataset, xr.DataArray)
    if is_da:
        orig_name = dataset.name or "var"
        ds = dataset.to_dataset(name=orig_name)
    else:
        if not isinstance(dataset, xr.Dataset):
            msg = "standardise expects an xarray.Dataset or xarray.DataArray"
            raise TypeError(msg)
        ds = dataset

    rename_map = {}
    for std, aliases in variable_dict.items():
        for alias in aliases:
            if alias in ds.variables:
                rename_map[alias] = std
                break
            if alias in ds.coords:
                rename_map[alias] = std
                break
            if alias in ds.dims:
                rename_map[alias] = std
                break

    ds = ds.rename(rename_map)

    if "x" in ds.dims:
        ds = ds.rename({"x": "nav_lon"})
    if "y" in ds.dims:
        ds = ds.rename({"y": "nav_lat"})

    if is_da:
        new_name = rename_map.get(orig_name, orig_name)
        return ds[new_name]

    return ds
