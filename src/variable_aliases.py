"""Lookup table for mapping standard variables with different variants."""

VARIABLE_ALIASES = {
    "temperature": [
        "toce",
        "tn",
        "temperature",
    ],  # Example: can be called 'toce' or 'tn' in various datasets
    "velocity_zonal": [
        "un",
        "u",
        "uoce",
        "zonal_velocity",
    ],  # Zonal velocity can be 'un', 'u', or 'zonal_velocity'
    "velocity_meridional": [
        "vn",
        "v",
        "voce",
        "meridional_velocity",
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


def standardize_variables(dataset, variable_dict):
    """
    Standardizes variable, coordinate, and dimension names across datasets.

    Parameters
    ----------
    dataset : xarray.Dataset
        The dataset containing variables, coordinates, and dimensions.
    variable_dict : dict
        A dictionary mapping standardized names to possible alternatives.

    Returns
    -------
    xarray.Dataset
        The dataset with standardized names.
    """
    rename_map = {}

    # Iterate over the variable aliases dictionary
    for standard_name, aliases in variable_dict.items():
        for alias in aliases:
            # Check for variables to rename
            if alias in dataset.variables:
                rename_map[alias] = standard_name
                break
            # Check for coordinates to rename
            if alias in dataset.coords:
                rename_map[alias] = standard_name
                break
            # Check for dimensions to rename (if the alias exists in dims)
            if alias in dataset.dims:
                rename_map[alias] = standard_name
                break

    # Rename variables, coordinates, and dimensions
    dataset = dataset.rename(rename_map)

    # Explicitly handle renaming of 'x' and 'y' dimensions to 'nav_lon' and 'nav_lat'
    if "x" in dataset.dims:
        dataset = dataset.rename({"x": "nav_lon"})
    if "y" in dataset.dims:
        dataset = dataset.rename({"y": "nav_lat"})

    return dataset
