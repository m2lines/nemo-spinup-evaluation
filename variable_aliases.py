##Dictionary for file to file variable mapping

VARIABLE_ALIASES = {
    "temperature": ["toce", "tn", "temp", "temperature"],
    "salinity": ["soce", "sn", "salt", "salinity"],
    "u_velocity": ["uoce", "un", "u"],
    "v_velocity": ["voce", "vn", "v"],
    "w_velocity": ["woce", "wn", "w"],
    "time": ["time_counter", "time"],
    "depth": ["deptht", "depthu", "depthv", "depthw", "nav_lev"],
    "latitude": ["y", "nav_lat"],
    "longitude": ["x", "nav_lon"],
}

def standardize_variables(dataset, alias_dict):
    rename_map = {}
    for standard_name, aliases in alias_dict.items():
        for alias in aliases:
            if alias in dataset.variables:
                rename_map[alias] = standard_name
                break  # Use the first match
    return dataset.rename(rename_map)
