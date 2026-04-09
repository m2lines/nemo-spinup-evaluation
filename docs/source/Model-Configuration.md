# Model Configuration

## DINO Config File

The YAML config (e.g., `configs/DINO-setup.yaml`) maps each canonical field name to the NetCDF file it lives in and the variable name used inside that file. Both must be specified explicitly — the loader does not infer variable names.

```yaml
output_variables:
  temperature:
    filename: grid_T_3D.nc
    term: toce
  salinity:
    filename: grid_T_3D.nc
    term: soce
  # ...
```

> **Note:** Support for specifying temporal granularities and resampling (e.g., daily, monthly, seasonal means) is under active development and will be available in a future release.

**Example config:**

```yaml
mesh_mask: mesh_mask.nc
restart_files: 'restart'
output_variables:
  temperature:
    filename: grid_T_3D.nc
    term: toce
  salinity:
    filename: grid_T_3D.nc
    term: soce
  density:
    filename: grid_T_3D.nc
    term: rhop
  ssh:
    filename: grid_T_2D.nc
    term: ssh
  velocity_u:
    filename: grid_U_3D.nc
    term: uoce
  velocity_v:
    filename: grid_V_3D.nc
    term: voce
```


## OM4 Config File

## ORCA2 Config File
