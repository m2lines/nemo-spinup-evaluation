# Model Configuration

## DINO Config File

The YAML config (e.g., `configs/DINO-setup.yaml`) maps each canonical field name to a NetCDF file along with the corresponding variable name. Both must be specified explicitly — the loader does not infer variable names.

```yaml
output_variables:
  temperature:
    file: grid_T_3D.nc
    var: toce
  salinity:
    file: grid_T_3D.nc
    var: soce
  # ...
```

**Example config:**

```yaml
mesh_mask: mesh_mask.nc
restart_files: 'restart'
output_variables:
  temperature:
    file: grid_T_3D.nc
    var: toce
  salinity:
    file: grid_T_3D.nc
    var: soce
  density:
    file: grid_T_3D.nc
    var: rhop
  ssh:
    file: grid_T_2D.nc
    var: ssh
  velocity_u:
    file: grid_U_3D.nc
    var: uoce
  velocity_v:
    file: grid_V_3D.nc
    var: voce
```


## OM4 Config File

## ORCA2 Config File
