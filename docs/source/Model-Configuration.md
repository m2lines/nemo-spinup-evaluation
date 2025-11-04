""" This page id dedicated to elaborate the configuration files for ocean models """


# Model Configuration

## DINO Config File

The YAML config (e.g., `configs/DINO-setup.yaml`) maps variable names to NetCDF files. You can specify variables in two ways:

### 1. Simple Form

```yaml
output_variables:
  temperature: grid_T_3D.nc
  salinity: grid_T_3D.nc
  # ...
```
- **Behavior:** The loader will try to infer the correct variable name (e.g., `toce` for temperature) from a list of likely candidates for each field.

### 2. Rich Form

```yaml
output_variables:
  temperature:
    file: grid_T_3D.nc
    var: toce
    time_from: density  # (optional) use time axis from another variable
  # ...
```
- **Behavior:** You can explicitly specify the file, the variable name within the file, and optionally a `time_from` field to use the time axis from another variable.

You can mix and match simple and rich forms in the same config. The loader will handle both.

> **Note:** Support for specifying temporal granularities and resampling (e.g., daily, monthly, seasonal means) is under active development and will be available in a future release.

**Example config:**

```yaml
mesh_mask: mesh_mask.nc
restart_files: 'restart'
output_variables:
  temperature: grid_T_3D.nc
  salinity:
    file: grid_T_3D.nc
    var: soce
  density: grid_T_3D.nc
  ssh: grid_T_2D.nc
  velocity_u: grid_U_3D.nc
  velocity_v: grid_V_3D.nc
```
