# Model Configuration

## DINO Config File

The YAML config (e.g., `configs/DINO-setup.yaml`) is used to specify the names of NetCDF files and variables to load along with a map for converting variable names used in the files to the canonical field names.

The canonical field names used internally:

```
time_counter
nav_lat
nav_lon
temperature
velocity_u
velocity_v
depth
ssh
salinity
density
```

If NetCDF input files use different variable names to these then use the `variable_map` to remap them. All loaded files are remapped using the specified map. If input files use differing names for the same variable, list all names against each canonical field name, e.g.

```yaml
variable_map:
  time_counter: [time]
  nav_lat: [latitude]
  nav_lon: [longitude]
  temperature: [toce, tn]
  velocity_u: [un, u, uoce, zonal_velocity]
  velocity_v: [vn, v, voce, meridional_velocity]
  depth: [nav_lev, deptht, depthu, depthv]
  ssh: [sshn]
  salinity: [so, sn, soce]
  density: [rhop]
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

variable_map:
  temperature: [tn]
  velocity_u: [un]
  velocity_v: [vn]
  depth: [nav_lev, deptht, depthu, depthv]
  ssh: [sshn]
  salinity: [sn]
  density: [rhop]
```


## OM4 Config File

## ORCA2 Config File
