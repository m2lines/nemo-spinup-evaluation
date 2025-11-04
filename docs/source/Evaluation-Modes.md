## Evaluation Modes

Spinup-Evaluation supports three modes, controlled by the `--mode` argument:

### `restart` mode (Instantaneous Output)
- **Purpose:** Evaluate a single model state (snapshot) from a NEMO/DINO `restart.nc` file.
- **Input:** `restart.nc` (and `mesh_mask.nc`)
- **Use case:** Assess the physical realism or convergence of a single model state, e.g., after a spin-up or forecast.
- **Output:** `results/metrics_results_restart.csv` (or your chosen prefix)
- **Reference:** If `--ref-sim-path` is provided, computes diffs/stats vs. a reference restart file.

### `output` mode (Time-Averaged State)
- **Purpose:** Evaluate time-averaged or multi-time-step model output, typically from files like `grid_T_3D.nc`, `grid_U_3D.nc`, `grid_V_3D.nc`, `grid_T_2D.nc`.
- **Input:** Grid files as mapped in the config YAML (see below).
- **Use case:** Assess the mean state or variability over a period, or compare time-averaged fields between runs.
- **Output:** `results/metrics_results_grid.csv` (or your chosen prefix)
- **Reference:** If `--ref-sim-path` is provided, computes diffs/stats vs. a reference output set.

### `both` mode
- **Purpose:** Run both `restart` and `output` metric suites in one command.
- **Output:** Both CSVs as above.
