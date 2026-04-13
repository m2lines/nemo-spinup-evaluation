# Command-Line Usage

Once the package is installed it exposes the `nemo-spinup-evaluation` console script. You can also run it as a module with `python -m nemo_spinup_evaluation`.

```sh
nemo-spinup-evaluation \
  [--sim-path <simulation_dir>]          # Path to the simulation directory
  [--ref-sim-path <reference_sim_dir>]   # Path to a reference simulation directory
  [--config configs/DINO-setup.yaml]     # YAML config file (default shown)
  [--mode output|restart|both]           # Which metric suite(s) to run (default: both)
  [--results-dir results]                # Output directory (default shown)
  [--result-file-prefix metrics_results] # Output file prefix (default shown)
  [--eager]                              # Load data eagerly with numpy (default: off)
```

**Arguments:**
- `--sim-path`: Path to the simulation directory.
- `--ref-sim-path`: Path to a reference simulation directory (enables comparison mode).
- `--config`: Path to the YAML config file (default: `configs/DINO-setup.yaml`).
- `--mode`: Which metric suite(s) to run: `output`, `restart`, or `both` (default: `both`).
- `--results-dir`: Directory to save output CSVs (default: `results`).
- `--result-file-prefix`: Prefix for output files (default: `metrics_results`).
- `--eager`: Load data eagerly using numpy instead of lazily with dask. Faster for small files.

At least one of `--sim-path` or `--ref-sim-path` must be supplied. Providing only `--sim-path` runs the metrics on that simulation in isolation; providing both runs a diff between the two.
