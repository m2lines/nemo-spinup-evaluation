"""This page is dedicated to providing the linux commands to run the data files for evaluation."""

## Command-Line Usage

The main entry point is `src/spinup_evaluation/cli.py` (or the installed `spinup-eval` script):


```sh
python -m spinup_evaluation.cli \
  --sim-path <simulation_dir>            # Required: path to simulation directory
  [--ref-sim-path <reference_sim_dir>]   # Optional: path to reference simulation
  [--config configs/DINO-setup.yaml]     # Optional: YAML config file (default shown)
  [--results-dir results]                # Optional: output directory (default shown)
  [--result-file-prefix metrics_results] # Optional: output file prefix (default shown)
  [--mode output|restart|both]           # Optional: which metric suite(s) to run
```

**Arguments:**
- `--sim-path`: Path to the simulation directory (required).
- `--ref-sim-path`: Path to a reference simulation directory (optional, enables comparison).
- `--config`: Path to the YAML config file (default: `configs/DINO-setup.yaml`).
- `--results-dir`: Directory to save output CSVs (default: `results`).
- `--result-file-prefix`: Prefix for output files (default: `metrics_results`).
- `--mode`: Which metric suite(s) to run: `output`, `restart`, or `both` (default: `both`).
