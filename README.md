[![Documentation Status](https://readthedocs.org/projects/nemo-spinup-evaluation/badge/?version=latest)](https://nemo-spinup-evaluation.readthedocs.io/en/latest/)[![CI](https://github.com/m2lines/Spinup-Evaluation/actions/workflows/ci-eval.yml/badge.svg)](https://github.com/m2lines/nemo-spinup-evaluation/actions/workflows/ci-eval.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

# Spinup-Evaluation

Spinup-Evaluation provides a command-line tool and Python API for benchmarking the spin-up and restart performance of NEMO/DINO ocean models and machine learning emulators. It supports both single-run and comparison (reference) evaluation, and outputs detailed metrics and difference statistics.

📖 Full documentation is available on [ReadTheDocs](https://nemo-spinup-evaluation.readthedocs.io/en/latest/).


## Features

- **Flexible CLI**: Evaluate restart and/or output files, with or without a reference simulation.
- **Configurable**: Uses a YAML config file (e.g., `configs/DINO-setup.yaml`) to map variables to files.
- **Comparison Mode**: Computes diffs, MAE, and RMSE between a simulation and a reference.
- **Modern Output**: Results are written as CSV files (one for restart, one for output).
- **Test Suite**: Integration and regression tests using real and subsampled NetCDF data.
- **Extensible**: Add new metrics by editing `src/spinup_evaluation/metrics.py`.



## Evaluation Flow

Spinup-Evaluation is designed to assess the quality and stability of ocean model spin-up and restart states, as well as time-averaged outputs. The evaluation workflow is flexible: you can analyse a single simulation, or compare a simulation against a reference (e.g., a previous spin-up, a control run, or a forecast). The tool supports both instantaneous (restart) and time-averaged (output) evaluation modes.

The diagram below (Figure 1) illustrates the typical evaluation procedure. Model output files (restart and/or time-averaged NetCDFs) are loaded and standardized according to the YAML config. Metrics are computed, and—if a reference is provided—differences, MAE, and RMSE are calculated.

Spinup-Evaluation is often used alongside [spinup-forecast](https://github.com/m2lines/nemo-spinup-forecast), which automates the generation of machine learned spin-up states for NEMO/DINO models. Together, these tools provide a robust workflow for accelerating ocean spin-up.

<p align="center">
<img src="diagram.png" alt="NEMO flow" width="500"/>
<figcaption>Fig 1. Evaluation flow diagram illustrating the coupling to spinup-forecast, but spinup-evaluation can in theory be used to evaluate any ocean model, be it ML data driven, numerical or otherwise. </figcaption>
</p>

## Repository Layout


```
.
├── pyproject.toml                  # Project metadata, dependencies, and build system
├── README.md                       # Main project documentation (this file)
├── configs/                        # Configuration files for variable/file mapping
│   └── DINO-setup.yaml             # Example YAML config for DINO/NEMO variables
├── src/
│   └── spinup_evaluation/          # Main Python package
│       ├── cli.py                  # Command-line interface (CLI) entry point
│       ├── loader.py               # Data loading and preprocessing utilities
│       ├── metrics_io.py           # Output helpers (CSV writing, formatting)
│       ├── metrics.py              # Metric calculation functions
│       ├── standardise_inputs.py   # Input standardization helpers
│       └── utils.py                # General utilities
├── tests/                          # Test suite, test data, and data download scripts
│   └── get-data.sh                 # Script to fetch test data from THREDDS
└── results/                        # Default output directory for metrics CSVs
```

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

## Modes: What Do They Mean?

Spinup-Evaluation supports three modes, controlled by the `--mode` argument:

### 1. `restart` mode (Instantaneous Output)
- **Purpose:** Evaluate a single model state (snapshot) from a NEMO/DINO `restart.nc` file.
- **Input:** `restart.nc` (and `mesh_mask.nc`)
- **Use case:** Assess the physical realism or convergence of a single model state, e.g., after a spin-up or forecast.
- **Output:** `results/metrics_results_restart.csv` (or your chosen prefix)
- **Reference:** If `--ref-sim-path` is provided, computes diffs/stats vs. a reference restart file.

### 2. `output` mode (Time-Averaged State)
- **Purpose:** Evaluate time-averaged or multi-time-step model output, typically from files like `grid_T_3D.nc`, `grid_U_3D.nc`, `grid_V_3D.nc`, `grid_T_2D.nc`.
- **Input:** Grid files as mapped in the config YAML (see below).
- **Use case:** Assess the mean state or variability over a period, or compare time-averaged fields between runs.
- **Output:** `results/metrics_results_grid.csv` (or your chosen prefix)
- **Reference:** If `--ref-sim-path` is provided, computes diffs/stats vs. a reference output set.

### 3. `both` mode
- **Purpose:** Run both `restart` and `output` metric suites in one command.
- **Output:** Both CSVs as above.

## Config File

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

## Output Files

- Results are written as CSV files in the results directory, e.g.:
  - `results/metrics_results_restart.csv`
  - `results/metrics_results_grid.csv`
- Each file contains metric values, and if a reference is provided, also includes:
  - Reference metric values (prefixed with `ref_`)
  - Differences (`diff_*`)
- A separate file with MAE and RMSE statistics is also generated if a reference directory is provided.

## Testing

Tests are in the `tests/` directory and use real subsampled NetCDF data. First download the dataset as follows:

```sh
sh tests/get-data.sh
```
To run all tests:

```sh
pytest tests/
```

## Development & Installation

Clone the repo and install in development mode:

```sh
git clone https://github.com/m2lines/Spinup-Evaluation.git
cd Spinup-Evaluation
python -m venv venv
source venv/bin/activate
pip install -e .[dev]
pre-commit install
```

## Adding New Metrics

Add new metric functions to `src/spinup_evaluation/metrics.py` and update the metric function lists in `cli.py` as needed.

## Acknowledgements

This work builds on significant contributions by [Etienne Meunier](https://github.com/Etienne-Meunier), whose efforts on the [Metrics-Ocean](https://github.com/Etienne-Meunier/Metrics-Ocean) repository laid the foundation for several components used here.
