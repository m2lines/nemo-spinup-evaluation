[![CI](https://github.com/m2lines/nemo-spinup-evaluation/actions/workflows/ci-eval.yml/badge.svg)](https://github.com/m2lines/nemo-spinup-evaluation/actions/workflows/ci-eval.yml)[![Documentation Status](https://readthedocs.org/projects/nemo-spinup-evaluation/badge/?version=latest)](https://nemo-spinup-evaluation.readthedocs.io/en/latest/) [![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

# NEMO Spinup-Evaluation

NEMO Spinup-Evaluation (`nemo-spinup-evaluation`) provides a command-line tool and Python API for benchmarking the spin-up and restart performance of NEMO/DINO ocean models and machine learning emulators. It supports both single-run and comparison (reference) evaluation, and outputs detailed metrics and difference statistics.

📖 Full documentation is available on [Read the Docs](https://nemo-spinup-evaluation.readthedocs.io/en/latest/).


## Features

- **Flexible CLI**: Evaluate restart and/or output files, with or without a reference simulation.
- **Configurable**: Uses a YAML config file (e.g., `configs/DINO-setup.yaml`) to map variables to files.
- **Comparison Mode**: Computes diffs, MAE, and RMSE between a simulation and a reference.
- **Modern Output**: Results are written as CSV files (one for restart, one for output).
- **Test Suite**: Integration and regression tests using real and subsampled NetCDF data.
- **Extensible**: Add new metrics by editing `src/nemo_spinup_evaluation/metrics.py`.


## Installation

Requires Python ≥ 3.10

1. **Clone the repository**

    ```sh
    git clone https://github.com/m2lines/nemo-spinup-evaluation.git
    cd nemo-spinup-evaluation
    ```

2. **Create and activate a virtual environment**

    ```sh
    python3 -m venv .venv
    source .venv/bin/activate
    ```

3. **Install the package and dependencies**

    ```sh
    pip install .
    ```

    For developer installs, include development dependencies and enable pre-commit hooks:

    ```sh
    pip install -e .[dev]
    pre-commit install
	```


## Quick Start

Download the example DINO dataset from Zenodo:

```sh
# Download the 50-year baseline dataset (551.6 MB)
wget https://zenodo.org/records/19474414/files/50.zip
unzip 50.zip

# For more comprehensive testing, you can also download:
# 200-year dataset (2.1 GB): wget https://zenodo.org/records/19474414/files/200.zip
# Restart files (6.5 GB): wget https://zenodo.org/records/19474414/files/restart.zip
```

Then run the evaluation against the example data:

```sh
nemo-spinup-evaluation --sim-path <simulation_dir> --config configs/DINO-setup.yaml
```

Or compare against a reference simulation:

```sh
nemo-spinup-evaluation --sim-path <simulation_dir> --ref-sim-path <reference_dir> --config configs/DINO-setup.yaml
```


## Command-Line Usage

| Argument               | Description                                                   | Default                   | Required |
|------------------------|---------------------------------------------------------------|---------------------------|----------|
| `--sim-path`           | Path to the simulation directory                              |                           | Yes      |
| `--ref-sim-path`       | Path to a reference simulation directory to enable comparison |                           | No       |
| `--config`             | Path to the YAML config file                                  | `configs/DINO-setup.yaml` | No       |
| `--results-dir`        | Directory to save output CSVs                                 | `results`                 | No       |
| `--result-file-prefix` | Prefix for output files                                       | `metrics_results`         | No       |
| `--mode`               | Which metric suite(s) to run: `output`, `restart`, or `both`  | `both`                    | No       |


## Tests

Use the test suite to check integrity with real subsampled DINO data. Download the dataset using the script:

```sh
./tests/get-data.sh
```

Run all tests:

```sh
pytest tests/
```


## Evaluation Flow

NEMO Spinup-Evaluation is designed to assess the quality and stability of ocean model spin-up and restart states, as well as time-averaged outputs. The evaluation workflow is flexible: you can analyse a single simulation, or compare a simulation against a reference (e.g., a previous spin-up, a control run, or a forecast). The tool supports both instantaneous (restart) and time-averaged (output) evaluation modes.

The diagram below (Figure 1) illustrates the typical evaluation procedure. Model output files (restart and/or time-averaged NetCDFs) are loaded and standardized according to the YAML config. Metrics are computed, and—if a reference is provided—differences, MAE, and RMSE are calculated.

NEMO Spinup-Evaluation is often used alongside [spinup-forecast](https://github.com/m2lines/nemo-spinup-forecast), which automates the generation of machine learned spin-up states for NEMO/DINO models. Together, these tools provide a robust workflow for accelerating ocean spin-up.

<p align="center">
<img src="diagram.png" alt="NEMO flow" width="500"/>
<figcaption>Fig 1. Evaluation flow diagram illustrating the coupling to spinup-forecast, but nemo-spinup-evaluation can in theory be used to evaluate any ocean model, be it ML data driven, numerical or otherwise. </figcaption>
</p>


## Repository Layout

```
.
├── pyproject.toml                  Project metadata, dependencies, and build system
├── README.md                       Main project documentation (this file)
├── configs/                        Configuration files for variable/file mapping
│   └── DINO-setup.yaml             Example YAML config for DINO/NEMO variables
├── src/
│   └── nemo_spinup_evaluation/     Main Python package
│       ├── cli.py                  Command-line interface (CLI) entry point
│       ├── loader.py               Data loading and preprocessing utilities
│       ├── metrics_io.py           Output helpers (CSV writing, formatting)
│       ├── metrics.py              Metric calculation functions
│       └── utils.py                General utilities
├── tests/                          Test suite, test data, and data download scripts
│   └── get-data.sh                 Script to fetch test data from THREDDS
└── results/                        Default output directory for metrics CSVs
```


## Modes

NEMO Spinup-Evaluation supports three modes, controlled by the `--mode` argument:

### `restart` (Instantaneous Output)

- **Purpose:** Evaluate a single model state (snapshot) from a NEMO/DINO `restart.nc` file.
- **Input:** `restart.nc` (and `mesh_mask.nc`)
- **Use case:** Assess the physical realism or convergence of a single model state, e.g., after a spin-up or forecast.
- **Output:** `results/metrics_results_restart.csv` (or your chosen prefix)
- **Reference:** If `--ref-sim-path` is provided, computes diffs/stats vs. a reference restart file.

### `output` (Time-Averaged State)

- **Purpose:** Evaluate time-averaged or multi-time-step model output, typically from files like `grid_T_3D.nc`, `grid_U_3D.nc`, `grid_V_3D.nc`, `grid_T_2D.nc`.
- **Input:** Grid files as mapped in the config YAML (see below).
- **Use case:** Assess the mean state or variability over a period, or compare time-averaged fields between runs.
- **Output:** `results/metrics_results_grid.csv` (or your chosen prefix)
- **Reference:** If `--ref-sim-path` is provided, computes diffs/stats vs. a reference output set.

### `both`

- **Purpose:** Run both `restart` and `output` metric suites in one command.
- **Output:** Both CSVs as above.


## Config File

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

If the NetCDF mesh mask, restart or grid files use different variable names to these then use the `variable_map` to remap them. All loaded files are remapped using the specified map. If files use differing names for the same variable, list all names against each canonical field name, e.g.

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


## Grid Files

All 2D and 3D grid files specified as `output_variables` must be temporally aligned. This can be easily done using the [CDO tools](https://code.mpimet.mpg.de/projects/cdo) (Climate Data Operators), e.g. resampling a SSH 2D grid to yearly cadence:

```bash
cdo yearmean DINO_1m_grid_T.nc DINO_1y_grid_T.nc
```


## Output Files

Results are written as CSV files in the results directory, e.g.:
- `results/metrics_results_restart.csv`
- `results/metrics_results_grid.csv`

Each file contains metric values, and if a reference is provided, also includes:
- Reference metric values (prefixed with `ref_`)
- Differences (`diff_*`)

A separate file with MAE and RMSE statistics is also generated if a reference directory is provided.


## Adding New Metrics

Add new metric functions to `src/nemo_spinup_evaluation/metrics.py` and update the metric function lists in `cli.py` as needed.


## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details on how to contribute to this project, including:

- How to report bugs and suggest enhancements
- Development setup and coding guidelines
- Testing procedures
- How to add new metrics
- Release process

## Code of Conduct

This project follows a [Code of Conduct](CODE_OF_CONDUCT.md) to ensure a welcoming and inclusive community. All contributors are expected to follow these guidelines.

## Acknowledgements

This work builds on significant contributions by [Etienne Meunier](https://github.com/Etienne-Meunier), whose efforts on the [Metrics-Ocean](https://github.com/Etienne-Meunier/Metrics-Ocean) repository laid the foundation for several components used here.
