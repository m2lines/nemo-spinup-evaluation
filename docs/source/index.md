# NEMO Spinup-Evaluation

`nemo-spinup-evaluation` provides a command-line tool and Python API for benchmarking
the spin-up and restart performance of NEMO/DINO ocean models and machine learning
emulators. It supports both single-run and comparison (reference) evaluation, and
outputs detailed metrics and difference statistics.

It is often used alongside
[nemo-spinup-forecast](https://github.com/m2lines/nemo-spinup-forecast), which
automates the generation of machine-learned spin-up states for NEMO/DINO models.

## Features

- **Flexible CLI**: Evaluate restart and/or output files, with or without a reference simulation.
- **Configurable**: Uses a YAML config file to map variables to files.
- **Comparison mode**: Computes diffs, MAE, and RMSE between a simulation and a reference.
- **Extensible**: Add new metrics by editing `src/nemo_spinup_evaluation/metrics.py`.
- **Test suite**: Integration and regression tests using real and subsampled NetCDF data.

## Installation

Requires Python ≥ 3.10.

```sh
git clone https://github.com/m2lines/nemo-spinup-evaluation.git
cd nemo-spinup-evaluation
pip install .
```

For development installs:

```sh
pip install -e .[dev]
pre-commit install
```

## Quick Start

```sh
nemo-spinup-evaluation --sim-path <simulation_dir> --config configs/DINO-setup.yaml
```

To compare against a reference simulation:

```sh
nemo-spinup-evaluation --sim-path <simulation_dir> \
                       --ref-sim-path <reference_dir> \
                       --config configs/DINO-setup.yaml
```

## Acknowledgements

This work builds on contributions by
[Etienne Meunier](https://github.com/Etienne-Meunier), whose work on
[Metrics-Ocean](https://github.com/Etienne-Meunier/Metrics-Ocean) laid the
foundation for several components used here.
