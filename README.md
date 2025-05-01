# Spinup-Evaluation

This repository contains code for benchmarking the machine learning spin-up of ocean models. It is designed to pair with the `Spinup-Forecast` repository, which provides the machine learning models for spin-up. The goal of this evaluation is to assess the performance of the spin-up process in terms of stability and convergence.

The evaluation is performed using the `main.py` script, which calls a set of metrics defined in the `metrics.py` file. 

The results are saved in a JSON file, which can be easily parsed and analyzed.

The API is as follows:

- `main.py`: The main script to run the evaluation.
analysis.
- `metrics.py`: Contains the definitions of the metrics used for evaluation.
- `results.json`: The output file containing the evaluation results.
analysis.
- `utils.py`: Contains utility functions for data processing and visualization.

`main.py` is the entry point for the evaluation process. It takes the following command-line arguments:
- `--input_dir`: The directory containing the spin-up data files.
- `--output_file`: The name of the output JSON file where the results will be saved.


## Installation
To install Spinup-Evaluation, clone the repository and create a virtual environment:
```bash
git clone
cd Spinup-Evaluation
python -m venv venv
source venv/bin/activate  # On Windows use `venv\Scripts\activate`
```

Then, install the required packages:

```bash
pip install -e .
```
For a development install, some further steps are recommended:

```sh
cd Spinup-Evaluation

# Install optional dev dependencies
pip install -e .[dev]

# Configure pre-commit hooks
pre-commit install
```

## Usage

These metrics were developed to assess the DINO configuration of NEMO, but they can be used for any spin-up ocean model, however, it will be necessary to add new metrics to the `metrics.py` file. Please see the section [Adding New Metrics](#adding-new-metrics) for more information.

Spinup-Evaluation expects that evaluation proceeds following the prediction of a new spin-up state using [Spinup-Forecast](https://github.com/m2lines/Spinup-Forecast). 
To run the evaluation on NEMO/DINO it is necessary to provide the following:

* `--predictions` : The path to directory containing the new `pred_[variable].npy` spin-up states from `Spinup-Forecast`.
    - pred_so.npy
    - pred_thetao.npy
    - pred_zos.npy
* `--mesh-mask` : Path to the `mesh_mask.nc` file. This file contains the grid information for the model.
* [Optional] The path to a reference spin-up state. This is used to compare the new spin-up state against a known good state. If not provided, the evaluation will only assess the new spin-up state.
* [Optional] The restart file from the NEMO/DINO model. 
    i.e. ``

![alt text](image.png)

## Adding New Metrics
To add new metrics to the evaluation, modify the `metrics.py` file. More to follow...

## Testing

Tests are provided in the `tests` directory. To run the tests, use the following command:

```bash
pytest tests/
```

## Restarting NEMO/DINO

To use the evaluation in a restart file, you can use the `--restart` flag when running the `main.py` script. This will save the evaluation results in a format that can be used for restarting the model.

See this Github Gist for more information on how to create a compatible restart file: https://gist.github.com/yourusername/your-gist-id. 
