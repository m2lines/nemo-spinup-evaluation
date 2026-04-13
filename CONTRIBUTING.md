# Contributing to NEMO Spinup-Evaluation

Thank you for your interest in contributing to NEMO Spinup-Evaluation! This project provides benchmarking tools for NEMO/DINO ocean model spin-up, and we welcome contributions from the community.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Data](#data)
- [Development Setup](#development-setup)
- [Making Changes](#making-changes)
- [Testing](#testing)
- [Adding New Metrics](#adding-new-metrics)
- [Documentation](#documentation)
- [Submitting a Pull Request](#submitting-a-pull-request)
- [Reporting Bugs](#reporting-bugs)
- [Suggesting Enhancements](#suggesting-enhancements)
- [Release Process](#release-process)
- [Getting Help](#getting-help)

## Code of Conduct

All participants are expected to follow our [Code of Conduct](CODE_OF_CONDUCT.md). Please read it before contributing.

## Getting Started

The project lives at [github.com/m2lines/nemo-spinup-evaluation](https://github.com/m2lines/nemo-spinup-evaluation). It is a Python package (Python >= 3.10) that provides a CLI and Python API for evaluating ocean model spin-up performance. The package is structured as:

```
src/nemo_spinup_evaluation/
    cli.py                  # Command-line entry point
    loader.py               # Data loading and preprocessing
    metrics.py              # Metric calculation functions
    metrics_io.py           # CSV output helpers
    standardise_inputs.py   # Input standardisation
    utils.py                # General utilities
configs/
    DINO-setup.yaml         # Variable-to-file mapping config
tests/                      # Test suite
docs/                       # Sphinx documentation (hosted on Read the Docs)
```

## Data

All evaluation datasets are hosted on Zenodo:

**[https://zenodo.org/records/19474414](https://zenodo.org/records/19474414)**

Available datasets:

| Dataset | File | Size | Description |
|---------|------|------|-------------|
| 50-year baseline | `50.zip` | ~552 MB | Standard evaluation dataset |
| 200-year extended | `200.zip` | ~2.1 GB | Longer spin-up for comprehensive testing |
| Restart files | `restart.zip` | ~6.5 GB | Restart snapshots for restart-mode evaluation |

To download for local development and testing:

```bash
# Download the 50-year baseline (recommended starting point)
wget https://zenodo.org/records/19474414/files/50.zip
unzip 50.zip
```

Test data (subsampled NetCDFs for the test suite) is fetched separately via the script described in [Testing](#testing).

## Development Setup

### Prerequisites

- Python >= 3.10
- Git

### Installation

```bash
# Fork and clone the repository
git clone https://github.com/<your-username>/nemo-spinup-evaluation.git
cd nemo-spinup-evaluation

# Create and activate a virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install in editable mode with dev dependencies
pip install -e .[dev]

# Set up pre-commit hooks
pre-commit install
```

### Development Tools

| Tool | Purpose | Command |
|------|---------|---------|
| `ruff` | Linting and formatting | `ruff check .` / `ruff format .` |
| `mypy` | Type checking (not enforced in CI) | `mypy src/` |
| `pytest` | Testing | `pytest tests/` |
| `pre-commit` | Git hooks (runs ruff automatically) | Runs on `git commit` |

## Making Changes

1. **Create a branch** from `main` with a descriptive name:
   ```bash
   git checkout -b add-kinetic-energy-metric
   ```

2. **Follow the coding style:**
   - [PEP 8](https://peps.python.org/pep-0008/) via `ruff`
   - [NumPy docstring format](https://numpydoc.readthedocs.io/en/latest/format.html)
   - Type hints are encouraged
   - Run `ruff check .` and `ruff format .` before committing (or let pre-commit handle it)
   - Type hints are encouraged when adding new code, but `mypy` is not currently enforced in CI and some existing type issues remain

3. **Write tests** for any new or changed functionality.

4. **Keep commits focused.** Use clear, descriptive commit messages and reference issue numbers where applicable (e.g., "Fixes #42"). We don't enforce [Conventional Commits](https://www.conventionalcommits.org/) format, but prefixing with a category (e.g., `fix:`, `feat:`, `docs:`, `test:`) is welcome if it comes naturally.

## Testing

### Fetching Test Data

The test suite uses subsampled DINO NetCDF files. Download them with:

```bash
./tests/get-data.sh
```

### Running Tests

```bash
# Run all tests
pytest tests/

# Run with coverage
pytest --cov=src tests/

# Run a specific test file
pytest tests/test_evaluation.py -v
```

### CI

Tests and linting run automatically on pull requests via GitHub Actions. The CI pipeline:
1. Runs `ruff format --diff .` and `ruff check --diff .`
2. Installs the package and downloads test data
3. Runs `pytest -v tests/`

Your PR must pass CI before it can be merged.

## Adding New Metrics

New metrics are added in `src/nemo_spinup_evaluation/metrics.py` and registered in `src/nemo_spinup_evaluation/cli.py`.

1. **Define your metric function** in `metrics.py`:

   ```python
   def kinetic_energy(data: xr.DataArray) -> float:
       """
       Compute mean kinetic energy.

       Parameters
       ----------
       data : xr.DataArray
           Velocity field data.

       Returns
       -------
       float
           Mean kinetic energy value.
       """
       return float((0.5 * data**2).mean())
   ```

2. **Register the metric** by adding it to the appropriate metric list in `cli.py`.

3. **Add tests** in `tests/` to verify the metric produces correct values.

4. **Update docs** if the metric introduces new concepts or config requirements.

## Documentation

Documentation is built with Sphinx and hosted on [Read the Docs](https://nemo-spinup-evaluation.readthedocs.io/).

### Building Locally

```bash
pip install .[docs]
cd docs
make html
# Open build/html/index.html in your browser
```

### Writing Docs

- Source files are in `docs/source/` (Markdown via MyST parser)
- Use NumPy docstring format in code -- Sphinx autodoc picks these up
- Update relevant doc pages when adding features

## Submitting a Pull Request

1. Push your branch to your fork
2. Open a PR against `main` on [m2lines/nemo-spinup-evaluation](https://github.com/m2lines/nemo-spinup-evaluation)
3. Fill in the PR description with:
   - What changed and why
   - How to test the changes
   - Any related issues (e.g., "Closes #15")
4. Ensure CI passes (linting + tests)
5. A maintainer will review your PR. Address any feedback, then it will be merged.

## Reporting Bugs

Before filing a bug, check [existing issues](https://github.com/m2lines/nemo-spinup-evaluation/issues) to avoid duplicates.

When reporting, please include:

- Python version and OS
- The exact command you ran
- The full error message / traceback
- The config file and dataset you used (if applicable)
- Steps to reproduce

## Suggesting Enhancements

Open an issue describing:

- The proposed feature and its use case
- Why it would benefit the project
- Any implementation ideas you have

## Release Process

The project follows [Semantic Versioning](https://semver.org/):

- **MAJOR**: breaking changes
- **MINOR**: new features (backwards-compatible)
- **PATCH**: bug fixes

To create a release:

1. Update the version in `pyproject.toml`
2. Update `CHANGELOG.md`
3. Create a GitHub release with notes

## Getting Help

- Read the [documentation](https://nemo-spinup-evaluation.readthedocs.io/)
- Open a [GitHub issue](https://github.com/m2lines/nemo-spinup-evaluation/issues)
- Contact the maintainers (see `pyproject.toml` for emails)

## License

By contributing, you agree that your contributions will be licensed under the [MIT License](LICENSE).
