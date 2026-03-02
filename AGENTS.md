# AGENTS.md for NEMO Spinup Evaluation

## Project Overview
Provides a CLI and Python API for benchmarking the spinup and restart performance of NEMO/DINO ocean models and machine learning emulators. Supports both single-run evaluation and optional comparison against a reference simulation, outputting metrics and difference statistics.

## Development Environment

### Setup
```bash
# Create virtual environment
python -m venv .venv
source .venv/bin/activate

# Install package with dev dependencies
pip install -e .[dev]
```

### Key Project Directories
- `./src/spinup_evaluation/`: Main Python package
- `./tests/`: Test suites
- `./configs/`: Configuration files (DINO-setup.yaml, gen-setup.yaml)

## Commands

### Testing
```bash
# Download test data (required before first run)
bash tests/get-data.sh

# Run all tests
python -m pytest tests/

# Run specific test
python -m pytest tests/test_evaluation.py
```

### Evaluation
```bash
# Run main evaluation (--sim-path and/or --ref-sim-path are required)
spinup-eval --config configs/DINO-setup.yaml --sim-path /path/to/simulation
spinup-eval --config configs/DINO-setup.yaml --sim-path /path/to/simulation --ref-sim-path /path/to/reference

# Equivalent forms
python -m spinup_evaluation --config configs/DINO-setup.yaml --sim-path /path/to/simulation
python -m spinup_evaluation.cli --config configs/DINO-setup.yaml --sim-path /path/to/simulation
```

### CLI Options
- `--mode`: `output`, `restart`, or `both` (default: `both`) — controls which file types are processed
- `--results-dir`: directory for CSV output (default: `results/`)
- `--result-file-prefix`: prefix for output CSV filenames (default: `metrics_results`)

## Project Structure

```
.
├── .pre-commit-config.yaml      # Pre-commit hooks
├── pyproject.toml               # Project configuration
├── src/
│   └── spinup_evaluation/
│       ├── __init__.py          # Package entry
│       ├── __main__.py          # Module entry point (python -m spinup_evaluation)
│       ├── cli.py               # CLI interface
│       ├── metrics.py           # Metric calculations
│       ├── loader.py            # Data loading
│       ├── metrics_io.py        # Metrics I/O
│       ├── standardise_inputs.py # Variable alias lookup and standardisation
│       └── utils.py             # Utilities
├── tests/
│   ├── test_evaluation.py      # Main tests
│   └── data/                   # Test data
└── configs/
    ├── DINO-setup.yaml         # Main config
    └── gen-setup.yaml          # Generation config
```

## Code Style & Quality

### Pre-commit Hooks
```bash
# Install pre-commit hooks
pre-commit install

# Run all hooks
pre-commit run --all-files

# Update hooks
pre-commit autoupdate
```

### Active Hooks
- `trailing-whitespace`: Remove trailing whitespace
- `end-of-file-fixer`: Ensure files end with newline
- `check-yaml`: Validate YAML files
- `check-added-large-files`: Prevent large file commits
- `ruff-format`: Auto-format Python code
- `ruff-check`: Lint Python code

### Formatting & Linting
```bash
# Format code with ruff
ruff format src/

# Lint code with ruff
ruff check src/

# Check imports
ruff check --select I src/
```

### Style Guide
- **Naming**: `snake_case` for functions/variables, `PascalCase` for classes
- **Type hints**: Required for all public functions (PEP 484)
- **Line length**: 88 characters maximum
- **Docstrings**: Google-style for public functions
- **Imports**: Grouped (stdlib, third-party, local) with blank lines
- **Quotes**: Single quotes for strings, double quotes for docstrings

### Configuration Files
- `.pre-commit-config.yaml`: Pre-commit hook configuration
- `pyproject.toml`: Ruff and project configuration
- `.ruff_cache/`: Ruff cache directory (auto-generated)

## Quality Checks
```bash
# Run all quality checks
pre-commit run --all-files

# Check specific file
ruff check src/spinup_evaluation/metrics.py

# Auto-fix formatting issues
ruff format --fix src/

# Run only ruff checks
pre-commit run ruff-format --all-files
```

## Git Workflow

### Commits
```bash
git commit -m "<type>(<scope>): <description>"
# feat(metrics): add RMSE calculation
# fix(loader): handle missing files
```

### PRs
- Target `main` for all changes
- Include tests for new functionality
- PR title: `[<component>] <description>`

## Testing
- All new features require tests
- Use pytest framework
- Cover edge cases

## Boundaries

### Agents CAN
- Read/analyze code
- Run tests/evaluations
- Generate metrics
- Create test data

### Agents CANNOT
- Modify production configs
- Delete data files
- Push to main branch
- Execute arbitrary system commands

## CI/CD Pipeline

### GitHub Actions
- Workflow file: `.github/workflows/ci-eval.yml`
- Runs on push to `main` branch and on pull requests targeting `main`
- Includes: linting and testing steps

### Pipeline Steps
```bash
# Check what CI runs
cat .github/workflows/ci-eval.yml

# Run CI checks locally
pre-commit run --all-files
python -m pytest tests/
```

### Quality Gates
- All pre-commit hooks must pass
- All tests must pass
- No linting errors allowed

## Debugging

### Common Issues
```bash
# Check imports
python -c "import spinup_evaluation; print('OK')"

# List test files
find tests/ -name "*.py"

# Check pre-commit status
pre-commit run --all-files
```

### Debug Commands
```bash
# Show ruff violations
ruff check src/ --show-fixes

# Show import issues
ruff check --select I src/

# Check specific test
python -m pytest tests/test_evaluation.py::test_specific_function -v
```
