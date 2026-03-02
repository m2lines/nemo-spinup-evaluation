# AGENTS.md for NEMO Spinup Evaluation

## Project Overview
Evaluates NEMO ocean model spinup simulations by comparing generated states with reference data using various metrics.

## Development Environment

### Setup
```bash
# Create virtual environment
python -m venv .venv
source .venv/bin/activate

# Install package
pip install -e .
```

### Key Project Directories
- `./src/spinup_evaluation/`: Main Python package
- `./tests/`: Test suites
- `./configs/`: Configuration files (DINO-setup.yaml, gen-setup.yaml)

## Commands

### Testing
```bash
# Run all tests
python -m pytest tests/

# Run specific test
python -m pytest tests/test_evaluation.py
```

### Evaluation
```bash
# Run main evaluation
python -m spinup_evaluation.cli --config configs/DINO-setup.yaml

# Generate metrics
python -m spinup_evaluation.metrics --input ./data/ --output ./results/
```

### Data Processing
```bash
# Standardize inputs
python -m spinup_evaluation.standardise_inputs --config configs/gen-setup.yaml

# Load data
python -m spinup_evaluation.loader --path ./data/ --output ./processed/
```

## Project Structure

```
.
├── .pre-commit-config.yaml      # Pre-commit hooks
├── pyproject.toml               # Project configuration
├── src/
│   └── spinup_evaluation/
│       ├── __init__.py          # Package entry
│       ├── cli.py               # CLI interface
│       ├── metrics.py           # Metric calculations
│       ├── loader.py            # Data loading
│       ├── metrics_io.py        # Metrics I/O
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
- Target `dev` for features, `main` for bugfixes
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
- Runs on push to `main` and `dev` branches
- Includes: linting, testing, and build steps

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
- 100% line coverage for new features

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
