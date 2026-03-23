"""Sphinx configuration for Spinup-Evaluation documentation."""

import os
import sys
from importlib.metadata import version as _get_version

sys.path.insert(0, os.path.abspath("../../src"))  # assumes your code is in root

# -- Project information -----------------------------------------------------
project = "nemo-spinup-evaluation"
# ruff: noqa: A001
copyright = "2026, Matt Archer, Surbhi Goel"
author = "Matt Archer, Surbhi Goel"
release = _get_version("nemo-spinup-evaluation")

# -- General configuration ---------------------------------------------------
extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
]

autosummary_generate = True
templates_path = ["_templates"]
exclude_patterns = []

language = "English"

# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
