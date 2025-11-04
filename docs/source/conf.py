"""Sphinx configuration for Spinup-Evaluation documentation."""

import os
import sys

sys.path.insert(0, os.path.abspath("../../src"))  # assumes your code is in root

# -- Project information -----------------------------------------------------
project = "Spinup-Evaluation"
# ruff: noqa: A001
copyright = "2025, Surbhi Goel, Matt Archer"
author = "Surbhi Goel, Matt Archer"
release = "2025"

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
