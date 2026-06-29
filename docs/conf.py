"""Sphinx configuration for BioCompiler documentation."""

# -- Path setup --------------------------------------------------------------
import os
import sys

# Add the source directory to sys.path so autodoc can find the package
sys.path.insert(0, os.path.abspath("../src"))

# -- Project information -----------------------------------------------------

project = "BioCompiler"
copyright = "2024-2026, BioCompiler Team"
author = "BioCompiler Team"

# Read version from the package dynamically
try:
    import biocompiler

    release = biocompiler.__version__
    version = ".".join(release.split(".")[:2])  # e.g. "0.9"
except ImportError:
    release = "0.9.0"
    version = "0.9"

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
]

templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "api",
    "adr",
]

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

# -- Options for autodoc ----------------------------------------------------

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
    "member-order": "bysource",
}

autodoc_typehints = "description"
autodoc_typehints_format = "short"

# -- Options for Napoleon ----------------------------------------------------

napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True

# -- Options for intersphinx -------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

# -- Options for todo --------------------------------------------------------

todo_include_todos = True

# -- Options for autosummary -------------------------------------------------

autosummary_generate = True
