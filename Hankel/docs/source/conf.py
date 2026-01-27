# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys


sys.path.insert(0, os.path.abspath("../.."))  # from docs/source to project root

project = 'Hankel_docs'
copyright = '2026, Jan Vábek'
author = 'Jan Vábek'
release = '0.9'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# extensions = []

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",   # supports Google/NumPy style docstrings
    "sphinx.ext.viewcode",   # adds source links
    "sphinx.ext.autosummary" # optional, but handy
]

autosummary_generate = True

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']
