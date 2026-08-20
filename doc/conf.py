# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0, os.path.abspath(".."))

project = 'PyComposites'
copyright = '2026, Joris Remmers'
author = 'Joris Remmers'
release = '2026.9'
version = '2026.9'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',    # Automatically include docstrings
    'sphinx.ext.viewcode',   # Add links to source code
    'sphinx.ext.napoleon',   # Support for NumPy/Google-style docstrings
    'sphinx_rtd_theme',      # Use the Read the Docs theme
    'sphinx.ext.mathjax',
    'myst_parser',           # Support Markdown documentation files
]

source_suffix = {
    '.md': 'markdown',
}

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

autodoc_default_flags = ['members', 'undoc-members', 'private-members']
add_module_names = False
autodoc_member_order = 'bysource'


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_logo = 'img/pycomposites_logo.png'
html_static_path = ['_static']
html_css_files = ['custom.css']

html_theme_options = {
    "navigation_depth": 3,   # controls depth in sidebar
}

html_context = {
    "default_mode": "sidebar",
}
