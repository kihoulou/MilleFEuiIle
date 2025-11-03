# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'MilleFEuiIle'
copyright = '2025, Martin Kihoulou'
author = 'Martin Kihoulou'

import os 
import sys 
from sphinxawesome_theme.postprocess import Icons

sys.path.insert(0, os.path.abspath('../'))
# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.autodoc"]

# templates_path = ['_templates']
# exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

master_doc = "index"

# # Prevents the TOC to order alphabetically
# autodoc_member_order = 'bysource'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinxawesome_theme'
html_permalinks_icon = Icons.permalinks_icon

# html_theme = 'conestack'
# html_theme = 'sphinx_book_theme'
# html_theme =  "sphinx_rtd_theme"
html_favicon = "_static/mf.png"

# For custom css changes
html_static_path = ["_static"]
html_css_files = ["custom.css"]

html_theme_options = {
    "logo_light": "_static/mf.png",
    "logo_dark": "_static/mf.png",
    "show_breadcrumbs": True
}
