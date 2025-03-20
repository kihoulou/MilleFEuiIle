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
sys.path.insert(0, os.path.abspath('../'))
# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.autodoc"]

# templates_path = ['_templates']
# exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

master_doc = "contents"


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_favicon = "_static/mf.png"
html_theme_options = {
    "home_page_in_toc": True,
    "toc_title": "On this page:",
    "logo": {
        # "text": "MilleFEuiIle Home",
        "image_light": "_static/mf.png",
        "image_dark": "_static/mf.png",
        "link": "index"
    }

}

# html_static_path = ['_static']
# html_logo = "_static/mf.png"
