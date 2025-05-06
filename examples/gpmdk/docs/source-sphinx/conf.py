# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://gitlab.lanl.gov/qmd-progress/exmaples/gpmdk

project = 'GPMDK'
copyright = 'See BSD3 licence'
author = 'GPMDK team'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.mathjax', 'sphinx.ext.autodoc', 'sphinx_mdinclude']
#extensions = ['sphinx.ext.mathjax', 'sphinx.ext.autodoc', 'sphinx_mdinclude', 'autoapi.extension', 'sphinxfortran.fortran_domain','sphinxfortran.fortran_autodoc']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
#fortran_ext = ['F90']
#fortran_src = ['../../src']

#autoapi_type = 'python'
#autoapi_dirs = ['../../src']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'press'
html_static_path = ['_static']

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown'
#    '.md': 'markdown',
}
