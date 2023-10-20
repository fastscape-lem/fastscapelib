# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import os
import subprocess

from matplotlib.font_manager import FontManager

# -- Doxygen build on RTD ----------------------------------------------------

on_rtd = os.environ.get("READTHEDOCS", None) == "True"

if on_rtd:
    subprocess.call("cd ..; doxygen", shell=True)

# -- Project information -----------------------------------------------------

project = "Fastscapelib"
copyright = "since 2018, Fastscapelib developers"
author = "Benoit Bovy"
# The short X.Y version
version = "0.2.1"
# The full version, including alpha/beta/rc tags
release = version

print("building doc for fastscapelib version ", version)

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.extlinks",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "breathe",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",
    "sphinx_remove_toctrees",
    "sphinxcontrib.bibtex",
    "sphinxcontrib.mermaid",
    "sphinx_design",
    "myst_nb",
]

extlinks = {
    "issue": ("https://github.com/fastscape-lem/fastscapelib/issues/%s", "#%s"),
    "pull": ("https://github.com/fastscape-lem/fastscapelib/pull/%s", "#%s"),
}

templates_path = ["_templates"]

source_suffix = [".rst", ".md"]

master_doc = "index"

language = "en"

exclude_patterns = [
    "**.ipynb_checkpoints",
    "build/**.ipynb",
    "Thumbs.db",
    ".DS_Store",
]

bibtex_bibfiles = ["references.bib"]

# pygments_style = "sphinx"
highlight_language = "c++"
myst_enable_extensions = ["colon_fence", "substitution"]

# not working yet (https://github.com/sphinx-doc/sphinx/issues/5379)
autodoc_member_order = "bysource"

napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_rtype = False
napoleon_preprocess_types = True
napoleon_type_aliases = {
    "array_like": ":term:`array_like`",
    "array-like": ":term:`array-like <array_like>`",
}
typehints_defaults = "comma"
typehints_use_rtype = False

remove_from_toctrees = ["api_python/_api_generated/*"]

# -- Myst-NB config ----------------------------------------------------------

# This is to mitigate errors on CI VMs, where you can get the message:
# Matplotlib is building the font cache" in output notebooks
FontManager()

nb_kernel_rgx_aliases = {".*fastscapelib.*": "python3"}

# -- Options for HTML output -------------------------------------------------

html_theme = "furo"
html_title = ""

html_theme_options = {
    # Furo Theme
    "sidebar_hide_name": True,
    "light_logo": "fastscapelib_logo.svg",
    "dark_logo": "fastscapelib_logo.svg",
    # "announcement": "Fastscapelib is in <em>active development</em> (refactoring) and will be released soon, stay tuned!",
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/fastscape-lem/fastscapelib",
            "html": """
                <svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 16 16">
                    <path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path>
                </svg>
            """,
            "class": "",
        },
    ],
}
html_static_path = ["_static"]
html_css_files = ["style.css"]
html_favicon = "_static/favicon.ico"
htmlhelp_basename = "fastscapelibdoc"

# -- Breathe configuration -------------------------------------------------

breathe_projects = {"fastscapelib": "../_xml"}
breathe_default_project = "fastscapelib"
breathe_show_include = True

# -- Options for intersphinx extension ---------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "xtensor": ("https://xtensor.readthedocs.io/en/latest/", None),
    "xtensor-python": ("https://xtensor-python.readthedocs.io/en/latest/", None),
}
