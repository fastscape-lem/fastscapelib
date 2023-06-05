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
from __future__ import print_function

import os
import platform
import re
import subprocess

# -- Doxygen build on RTD ----------------------------------------------------

on_rtd = os.environ.get("READTHEDOCS", None) == "True"

if on_rtd:
    subprocess.call("cd ..; doxygen", shell=True)


# -- Project information -----------------------------------------------------

project = "Fastscapelib"
copyright = "2018, Benoit Bovy"
author = "Benoit Bovy"


# -- Version -----------------------------------------------------------------


def check_cmake_version():
    try:
        out = subprocess.check_output(["cmake", "--version"])
    except OSError:
        raise RuntimeError("CMake must be installed to build this project")

    if platform.system() == "Windows":
        cmake_version = LooseVersion(
            re.search(r"version\s*([\d.]+)", out.decode()).group(1)
        )
        if cmake_version < "3.1.0":
            raise RuntimeError("CMake >= 3.1.0 is required on Windows")


def get_version():
    """Get version string using git, formatted according to pep440
    (call cmake script).

    """
    check_cmake_version()

    cmake_script = os.path.join("cmake", "PrintVersion.cmake")
    out = subprocess.check_output(
        ["cmake", "-DINCLUDEDIR=include", "-P", cmake_script],
        cwd=os.path.join(os.pardir, os.pardir),
        stderr=subprocess.STDOUT,
    )
    version_str = out.decode().strip()
    return version_str


_version_str = get_version()
print("building doc for fastscapelib version ", _version_str)

# The short X.Y version
version = _version_str
# The full version, including alpha/beta/rc tags
release = _version_str


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
    "myst_nb",
]

extlinks = {
    "issue": ("https://github.com/fastscape-lem/fastscapelib/issues/%s", "GH"),
    "pull": ("https://github.com/fastscape-lem/fastscapelib/pull/%s", "PR"),
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

pygments_style = "sphinx"
highlight_language = "c++"

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

# -- Options for HTML output -------------------------------------------------

html_theme = "furo"
html_title = ""

html_theme_options = {
    # Sphinx book theme
    # repository_url="https://github.com/fastscape-lem/fastscapelib",
    # repository_branch="main",
    # path_to_docs="doc",
    # use_edit_page_button=True,
    # use_repository_button=True,
    # use_issues_button=True,
    # home_page_in_toc=False,
    # extra_navbar="",
    # navbar_footer_text="",
    # Furo Theme
    "sidebar_hide_name": True,
    "light_logo": "fastscapelib_logo.svg",
    "dark_logo": "fastscapelib_logo.svg",
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
