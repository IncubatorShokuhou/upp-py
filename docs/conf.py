"""Sphinx configuration for upp-py.

The layout follows the scikit-learn documentation: a short getting-started
section, a narrative user guide, an autosummary API, a sphinx-gallery, and
a validation page. The user guide is written in Chinese; docstrings stay
in English.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.abspath("../src"))

project = "upp-py"
author = "upp-py contributors"
release = "0.2.0"
version = "0.2"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "numpydoc",
    "sphinx_gallery.gen_gallery",
    "sphinx_copybutton",
]

autosummary_generate = True
autosummary_imported_members = False
numpydoc_show_class_members = False
numpydoc_class_members_toctree = False
numpydoc_xref_param_type = False
autodoc_typehints = "description"
autodoc_default_options = {"members": True, "undoc-members": False, "show-inheritance": False}

sphinx_gallery_conf = {
    "examples_dirs": "../examples",
    "gallery_dirs": "auto_examples",
    "filename_pattern": r"/plot_",
    "within_subsection_order": "FileNameSortKey",
    "download_all_examples": False,
    "plot_gallery": os.environ.get("UPP_PY_PLOT_GALLERY", "True"),
}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "xarray": ("https://docs.xarray.dev/en/stable/", None),
}

templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "validation.md",
    "validation_generated.rst",
]
html_theme = "pydata_sphinx_theme"
html_theme_options = {
    "navigation_depth": 3,
    "show_nav_level": 2,
    "collapse_navigation": False,
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
}
html_title = "upp-py"
html_static_path = ["_static"]
master_doc = "index"
pygments_style = "sphinx"
nitpicky = False
