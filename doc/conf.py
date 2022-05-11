from sphinx_mdolab_theme.config import *

# -- Path setup --------------------------------------------------------------
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath("../"))

extensions.extend(["sphinxarg.ext"])

# -- Project information -----------------------------------------------------

project = "CGNS Utilities"
