from sphinx_mdolab_theme.config import *

# -- Path setup --------------------------------------------------------------
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
from unittest.mock import MagicMock

sys.path.insert(0, os.path.abspath("../"))

extensions.extend(["sphinxarg.ext"])

# -- Project information -----------------------------------------------------

project = "CGNS Utilities"


# We have to mock some stuff
# and sphinxarg.ext does NOT use the autodoc mock feature

for mod in ["cgnsutilities.libcgns_utils", "numpy"]:
    sys.modules[mod] = MagicMock()
