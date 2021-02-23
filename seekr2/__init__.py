"""
seekr2
Simulation-Enabled Estimation of Kinetic Rates - Version 2
"""

# Add imports here
from .seekr2 import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
