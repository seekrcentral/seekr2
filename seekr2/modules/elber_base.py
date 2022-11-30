"""
elber_base.py

Base classes, objects, and constants used in multiple stages of the
Elber milestoning calculations.
"""

import math

import numpy as np
from parmed import unit
import mdtraj

try:
    import openmm
except ImportError:
    import simtk.openmm as openmm
    
try:
    import openmm.app as openmm_app
except ImportError:
    import simtk.openmm.app as openmm_app

from abserdes import Serializer

import seekr2.modules.common_base as base






# TODO: remove this module