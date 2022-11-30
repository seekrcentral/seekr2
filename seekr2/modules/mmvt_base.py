"""
mmvt_base.py

Base classes, objects, and constants used in multiple stages of the
MMVT calculations.
"""
import math
import os
import re

import numpy as np
from scipy.spatial import transform
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
# TEMPORARY until restructure finished
from seekr2.modules.mmvt_cvs.mmvt_cv_base import MMVT_collective_variable


            

# TODO: remove module





