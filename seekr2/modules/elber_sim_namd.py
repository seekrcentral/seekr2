"""
elber/sim_namd.py

Base objects and routines for preparing NAMD simulations
to run for Elber milestoning calculations.
"""

import os
import datetime

from parmed import unit
import parmed

import seekr2.modules.common_base as base
import seekr2.modules.common_sim_namd as common_sim_namd
from seekr2.modules.common_sim_namd import add_string_buffer

dt_date = datetime.datetime.now()
SEEKR_ELBER_PROD_HEADER = "#" + "*"*64 + "\n" \
    + "#  Name SEEKR Elber Milestoning - Production Stage\n" \
    + "#  Date " + dt_date.strftime("%Y-%b-%d") + "\n" \
    + "#  Generated with the SEEKR2 automatic NAMD input file generator\n" \
    + "#" + "*"*64 + "\n"