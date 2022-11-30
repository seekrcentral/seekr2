"""
mmvt_cv.py

Define any type of collective variable (or milestone shape) that might
be used in an MMVT calculation.
"""
import os
import shutil

import seekr2.modules.common_base as base
import seekr2.modules.mmvt_base as mmvt_base
import seekr2.modules.common_cv as common_cv

""" # MOVED to mmvt_cvs.mmvt_spherical_cv
def make_mmvt_spherical_cv_object(spherical_cv_input, index):
    ""
    Create a SphericalCV object to be placed into the Model.
    ""
    group1 = base.parse_xml_list(spherical_cv_input.group1)
    group2 = base.parse_xml_list(spherical_cv_input.group2)
    groups = [group1, group2]
    cv = mmvt_base.MMVT_spherical_CV(index, groups)
    return cv
"""

# TODO: remove module