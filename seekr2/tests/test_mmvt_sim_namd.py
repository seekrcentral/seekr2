"""
test_mmvt_sim_namd.py
"""

import pytest
import os

from seekr2.modules import mmvt_sim_namd
from seekr2.tests import make_test_model

def test_mmvt_sim_namd_amber(tmp_path):
    
    mymodel = make_test_model.make_test_model(tmp_path, engine="namd")
    myanchor = mymodel.anchors[1]
    output_file = os.path.join(tmp_path, "output.txt")
    my_sim_openmm = mmvt_sim_namd.create_sim_namd(mymodel, myanchor, 
                                                      output_file)
    return
