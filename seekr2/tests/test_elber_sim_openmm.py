"""
test_elber_sim_openmm.py
"""

import pytest
import os

from seekr2.modules import elber_sim_openmm
from seekr2.tests import make_test_model
""" # TODO: re-establish when Elber plugin is working
def test_mmvt_sim_openmm_amber(tmp_path):
    
    mymodel = make_test_model.make_test_model(tmp_path)
    myanchor = mymodel.anchors[1]
    output_file = os.path.join(tmp_path, "output.txt")
    my_sim_openmm = elber_sim_openmm.create_sim_openmm(mymodel, myanchor, 
                                                      output_file)
    return

def test_mmvt_sim_openmm_forcefield(tmp_path):
    
    mymodel = make_test_model.make_test_model(tmp_path, mode='forcefield')
    myanchor = mymodel.anchors[1]
    output_file = os.path.join(tmp_path, "output.txt")
    my_sim_openmm = elber_sim_openmm.create_sim_openmm(mymodel, myanchor, 
                                                      output_file)
    return
"""