"""
test_mmvt_sim_openmm.py
"""

import pytest
import os

from seekr2.modules import mmvt_sim_openmm
from seekr2.tests import make_test_model

def test_mmvt_sim_openmm_amber(tmp_path):
    
    mymodel = make_test_model.make_test_model(tmp_path)
    sim_openmm_factory = mmvt_sim_openmm.MMVT_sim_openmm_factory()
    myanchor = mymodel.anchors[1]
    output_file = os.path.join(tmp_path, "output.txt")
    my_sim_openmm = sim_openmm_factory.create_sim_openmm(mymodel, myanchor, 
                                                      output_file)
    return

def test_mmvt_sim_openmm_forcefield(tmp_path):
    
    mymodel = make_test_model.make_test_model(tmp_path, mode='forcefield')
    sim_openmm_factory = mmvt_sim_openmm.MMVT_sim_openmm_factory()
    myanchor = mymodel.anchors[1]
    output_file = os.path.join(tmp_path, "output.txt")
    my_sim_openmm = sim_openmm_factory.create_sim_openmm(mymodel, myanchor, 
                                                      output_file)
    return