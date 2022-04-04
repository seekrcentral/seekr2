"""
test_mmvt_sim_openmm.py
"""

import pytest
import os

from seekr2.modules import mmvt_sim_openmm

def test_mmvt_sim_openmm_amber(tmp_path, host_guest_mmvt_model):
    
    myanchor = host_guest_mmvt_model.anchors[1]
    output_file = os.path.join(tmp_path, "output.txt")
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(
        host_guest_mmvt_model, myanchor, output_file)
    assert my_sim_openmm.system is not None
    assert my_sim_openmm.integrator is not None
    assert my_sim_openmm.simulation is not None
    return

# No forcefield inputs currently working correctly
@pytest.mark.skip()
def test_mmvt_sim_openmm_forcefield(tmp_path, host_guest_mmvt_model_forcefield):
    
    myanchor = host_guest_mmvt_model_forcefield.anchors[1]
    output_file = os.path.join(tmp_path, "output.txt")
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(
        host_guest_mmvt_model_forcefield, myanchor, output_file)
    assert my_sim_openmm.system is not None
    assert my_sim_openmm.integrator is not None
    assert my_sim_openmm.simulation is not None
    return