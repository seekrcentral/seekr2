"""
test_elber_sim_openmm.py
"""

import pytest
import os

from seekr2.modules import elber_sim_openmm

def test_elber_sim_openmm_amber(tmp_path, host_guest_elber_model):
    myanchor = host_guest_elber_model.anchors[1]
    host_guest_elber_model.openmm_settings.cuda_platform_settings = None
    host_guest_elber_model.openmm_settings.reference_platform = True
    output_file = os.path.join(tmp_path, "output.txt")
    my_sim_openmm = elber_sim_openmm.create_sim_openmm(
        host_guest_elber_model, myanchor, output_file)
    return

@pytest.mark.skip
def test_elber_sim_openmm_forcefield(tmp_path, host_guest_elber_model):
    # TODO: un-skip this test if forcefield inputs ever become widely-used
    myanchor = host_guest_elber_model.anchors[1]
    host_guest_elber_model.openmm_settings.cuda_platform_settings = None
    host_guest_elber_model.openmm_settings.reference_platform = True
    output_file = os.path.join(tmp_path, "output.txt")
    my_sim_openmm = elber_sim_openmm.create_sim_openmm(
        host_guest_elber_model, myanchor, output_file)
    return