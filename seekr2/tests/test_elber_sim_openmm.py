"""
test_elber_sim_openmm.py
"""

import pytest
import os

from seekr2.modules import elber_sim_openmm

def test_elber_sim_openmm_amber(tmp_path, tryp_ben_elber_model):
    myanchor = tryp_ben_elber_model.anchors[1]
    output_file = os.path.join(tmp_path, "output.txt")
    my_sim_openmm = elber_sim_openmm.create_sim_openmm(
        tryp_ben_elber_model, myanchor, output_file)
    return

@pytest.mark.skip
def test_elber_sim_openmm_forcefield(tmp_path, tryp_ben_elber_model):
    # TODO: un-skip this test if forcefield inputs ever become widely-used
    myanchor = tryp_ben_elber_model.anchors[1]
    output_file = os.path.join(tmp_path, "output.txt")
    my_sim_openmm = elber_sim_openmm.create_sim_openmm(
        tryp_ben_elber_model, myanchor, output_file)
    return