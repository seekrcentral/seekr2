"""
test_mmvt_sim_namd.py
"""

import os

from seekr2.modules import mmvt_sim_namd

def test_mmvt_sim_namd_amber(tmp_path, host_guest_mmvt_model_namd):
    myanchor = host_guest_mmvt_model_namd.anchors[1]
    output_file = os.path.join(tmp_path, "output.txt")
    my_sim_openmm = mmvt_sim_namd.create_sim_namd(host_guest_mmvt_model_namd, 
                                                  myanchor, output_file)
    assert my_sim_openmm.colvars_config is not None
    assert my_sim_openmm.seekr_namd_settings is not None
    return
