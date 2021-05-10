"""
test_runner_openmm.py

test runner_openmm.py script(s)
"""

import pytest
import os
import glob

import seekr2.modules.common_base as base
import seekr2.modules.mmvt_base as mmvt_base
import seekr2.modules.mmvt_sim_openmm as mmvt_sim_openmm
import seekr2.modules.runner_openmm as runner_openmm
import seekr2.tests.make_test_model as make_test_model

def test_Runner_openmm_default(tmp_path):
    if not os.path.exists(tmp_path):
        os.mkdir(tmp_path)
    mymodel = make_test_model.make_test_model(tmp_path)
    #sim_openmm_factory = sim_openmm.Sim_openmm_factory()
    myanchor = mymodel.anchors[1]
    mmvt_output_filename = os.path.join(
                    tmp_path, myanchor.name, "prod", 
                    "%s%d.%s" % (mmvt_base.OPENMMVT_BASENAME, 1, 
                                         mmvt_base.OPENMMVT_EXTENSION))
    runner = runner_openmm.Runner_openmm(mymodel, myanchor)
    default_output_file, state_file_prefix, restart_index = runner.prepare()
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(mymodel, myanchor, 
                                                      mmvt_output_filename)
    
    runner.run(my_sim_openmm, False)
    
    assert os.path.exists(mmvt_output_filename)
    lines = 0
    with open(mmvt_output_filename,"r") as f:
        for line in f:
            lines += 1
    #assert lines > 0

def test_Runner_openmm_other_settings(tmp_path):
    if not os.path.exists(tmp_path):
        os.mkdir(tmp_path)
    mymodel = make_test_model.make_test_model(tmp_path)
    #sim_openmm_factory = sim_openmm.Sim_openmm_factory()
    myanchor = mymodel.anchors[1]
    mmvt_output_filename = os.path.join(
                    tmp_path, myanchor.name, "prod", 
                    "%s%d.%s" % (mmvt_base.OPENMMVT_BASENAME, 1, 
                                         mmvt_base.OPENMMVT_EXTENSION))
    loading_state_filename = os.path.join(tmp_path, "start.state")
    runner = runner_openmm.Runner_openmm(mymodel, myanchor)
    default_output_file, state_file_prefix, restart_index = runner.prepare(
        False, True)
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(mymodel, myanchor, 
                                                      mmvt_output_filename,
                                                      state_file_prefix)
    
    my_sim_openmm.simulation.saveState(loading_state_filename)
    
    runner.run(my_sim_openmm, False, load_state_file=loading_state_filename)
    assert my_sim_openmm.integrator.getOutputFileName() == mmvt_output_filename
    lines = 0
    assert os.path.exists(mmvt_output_filename)
    assert os.path.exists(loading_state_filename)
    #print("glob.glob(state_file_prefix+'*')", glob.glob(state_file_prefix+'*'))
    #assert len(glob.glob(state_file_prefix+'*')) > 0
    with open(mmvt_output_filename,"r") as f:
        for line in f:
            lines += 1
    #assert lines > 0

def test_Runner_openmm_forcefield(tmp_path):
    if not os.path.exists(tmp_path):
        os.mkdir(tmp_path)
    mymodel = make_test_model.make_test_model(tmp_path, mode="forcefield")
    #sim_openmm_factory = sim_openmm.Sim_openmm_factory()
    myanchor = mymodel.anchors[1]
    mmvt_output_filename = os.path.join(
                    tmp_path, myanchor.name, "prod", 
                    "%s%d.%s" % (mmvt_base.OPENMMVT_BASENAME, 1, 
                                         mmvt_base.OPENMMVT_EXTENSION))
    runner = runner_openmm.Runner_openmm(mymodel, myanchor)
    default_output_file, state_file_prefix, restart_index = runner.prepare()
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(mymodel, myanchor, 
                                                      mmvt_output_filename)
    
    runner.run(my_sim_openmm, False)
    
    assert os.path.exists(mmvt_output_filename)
    lines = 0
    with open(mmvt_output_filename,"r") as f:
        for line in f:
            lines += 1
    #assert lines > 0

if __name__ == "__main__":
    tmp_path = "/tmp"
    test_Runner_openmm_default(tmp_path)