"""
test_runner_openmm.py

test runner_openmm.py script(s)
"""

import pytest
import os
import glob
import shutil

import seekr2.modules.common_base as base
import seekr2.modules.mmvt_base as mmvt_base
import seekr2.modules.elber_base as elber_base
import seekr2.modules.mmvt_sim_openmm as mmvt_sim_openmm
import seekr2.modules.elber_sim_openmm as elber_sim_openmm
import seekr2.modules.runner_openmm as runner_openmm
import seekr2.tests.make_test_model as make_test_model

TEST_DIRECTORY = os.path.dirname(__file__)
TEST_DATA_DIRECTORY = os.path.join(TEST_DIRECTORY, "data")

def test_Runner_openmm_default(host_guest_mmvt_model):
    host_guest_mmvt_model.calculation_settings.num_production_steps = 10
    host_guest_mmvt_model.openmm_settings.cuda_platform_settings = None
    host_guest_mmvt_model.openmm_settings.reference_platform = True
    myanchor = host_guest_mmvt_model.anchors[1]
    mmvt_output_filename = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, myanchor.name, "prod", 
        "%s%d.%s" % (mmvt_base.OPENMMVT_BASENAME, 1, 
                             mmvt_base.OPENMMVT_EXTENSION))
    runner = runner_openmm.Runner_openmm(host_guest_mmvt_model, myanchor)
    default_output_file, state_file_prefix, restart_index = runner.prepare(
        force_overwrite=True)
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(
        host_guest_mmvt_model, myanchor, mmvt_output_filename)
    runner.run(my_sim_openmm, False)
    assert os.path.exists(mmvt_output_filename)
    
    # restart
    host_guest_mmvt_model.calculation_settings.num_production_steps = 20
    mmvt_output_filename = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, myanchor.name, "prod", 
        "%s%d.%s" % (mmvt_base.OPENMMVT_BASENAME, 2, 
                             mmvt_base.OPENMMVT_EXTENSION))
    runner = runner_openmm.Runner_openmm(host_guest_mmvt_model, myanchor)
    default_output_file, state_file_prefix, restart_index = runner.prepare(
        restart=True)
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(
        host_guest_mmvt_model, myanchor, mmvt_output_filename)
    runner.run(my_sim_openmm, True, restart_index=restart_index)
    return

@pytest.mark.needs_cuda
def test_Runner_openmm_default_cuda(host_guest_mmvt_model):
    host_guest_mmvt_model.calculation_settings.num_production_steps = 1000
    myanchor = host_guest_mmvt_model.anchors[1]
    mmvt_output_filename = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, myanchor.name, "prod", 
        "%s%d.%s" % (mmvt_base.OPENMMVT_BASENAME, 1, 
                             mmvt_base.OPENMMVT_EXTENSION))
    runner = runner_openmm.Runner_openmm(host_guest_mmvt_model, myanchor)
    default_output_file, state_file_prefix, restart_index = runner.prepare(
        force_overwrite=True)
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(
        host_guest_mmvt_model, myanchor, mmvt_output_filename)
    runner.run(my_sim_openmm, False)
    assert os.path.exists(mmvt_output_filename)
    return

def test_Runner_openmm_load_state(host_guest_mmvt_model):
    host_guest_mmvt_model.calculation_settings.num_production_steps = 10
    host_guest_mmvt_model.openmm_settings.cuda_platform_settings = None
    host_guest_mmvt_model.openmm_settings.reference_platform = True
    myanchor = host_guest_mmvt_model.anchors[1]
    mmvt_output_filename = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, myanchor.name, "prod", 
        "%s%d.%s" % (mmvt_base.OPENMMVT_BASENAME, 1, 
                             mmvt_base.OPENMMVT_EXTENSION))
    loading_state_filename = os.path.join(host_guest_mmvt_model.anchor_rootdir, 
                                          "start.state")
    runner = runner_openmm.Runner_openmm(host_guest_mmvt_model, myanchor)
    default_output_file, state_file_prefix, restart_index = runner.prepare(
        force_overwrite=True)
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(
        host_guest_mmvt_model, myanchor, mmvt_output_filename)
    my_sim_openmm.simulation.saveState(loading_state_filename)
    runner.run(my_sim_openmm, False, load_state_file=loading_state_filename)
    assert os.path.exists(mmvt_output_filename)
    return

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
    default_output_file, state_file_prefix, restart_index = runner.prepare(
        force_overwrite=True)
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(mymodel, myanchor, 
                                                      mmvt_output_filename)
    runner.run(my_sim_openmm, False)
    assert os.path.exists(mmvt_output_filename)

if __name__ == "__main__":
    tmp_path = "/tmp"
    test_Runner_openmm_default(tmp_path)
    
def test_runner_openmm_elber(host_guest_elber_model):
    runner_openmm.REV_STAGE_STEPS_PER_BLOCK = 10
    runner_openmm.FWD_STAGE_STEPS_PER_BLOCK = 10
    host_guest_elber_model.calculation_settings.num_umbrella_stage_steps = 2
    host_guest_elber_model.calculation_settings.umbrella_energy_reporter_interval = 1
    host_guest_elber_model.calculation_settings.umbrella_trajectory_reporter_interval = 1
    host_guest_elber_model.calculation_settings.fwd_rev_interval = 1
    host_guest_elber_model.openmm_settings.cuda_platform_settings = None
    host_guest_elber_model.openmm_settings.reference_platform = True
    myanchor = host_guest_elber_model.anchors[1]
    myanchor.milestones[0].variables["radius"] = 0.140
    myanchor.milestones[2].variables["radius"] = 0.160
    runner = runner_openmm.Runner_openmm(host_guest_elber_model, myanchor)
    default_output_file, state_file_prefix, restart_index = runner.prepare(
        force_overwrite=True)
    elber_output_filename = os.path.join(
        host_guest_elber_model.anchor_rootdir, myanchor.name, "prod", 
        "elber1.out")
    my_sim_openmm = elber_sim_openmm.create_sim_openmm(
        host_guest_elber_model, myanchor, elber_output_filename)
    runner.run(my_sim_openmm, False)
    return
    
def test_runner_openmm_save_states(host_guest_mmvt_model):
    host_guest_mmvt_model.calculation_settings.num_production_steps = 10
    host_guest_mmvt_model.openmm_settings.cuda_platform_settings = None
    host_guest_mmvt_model.openmm_settings.reference_platform = True
    myanchor = host_guest_mmvt_model.anchors[1]
    mmvt_output_filename = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, myanchor.name, "prod", 
        "%s%d.%s" % (mmvt_base.OPENMMVT_BASENAME, 1, 
                             mmvt_base.OPENMMVT_EXTENSION))
    runner = runner_openmm.Runner_openmm(host_guest_mmvt_model, myanchor)
    default_output_file, state_file_prefix, restart_index = runner.prepare(
        force_overwrite=True, save_state_file=True)
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(
        host_guest_mmvt_model, myanchor, mmvt_output_filename)
    runner.run(my_sim_openmm, False)
    assert os.path.exists(mmvt_output_filename)
    return
    
def test_runner_openmm_save_states_until_all_bounds(host_guest_mmvt_model):
    host_guest_mmvt_model.calculation_settings.num_production_steps = 10
    host_guest_mmvt_model.openmm_settings.cuda_platform_settings = None
    host_guest_mmvt_model.openmm_settings.reference_platform = True
    myanchor = host_guest_mmvt_model.anchors[1]
    mmvt_output_filename = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, myanchor.name, "prod", 
        "%s%d.%s" % (mmvt_base.OPENMMVT_BASENAME, 1, 
                             mmvt_base.OPENMMVT_EXTENSION))
    runner = runner_openmm.Runner_openmm(host_guest_mmvt_model, myanchor)
    default_output_file, state_file_prefix, restart_index = runner.prepare(
        force_overwrite=True, save_state_boundaries=True)
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(
        host_guest_mmvt_model, myanchor, mmvt_output_filename)
    runner.run(my_sim_openmm, False)
    assert os.path.exists(mmvt_output_filename)
    return

def test_mmvt_swarm(host_guest_mmvt_model):
    """
    If a multi-frame trajectory is provided as the input PDB, then an MMVT
    swarm should be started.
    """
    swarm_file_name = "hostguest_at0.5_swarm.pdb"
    host_guest_mmvt_model.anchors[0].amber_params.pdb_coordinates_filename \
        = swarm_file_name
    anchor_building_dir = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, 
        host_guest_mmvt_model.anchors[0].directory, 
        host_guest_mmvt_model.anchors[0].building_directory)
    assert os.path.exists(anchor_building_dir)
    src_filename = os.path.join(TEST_DATA_DIRECTORY, swarm_file_name)
    dest_filename = os.path.join(anchor_building_dir, swarm_file_name)
    shutil.copyfile(src_filename, dest_filename)
    host_guest_mmvt_model.calculation_settings.num_production_steps = 10
    host_guest_mmvt_model.openmm_settings.cuda_platform_settings = None
    host_guest_mmvt_model.openmm_settings.reference_platform = True
    myanchor = host_guest_mmvt_model.anchors[0]
    for i in range(10):
        runner = runner_openmm.Runner_openmm(host_guest_mmvt_model, myanchor)
        default_output_file, state_file_prefix, restart_index = runner.prepare(
            force_overwrite=True, swarm_index=i)
        my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(
            host_guest_mmvt_model, myanchor, default_output_file)
        runner.run(my_sim_openmm, False)
        assert os.path.exists(default_output_file)
    
    return
    