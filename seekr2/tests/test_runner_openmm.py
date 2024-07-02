"""
test_runner_openmm.py

test runner_openmm.py script(s)
"""

import pytest
import os
import shutil
import pathlib
import glob

import seekr2.modules.mmvt_cvs.mmvt_cv_base as mmvt_cv_base
import seekr2.modules.mmvt_sim_openmm as mmvt_sim_openmm
import seekr2.modules.elber_sim_openmm as elber_sim_openmm
import seekr2.modules.runner_openmm as runner_openmm
import seekr2.run as run

TEST_DIRECTORY = os.path.dirname(__file__)
TEST_DATA_DIRECTORY = os.path.join(TEST_DIRECTORY, "data")

def test_all_boundaries_have_state(toy_mmvt_model):
    anchor = toy_mmvt_model.anchors[1]
    state_dir = os.path.join(
        toy_mmvt_model.anchor_rootdir, anchor.directory,
        anchor.production_directory, "states")
    if not os.path.exists(state_dir):
        os.mkdir(state_dir)
    state_glob = os.path.join(state_dir, "openmm_*")
    state1 = os.path.join(state_dir, "openmm_0_1")
    pathlib.Path(state1).touch()
    assert not runner_openmm.all_boundaries_have_state(state_glob, anchor)
    state2 = os.path.join(state_dir, "openmm_1_1")
    pathlib.Path(state2).touch()
    assert not runner_openmm.all_boundaries_have_state(state_glob, anchor)
    state3 = os.path.join(state_dir, "openmm_2_2")
    pathlib.Path(state3).touch()
    assert runner_openmm.all_boundaries_have_state(state_glob, anchor)
    return

def test_elber_anchor_has_umbrella_files(toy_elber_model):
    anchor = toy_elber_model.anchors[1]
    runner_openmm.cleanse_anchor_outputs(
        toy_elber_model, anchor, skip_umbrella_files=False)
    assert not runner_openmm.elber_anchor_has_umbrella_files(
        toy_elber_model, anchor)
    prod_dir = os.path.join(
        toy_elber_model.anchor_rootdir, anchor.directory,
        anchor.production_directory)
    umbrella1 = os.path.join(prod_dir, "umbrella1.dcd")
    pathlib.Path(umbrella1).touch()
    assert runner_openmm.elber_anchor_has_umbrella_files(
        toy_elber_model, anchor)
    return

def test_search_for_state_to_load(toy_mmvt_model):
    anchor0 = toy_mmvt_model.anchors[1]
    state_dir = os.path.join(
        toy_mmvt_model.anchor_rootdir, anchor0.directory,
        anchor0.production_directory, "states")
    if not os.path.exists(state_dir):
        os.mkdir(state_dir)
    state0 = os.path.join(state_dir, "openmm_2_2")
    pathlib.Path(state0).touch()
    anchor1 = toy_mmvt_model.anchors[2]
    state1 = runner_openmm.search_for_state_to_load(toy_mmvt_model, anchor1)
    assert state1 == os.path.join(state_dir, "openmm_2_2")
    anchor2 = toy_mmvt_model.anchors[3]
    state2 = runner_openmm.search_for_state_to_load(toy_mmvt_model, anchor2)
    assert state2 is None
    return
    
def test_get_data_file_length(tmp_path):
    assert runner_openmm.get_data_file_length("doesnt_exist.txt") == 0
    num_lines1 = 13
    data_file1 = os.path.join(tmp_path, "data_file1.txt")
    with open(data_file1, "w") as f:
        for i in range(num_lines1):
            f.write("line {}\n".format(i))
            
    assert runner_openmm.get_data_file_length(data_file1) == num_lines1
    return

def test_read_reversal_data_file_last(tmp_path):
    data_file1 = os.path.join(tmp_path, "data_file1.txt")
    with open(data_file1, "w") as f:
        f.write("1,0,123\n")
    assert runner_openmm.read_reversal_data_file_last(data_file1)
    with open(data_file1, "a") as f:
        f.write("2,1,234\n")
    assert not runner_openmm.read_reversal_data_file_last(data_file1)
    with open(data_file1, "a") as f:
        f.write("3,2,345\n")
    assert runner_openmm.read_reversal_data_file_last(data_file1)
    return

def test_cleanse_anchor_outputs(toy_mmvt_model, toy_elber_model):
    num_steps = 10000
    toy_mmvt_model.calculation_settings.num_production_steps = num_steps
    toy_mmvt_model.openmm_settings.cuda_platform_settings = None
    toy_mmvt_model.openmm_settings.reference_platform = True
    anchor1 = toy_mmvt_model.anchors[0]
    prod_directory = os.path.join(
        toy_mmvt_model.anchor_rootdir, anchor1.directory,
        anchor1.production_directory)
    run.run(toy_mmvt_model, "0", force_overwrite=True, save_state_file=True)
    runner_openmm.cleanse_anchor_outputs(
        toy_mmvt_model, anchor1, skip_umbrella_files=False)
    assert len(os.listdir(prod_directory)) == 0
    
    num_steps = 10000
    toy_elber_model.calculation_settings.num_umbrella_stage_steps = num_steps
    toy_elber_model.openmm_settings.cuda_platform_settings = None
    toy_elber_model.openmm_settings.reference_platform = True
    anchor1 = toy_mmvt_model.anchors[1]
    prod_directory = os.path.join(
        toy_elber_model.anchor_rootdir, anchor1.directory,
        anchor1.production_directory)
    run.run(toy_elber_model, "1", force_overwrite=True, save_state_file=True)
    runner_openmm.cleanse_anchor_outputs(
        toy_elber_model, anchor1, skip_umbrella_files=False)
    assert len(os.listdir(prod_directory)) == 0
    
    run.run(toy_elber_model, "1", force_overwrite=True, save_state_file=True)
    runner_openmm.cleanse_anchor_outputs(
        toy_elber_model, anchor1, skip_umbrella_files=True)
    assert os.listdir(prod_directory)[0].endswith("umbrella1.dcd")

def test_get_last_bounce(tmp_path):
    num_lines1 = 13
    data_file1 = os.path.join(tmp_path, "data_file1.txt")
    with open(data_file1, "w") as f:
        for i in range(num_lines1):
            f.write("1,{},{}\n".format(i, i*0.002))
            
    assert runner_openmm.get_last_bounce(data_file1) == 13
    return

def test_saveCheckpoint(host_guest_mmvt_model):
    host_guest_mmvt_model.calculation_settings.num_production_steps = 10
    host_guest_mmvt_model.openmm_settings.cuda_platform_settings = None
    host_guest_mmvt_model.openmm_settings.reference_platform = True
    myanchor = host_guest_mmvt_model.anchors[1]
    mmvt_output_filename = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, myanchor.directory, "prod", 
        "%s%d.%s" % (mmvt_cv_base.OPENMMVT_BASENAME, 1, 
                             mmvt_cv_base.OPENMMVT_EXTENSION))
    runner = runner_openmm.Runner_openmm(host_guest_mmvt_model, myanchor)
    default_output_file, state_file_prefix, restart_index = runner.prepare(
        force_overwrite=True)
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(
        host_guest_mmvt_model, myanchor, mmvt_output_filename)
    #assert not os.path.exists(runner.restart_checkpoint_filename)
    if os.path.exists(runner.restart_checkpoint_filename):
        os.remove(runner.restart_checkpoint_filename)
    runner_openmm.saveCheckpoint(my_sim_openmm, 
                                 runner.restart_checkpoint_filename)
    assert os.path.exists(runner.restart_checkpoint_filename)
    lastline = ""
    with open(my_sim_openmm.output_filename, "r") as f:
        for line in f.readlines():
            lastline = line
            print("line:", line)
    assert lastline.startswith("CHECKPOINT")

def test_Runner_openmm_default(host_guest_mmvt_model):
    host_guest_mmvt_model.calculation_settings.num_production_steps = 10
    host_guest_mmvt_model.calculation_settings.restart_checkpoint_interval = 10
    host_guest_mmvt_model.openmm_settings.cuda_platform_settings = None
    host_guest_mmvt_model.openmm_settings.reference_platform = True
    myanchor = host_guest_mmvt_model.anchors[1]
    mmvt_output_filename = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, myanchor.directory, "prod", 
        "%s%d.%s" % (mmvt_cv_base.OPENMMVT_BASENAME, 1, 
                             mmvt_cv_base.OPENMMVT_EXTENSION))
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
        host_guest_mmvt_model.anchor_rootdir, myanchor.directory, "prod", 
        "%s%d.%s" % (mmvt_cv_base.OPENMMVT_BASENAME, 2, 
                             mmvt_cv_base.OPENMMVT_EXTENSION))
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
        host_guest_mmvt_model.anchor_rootdir, myanchor.directory, "prod", 
        "%s%d.%s" % (mmvt_cv_base.OPENMMVT_BASENAME, 1, 
                             mmvt_cv_base.OPENMMVT_EXTENSION))
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
    host_guest_mmvt_model.calculation_settings.restart_checkpoint_interval = 10
    host_guest_mmvt_model.openmm_settings.cuda_platform_settings = None
    host_guest_mmvt_model.openmm_settings.reference_platform = True
    myanchor = host_guest_mmvt_model.anchors[1]
    mmvt_output_filename = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, myanchor.directory, "prod", 
        "%s%d.%s" % (mmvt_cv_base.OPENMMVT_BASENAME, 1, 
                             mmvt_cv_base.OPENMMVT_EXTENSION))
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

# No working XML inputs for hostguest
@pytest.mark.skip()
def test_Runner_openmm_forcefield(tmp_path, host_guest_mmvt_model_forcefield):
    if not os.path.exists(tmp_path):
        os.mkdir(tmp_path)
    mymodel = host_guest_mmvt_model_forcefield
    myanchor = mymodel.anchors[1]
    mmvt_output_filename = os.path.join(
                    mymodel.anchor_rootdir, myanchor.directory, "prod", 
                    "%s%d.%s" % (mmvt_cv_base.OPENMMVT_BASENAME, 1, 
                                         mmvt_cv_base.OPENMMVT_EXTENSION))
    runner = runner_openmm.Runner_openmm(mymodel, myanchor)
    default_output_file, state_file_prefix, restart_index = runner.prepare(
        force_overwrite=True)
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(mymodel, myanchor, 
                                                      mmvt_output_filename)
    runner.run(my_sim_openmm, False)
    assert os.path.exists(mmvt_output_filename)
    
def test_Runner_openmm_system(tmp_path, host_guest_mmvt_model_system):
    if not os.path.exists(tmp_path):
        os.mkdir(tmp_path)
    mymodel = host_guest_mmvt_model_system
    myanchor = mymodel.anchors[1]
    mmvt_output_filename = os.path.join(
                    mymodel.anchor_rootdir, myanchor.directory, "prod", 
                    "%s%d.%s" % (mmvt_cv_base.OPENMMVT_BASENAME, 1, 
                                         mmvt_cv_base.OPENMMVT_EXTENSION))
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
        host_guest_elber_model.anchor_rootdir, myanchor.directory, "prod", 
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
        host_guest_mmvt_model.anchor_rootdir, myanchor.directory, "prod", 
        "%s%d.%s" % (mmvt_cv_base.OPENMMVT_BASENAME, 1, 
                             mmvt_cv_base.OPENMMVT_EXTENSION))
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
        host_guest_mmvt_model.anchor_rootdir, myanchor.directory, "prod", 
        "%s%d.%s" % (mmvt_cv_base.OPENMMVT_BASENAME, 1, 
                             mmvt_cv_base.OPENMMVT_EXTENSION))
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
