"""
test_runner_namd.py

test runner_namd.py script(s)
"""

import os
import glob
import pathlib

import seekr2.modules.common_base as base
import seekr2.modules.mmvt_cvs.mmvt_cv_base as mmvt_cv_base
import seekr2.modules.mmvt_sim_namd as mmvt_sim_namd
import seekr2.modules.runner_namd as runner_namd
import seekr2.run as run

TEST_DIRECTORY = os.path.dirname(__file__)
TEST_DATA_DIRECTORY = os.path.join(TEST_DIRECTORY, "data")

def test_search_for_state_to_load(host_guest_mmvt_model):
    anchor0 = host_guest_mmvt_model.anchors[1]
    state_dir = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, anchor0.directory,
        anchor0.production_directory, "states")
    if not os.path.exists(state_dir):
        os.mkdir(state_dir)
    state0 = os.path.join(state_dir, "namdmmvt_2_2")
    pathlib.Path(state0).touch()
    anchor1 = host_guest_mmvt_model.anchors[2]
    state1 = runner_namd.search_for_state_to_load(host_guest_mmvt_model, anchor1)
    assert state1 == os.path.join(state_dir, "namdmmvt_2_2")
    anchor2 = host_guest_mmvt_model.anchors[3]
    state2 = runner_namd.search_for_state_to_load(host_guest_mmvt_model, anchor2)
    assert state2 is None
    return

def test_read_xsc_step_number():
    xsc_filename = os.path.join(TEST_DATA_DIRECTORY, "test_namd.xsc")
    assert runner_namd.read_xsc_step_number(xsc_filename) == 1000

def test_cleanse_anchor_outputs(host_guest_mmvt_model_namd):
    num_steps = 10
    dcd_interval = 1
    host_guest_mmvt_model_namd.calculation_settings.num_production_steps \
        = num_steps
    host_guest_mmvt_model_namd.calculation_settings\
        .restart_checkpoint_interval = 1
    host_guest_mmvt_model_namd.calculation_settings\
        .trajectory_reporter_interval = dcd_interval
    anchor1 = host_guest_mmvt_model_namd.anchors[0]
    prod_directory = os.path.join(
        host_guest_mmvt_model_namd.anchor_rootdir, anchor1.directory,
        anchor1.production_directory)
    run.run(host_guest_mmvt_model_namd, "0", force_overwrite=True, save_state_file=True)
    runner_namd.cleanse_anchor_outputs(
        host_guest_mmvt_model_namd, anchor1)
    #print("os.listdir(prod_directory):", os.listdir(prod_directory))
    assert len(os.listdir(prod_directory)) == 0
    
def test_Runner_namd_default(host_guest_mmvt_model):
    host_guest_mmvt_model.calculation_settings.num_production_steps = 10
    host_guest_mmvt_model.namd_settings = base.Namd_settings()
    myanchor = host_guest_mmvt_model.anchors[1]
    namd_command = "namd2"
    namd_arguments = ""
    runner = runner_namd.Runner_namd(
        host_guest_mmvt_model, myanchor, namd_command, namd_arguments)
    default_output_file, output_basename, state_file_prefix, restart_index = \
        runner.prepare(force_overwrite=True)
    mmvt_output_filename = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, myanchor.directory, "prod", 
        "%s%d" % (mmvt_cv_base.NAMDMMVT_BASENAME, 1))
    my_sim_namd = mmvt_sim_namd.create_sim_namd(
        host_guest_mmvt_model, myanchor, mmvt_output_filename)
    my_sim_namd.namd_root.periodic_boundary_conditions.PMEGridSpacing = None
    runner.run(my_sim_namd, mmvt_output_filename+".out")
    myglob = mmvt_output_filename+".out*"
    outfiles = glob.glob(myglob)
    outfile = outfiles[0]
    assert os.path.exists(outfile)
    return

def test_Runner_namd_load_state(host_guest_mmvt_model):
    host_guest_mmvt_model.calculation_settings.num_production_steps = 10
    host_guest_mmvt_model.namd_settings = base.Namd_settings()
    myanchor = host_guest_mmvt_model.anchors[1]
    namd_command = "namd2"
    namd_arguments = ""
    runner = runner_namd.Runner_namd(
        host_guest_mmvt_model, myanchor, namd_command, namd_arguments)
    default_output_file, output_basename, state_file_prefix, restart_index = \
        runner.prepare(force_overwrite=True)
    mmvt_output_filename = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, myanchor.directory, "prod", 
        "%s%d" % (mmvt_cv_base.NAMDMMVT_BASENAME, 1))
    loading_state_filename = os.path.join(host_guest_mmvt_model.anchor_rootdir, 
                                          "start.state")
    my_sim_namd = mmvt_sim_namd.create_sim_namd(
        host_guest_mmvt_model, myanchor, mmvt_output_filename)
    my_sim_namd.namd_root.periodic_boundary_conditions.PMEGridSpacing = None
    runner.run(my_sim_namd, mmvt_output_filename+".out")
    myglob = mmvt_output_filename+".out*"
    outfiles = glob.glob(myglob)
    outfile = outfiles[0]
    assert os.path.exists(outfile)
    return

def test_Runner_namd_save_states(host_guest_mmvt_model):
    host_guest_mmvt_model.calculation_settings.num_production_steps = 10
    host_guest_mmvt_model.namd_settings = base.Namd_settings()
    myanchor = host_guest_mmvt_model.anchors[1]
    namd_command = "namd2"
    namd_arguments = ""
    runner = runner_namd.Runner_namd(
        host_guest_mmvt_model, myanchor, namd_command, namd_arguments)
    default_output_file, output_basename, state_file_prefix, restart_index = \
        runner.prepare(force_overwrite=True, save_state_file=True)
    mmvt_output_filename = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, myanchor.directory, "prod", 
        "%s%d" % (mmvt_cv_base.NAMDMMVT_BASENAME, 1))
    my_sim_namd = mmvt_sim_namd.create_sim_namd(
        host_guest_mmvt_model, myanchor, mmvt_output_filename)
    my_sim_namd.namd_root.periodic_boundary_conditions.PMEGridSpacing = None
    runner.run(my_sim_namd, mmvt_output_filename+".out")
    myglob = mmvt_output_filename+".out*"
    outfiles = glob.glob(myglob)
    outfile = outfiles[0]
    assert os.path.exists(outfile)
    return

def test_Runner_namd_save_states_until_all_bounds(host_guest_mmvt_model):
    host_guest_mmvt_model.calculation_settings.num_production_steps = 10
    host_guest_mmvt_model.namd_settings = base.Namd_settings()
    myanchor = host_guest_mmvt_model.anchors[1]
    namd_command = "namd2"
    namd_arguments = ""
    runner = runner_namd.Runner_namd(
        host_guest_mmvt_model, myanchor, namd_command, namd_arguments)
    default_output_file, output_basename, state_file_prefix, restart_index = \
        runner.prepare(force_overwrite=True, save_state_boundaries=True)
    mmvt_output_filename = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, myanchor.directory, "prod", 
        "%s%d" % (mmvt_cv_base.NAMDMMVT_BASENAME, 1))
    my_sim_namd = mmvt_sim_namd.create_sim_namd(
        host_guest_mmvt_model, myanchor, mmvt_output_filename)
    my_sim_namd.namd_root.periodic_boundary_conditions.PMEGridSpacing = None
    runner.run(my_sim_namd, mmvt_output_filename+".out")
    myglob = mmvt_output_filename+".out*"
    outfiles = glob.glob(myglob)
    outfile = outfiles[0]
    assert os.path.exists(outfile)
    return