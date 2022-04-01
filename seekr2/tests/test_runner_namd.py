"""
test_runner_namd.py

test runner_namd.py script(s)
"""

import os
import pytest
import glob
import re

import seekr2.modules.common_base as base
import seekr2.modules.mmvt_base as mmvt_base
import seekr2.modules.common_sim_namd as sim_namd
import seekr2.modules.mmvt_sim_namd as mmvt_sim_namd
import seekr2.modules.runner_namd as runner_namd
import seekr2.tests.make_test_model as make_test_model

def test_Runner_namd_default(host_guest_mmvt_model):
    host_guest_mmvt_model.calculation_settings.num_production_steps = 10
    host_guest_mmvt_model.namd_settings = base.Namd_settings()
    myanchor = host_guest_mmvt_model.anchors[1]
    namd_command = "namd2"
    namd_arguments = ""
    runner = runner_namd.Runner_namd(
        host_guest_mmvt_model, myanchor, namd_command, namd_arguments)
    default_output_file, output_basename, state_file_prefix, restart_index = \
        runner.prepare()
    mmvt_output_filename = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, myanchor.directory, "prod", 
        "%s%d" % (mmvt_base.NAMDMMVT_BASENAME, 1))
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
        runner.prepare()
    mmvt_output_filename = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, myanchor.directory, "prod", 
        "%s%d" % (mmvt_base.NAMDMMVT_BASENAME, 1))
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
        "%s%d" % (mmvt_base.NAMDMMVT_BASENAME, 1))
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
        "%s%d" % (mmvt_base.NAMDMMVT_BASENAME, 1))
    my_sim_namd = mmvt_sim_namd.create_sim_namd(
        host_guest_mmvt_model, myanchor, mmvt_output_filename)
    my_sim_namd.namd_root.periodic_boundary_conditions.PMEGridSpacing = None
    runner.run(my_sim_namd, mmvt_output_filename+".out")
    myglob = mmvt_output_filename+".out*"
    outfiles = glob.glob(myglob)
    outfile = outfiles[0]
    assert os.path.exists(outfile)
    return