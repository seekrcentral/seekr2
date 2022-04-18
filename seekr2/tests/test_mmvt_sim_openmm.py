"""
test_mmvt_sim_openmm.py
"""

import pytest
import os

try:
    import openmm.unit as unit
except ImportError:
    import simtk.openmm.unit as unit
from seekr2.modules import mmvt_sim_openmm

TEST_DIRECTORY = os.path.dirname(__file__)

def test_add_integrator(toy_mmvt_model):
    my_sim_openmm = mmvt_sim_openmm.MMVT_sim_openmm()
    state_prefix = "test_state_prefix"
    mmvt_sim_openmm.add_integrator(my_sim_openmm, toy_mmvt_model, 
                                   state_prefix=state_prefix)
    assert my_sim_openmm.integrator.getStepSize() \
        == toy_mmvt_model.get_timestep()
    assert my_sim_openmm.integrator.getTemperature() \
        == toy_mmvt_model.temperature
    assert my_sim_openmm.integrator.getFriction() \
        == toy_mmvt_model.openmm_settings.langevin_integrator\
            .friction_coefficient
    assert my_sim_openmm.integrator.getRandomNumberSeed() \
        == toy_mmvt_model.openmm_settings.langevin_integrator.random_seed
    assert my_sim_openmm.integrator.getOutputFileName() \
        == my_sim_openmm.output_filename
    assert my_sim_openmm.integrator.getSaveStateFileName() \
        == state_prefix
    return

def test_add_forces(tmp_path, toy_mmvt_model):
    anchor = toy_mmvt_model.anchors[0]
    output_file = os.path.join(tmp_path, "output.txt")
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(
        toy_mmvt_model, anchor, output_file)
    forces = my_sim_openmm.system.getForces()
    old_force_len = len(forces)
    assert my_sim_openmm.integrator.getNumMilestoneGroups() \
        == len(anchor.milestones)
    mmvt_sim_openmm.add_forces(my_sim_openmm, toy_mmvt_model, anchor)
    forces2 = my_sim_openmm.system.getForces()
    assert len(forces2) == old_force_len + len(anchor.milestones)
    return

def test_add_simulation(tmp_path, toy_mmvt_model):
    anchor = toy_mmvt_model.anchors[0]
    output_file = os.path.join(tmp_path, "output.txt")
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(
        toy_mmvt_model, anchor, output_file)
    assert my_sim_openmm.simulation is not None

@pytest.mark.needs_cuda
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

def test_make_mmvt_boundary_definitions(toy_mmvt_model):
    myanchor = toy_mmvt_model.anchors[1]
    mymilestone = myanchor.milestones[0]
    cv_index = mymilestone.cv_index
    mycv = toy_mmvt_model.collective_variables[cv_index]
    myforce = mmvt_sim_openmm.make_mmvt_boundary_definitions(
            mycv, mymilestone)
    assert myforce.getPerBondParameterName(0) == "bitcode"
    assert myforce.getPerBondParameterName(1) == "k"
    assert myforce.getPerBondParameterName(2) == "value"
    return

def test_get_starting_structure_num_frames(tmp_path, toy_mmvt_model):
    myanchor = toy_mmvt_model.anchors[0]
    output_file = os.path.join(tmp_path, "output.txt")
    num_frames = mmvt_sim_openmm.get_starting_structure_num_frames(
        toy_mmvt_model, myanchor, output_file, load_state_file=None)
    assert num_frames == 2
    num_frames2 = mmvt_sim_openmm.get_starting_structure_num_frames(
        toy_mmvt_model, myanchor, output_file, load_state_file="dummy")
    assert num_frames2 == 1

@pytest.mark.needs_cuda
def test_check_if_state_in_anchor(host_guest_mmvt_model):
    state_file = os.path.join(TEST_DIRECTORY, "data/hostguest_state_a1.xml")
    myanchor1 = host_guest_mmvt_model.anchors[1]
    result1 = mmvt_sim_openmm.check_if_state_in_anchor(
        host_guest_mmvt_model, myanchor1, state_file)
    assert result1
    myanchor2 = host_guest_mmvt_model.anchors[0]
    result2 = mmvt_sim_openmm.check_if_state_in_anchor(
        host_guest_mmvt_model, myanchor2, state_file)
    assert not result2
    myanchor3 = host_guest_mmvt_model.anchors[2]
    result3 = mmvt_sim_openmm.check_if_state_in_anchor(
        host_guest_mmvt_model, myanchor3, state_file)
    assert not result3