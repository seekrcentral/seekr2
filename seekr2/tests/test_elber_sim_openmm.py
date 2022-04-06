"""
test_elber_sim_openmm.py
"""

import pytest
import os

try:
    import openmm.unit as unit
except ImportError:
    import simtk.openmm.unit as unit

from seekr2.modules import elber_sim_openmm

def test_add_integrators(toy_elber_model):
    my_sim_openmm = elber_sim_openmm.Elber_sim_openmm()
    state_prefix = "test_state_prefix"
    elber_sim_openmm.add_integrators(my_sim_openmm, toy_elber_model, 
                                   state_prefix=state_prefix)
    assert my_sim_openmm.umbrella_integrator.getStepSize() \
        == toy_elber_model.get_timestep() * unit.picoseconds
    assert my_sim_openmm.umbrella_integrator.getTemperature() \
        == toy_elber_model.temperature * unit.kelvin
    assert my_sim_openmm.umbrella_integrator.getFriction() \
        == toy_elber_model.openmm_settings.langevin_integrator\
            .friction_coefficient / unit.picoseconds
    assert my_sim_openmm.umbrella_integrator.getRandomNumberSeed() \
        == toy_elber_model.openmm_settings.langevin_integrator.random_seed
    
    assert my_sim_openmm.rev_integrator.getOutputFileName() \
        == my_sim_openmm.rev_output_filename
    
    assert my_sim_openmm.fwd_integrator.getOutputFileName() \
        == my_sim_openmm.fwd_output_filename
    assert my_sim_openmm.fwd_integrator.getSaveStateFileName() \
        == state_prefix
    return

def test_add_forces(tmp_path, toy_elber_model):
    anchor = toy_elber_model.anchors[1]
    output_file = os.path.join(tmp_path, "output.txt")
    my_sim_openmm = elber_sim_openmm.create_sim_openmm(
        toy_elber_model, anchor, output_file)
    umbrella_forces = my_sim_openmm.umbrella_system.getForces()
    old_force_len_umbrella = len(umbrella_forces)
    rev_forces = my_sim_openmm.rev_system.getForces()
    old_force_len_rev = len(rev_forces)
    fwd_forces = my_sim_openmm.fwd_system.getForces()
    old_force_len_fwd = len(fwd_forces)
    assert my_sim_openmm.rev_integrator.getNumSrcMilestoneGroups() \
        == 1
    assert my_sim_openmm.rev_integrator.getNumDestMilestoneGroups() \
        == len(anchor.milestones)-1
    assert my_sim_openmm.fwd_integrator.getNumSrcMilestoneGroups() \
        == 1
    assert my_sim_openmm.fwd_integrator.getNumDestMilestoneGroups() \
        == len(anchor.milestones)-1
    elber_sim_openmm.add_forces(my_sim_openmm, toy_elber_model, anchor)
    
    umbrella_forces2 = my_sim_openmm.umbrella_system.getForces()
    assert len(umbrella_forces2) == old_force_len_umbrella + 1
    rev_forces2 = my_sim_openmm.rev_system.getForces()
    assert len(rev_forces2) == old_force_len_rev + len(anchor.milestones)
    fwd_forces2 = my_sim_openmm.fwd_system.getForces()
    assert len(fwd_forces2) == old_force_len_fwd + len(anchor.milestones)
    return

def test_add_simulation(tmp_path, toy_elber_model):
    anchor = toy_elber_model.anchors[0]
    output_file = os.path.join(tmp_path, "output.txt")
    my_sim_openmm = elber_sim_openmm.create_sim_openmm(
        toy_elber_model, anchor, output_file)
    assert my_sim_openmm.umbrella_simulation is not None
    assert my_sim_openmm.rev_simulation is not None
    assert my_sim_openmm.fwd_simulation is not None

def test_elber_sim_openmm_toy(tmp_path, toy_elber_model):
    myanchor = toy_elber_model.anchors[1]
    toy_elber_model.openmm_settings.cuda_platform_settings = None
    toy_elber_model.openmm_settings.reference_platform = True
    output_file = os.path.join(tmp_path, "output.txt")
    my_sim_openmm = elber_sim_openmm.create_sim_openmm(
        toy_elber_model, myanchor, output_file)
    return

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

def test_make_elber_umbrella_force(toy_elber_model):
    myanchor = toy_elber_model.anchors[1]
    mymilestone = myanchor.milestones[1]
    cv_index = mymilestone.cv_index
    mycv = toy_elber_model.collective_variables[cv_index]
    myforce = elber_sim_openmm.make_elber_umbrella_force(
            toy_elber_model, myanchor)
    assert myforce.getPerBondParameterName(0) == "k"
    assert myforce.getPerBondParameterName(1) == "value"
    return

def test_make_elber_boundary_definitions(toy_elber_model):
    myanchor = toy_elber_model.anchors[1]
    mymilestone = myanchor.milestones[1]
    cv_index = mymilestone.cv_index
    mycv = toy_elber_model.collective_variables[cv_index]
    myforce = elber_sim_openmm.make_elber_boundary_definitions(
            mycv, mymilestone)
    assert myforce.getPerBondParameterName(0) == "k"
    assert myforce.getPerBondParameterName(1) == "value"
    return
