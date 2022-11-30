"""
test_mmvt_closest_pair_cv.py
"""

import os
import shutil

import numpy as np

try:
    import openmm.unit as unit
except ImportError:
    import simtk.unit as unit

try:
    import openmm
except ImportError:
    import simtk.openmm as openmm
    
try:
    import openmm.app as openmm_app
except ImportError:
    import simtk.openmm.app as openmm_app
    
import mdtraj

import seekr2.modules.common_base as base
import seekr2.modules.common_sim_openmm as common_sim_openmm
import seekr2.modules.mmvt_cvs.mmvt_closest_pair_cv as mmvt_closest_pair_cv
import seekr2.tests.create_toy_system as create_toy_system


def test_mmvt_closest_pair_boundary(tmp_path):
    # Set initial variables
    milestone_variables = {"k":-1.0, "value":0.1, "exponent":50.0}
    my_variables = [1.0, -1.0, 0.1, 50.0]
    positions1 = [[0.0, 0.0, 0.0], [0.1, 0.0, 0.1],
                  [0.3, 0.0, 0.1], [0.4, 0.0, 0.0]]
    
    # Initialize the toy system with boundary force
    toy_system, toy_topology = create_toy_system.make_toy_system_and_topology(4)
    my_cv = mmvt_closest_pair_cv.MMVT_closest_pair_CV(
        index=0, groups=[[0,1],[2,3]])
    
    box_vectors = base.Box_vectors()
    
    box_vectors.from_quantity(np.array(
        [[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]]) \
        * unit.nanometers)
    
    toy_system.setDefaultPeriodicBoxVectors(*box_vectors.to_quantity())
    toy_topology.setPeriodicBoxVectors(box_vectors.to_quantity())
    
    create_toy_system.assign_nonbonded_cv_info(my_cv, toy_system, box_vectors)
    
    boundary_force1 = my_cv.make_boundary_force(alias_id=1)
    boundary_force1.setForceGroup(1)
    my_cv.add_parameters(boundary_force1)
    my_cv.add_groups_and_variables(
        force=boundary_force1, variables=my_variables, alias_id=1)
    toy_system.addForce(boundary_force1)
    toy_simulation = create_toy_system.make_toy_simulation_and_context(
        toy_system, toy_topology, positions1)
    
    # Check that the particle is within the boundary
    toy_state = toy_simulation.context.getState(getEnergy=True, groups={1})
    boundary_val1 = toy_state.getPotentialEnergy()
    assert boundary_val1.value_in_unit(unit.kilojoules/unit.mole) == 0.0
    
    # Ensure that check_mdtraj_within_boundary works for this cv
    pdb_filename1 = os.path.join(tmp_path, "closest_mmvt_toy1.pdb")
    common_sim_openmm.write_toy_pdb_file(toy_topology, positions1, pdb_filename1)
    pdbtraj1 = mdtraj.load(pdb_filename1)
    assert my_cv.check_mdtraj_within_boundary(
        pdbtraj1, milestone_variables, verbose=False, TOL=0.0)
    
    # Ensure that check_openmm_context_within_boundary works for this cv
    assert my_cv.check_openmm_context_within_boundary(
        toy_simulation.context, milestone_variables)
    
    # Move the particle outside the boundary
    positions2 = [[0.0, 0.0, 0.0], [0.1, 0.0, 0.1],
                  [0.15, 0.0, 0.1], [0.25, 0.0, 0.0]]
    toy_simulation.context.setPositions(positions2)
    
    # Check that the particle is outside the boundary
    toy_state = toy_simulation.context.getState(getEnergy=True, groups={1})
    boundary_val2 = toy_state.getPotentialEnergy()
    assert boundary_val2.value_in_unit(unit.kilojoules/unit.mole) == 1.0
    
    # Move the boundary to inside the interparticle distance
    milestone_variables = {"k":-1.0, "value":0.01, "exponent":50.0}
    new_variables = [1.0, -1.0, 0.01, 50.0]
    
    my_cv.update_groups_and_variables(
        force=boundary_force1, variables=new_variables, alias_id=1,
        context=toy_simulation.context)
    toy_simulation.context.reinitialize(preserveState=True)
    
    # Check that the particle is within the new boundary
    toy_state = toy_simulation.context.getState(getEnergy=True, groups={1})
    boundary_val3 = toy_state.getPotentialEnergy()
    assert boundary_val3.value_in_unit(unit.kilojoules/unit.mole) == 0.0
    
    # Ensure that check_mdtraj_within_boundary still works
    pdb_filename2 = os.path.join(tmp_path, "closest_mmvt_toy2.pdb")
    common_sim_openmm.write_toy_pdb_file(toy_topology, positions2, pdb_filename2)
    pdbtraj2 = mdtraj.load(pdb_filename2)
    assert my_cv.check_mdtraj_within_boundary(
        pdbtraj2, milestone_variables, verbose=False, TOL=0.0)
    
    # Ensure that check_openmm_context_within_boundary still works
    assert my_cv.check_openmm_context_within_boundary(
        toy_simulation.context, milestone_variables)
    
    return

def test_mmvt_count_contacts_voronoi_boundary():
    # Set initial variables
    me_val = 0.2
    neighbor_val = 0.4
    positions1 = [[0.0, 0.0, 0.0], [0.1, 0.0, 0.1],
                  [0.3, 0.0, 0.1], [0.4, 0.0, 0.0]]
    toy_system, toy_topology = create_toy_system.make_toy_system_and_topology(4)
    # Initialize the toy system with boundary force
    box_vectors = base.Box_vectors()
    box_vectors.from_quantity(np.array(
        [[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]]) \
        * unit.nanometers)
    toy_system.setDefaultPeriodicBoxVectors(*box_vectors.to_quantity())
    toy_topology.setPeriodicBoxVectors(box_vectors.to_quantity())
    my_cv = mmvt_closest_pair_cv.MMVT_closest_pair_CV(
        index=0, groups=[[0,1],[2,3]])
    
    create_toy_system.assign_nonbonded_cv_info(my_cv, toy_system, box_vectors)
    
    me_force, neighbor_force = my_cv.make_voronoi_cv_boundary_forces(
        me_val, neighbor_val, 1)
    me_force.setForceGroup(1)
    neighbor_force.setForceGroup(2)
    toy_system.addForce(me_force)
    toy_system.addForce(neighbor_force)
    toy_simulation = create_toy_system.make_toy_simulation_and_context(
        toy_system, toy_topology, positions1)
    
    # Test if the Voronoi boundary returns the correct number
    toy_state_me = toy_simulation.context.getState(getEnergy=True, groups={1})
    toy_state_neighbor = toy_simulation.context.getState(getEnergy=True, groups={2})
    boundary_val_me = toy_state_me.getPotentialEnergy()
    boundary_val_neighbor = toy_state_neighbor.getPotentialEnergy()
    assert np.isclose(boundary_val_me.value_in_unit(unit.kilojoules/unit.mole), 0.0)
    assert not np.isclose(boundary_val_neighbor.value_in_unit(unit.kilojoules/unit.mole), 0.0)
    
    # Move the particle outside the boundary and test Voronoi again
    positions2 = [[0.0, 0.0, 0.0], [0.1, 0.0, 0.1],
                  [0.5, 0.0, 0.1], [0.6, 0.0, 0.0]]
    toy_simulation.context.setPositions(positions2)
    toy_state_me = toy_simulation.context.getState(getEnergy=True, groups={1})
    toy_state_neighbor = toy_simulation.context.getState(getEnergy=True, groups={2})
    boundary_val_me = toy_state_me.getPotentialEnergy()
    boundary_val_neighbor = toy_state_neighbor.getPotentialEnergy()
    assert not np.isclose(boundary_val_me.value_in_unit(unit.kilojoules/unit.mole), 0.0)
    assert np.isclose(boundary_val_neighbor.value_in_unit(unit.kilojoules/unit.mole), 0.0)
    
    # Move the boundary inside the particle and test
    new_me_val = 0.4
    new_neighbor_val = 0.6
    my_cv.update_voronoi_cv_boundary_forces(
        me_force, new_me_val, neighbor_force, new_neighbor_val, 1, 
        toy_simulation.context)
    toy_simulation.context.reinitialize(preserveState=True)
    toy_state_me = toy_simulation.context.getState(getEnergy=True, groups={1})
    toy_state_neighbor = toy_simulation.context.getState(getEnergy=True, groups={2})
    boundary_val_me = toy_state_me.getPotentialEnergy()
    boundary_val_neighbor = toy_state_neighbor.getPotentialEnergy()
    assert np.isclose(boundary_val_me.value_in_unit(unit.kilojoules/unit.mole), 0.0)
    assert not np.isclose(boundary_val_neighbor.value_in_unit(unit.kilojoules/unit.mole), 0.0)
    
    return
    