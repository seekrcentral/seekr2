"""
test_mmvt_voronoi.py
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

import seekr2.modules.common_sim_openmm as common_sim_openmm
import seekr2.modules.mmvt_cvs.mmvt_spherical_cv as mmvt_spherical_cv
import seekr2.modules.mmvt_cvs.mmvt_voronoi_cv as mmvt_voronoi_cv
import seekr2.tests.create_toy_system as create_toy_system

def test_mmvt_voronoi_boundary(tmp_path):
    # Set initial variables
    milestone_variables = {"k":1.0, "me_0":0.4, "neighbor_0":0.7}
    my_variables = [1.0, 1.0, 0.4, 0.7]
    positions1 = [[0.0, 0.0, 0.0], [0.3, 0.4, 0.0]]
    
    # Initialize the toy system with boundary force
    toy_system, toy_topology = create_toy_system.make_toy_system_and_topology(2)
    my_cv = mmvt_voronoi_cv.MMVT_Voronoi_CV(index=0)
    child_cv = mmvt_spherical_cv.MMVT_spherical_CV(index=0, groups=[[0],[1]])
    my_cv.child_cvs.append(child_cv)
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
    pdb_filename1 = os.path.join(tmp_path, "voronoi_mmvt_toy1.pdb")
    common_sim_openmm.write_toy_pdb_file(toy_topology, positions1, pdb_filename1)
    pdbtraj1 = mdtraj.load(pdb_filename1)
    assert my_cv.check_mdtraj_within_boundary(
        pdbtraj1, milestone_variables, verbose=False, TOL=0.0)
    
    # Ensure that check_openmm_context_within_boundary works for this cv
    assert my_cv.check_openmm_context_within_boundary(
        toy_simulation.context, milestone_variables)
    
    # Move the particle outside the boundary
    positions2 = [[0.0, 0.0, 0.0], [0.0, 0.6, 0.0]]
    toy_simulation.context.setPositions(positions2)
    
    # Check that the particle is outside the boundary
    toy_state = toy_simulation.context.getState(getEnergy=True, groups={1})
    boundary_val2 = toy_state.getPotentialEnergy()
    assert boundary_val2.value_in_unit(unit.kilojoules/unit.mole) == 1.0
    
    return
