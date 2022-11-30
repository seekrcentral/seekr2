"""
test_elber_external_cv.py
"""

try:
    import openmm.unit as unit
except ImportError:
    import simtk.unit as unit

import seekr2.modules.common_sim_openmm as common_sim_openmm
import seekr2.modules.elber_cvs.elber_external_cv as elber_external_cv
import seekr2.tests.create_toy_system as create_toy_system

def test_elber_external_boundary(tmp_path):
    # Set initial variables
    milestone_variables = {"k":-1.0, "value":0.5}
    my_variables = [-1.0, 0.5]
    positions1 = [[0.7, 0.0, 0.0], [0.7, 0.4, 0.0]]
    
    # Initialize the toy system with boundary force
    toy_system, toy_topology = create_toy_system.make_toy_system_and_topology(2)
    my_cv = elber_external_cv.Elber_external_CV(index=0, groups=[[0],[1]])
    my_cv.openmm_fwd_rev_expression = "step(k*(x1 - value))"
    boundary_force1 = my_cv.make_fwd_rev_force_object()
    boundary_force1.setForceGroup(1)
    my_cv.add_fwd_rev_parameters(boundary_force1)
    my_cv.add_groups_and_variables(
        force=boundary_force1, group_list=[0,1], variables=my_variables)
    toy_system.addForce(boundary_force1)
    toy_simulation = create_toy_system.make_toy_simulation_and_context(
        toy_system, toy_topology, positions1)
    
    # Check that the particle is within the boundary
    toy_state = toy_simulation.context.getState(getEnergy=True, groups={1})
    boundary_val1 = toy_state.getPotentialEnergy()
    assert boundary_val1.value_in_unit(unit.kilojoules/unit.mole) == 0.0
    
    # Move the particle outside the boundary
    positions2 = [[0.4, 0.0, 0.0], [0.4, 0.4, 0.0]]
    toy_simulation.context.setPositions(positions2)
    
    # Check that the particle is outside the boundary
    toy_state = toy_simulation.context.getState(getEnergy=True, groups={1})
    boundary_val2 = toy_state.getPotentialEnergy()
    assert boundary_val2.value_in_unit(unit.kilojoules/unit.mole) == 1.0
    