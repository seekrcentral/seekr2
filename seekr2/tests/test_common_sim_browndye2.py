"""
test_sim_browndye2.py
"""

import pytest
import os
import re

import parmed
import numpy as np

import seekr2.modules.common_sim_browndye2 as sim_browndye2

def verify_xml_text(filename, tag, text):
    with open(filename, 'r') as f:
        for line in f.readlines():
            #print("line:", line)
            result = re.search("<%s>%s</%s>" % (tag, text, tag), line)
            if result:
                return
            """
            result = re.search("<%s>(.+?)</%s>" % (tag, tag), line)
            if result:
                xml_text = result.group(1)
                assert xml_text == text
                return
            """
    
    raise Exception("tag not found in file: %s" % tag)
    return

def make_test_browndye_input(tmp_path):
    test_filename = os.path.join(tmp_path, "test_input.xml")
    root = sim_browndye2.Root()
    root.n_threads = 3
    root.seed = 4567
    root.output = "test_output.xml"
    root.n_trajectories = 65535
    root.n_trajectories_per_output = 5
    root.max_n_steps = 123456
    root.trajectory_file = "test_traj"
    root.n_steps_per_output = 123
    root.system.start_at_site = "true"
    root.system.reaction_file = "my_rxns.xml"
    root.system.hydrodynamic_interactions = "false"
    root.system.solvent.dielectric = 77.4
    root.system.solvent.relative_viscosity = 0.4
    root.system.solvent.debye_length = 9.8
    root.system.solvent.kT = 1.456
    root.system.solvent.desolvation_parameter = 0.8
    ion1 = sim_browndye2.Ion()
    ion1.radius = 1.34
    ion1.charge = -2.0
    ion1.conc = 0.2
    ion2 = sim_browndye2.Ion()
    ion2.radius = 0.85
    ion2.charge = 1.0
    ion2.conc = 0.4
    root.system.solvent.ions = [ion1, ion2]
    root.system.time_step_tolerances = sim_browndye2.Time_step_tolerances()
    root.system.time_step_tolerances.minimum_core_dt = 0.1
    root.system.time_step_tolerances.minimum_core_reaction_dt = 0.01
    group1 = sim_browndye2.Group()
    group1.name = "receptor"
    core1 = sim_browndye2.Core()
    core1.name = "receptor"
    core1.atoms = "receptor.xml"
    core1.all_in_surface = "true"
    core1.is_protein = "true"
    core1.dielectric = 5.0
    core1.grid_spacing = 0.75
    grid1 = "receptor.dx"
    core1.electric_field.grid_list.append(grid1)
    group1.core_list.append(core1)
    
    group2 = sim_browndye2.Group()
    group2.name = "ligand"
    core2 = sim_browndye2.Core()
    core2.name = "ligand"
    core2.atoms = "ligand.xml"
    core1.all_in_surface = "true"
    core1.is_protein = "true"
    core1.dielectric = 5.0
    core2.grid_spacing = 0.5
    grid2 = "ligand.dx"
    core2.electric_field.grid_list.append(grid2)
    group2.core_list.append(core2)
    
    root.system.group_list.append(group1)
    root.system.group_list.append(group2)
    
    root.write(test_filename)
    return test_filename

def make_test_browndye_rxn(tmp_path):
    test_filename = os.path.join(tmp_path, "test_rxn.xml")
    rxnroot = sim_browndye2.Reaction_root()
    rxnroot.first_state = "something1"
    rxn = sim_browndye2.Reaction()
    rxn.name = "test_rxn1"
    rxn.state_before = "something1"
    rxn.state_after = "something2"
    rxn.molecule0_group = "molec1"
    rxn.molecule0_core = "molec1"
    rxn.molecule1_group = "molec2"
    rxn.molecule1_core = "molec2"
    rxn.n_needed = 4
    pair = sim_browndye2.Pair()
    pair.atom1_index = 5432
    pair.atom2_index = 789
    pair.distance = 8.45
    rxn.pair_list.append(pair)
    rxnroot.reaction_list.append(rxn)
    rxnroot.write(test_filename)
    return test_filename

def test_add_ghost_atom_to_pqr_from_atoms_center_of_mass(tmp_path):
    input_pqr_filename = \
        os.path.join(os.path.dirname(__file__), 
                     "../data/hostguest_files/hostguest_ligand.pqr")
    #print("input_pqr_filename:", input_pqr_filename)
    atom_index_list = list(range(15))
    output_pqr_filename = os.path.join(tmp_path, "test.pqr")
    ghost_index = sim_browndye2.add_ghost_atom_to_pqr_from_atoms_center_of_mass(
        input_pqr_filename, atom_index_list, output_pqr_filename)
    pqr_struct = parmed.load_file(output_pqr_filename, skip_bonds=True)
    ghost_atom = pqr_struct.atoms[ghost_index-1]
    assert(ghost_atom.name == "GHO")
    expected_ghost_location = np.array([[-0.056, -0.323, 2.439]])
    ghost_location = pqr_struct.coordinates[ghost_index-1]
    difference = np.linalg.norm(expected_ghost_location - ghost_location)
    assert(difference == 0.0)
    return

def test_browndye_input_serializer_apbs_mode(tmp_path):
    test_filename = make_test_browndye_input(tmp_path)
    assert os.path.exists(test_filename)
    
    verify_xml_text(test_filename, "n_trajectories", "65535")
    verify_xml_text(test_filename, "n_threads", "3")
    verify_xml_text(test_filename, "seed", "4567")
    verify_xml_text(test_filename, "output", "test_output.xml")
    verify_xml_text(test_filename, "n_trajectories_per_output", "5")
    verify_xml_text(test_filename, "max_n_steps", "123456")
    verify_xml_text(test_filename, "trajectory_file", "test_traj")
    verify_xml_text(test_filename, "n_steps_per_output", "123")
    verify_xml_text(test_filename, "start_at_site", "true")
    verify_xml_text(test_filename, "reaction_file", "my_rxns.xml")
    verify_xml_text(test_filename, "hydrodynamic_interactions", "false")
    verify_xml_text(test_filename, "dielectric", "77.4")
    verify_xml_text(test_filename, "relative_viscosity", "0.4")
    verify_xml_text(test_filename, "kT", "1.456")
    verify_xml_text(test_filename, "desolvation_parameter", "0.8")
    verify_xml_text(test_filename, "radius", "1.34")
    verify_xml_text(test_filename, "charge", "-2.0")
    verify_xml_text(test_filename, "conc", "0.2")
    verify_xml_text(test_filename, "radius", "0.85")
    verify_xml_text(test_filename, "charge", "1.0")
    verify_xml_text(test_filename, "conc", "0.4")
    verify_xml_text(test_filename, "minimum_core_dt", "0.1")
    verify_xml_text(test_filename, "minimum_core_reaction_dt", "0.01")
    verify_xml_text(test_filename, "name", "receptor")
    verify_xml_text(test_filename, "atoms", "receptor.xml")
    verify_xml_text(test_filename, "all_in_surface", "true")
    verify_xml_text(test_filename, "is_protein", "true")
    verify_xml_text(test_filename, "dielectric", "5.0")
    verify_xml_text(test_filename, "grid_spacing", "0.75")
    verify_xml_text(test_filename, "name", "ligand")
    verify_xml_text(test_filename, "atoms", "ligand.xml")
    verify_xml_text(test_filename, "all_in_surface", "true")
    verify_xml_text(test_filename, "is_protein", "true")
    verify_xml_text(test_filename, "dielectric", "5.0")
    verify_xml_text(test_filename, "grid_spacing", "0.5")
    
    pytest.raises(Exception, verify_xml_text, test_filename, "gobblygook", "18")
    
    return

def test_browndye_input_serializer_apbs_mode_rxn(tmp_path):
    test_filename = make_test_browndye_rxn(tmp_path)
    assert os.path.exists(test_filename)
    
    verify_xml_text(test_filename, "first_state", "something1")
    verify_xml_text(test_filename, "name", "test_rxn1")
    verify_xml_text(test_filename, "state_before", "something1")
    verify_xml_text(test_filename, "state_after", "something2")
    verify_xml_text(test_filename, "molecule0", "molec1 molec1")
    verify_xml_text(test_filename, "molecule1", "molec2 molec2")
    verify_xml_text(test_filename, "n_needed", "4")
    verify_xml_text(test_filename, "atoms", "5432 789")
    verify_xml_text(test_filename, "distance", "8.45")
    pytest.raises(Exception, verify_xml_text, test_filename, "gobblygook", "18")
