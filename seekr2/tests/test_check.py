"""
test_check.py

Testing check.py
"""

import os
import tempfile
import pytest
import shutil

import parmed

import seekr2.modules.check as check

TEST_DIRECTORY = os.path.dirname(__file__)

def test_is_ion():
    """
    Validate the is_ion() function by verifying it gives the correct output
    on a number of atoms in a sample pdb structure.
    """
    
    test_pdb_filename = os.path.join(TEST_DIRECTORY, "data/tryp_ben_at0.pdb")
    structure = parmed.load_file(test_pdb_filename)
    ion_indices = [3239, 3245, 3220, 3246]
    non_ion_indices = [3233, 9068, 12116, 2029, 2051]
    for ion_index in ion_indices:
        ion_atom = structure.atoms[ion_index]
        assert check.is_ion(ion_atom)
    for non_ion_index in non_ion_indices:
        non_ion_atom = structure.atoms[non_ion_index]
        assert not check.is_ion(non_ion_atom)
    return

# TODO: uncomment when other tests complete
#def test_check_pre_sim_bubbles(tryp_ben_mmvt_model):
#    """
#    Test whether the check can properly detect bubbles in the system.
#    """
#    # Test check success: no bubbles
#    no_bubbles_output = check.check_pre_sim_bubbles(tryp_ben_mmvt_model)
#    assert no_bubbles_output
#    
#    # Test check failure: bubbles in structure
#    # copy over the PDB file with the bubble
#    bubbles_files = ["tryp_ben_at0_bubble1.pdb",
#                     "tryp_ben_at0_bubble2.pdb"]
#    for bubbles_file in bubbles_files:
#        src_filename = os.path.join(TEST_DIRECTORY, "data", bubbles_file)
#        dest_filename = os.path.join(
#            tryp_ben_mmvt_model.anchor_rootdir, 
#            tryp_ben_mmvt_model.anchors[0].directory,
#            tryp_ben_mmvt_model.anchors[0].building_directory, 
#            bubbles_file)
#        shutil.copyfile(src_filename, dest_filename)
#        tryp_ben_mmvt_model.anchors[0].amber_params.pdb_coordinates_filename \
#            = bubbles_file
#        bubbles_output = check.check_pre_sim_bubbles(tryp_ben_mmvt_model)
#        #print("testing check of structure with bubbles: ", bubbles_file)
#        assert not bubbles_output
#    
#    return

def test_check_pre_sim_MD_and_BD_salt_concentration(tryp_ben_mmvt_model):
    """
    Test whether the check can find if the MD and BD stages have differring
    salt concentrations.
    """
    # Test check success: equal MD and BD salt concentrations
    tryp_ben_mmvt_model.k_on_info.ions[0].conc = 0.0
    tryp_ben_mmvt_model.k_on_info.ions[1].conc = 0.1
    tryp_ben_mmvt_model.k_on_info.ions[2].conc = 0.0
    equal_conc_output = check.check_pre_sim_MD_and_BD_salt_concentration(
        tryp_ben_mmvt_model)
    assert equal_conc_output
    
    # Test check failure: unequal MD and BD salt concentrations
    tryp_ben_mmvt_model.k_on_info.ions[0].conc = 0.04
    tryp_ben_mmvt_model.k_on_info.ions[1].conc = 0.14
    unequal_conc_output = check.check_pre_sim_MD_and_BD_salt_concentration(
        tryp_ben_mmvt_model)
    assert not unequal_conc_output
    
    for ion in tryp_ben_mmvt_model.k_on_info.ions:
        ion.conc = 0.0
    
    unequal_conc_output = check.check_pre_sim_MD_and_BD_salt_concentration(
        tryp_ben_mmvt_model)
    assert not unequal_conc_output
    return

def test_check_systems_within_Voronoi_cells(tryp_ben_mmvt_model, tmp_path):
    """
    Test whether the check can find systems that exist outside the proper
    Voronoi cells.
    """
    # Test check success: system starting within Voronoi cells
    in_cell_output = check.check_systems_within_Voronoi_cells(
        tryp_ben_mmvt_model)
    assert in_cell_output
    
    # Test check failure(s): system outside Voronoi cells
    anchor_pdb1 \
        = tryp_ben_mmvt_model.anchors[0].amber_params.pdb_coordinates_filename
    anchor_pdb2 \
        = tryp_ben_mmvt_model.anchors[1].amber_params.pdb_coordinates_filename
    anchor_pdb_src_path1 = os.path.join(
        tryp_ben_mmvt_model.anchor_rootdir, 
        tryp_ben_mmvt_model.anchors[0].directory,
        tryp_ben_mmvt_model.anchors[0].building_directory, 
        anchor_pdb1)
    anchor_pdb_src_path2 = os.path.join(
        tryp_ben_mmvt_model.anchor_rootdir, 
        tryp_ben_mmvt_model.anchors[1].directory,
        tryp_ben_mmvt_model.anchors[1].building_directory, 
        anchor_pdb2)
    anchor_pdb_dest_path1 = os.path.join(
        tryp_ben_mmvt_model.anchor_rootdir, 
        tryp_ben_mmvt_model.anchors[0].directory,
        tryp_ben_mmvt_model.anchors[0].building_directory, 
        anchor_pdb2)
    anchor_pdb_dest_path2 = os.path.join(
        tryp_ben_mmvt_model.anchor_rootdir, 
        tryp_ben_mmvt_model.anchors[1].directory,
        tryp_ben_mmvt_model.anchors[1].building_directory, 
        anchor_pdb1)
    shutil.copyfile(anchor_pdb_src_path1, anchor_pdb_dest_path2)
    shutil.copyfile(anchor_pdb_src_path2, anchor_pdb_dest_path1)
    tryp_ben_mmvt_model.anchors[0].amber_params.pdb_coordinates_filename \
        = anchor_pdb2
    tryp_ben_mmvt_model.anchors[1].amber_params.pdb_coordinates_filename \
        = anchor_pdb1
    out_of_cell_output = check.check_systems_within_Voronoi_cells(
        tryp_ben_mmvt_model)
    assert not out_of_cell_output
    return

def test_check_atom_selections_on_same_molecule(tryp_ben_mmvt_model):
    """
    Test whether the check can find atom selections that may exist on
    different molecules.
    """
    # Check success: selections are on same molecule
    same_mol_output = check.check_atom_selections_on_same_molecule(
        tryp_ben_mmvt_model)
    assert same_mol_output
    
    # Check failure: selections are on different molecules
    tryp_ben_mmvt_model.collective_variables[0].group1 = [
        2478, 2489, 2499, 2535, 2718, 2745, 2769, 2787, 2794, 2867, 3221]
    tryp_ben_mmvt_model.collective_variables[0].group2 = [
        2926, 3222, 3223, 3224, 3225, 3226, 3227, 3228, 3229]
    diff_mol_output = check.check_atom_selections_on_same_molecule(
        tryp_ben_mmvt_model)
    assert not diff_mol_output
    
    tryp_ben_mmvt_model.collective_variables[0].group1 = [
        2478, 2489, 2499, 2535, 2718, 2745, 2769, 2787, 2794, 2867, 2926]
    tryp_ben_mmvt_model.collective_variables[0].group2 = [
        3221, 3222, 3223, 3224, 3225, 3226, 3227, 3228, 3250]
    diff_mol_output = check.check_atom_selections_on_same_molecule(
        tryp_ben_mmvt_model)
    assert not diff_mol_output
    
    tryp_ben_mmvt_model.collective_variables[0].group1 = [
        2478, 2489, 2499, 2535, 2718, 2745, 2769, 2787, 2794, 2867, 2926]
    tryp_ben_mmvt_model.collective_variables[0].group2 = [
        3220, 3222, 3223, 3224, 3225, 3226, 3227, 3228, 3229]
    diff_mol_output = check.check_atom_selections_on_same_molecule(
        tryp_ben_mmvt_model)
    assert not diff_mol_output
    
    tryp_ben_mmvt_model.collective_variables[0].group1 = [
        2478, 2489, 2499, 2535, 2718, 2745, 2769, 2787, 2794, 2867, 2926]
    tryp_ben_mmvt_model.collective_variables[0].group2 = [
        3221, 3222, 3223, 3224, 3225, 3226, 3227, 3228, 40000]
    diff_mol_output = check.check_atom_selections_on_same_molecule(
        tryp_ben_mmvt_model)
    assert not diff_mol_output
    return

def test_check_atom_selections_MD_BD(tryp_ben_mmvt_model):
    """
    Test whether the check can identify a mismatch between MD and BD atom
    selections
    """
    # Check success: selections are on same molecule
    same_sel_output = check.check_atom_selections_MD_BD(
        tryp_ben_mmvt_model)
    assert same_sel_output
    
    # Check failure: selections are on different molecules
    tryp_ben_mmvt_model.k_on_info.bd_milestones[0].receptor_indices = [
        2478, 2489, 2499, 2535, 2718, 2745, 2769, 2787, 2794, 2867, 2927]
    diff_sel_output = check.check_atom_selections_MD_BD(
        tryp_ben_mmvt_model)
    assert not diff_sel_output
    
    tryp_ben_mmvt_model.k_on_info.bd_milestones[0].receptor_indices = [
        2478, 2489, 2499, 2535, 2718, 2745, 2769, 2787, 2794, 2867, 2926]
    tryp_ben_mmvt_model.k_on_info.bd_milestones[0].ligand_indices = [
        0, 1, 2, 3, 4, 5, 6, 7, 9]
    diff_sel_output = check.check_atom_selections_MD_BD(
        tryp_ben_mmvt_model)
    assert not diff_sel_output

def test_check_pqr_residues(tryp_ben_mmvt_model):
    """
    Test whether the check can correctly identify if only one resid exists
    in a PQR file.
    """
    pqr_resid_output = check.check_pqr_residues(
        tryp_ben_mmvt_model)
    assert pqr_resid_output
    
    one_resid_pqr_files = ["benzamidine_same_resid.pqr"]
    for one_resid_pqr_file in one_resid_pqr_files:
        src_filename = os.path.join(TEST_DIRECTORY, "data", 
                                    one_resid_pqr_file)
        dest_filename = os.path.join(
            tryp_ben_mmvt_model.anchor_rootdir, 
            tryp_ben_mmvt_model.k_on_info.b_surface_directory,
            one_resid_pqr_file)
        shutil.copyfile(src_filename, dest_filename)
        tryp_ben_mmvt_model.browndye_settings.ligand_pqr_filename \
            = one_resid_pqr_file
        pqr_resid_output = check.check_pqr_residues(tryp_ben_mmvt_model)
        print("testing check of structure with bubbles: ", one_resid_pqr_file)
        assert not pqr_resid_output

def test_check_elber_umbrella_stage(tryp_ben_elber_model):
    """
    Test the check which ensures that Elber umbrella sims remain
    close to the central milestone.
    """
    # Test check success: system starting within Voronoi cells
    small_umbrella_anchor_0_filename = "small_umbrella.dcd"
    src_filename = os.path.join(TEST_DIRECTORY, "data", 
                                small_umbrella_anchor_0_filename)
    dest_filename = os.path.join(
        tryp_ben_elber_model.anchor_rootdir, 
        tryp_ben_elber_model.anchors[0].directory,
        tryp_ben_elber_model.anchors[0].production_directory, 
        "umbrella1.dcd")
    shutil.copyfile(src_filename, dest_filename)
    umbrella_output = check.check_elber_umbrella_stage(
        tryp_ben_elber_model)
    assert umbrella_output
    
    # Test check failure: umbrella traj deviates from milestone
    dest_filename = os.path.join(
        tryp_ben_elber_model.anchor_rootdir, 
        tryp_ben_elber_model.anchors[1].directory,
        tryp_ben_elber_model.anchors[1].production_directory, 
        "umbrella1.dcd")
    shutil.copyfile(src_filename, dest_filename)
    umbrella_output = check.check_elber_umbrella_stage(
        tryp_ben_elber_model)
    assert not umbrella_output

def test_check_xml_boundary_states(tryp_ben_mmvt_model):
    """
    Test the check that determines the location of state files to ensure
    that they exist along the milestone boundaries.
    """
    # Test check success: system starting within Voronoi cells
    state_anchor_1_filename = "tryp_ben_anchor_1_state.xml"
    src_filename = os.path.join(TEST_DIRECTORY, "data", 
                                state_anchor_1_filename)
    state_dir = os.path.join(
        tryp_ben_mmvt_model.anchor_rootdir, 
        tryp_ben_mmvt_model.anchors[1].directory,
        tryp_ben_mmvt_model.anchors[1].production_directory, 
        "states")
    dest_filename = os.path.join(state_dir, "openmm_0_2")
    if not os.path.exists(state_dir):
        os.mkdir(state_dir)
    shutil.copyfile(src_filename, dest_filename)
    xml_state_output = check.check_xml_boundary_states(
        tryp_ben_mmvt_model)
    assert xml_state_output
    
    # Test check failure: umbrella traj deviates from milestone
    state_dir = os.path.join(
        tryp_ben_mmvt_model.anchor_rootdir, 
        tryp_ben_mmvt_model.anchors[3].directory,
        tryp_ben_mmvt_model.anchors[3].production_directory, 
        "states")
    dest_filename = os.path.join(state_dir, "openmm_0_2")
    if not os.path.exists(state_dir):
        os.mkdir(state_dir)
    shutil.copyfile(src_filename, dest_filename)
    xml_state_output = check.check_xml_boundary_states(
        tryp_ben_mmvt_model)
    assert not xml_state_output

def test_check_mmvt_in_Voronoi_cell(tryp_ben_mmvt_model):
    """
    Test the check that makes sure whether MMVT trajectories
    remain within the correct VT cell.
    """
    # Test check success: system starting within Voronoi cells
    state_anchor_1_filename = "tryp_ben_anchor_0_mmvt_traj.dcd"
    src_filename = os.path.join(TEST_DIRECTORY, "data", 
                                state_anchor_1_filename)
    dest_filename = os.path.join(
        tryp_ben_mmvt_model.anchor_rootdir, 
        tryp_ben_mmvt_model.anchors[0].directory,
        tryp_ben_mmvt_model.anchors[0].production_directory, 
        "mmvt1.dcd")
    shutil.copyfile(src_filename, dest_filename)
    mmvt_in_cell_output = check.check_mmvt_in_Voronoi_cell(
        tryp_ben_mmvt_model)
    assert mmvt_in_cell_output
    
    # Test check failure: umbrella traj deviates from milestone
    dest_filename = os.path.join(
        tryp_ben_mmvt_model.anchor_rootdir, 
        tryp_ben_mmvt_model.anchors[1].directory,
        tryp_ben_mmvt_model.anchors[1].production_directory, 
        "mmvt1.dcd")
    shutil.copyfile(src_filename, dest_filename)
    mmvt_out_of_cell_output = check.check_mmvt_in_Voronoi_cell(
        tryp_ben_mmvt_model)
    assert not mmvt_out_of_cell_output

def test_find_parmed_structure_com():
    # Test should not be necessary
    pass

def test_check_bd_simulation_end_state(tryp_ben_mmvt_model):
    """
    Make sure this check can distinguish when BD simulations don't
    end on the proper encounter distance.
    """
    # Test check success: system starting within Voronoi cells
    encounter_rec = "tryp_ben_encounter_rec.pqr"
    encounter_lig = "tryp_ben_encounter_lig.pqr"
    
    bd_milestone_dir = os.path.join(
        tryp_ben_mmvt_model.anchor_rootdir,
        tryp_ben_mmvt_model.k_on_info.bd_milestones[0].directory)
    fhpd_dir = os.path.join(
        bd_milestone_dir,
        tryp_ben_mmvt_model.k_on_info.bd_milestones[0].fhpd_directory)
    if not os.path.exists(fhpd_dir):
        os.mkdir(fhpd_dir)
    lig1_0_0_dir = os.path.join(fhpd_dir, "lig1_0_0")
    if not os.path.exists(lig1_0_0_dir):
        os.mkdir(lig1_0_0_dir)
        
    src_filename = os.path.join(TEST_DIRECTORY, "data", 
                                encounter_rec)
    dest_filename = os.path.join(lig1_0_0_dir, "receptor.pqr")
    shutil.copyfile(src_filename, dest_filename)
    src_filename = os.path.join(TEST_DIRECTORY, "data", 
                                encounter_lig)
    dest_filename = os.path.join(lig1_0_0_dir, "ligand.pqr")
    shutil.copyfile(src_filename, dest_filename)
    bd_on_milestone_output = check.check_bd_simulation_end_state(
        tryp_ben_mmvt_model)
    assert bd_on_milestone_output
    
    # Test check failure
    tryp_ben_mmvt_model.k_on_info.bd_milestones[0].outer_milestone.variables[
        "radius"] = 1.6
    bd_on_milestone_output = check.check_bd_simulation_end_state(
        tryp_ben_mmvt_model)
    assert not bd_on_milestone_output

