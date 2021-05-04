"""
test_runner_browndye2.py

test runner_browndye2.py script(s)
"""

import os
import glob

import seekr2.tests.make_test_model as make_test_model
import seekr2.modules.common_sim_browndye2 as sim_browndye2
import seekr2.modules.runner_browndye2 as runner_browndye2

def test_runner_browndye2_b_surface_default(tmp_path):
    model = make_test_model.make_test_model(tmp_path)
    b_surface_abs_path = os.path.join(tmp_path, "b_surface")
    receptor_pqr_filename = os.path.join(
        b_surface_abs_path, model.browndye_settings.receptor_pqr_filename)
    ligand_pqr_filename = os.path.join(
        b_surface_abs_path, model.browndye_settings.ligand_pqr_filename)
    bd_milestone = model.k_on_info.bd_milestones[0]
    
    ghost_index_rec = \
                sim_browndye2.add_ghost_atom_to_pqr_from_atoms_center_of_mass(
                    receptor_pqr_filename, bd_milestone.receptor_indices)
    ghost_index_lig = \
                sim_browndye2.add_ghost_atom_to_pqr_from_atoms_center_of_mass(
                    ligand_pqr_filename, bd_milestone.ligand_indices)
    assert ghost_index_rec == 148
    assert ghost_index_lig == 16
    
    receptor_xml_filename = sim_browndye2.make_pqrxml(receptor_pqr_filename)
    ligand_xml_filename = sim_browndye2.make_pqrxml(ligand_pqr_filename)
    
    debye_length, reaction_filename = \
        runner_browndye2.make_browndye_input_xml(
        model, tmp_path, receptor_xml_filename, ligand_xml_filename,
        model.k_on_info.b_surface_num_steps)
    model.browndye_settings.debye_length = debye_length
    assert os.path.exists(os.path.join(b_surface_abs_path, "apbs_input.xml"))
    assert os.path.exists(os.path.join(b_surface_abs_path, "input.xml"))
    abs_reaction_path = os.path.join(b_surface_abs_path, 
                                     reaction_filename)
    runner_browndye2.make_browndye_reaction_xml(model, abs_reaction_path)
    assert os.path.exists(abs_reaction_path)
    bd_directory = b_surface_abs_path
    runner_browndye2.run_bd_top(model.browndye_settings.browndye_bin_dir, 
               bd_directory)
    #runner_browndye2.modify_variables(bd_directory, 10000)
    runner_browndye2.run_nam_simulation(
        model.browndye_settings.browndye_bin_dir, bd_directory, 
        model.k_on_info.bd_output_glob)
    #assert os.path.exists(os.path.join(b_surface_abs_path, "results.xml"))
    assert len(glob.glob(os.path.join(b_surface_abs_path, "results*.xml"))) > 0
    """
    bd_milestone_abs_path = os.path.join(tmp_path, bd_milestone.directory)
    assert os.path.exists(bd_milestone_abs_path)
    runner_browndye2.extract_bd_surface(model, bd_milestone, 10)
    bd_directory_list = runner_browndye2.make_fhpd_directories(
        model, bd_milestone)
    
    for bd_directory in bd_directory_list:
        runner_browndye2.run_bd_top(model.browndye_settings.browndye_bin_dir, 
                           bd_directory)
        runner_browndye2.modify_variables(bd_directory, n_trajectories=1000)
        runner_browndye2.run_nam_simulation(
            model.browndye_settings.browndye_bin_dir, bd_directory, 
            model.k_on_info.bd_output_glob)
    
    runner_browndye2.combine_fhpd_results(
        model, bd_milestone, bd_directory_list)
    assert os.path.exists(os.path.join(bd_milestone_abs_path, "results.xml"))
    """
    return

# TODO: bd_milestone test runs