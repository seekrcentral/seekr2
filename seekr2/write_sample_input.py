"""
write_sample_input.py

Writes sample input files which can be modified and used for custom
input.
"""

import sys
import os

import seekr2.common.base as base
import seekr2.common.collective_variables as common_cv
import seekr2.common.prepare as common_prepare


def create_model_input(root_dir, calculation_type="mmvt"):
    """
    TEMPORARY
    """
    model_input = common_prepare.Model_input()
    model_input.root_directory = root_dir
    if calculation_type == "mmvt":
        model_input.calculation_type = "mmvt"
        model_input.calculation_settings = common_prepare.MMVT_input_settings()
        model_input.calculation_settings.md_output_interval = 10000
        model_input.calculation_settings.md_steps_per_anchor = 100000
    elif calculation_type == "elber":
        model_input.calculation_type = "elber"
        model_input.calculation_settings = common_prepare.Elber_input_settings()
        model_input.calculation_settings.temperature_equil_progression = [310.0, 320.0]
        model_input.calculation_settings.num_temperature_equil_steps = 234
        model_input.calculation_settings.num_umbrella_stage_steps = 40000
        model_input.calculation_settings.fwd_rev_interval = 400
        model_input.calculation_settings.rev_output_interval = 100
        model_input.calculation_settings.fwd_output_interval = 100
        
    else:
        raise Exception("Option not available: %s", calculation_type)
    model_input.temperature = 277.8
    
    model_input.run_minimization = False
    cv_input1 = common_cv.Spherical_cv_input()
    cv_input1.group1 = list(range(147))
    cv_input1.group2 = list(range(147, 162))
    cv_input1.input_anchors = []
    input_anchor1 = common_cv.Spherical_cv_anchor()
    input_anchor1.radius = 0.1
    input_anchor1.starting_amber_params = base.Amber_params()
    input_anchor1.starting_amber_params.prmtop_filename = "data/hostguest.parm7"
    input_anchor1.starting_amber_params.inpcrd_filename = "data/hostguest.rst7"
    input_anchor1.starting_amber_params.pdb_coordinates_filename = ""
    input_anchor1.bound_state = True
    cv_input1.input_anchors.append(input_anchor1)
    
    input_anchor2 = common_cv.Spherical_cv_anchor()
    input_anchor2.radius = 0.2
    input_anchor2.starting_forcefield_params = base.Forcefield_params()
    input_anchor2.starting_forcefield_params.built_in_forcefield_filenames = \
        ["amber14/tip3pfb.xml"]
    input_anchor2.starting_forcefield_params.custom_forcefield_filenames = \
        ["data/hostguest.xml"]
    input_anchor2.starting_forcefield_params.pdb_filename = \
        "data/hostguest_for_xml.pdb"
    input_anchor2.bound_state = False
    cv_input1.input_anchors.append(input_anchor2)
    
    input_anchor3 = common_cv.Spherical_cv_anchor()
    input_anchor3.radius = 0.4
    input_anchor3.starting_amber_params = base.Amber_params()
    input_anchor3.starting_amber_params.prmtop_filename = "data/hostguest.parm7"
    input_anchor3.starting_amber_params.inpcrd_filename = "data/hostguest.rst7"
    input_anchor3.starting_amber_params.pdb_coordinates_filename = "data/hostguest_at3.8.pdb"
    input_anchor3.bound_state = False
    cv_input1.input_anchors.append(input_anchor3)
    
    input_anchor4 = common_cv.Spherical_cv_anchor()
    input_anchor4.radius = 0.5
    input_anchor4.starting_amber_params = base.Amber_params()
    input_anchor4.starting_amber_params.prmtop_filename = "data/hostguest.parm7"
    input_anchor4.starting_amber_params.inpcrd_filename = "data/hostguest.rst7"
    input_anchor4.starting_amber_params.pdb_coordinates_filename = "data/hostguest_at5.4.pdb"
    input_anchor4.bound_state = False
    cv_input1.input_anchors.append(input_anchor4)
    
    input_anchor5 = common_cv.Spherical_cv_anchor()
    input_anchor5.radius = 0.7
    input_anchor5.starting_amber_params = base.Amber_params()
    input_anchor5.starting_amber_params.prmtop_filename = "data/hostguest.parm7"
    input_anchor5.starting_amber_params.inpcrd_filename = "data/hostguest.rst7"
    input_anchor5.starting_amber_params.pdb_coordinates_filename = "data/hostguest_at7.0.pdb"
    input_anchor5.bound_state = False
    cv_input1.input_anchors.append(input_anchor5)
    
    input_anchor6 = common_cv.Spherical_cv_anchor()
    input_anchor6.radius = 0.8
    input_anchor6.starting_amber_params = base.Amber_params()
    input_anchor6.starting_amber_params.prmtop_filename = "data/hostguest.parm7"
    input_anchor6.starting_amber_params.inpcrd_filename = "data/hostguest.rst7"
    input_anchor6.starting_amber_params.pdb_coordinates_filename = "data/hostguest_at8.3.pdb"
    input_anchor6.bound_state = False
    cv_input1.input_anchors.append(input_anchor6)
    
    input_anchor7 = common_cv.Spherical_cv_anchor()
    input_anchor7.radius = 1.0
    input_anchor7.starting_amber_params = base.Amber_params()
    input_anchor7.starting_amber_params.prmtop_filename = "data/hostguest.parm7"
    input_anchor7.starting_amber_params.inpcrd_filename = "data/hostguest.rst7"
    input_anchor7.starting_amber_params.pdb_coordinates_filename = "data/hostguest_at9.8.pdb"
    input_anchor7.bound_state = False
    cv_input1.input_anchors.append(input_anchor7)
    
    input_anchor8 = common_cv.Spherical_cv_anchor()
    input_anchor8.radius = 1.1
    input_anchor8.starting_amber_params = base.Amber_params()
    input_anchor8.starting_amber_params.prmtop_filename = "data/hostguest.parm7"
    input_anchor8.starting_amber_params.inpcrd_filename = "data/hostguest.rst7"
    input_anchor8.starting_amber_params.pdb_coordinates_filename = "data/hostguest_at11.2.pdb"
    input_anchor8.bound_state = False
    cv_input1.input_anchors.append(input_anchor8)
    
    input_anchor9 = common_cv.Spherical_cv_anchor()
    input_anchor9.radius = 1.3
    input_anchor9.starting_amber_params = base.Amber_params()
    input_anchor9.starting_amber_params.prmtop_filename = "data/hostguest.parm7"
    input_anchor9.starting_amber_params.inpcrd_filename = "data/hostguest.rst7"
    input_anchor9.starting_amber_params.pdb_coordinates_filename = "data/hostguest_at12.8.pdb"
    input_anchor9.bound_state = False
    cv_input1.input_anchors.append(input_anchor9)
    
    input_anchor10 = common_cv.Spherical_cv_anchor()
    input_anchor10.radius = 1.4
    input_anchor10.starting_amber_params = base.Amber_params()
    input_anchor10.starting_amber_params.prmtop_filename = ""
    input_anchor10.starting_amber_params.inpcrd_filename = ""
    input_anchor10.starting_amber_params.pdb_coordinates_filename = ""
    input_anchor10.bound_state = False
    input_anchor10.bulk_anchor = True
    cv_input1.input_anchors.append(input_anchor10)
    
    model_input.cv_inputs = [cv_input1]
    
    model_input.browndye_settings = common_prepare.Browndye_settings_input()
    model_input.browndye_settings.receptor_pqr_filename = "data/hostguest_receptor.pqr"
    model_input.browndye_settings.ligand_pqr_filename = "data/hostguest_ligand.pqr"
    model_input.browndye_settings.apbs_grid_spacing = 0.5
    model_input.browndye_settings.receptor_indices = list(range(147))
    model_input.browndye_settings.ligand_indices = list(range(15))
    model_input.browndye_settings.binary_directory = ""
    
    ion1 = base.Ion()
    ion1.radius = 1.2
    ion1.charge = -1.0
    ion1.conc = 0.15
    ion2 = base.Ion()
    ion2.radius = 0.9
    ion2.charge = 1.0
    ion2.conc = 0.15
    model_input.browndye_settings.ions = [ion1, ion2]
    model_input.browndye_settings.num_bd_milestone_steps = 1000
    model_input.browndye_settings.num_b_surface_steps = 10000
    model_input.browndye_settings.n_threads = 1
    
    return model_input

def create_sample_input_model_with_spherical_cv():
    """
    Create an input_model object for a concentric spherical CV system.
    """
    model_input = common_prepare.Model_input()
    model_input.root_directory = "/path/to/directory/"
    cv_input1 = common_cv.Spherical_cv_input()
    cv_input1.index = 0
    cv_input1.group1 = [2,5,8]
    cv_input1.group2 = [3,6,9]
    cv_input1.radii = [1.5, 3.0, 4.5, 6.0]
    cv_input1.end_state_radii = [1.5, 6.0]
    cv_input1.last_radius_is_bulk = True
    
    model_input.collective_variable_inputs = [cv_input1]
    return model_input

def write_input_xml(filename):
    """
    Write the XML file of this sample model
    """
    
    model_input = create_sample_input_model_with_spherical_cv()
    model_input.serialize(filename)
    return

if __name__ == "__main__":
    assert len(sys.argv) > 1, "Please enter a filename for the input xml."
    dest_filename = sys.argv[1]
    root_dir = os.path.abspath(os.path.dirname(dest_filename))
    #model_input = create_model_input(root_dir)
    model_input = create_model_input(root_dir, calculation_type="elber")
    model_input.serialize(dest_filename)
