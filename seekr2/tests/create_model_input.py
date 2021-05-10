"""
create_model_input.py

Create generic Model input objects
"""

import seekr2.modules.common_base as base
import seekr2.modules.common_prepare as common_prepare
import seekr2.modules.common_cv as common_cv

def create_host_guest_mmvt_model_input(root_dir):
    """
    Create a generic host-guest model input object.
    """
    model_input = common_prepare.Model_input()
    model_input.calculation_type = "mmvt"
    model_input.calculation_settings = common_prepare.MMVT_input_settings()
    model_input.calculation_settings.md_output_interval = 10000
    model_input.calculation_settings.md_steps_per_anchor = 100000
    model_input.temperature = 277.8
    model_input.pressure = 1.0
    model_input.ensemble = "nvt"
    model_input.root_directory = root_dir
    model_input.md_program = "openmm"
    model_input.constraints = "HBonds"
    model_input.rigidWater = True
    model_input.hydrogenMass = None
    model_input.timestep = 0.002
    model_input.nonbonded_cutoff = 0.9
    cv_input1 = common_cv.Spherical_cv_input()
    cv_input1.group1 = list(range(147))
    cv_input1.group2 = list(range(147, 162))
    cv_input1.input_anchors = []
    
    input_anchor1 = common_cv.Spherical_cv_anchor()
    input_anchor1.radius = 0.1
    input_anchor1.starting_amber_params = base.Amber_params()
    input_anchor1.starting_amber_params.prmtop_filename \
        = "../data/hostguest_files/hostguest.parm7"
    input_anchor1.starting_amber_params.inpcrd_filename = \
        "../data/hostguest_files/hostguest.rst7"
    input_anchor1.starting_amber_params.pdb_coordinates_filename = \
        "../data/hostguest_files/hostguest_at0.4.pdb"
    input_anchor1.bound_state = True
    input_anchor1.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor1)
    
    input_anchor2 = common_cv.Spherical_cv_anchor()
    input_anchor2.radius = 0.2
    input_anchor2.starting_amber_params = base.Amber_params()
    input_anchor2.starting_amber_params.prmtop_filename = \
        "../data/hostguest_files/hostguest.parm7"
    input_anchor2.starting_amber_params.inpcrd_filename = \
        "../data/hostguest_files/hostguest.rst7"
    input_anchor2.starting_amber_params.pdb_coordinates_filename = \
        "../data/hostguest_files/hostguest_at1.9.pdb"
    input_anchor2.bound_state = False
    input_anchor2.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor2)
    
    input_anchor3 = common_cv.Spherical_cv_anchor()
    input_anchor3.radius = 0.4
    input_anchor3.starting_amber_params = base.Amber_params()
    input_anchor3.starting_amber_params.prmtop_filename = \
        "../data/hostguest_files/hostguest.parm7"
    input_anchor3.starting_amber_params.inpcrd_filename = \
        "../data/hostguest_files/hostguest.rst7"
    input_anchor3.starting_amber_params.pdb_coordinates_filename = \
        "../data/hostguest_files/hostguest_at3.8.pdb"
    input_anchor3.bound_state = False
    input_anchor3.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor3)
    
    input_anchor4 = common_cv.Spherical_cv_anchor()
    input_anchor4.radius = 0.5
    input_anchor4.starting_amber_params = base.Amber_params()
    input_anchor4.starting_amber_params.prmtop_filename = \
        "../data/hostguest_files/hostguest.parm7"
    input_anchor4.starting_amber_params.inpcrd_filename = \
        "../data/hostguest_files/hostguest.rst7"
    input_anchor4.starting_amber_params.pdb_coordinates_filename \
        = "../data/hostguest_files/hostguest_at5.4.pdb"
    input_anchor4.bound_state = False
    input_anchor4.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor4)
    
    input_anchor5 = common_cv.Spherical_cv_anchor()
    input_anchor5.radius = 0.7
    input_anchor5.starting_amber_params = base.Amber_params()
    input_anchor5.starting_amber_params.prmtop_filename \
        = "../data/hostguest_files/hostguest.parm7"
    input_anchor5.starting_amber_params.inpcrd_filename \
        = "../data/hostguest_files/hostguest.rst7"
    input_anchor5.starting_amber_params.pdb_coordinates_filename \
        = "../data/hostguest_files/hostguest_at7.0.pdb"
    input_anchor5.bound_state = False
    input_anchor5.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor5)
    
    input_anchor6 = common_cv.Spherical_cv_anchor()
    input_anchor6.radius = 0.8
    input_anchor6.starting_amber_params = base.Amber_params()
    input_anchor6.starting_amber_params.prmtop_filename \
        = "../data/hostguest_files/hostguest.parm7"
    input_anchor6.starting_amber_params.inpcrd_filename \
        = "../data/hostguest_files/hostguest.rst7"
    input_anchor6.starting_amber_params.pdb_coordinates_filename \
        = "../data/hostguest_files/hostguest_at8.3.pdb"
    input_anchor6.bound_state = False
    input_anchor6.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor6)
    
    input_anchor7 = common_cv.Spherical_cv_anchor()
    input_anchor7.radius = 1.0
    input_anchor7.starting_amber_params = base.Amber_params()
    input_anchor7.starting_amber_params.prmtop_filename \
        = "../data/hostguest_files/hostguest.parm7"
    input_anchor7.starting_amber_params.inpcrd_filename \
        = "../data/hostguest_files/hostguest.rst7"
    input_anchor7.starting_amber_params.pdb_coordinates_filename \
        = "../data/hostguest_files/hostguest_at9.8.pdb"
    input_anchor7.bound_state = False
    input_anchor7.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor7)
    
    input_anchor8 = common_cv.Spherical_cv_anchor()
    input_anchor8.radius = 1.1
    input_anchor8.starting_amber_params = base.Amber_params()
    input_anchor8.starting_amber_params.prmtop_filename \
        = "../data/hostguest_files/hostguest.parm7"
    input_anchor8.starting_amber_params.inpcrd_filename \
        = "../data/hostguest_files/hostguest.rst7"
    input_anchor8.starting_amber_params.pdb_coordinates_filename \
        = "../data/hostguest_files/hostguest_at11.2.pdb"
    input_anchor8.bound_state = False
    input_anchor8.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor8)
    
    input_anchor9 = common_cv.Spherical_cv_anchor()
    input_anchor9.radius = 1.3
    input_anchor9.starting_amber_params = base.Amber_params()
    input_anchor9.starting_amber_params.prmtop_filename \
        = "../data/hostguest_files/hostguest.parm7"
    input_anchor9.starting_amber_params.inpcrd_filename \
        = "../data/hostguest_files/hostguest.rst7"
    input_anchor9.starting_amber_params.pdb_coordinates_filename \
        = "../data/hostguest_files/hostguest_at12.8.pdb"
    input_anchor9.bound_state = False
    input_anchor9.bulk_anchor = False
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
    model_input.browndye_settings_input \
        = common_prepare.Browndye_settings_input()
    model_input.browndye_settings_input.binary_directory = ""
    model_input.browndye_settings_input.receptor_pqr_filename \
        = "../data/hostguest_files/hostguest_receptor.pqr"
    model_input.browndye_settings_input.ligand_pqr_filename \
        = "../data/hostguest_files/hostguest_ligand.pqr"
    model_input.browndye_settings_input.apbs_grid_spacing = 0.5
    model_input.browndye_settings_input.receptor_indices = list(range(147))
    model_input.browndye_settings_input.ligand_indices = list(range(15))
    
    ion1 = base.Ion()
    ion1.radius = 1.2
    ion1.charge = -1.0
    ion1.conc = 0.0
    ion2 = base.Ion()
    ion2.radius = 0.9
    ion2.charge = 1.0
    ion2.conc = 0.0
    model_input.browndye_settings_input.ions = [ion1, ion2]
    model_input.browndye_settings_input.num_bd_milestone_trajectories = 1000
    model_input.browndye_settings_input.num_b_surface_trajectories = 10000
    model_input.browndye_settings_input.max_b_surface_trajs_to_extract = 1000
    model_input.browndye_settings_input.n_threads = 1
    
    return model_input

def create_host_guest_elber_model_input(root_dir):
    """
    Create a generic host-guest model input object.
    """
    model_input = create_host_guest_mmvt_model_input(root_dir)
    model_input.calculation_type = "elber"
    model_input.calculation_settings = common_prepare.Elber_input_settings()
    model_input.calculation_settings.num_umbrella_stage_steps = 100000
    model_input.calculation_settings.umbrella_force_constant = 9000.0
    model_input.calculation_settings.fwd_rev_interval = 500
    model_input.calculation_settings.rev_output_interval = None
    model_input.calculation_settings.fwd_output_interval = None
    return model_input