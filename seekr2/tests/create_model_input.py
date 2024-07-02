"""
create_model_input.py

Create Model input objects for use in tests.
"""

import os

import numpy as np

import seekr2.modules.common_base as base
import seekr2.modules.common_prepare as common_prepare
import seekr2.modules.common_cv as common_cv

TEST_DIRECTORY = os.path.dirname(__file__)

def assign_amber_params(input_anchor, prmtop_filename, pdb_filename):
    input_anchor.starting_amber_params = base.Amber_params()
    input_anchor.starting_amber_params.prmtop_filename = prmtop_filename
    input_anchor.starting_amber_params.pdb_coordinates_filename = pdb_filename
    return

def assign_forcefield_params(input_anchor, built_in_ff_list, custom_ff_list, 
                             pdb_filename):
    input_anchor.starting_forcefield_params = base.Forcefield_params()
    input_anchor.starting_forcefield_params.built_in_forcefield_filenames \
        = built_in_ff_list
    input_anchor.starting_forcefield_params.custom_forcefield_filenames \
        = custom_ff_list
    input_anchor.starting_forcefield_params.pdb_coordinates_filename \
        = pdb_filename
    return

def assign_system_params(input_anchor, system_filename, pdb_filename):
    input_anchor.starting_forcefield_params = base.Forcefield_params()
    input_anchor.starting_forcefield_params.system_filename \
        = system_filename
    input_anchor.starting_forcefield_params.pdb_coordinates_filename \
        = pdb_filename
    return

def assign_charmm_params(input_anchor, psf_filename, charmm_ff_filenames, 
                         pdb_filename):
    input_anchor.starting_charmm_params = base.Charmm_params()
    input_anchor.starting_charmm_params.psf_filename = psf_filename
    input_anchor.starting_charmm_params.charmm_ff_files = charmm_ff_filenames
    input_anchor.starting_charmm_params.pdb_coordinates_filename = pdb_filename
    return

def create_host_guest_mmvt_model_input(root_dir, bd=True, ff="amber"):
    """
    Create a generic host-guest model input object.
    """
    os.chdir(TEST_DIRECTORY)
    model_input = common_prepare.Model_input()
    model_input.calculation_type = "mmvt"
    model_input.calculation_settings = common_prepare.MMVT_input_settings()
    model_input.calculation_settings.md_output_interval = 10000
    model_input.calculation_settings.md_steps_per_anchor = 100000 #1000000
    model_input.temperature = 298.15
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
    
    radius_list = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95,
                   1.05, 1.15, 1.25, 1.35]
    amber_prmtop_filename = os.path.abspath(
        "../data/hostguest_files/hostguest.parm7")
    forcefield_built_in_ff_list = ["amber14/tip3pfb.xml"]
    forcefield_custom_ff_list = [os.path.abspath(
        "../data/hostguest_files/hostguest.xml")]
    system_filename = os.path.abspath(
        "../data/hostguest_files/hostguest_system.xml")
    pdb_filenames = [os.path.abspath(
        "../data/hostguest_files/hostguest_at0.5.pdb"),
                     os.path.abspath(
                         "../data/hostguest_files/hostguest_at1.5.pdb"),
                     os.path.abspath(
                         "../data/hostguest_files/hostguest_at2.5.pdb"),
                     os.path.abspath(
                         "../data/hostguest_files/hostguest_at3.5.pdb"),
                     os.path.abspath(
                         "../data/hostguest_files/hostguest_at4.5.pdb"),
                     os.path.abspath(
                         "../data/hostguest_files/hostguest_at5.5.pdb"),
                     os.path.abspath(
                         "../data/hostguest_files/hostguest_at6.5.pdb"),
                     os.path.abspath(
                         "../data/hostguest_files/hostguest_at7.5.pdb"),
                     os.path.abspath(
                         "../data/hostguest_files/hostguest_at8.5.pdb"),
                     os.path.abspath(
                         "../data/hostguest_files/hostguest_at9.5.pdb"),
                     os.path.abspath(
                         "../data/hostguest_files/hostguest_at10.5.pdb"),
                     os.path.abspath(
                         "../data/hostguest_files/hostguest_at11.5.pdb"),
                     os.path.abspath(
                         "../data/hostguest_files/hostguest_at12.5.pdb"),
                     ""]
    for i, (radius, pdb_filename) in enumerate(zip(radius_list, pdb_filenames)):
        input_anchor = common_cv.Spherical_cv_anchor()
        input_anchor.radius = radius
        if ff == "amber":
            assign_amber_params(input_anchor, amber_prmtop_filename, 
                                pdb_filename)
        elif ff == "forcefield":
            print("test pdb_filename:", pdb_filename)
            assign_forcefield_params(input_anchor, forcefield_built_in_ff_list, 
                                     forcefield_custom_ff_list, pdb_filename)
        elif ff == "system":
            print("test pdb_filename:", pdb_filename)
            assign_system_params(input_anchor, system_filename, pdb_filename)
        else:
            raise Exception("ff type not supported: {}".format(ff))
        
        if i == 0:
            input_anchor.bound_state = True
        else:
            input_anchor.bound_state = False
            
        if i == len(radius_list)-1:
            input_anchor.bulk_anchor = True
        else:
            input_anchor.bulk_anchor = False
    
        cv_input1.input_anchors.append(input_anchor)
    
    model_input.cv_inputs = [cv_input1]
    
    if bd:
        model_input.browndye_settings_input \
            = common_prepare.Browndye_settings_input()
        model_input.browndye_settings_input.binary_directory = ""
        model_input.browndye_settings_input.receptor_pqr_filename \
            = os.path.abspath("../data/hostguest_files/hostguest_receptor.pqr")
        model_input.browndye_settings_input.ligand_pqr_filename \
            = os.path.abspath("../data/hostguest_files/hostguest_ligand.pqr")
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
        #model_input.browndye_settings_input.num_bd_milestone_trajectories = 100
        model_input.browndye_settings_input.num_b_surface_trajectories = 100
        #model_input.browndye_settings_input.max_b_surface_trajs_to_extract = 100
        model_input.browndye_settings_input.n_threads = 1
    else:
        model_input.browndye_settings_input = None
    
    return model_input

def create_host_guest_elber_model_input(root_dir, bd=True):
    """
    Create a generic host-guest model input object.
    """
    model_input = create_host_guest_mmvt_model_input(root_dir, bd=bd)
    model_input.calculation_type = "elber"
    model_input.calculation_settings = common_prepare.Elber_input_settings()
    model_input.calculation_settings.num_umbrella_stage_steps = 100000
    model_input.calculation_settings.umbrella_force_constant = 9000.0
    model_input.calculation_settings.fwd_rev_interval = 500
    model_input.calculation_settings.rev_output_interval = None
    model_input.calculation_settings.fwd_output_interval = None
    return model_input

def create_tiwary_mmvt_model_input(root_dir):
    """
    Create a generic host-guest model input object.
    """
    os.chdir(TEST_DIRECTORY)
    model_input = common_prepare.Model_input()
    model_input.calculation_type = "mmvt"
    model_input.calculation_settings = common_prepare.MMVT_input_settings()
    model_input.calculation_settings.md_output_interval = 10000
    model_input.calculation_settings.md_steps_per_anchor = 100000 #1000000
    model_input.temperature = 298.15
    model_input.pressure = 1.0
    model_input.ensemble = "nvt"
    model_input.root_directory = root_dir
    model_input.md_program = "openmm"
    model_input.constraints = "HBonds"
    model_input.rigidWater = True
    model_input.hydrogenMass = None
    model_input.timestep = 0.002
    model_input.nonbonded_cutoff = 0.9
    
    order_parameter1 = common_cv.Tiwary_cv_distance_order_parameter()
    order_parameter1.group1 = [
                    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 
                    17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 
                    31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 
                    45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 
                    59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 
                    73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 
                    87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 
                    101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 
                    112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 
                    123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 
                    134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 
                    145, 146]
    order_parameter1.group2 = [
                    147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 
                    158, 159, 160, 161]
    order_parameter2 = common_cv.Tiwary_cv_angle_order_parameter()
    order_parameter2.group1 = [0]
    order_parameter2.group2 = [147]
    order_parameter2.group3 = [161]
    order_parameter3 = common_cv.Tiwary_cv_angle_order_parameter()
    order_parameter3.group1 = [0]
    order_parameter3.group2 = [1]
    order_parameter3.group3 = [147]
    order_parameter3.group4 = [161]
    order_parameters = [order_parameter1, order_parameter2, order_parameter3]
    order_parameter_weights1 = [0.99, 0.005, 0.005]
    order_parameter_weights2 = [0.97, 0.015, 0.015]
    
    cv_input1 = common_cv.Tiwary_cv_input()
    cv_input1.order_parameters = order_parameters
    cv_input1.order_parameter_weights = order_parameter_weights1
    cv_input1.input_anchors = []
    
    cv_input2 = common_cv.Tiwary_cv_input()
    cv_input2.order_parameters = order_parameters
    cv_input2.order_parameter_weights = order_parameter_weights2
    cv_input2.input_anchors = []
    
    cv_input3 = common_cv.Spherical_cv_input()
    cv_input3.group1 = list(range(147))
    cv_input3.group2 = list(range(147, 162))
    cv_input3.input_anchors = []
    
    values_list1 = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55]
    values_list2 = [0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15]
    spherical_radii = [1.15, 1.25, 1.35]
    
    amber_prmtop_filename = os.path.abspath(
        "../data/hostguest_files/hostguest.parm7")
    pdb_filenames1 = [os.path.abspath(
        "../data/hostguest_files/hostguest_at0.5.pdb"),
                     os.path.abspath(
                         "../data/hostguest_files/hostguest_at1.5.pdb"),
                     os.path.abspath(
                         "../data/hostguest_files/hostguest_at2.5.pdb"),
                     os.path.abspath(
                         "../data/hostguest_files/hostguest_at3.5.pdb"),
                     os.path.abspath(
                         "../data/hostguest_files/hostguest_at4.5.pdb"),
                     os.path.abspath(
                         "../data/hostguest_files/hostguest_at5.5.pdb"),]
                     
    pdb_filenames2 = ["",
                      os.path.abspath(
                         "../data/hostguest_files/hostguest_at6.5.pdb"),
                     os.path.abspath(
                         "../data/hostguest_files/hostguest_at7.5.pdb"),
                     os.path.abspath(
                         "../data/hostguest_files/hostguest_at8.5.pdb"),
                     os.path.abspath(
                         "../data/hostguest_files/hostguest_at9.5.pdb"),
                     os.path.abspath(
                         "../data/hostguest_files/hostguest_at10.5.pdb"),
                     os.path.abspath(
                         "../data/hostguest_files/hostguest_at11.5.pdb"),]

    pdb_filenames3 = ["",
                      os.path.abspath(
                         "../data/hostguest_files/hostguest_at12.5.pdb"),
                     ""]
    for i, (value, pdb_filename) in enumerate(zip(values_list1, pdb_filenames1)):
        input_anchor = common_cv.Tiwary_cv_anchor()
        input_anchor.value = value
        assign_amber_params(input_anchor, amber_prmtop_filename, 
                            pdb_filename)
        if i == 0:
            input_anchor.bound_state = True
        else:
            input_anchor.bound_state = False
            
        
        input_anchor.bulk_anchor = False
    
        cv_input1.input_anchors.append(input_anchor)
    
    cv_input1.input_anchors[-1].connection_flags = [1]
    
    for i, (value, pdb_filename) in enumerate(zip(values_list2, pdb_filenames2)):
        input_anchor = common_cv.Tiwary_cv_anchor()
        input_anchor.value = value
        assign_amber_params(input_anchor, amber_prmtop_filename, 
                            pdb_filename)
        input_anchor.bound_state = False
        input_anchor.bulk_anchor = False
        cv_input2.input_anchors.append(input_anchor)
    
    cv_input2.input_anchors[0].connection_flags = [1]
    cv_input2.input_anchors[-1].connection_flags = [2]
        
    for i, (radius, pdb_filename) in enumerate(zip(spherical_radii, pdb_filenames3)):
        input_anchor = common_cv.Spherical_cv_anchor()
        input_anchor.radius = radius
        assign_amber_params(input_anchor, amber_prmtop_filename, 
                            pdb_filename)
        input_anchor.bound_state = False
            
        if i == len(pdb_filenames3)-1:
            input_anchor.bulk_anchor = True
        else:
            input_anchor.bulk_anchor = False
    
        cv_input3.input_anchors.append(input_anchor)
    
    cv_input3.input_anchors[0].connection_flags = [2]
    
    model_input.cv_inputs = [cv_input1, cv_input2, cv_input3]
    
    return model_input

def create_planar_mmvt_model_input(root_dir):
    """
    Create a generic host-guest model input object.
    """
    os.chdir(TEST_DIRECTORY)
    model_input = common_prepare.Model_input()
    model_input.calculation_type = "mmvt"
    model_input.calculation_settings = common_prepare.MMVT_input_settings()
    model_input.calculation_settings.md_output_interval = 10000
    model_input.calculation_settings.md_steps_per_anchor = 100000 #1000000
    model_input.temperature = 298.15
    model_input.pressure = 1.0
    model_input.ensemble = "nvt"
    model_input.root_directory = root_dir
    model_input.md_program = "openmm"
    model_input.constraints = "HBonds"
    model_input.rigidWater = True
    model_input.hydrogenMass = None
    model_input.timestep = 0.002
    model_input.nonbonded_cutoff = 0.9
    cv_input1 = common_cv.Planar_cv_input()
    cv_input1.start_group = [23]
    cv_input1.end_group = [112]
    cv_input1.mobile_group = list(range(147, 162))
    cv_input1.input_anchors = []
    
    planar_list = [0.0, 0.5, 1.5]
    amber_prmtop_filename = os.path.abspath(
        "../data/hostguest_files/hostguest.parm7")
    pdb_filename = os.path.abspath(
        "../data/hostguest_files/hostguest_at0.5.pdb")
    for i, value in enumerate(planar_list):
        input_anchor = common_cv.Planar_cv_anchor()
        input_anchor.value = value
        if i == 1:
            assign_amber_params(input_anchor, amber_prmtop_filename, 
                                pdb_filename)
        input_anchor.bound_state = False
        input_anchor.bulk_anchor = False
        cv_input1.input_anchors.append(input_anchor)
    
    model_input.cv_inputs = [cv_input1]
    model_input.browndye_settings_input = None
    return model_input

def create_rmsd_mmvt_model_input(root_dir):
    """
    Create a generic host-guest model input object.
    """
    os.chdir(TEST_DIRECTORY)
    model_input = common_prepare.Model_input()
    model_input.calculation_type = "mmvt"
    model_input.calculation_settings = common_prepare.MMVT_input_settings()
    model_input.calculation_settings.md_output_interval = 10000
    model_input.calculation_settings.md_steps_per_anchor = 100000 #1000000
    model_input.temperature = 298.15
    model_input.pressure = 1.0
    model_input.ensemble = "nvt"
    model_input.root_directory = root_dir
    model_input.md_program = "openmm"
    model_input.constraints = "HBonds"
    model_input.rigidWater = True
    model_input.hydrogenMass = None
    model_input.timestep = 0.002
    model_input.nonbonded_cutoff = 0.9
    cv_input1 = common_cv.RMSD_cv_input()
    cv_input1.group = list(range(147, 162))
    pdb_filename = os.path.abspath(
        "../data/hostguest_files/hostguest_at0.5.pdb")
    cv_input1.ref_structure = pdb_filename
    cv_input1.input_anchors = []
    
    rmsd_list = [0.0, 0.5, 1.0]
    amber_prmtop_filename = os.path.abspath(
        "../data/hostguest_files/hostguest.parm7")
    for i, value in enumerate(rmsd_list):
        input_anchor = common_cv.RMSD_cv_anchor()
        input_anchor.value = value
        if i == 0:
            assign_amber_params(input_anchor, amber_prmtop_filename, 
                                pdb_filename)
        input_anchor.bound_state = False
        input_anchor.bulk_anchor = False
    
        cv_input1.input_anchors.append(input_anchor)
    
    model_input.cv_inputs = [cv_input1]
    model_input.browndye_settings_input = None
    return model_input

def create_closest_pair_mmvt_model_input(root_dir):
    """
    Create a generic host-guest model input object.
    """
    os.chdir(TEST_DIRECTORY)
    model_input = common_prepare.Model_input()
    model_input.calculation_type = "mmvt"
    model_input.calculation_settings = common_prepare.MMVT_input_settings()
    model_input.calculation_settings.md_output_interval = 10000
    model_input.calculation_settings.md_steps_per_anchor = 100000 #1000000
    model_input.temperature = 298.15
    model_input.pressure = 1.0
    model_input.ensemble = "nvt"
    model_input.root_directory = root_dir
    model_input.md_program = "openmm"
    model_input.constraints = "HBonds"
    model_input.rigidWater = True
    model_input.hydrogenMass = None
    model_input.timestep = 0.002
    model_input.nonbonded_cutoff = 0.9
    cv_input1 = common_cv.Closest_pair_cv_input()
    cv_input1.group1 = [150]
    cv_input1.group2 = list(range(162, 5358, 3))
    cv_input1.exponent = 50
    pdb_filename = os.path.abspath(
        "../data/hostguest_files/hostguest_at0.5.pdb")
    cv_input1.input_anchors = []
    
    closest_list = [0.0, 1.0, 2.0]
    amber_prmtop_filename = os.path.abspath(
        "../data/hostguest_files/hostguest.parm7")
    for i, value in enumerate(closest_list):
        input_anchor = common_cv.Closest_pair_cv_anchor()
        input_anchor.value = value
        if i == 0:
            assign_amber_params(input_anchor, amber_prmtop_filename, 
                                pdb_filename)
        input_anchor.bound_state = False
        input_anchor.bulk_anchor = False
    
        cv_input1.input_anchors.append(input_anchor)
    
    model_input.cv_inputs = [cv_input1]
    model_input.browndye_settings_input = None
    return model_input

def create_count_contacts_mmvt_model_input(root_dir):
    """
    Create a generic host-guest model input object.
    """
    os.chdir(TEST_DIRECTORY)
    model_input = common_prepare.Model_input()
    model_input.calculation_type = "mmvt"
    model_input.calculation_settings = common_prepare.MMVT_input_settings()
    model_input.calculation_settings.md_output_interval = 10000
    model_input.calculation_settings.md_steps_per_anchor = 100000 #1000000
    model_input.temperature = 298.15
    model_input.pressure = 1.0
    model_input.ensemble = "nvt"
    model_input.root_directory = root_dir
    model_input.md_program = "openmm"
    model_input.constraints = "HBonds"
    model_input.rigidWater = True
    model_input.hydrogenMass = None
    model_input.timestep = 0.002
    model_input.nonbonded_cutoff = 0.9
    cv_input1 = common_cv.Count_contacts_cv_input()
    cv_input1.group1 = [150]
    cv_input1.group2 = list(range(162, 5358, 3))
    cv_input1.cutoff_distance = 0.3
    pdb_filename = os.path.abspath(
        "../data/hostguest_files/hostguest_at0.5.pdb")
    cv_input1.input_anchors = []
    
    closest_list = [2, 3, 4]
    amber_prmtop_filename = os.path.abspath(
        "../data/hostguest_files/hostguest.parm7")
    for i, value in enumerate(closest_list):
        input_anchor = common_cv.Count_contacts_cv_anchor()
        input_anchor.value = value
        if i == 0:
            assign_amber_params(input_anchor, amber_prmtop_filename, 
                                pdb_filename)
        input_anchor.bound_state = False
        input_anchor.bulk_anchor = False
    
        cv_input1.input_anchors.append(input_anchor)
    
    model_input.cv_inputs = [cv_input1]
    model_input.browndye_settings_input = None
    return model_input

def create_toy_mmvt_model_input(root_dir):
    """
    Create a toy mmvt model input.
    
    This is the Entropy barrier system.
    """
    os.chdir(TEST_DIRECTORY)
    model_input = common_prepare.Model_input()
    model_input.calculation_type = "mmvt"
    model_input.calculation_settings = common_prepare.MMVT_input_settings()
    model_input.calculation_settings.md_output_interval = 1000
    model_input.calculation_settings.md_steps_per_anchor = 100000
    model_input.temperature = 300.0
    model_input.pressure = 1.0
    model_input.ensemble = "nvt"
    model_input.root_directory = root_dir
    model_input.md_program = "openmm"
    model_input.constraints = "none"
    model_input.rigidWater = True
    model_input.hydrogenMass = None
    model_input.integrator_type = "langevin"
    model_input.timestep = 0.002
    model_input.nonbonded_cutoff = None
    
    cv_input1 = common_cv.Toy_cv_input()
    cv_input1.groups = [[0]]
    cv_input1.variable_name = "value"
    cv_input1.cv_expression = "y1"
    #cv_input1.openmm_expression = "step(k*(y1-value))"
    cv_input1.input_anchors = []
    
    values_list = [-0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7]
    positions_list = [
        [[[0.0, -0.7, 0.0]], [[0.3, -0.7, 0.0]]],
        [[[0.0, -0.5, 0.0]]],
        [[[0.0, -0.3, 0.0]]],
        [[[0.0, -0.1, 0.0]]],
        [[[0.0, 0.1, 0.0]]],
        [[[0.0, 0.3, 0.0]]],
        [[[0.0, 0.5, 0.0]]],
        [[[0.0, 0.7, 0.0]]]]
        
    
    for i, (value, position) in enumerate(zip(values_list, positions_list)):
        input_anchor = common_cv.Toy_cv_anchor()
        input_anchor.value = value
        input_anchor.starting_positions = np.array(position)
        
        if i == 0:
            input_anchor.bound_state = True
        else:
            input_anchor.bound_state = False
            
        if i == len(values_list)-1:
            input_anchor.bulk_anchor = True
        else:
            input_anchor.bulk_anchor = False
        #input_anchor.bulk_anchor = False
    
        cv_input1.input_anchors.append(input_anchor)
    
    model_input.cv_inputs = [cv_input1]
    model_input.browndye_settings_input = None
    model_input.toy_settings_input = common_prepare.Toy_settings_input()
    model_input.toy_settings_input.potential_energy_expression \
        = "5*(x1^6+y1^6+exp(-(10*y1)^2)*(1-exp(-(10*x1)^2)))"
    model_input.toy_settings_input.num_particles = 1
    model_input.toy_settings_input.masses = np.array([10.0])
    
    return model_input

def create_toy_elber_model_input(root_dir):
    """
    Create a generic host-guest model input object.
    """
    model_input = create_toy_mmvt_model_input(root_dir)
    model_input.calculation_type = "elber"
    model_input.calculation_settings = common_prepare.Elber_input_settings()
    model_input.calculation_settings.num_umbrella_stage_steps = 100000
    model_input.calculation_settings.umbrella_force_constant = 9000.0
    model_input.calculation_settings.fwd_rev_interval = 500
    model_input.calculation_settings.rev_output_interval = None
    model_input.calculation_settings.fwd_output_interval = None
    #model_input.cv_inputs[0].restraining_expression = "0.5*k*(y1-value)^2"
    model_input.cv_inputs[0].input_anchors[0].starting_positions = np.array(
        [[[0.0, -0.7, 0.0]]])
    return model_input

def create_toy_multi_model_input(root_dir):
    """
    Create a toy mmvt multidimensional model input.
    
    This is the Entropy barrier system.
    """
    os.chdir(TEST_DIRECTORY)
    model_input = common_prepare.Model_input()
    model_input.calculation_type = "mmvt"
    model_input.calculation_settings = common_prepare.MMVT_input_settings()
    model_input.calculation_settings.md_output_interval = 1000
    model_input.calculation_settings.md_steps_per_anchor = 100000
    model_input.temperature = 300.0
    model_input.pressure = 1.0
    model_input.ensemble = "nvt"
    model_input.root_directory = root_dir
    model_input.md_program = "openmm"
    model_input.constraints = "none"
    model_input.rigidWater = True
    model_input.hydrogenMass = None
    model_input.integrator_type = "langevin"
    model_input.timestep = 0.002
    model_input.nonbonded_cutoff = None
    
    combo = common_cv.Grid_combo()
    
    cv_input1 = common_cv.Toy_cv_input()
    cv_input1.groups = [[0]]
    cv_input1.variable_name = "value"
    cv_input1.cv_expression = "x1"
    cv_input1.openmm_expression = "step(k*(x1-value))"
    cv_input1.input_anchors = []
    
    cv_input2 = common_cv.Toy_cv_input()
    cv_input2.groups = [[0]]
    cv_input2.variable_name = "value"
    cv_input2.cv_expression = "y1"
    cv_input2.openmm_expression = "step(k*(y1-value))"
    cv_input2.input_anchors = []
    
    values_list = [-0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7]
    
    for i, value in enumerate(values_list):
        input_anchor = common_cv.Toy_cv_anchor()
        input_anchor.value = value
        input_anchor.starting_positions = None
        input_anchor.bound_state = False
        input_anchor.bulk_anchor = False    
        cv_input1.input_anchors.append(input_anchor)
        
    for i, value in enumerate(values_list):
        input_anchor = common_cv.Toy_cv_anchor()
        input_anchor.value = value
        input_anchor.starting_positions = None
        input_anchor.bound_state = False
        input_anchor.bulk_anchor = False    
        cv_input2.input_anchors.append(input_anchor)
    
    combo.cv_inputs = [cv_input1, cv_input2]
    stateA = common_cv.State_point()
    stateA.name = "stateA"
    stateA.location = [[0.1, -0.7, 0.0]]
    stateB = common_cv.State_point()
    stateB.name = "stateB"
    stateB.location = [[0.1, 0.7, 0.0]]
    combo.state_points = [stateA, stateB]
    
    model_input.cv_inputs = [combo]
    model_input.browndye_settings_input = None
    model_input.toy_settings_input = common_prepare.Toy_settings_input()
    model_input.toy_settings_input.potential_energy_expression \
        = "5*(x1^6+y1^6+exp(-(10*y1)^2)*(1-exp(-(10*x1)^2)))"
    model_input.toy_settings_input.num_particles = 1
    model_input.toy_settings_input.masses = np.array([10.0])
    
    return model_input

def create_toy_voronoi_model_input(root_dir):
    """
    Create a toy mmvt voronoi model input.
    
    This is the Entropy barrier system.
    """
    os.chdir(TEST_DIRECTORY)
    model_input = common_prepare.Model_input()
    model_input.calculation_type = "mmvt"
    model_input.calculation_settings = common_prepare.MMVT_input_settings()
    model_input.calculation_settings.md_output_interval = 1000
    model_input.calculation_settings.md_steps_per_anchor = 100000
    model_input.temperature = 300.0
    model_input.pressure = 1.0
    model_input.ensemble = "nvt"
    model_input.root_directory = root_dir
    model_input.md_program = "openmm"
    model_input.constraints = "none"
    model_input.rigidWater = True
    model_input.hydrogenMass = None
    model_input.integrator_type = "langevin"
    model_input.timestep = 0.002
    model_input.nonbonded_cutoff = None
    
    voronoi_cv = common_cv.Voronoi_cv_input()
    
    cv_input1 = common_cv.Toy_cv_input()
    cv_input1.groups = [[0]]
    cv_input1.variable_name = "value"
    cv_input1.cv_expression = "x1"
    cv_input1.openmm_expression = "step(k*(x1-value))"
    cv_input1.input_anchors = []
    
    cv_input2 = common_cv.Toy_cv_input()
    cv_input2.groups = [[0]]
    cv_input2.variable_name = "value"
    cv_input2.cv_expression = "y1"
    cv_input2.openmm_expression = "step(k*(y1-value))"
    cv_input2.input_anchors = []
    
    voronoi_cv.cv_inputs = [cv_input1, cv_input2]
    
    input_anchor1 = common_cv.Voronoi_cv_toy_anchor()
    input_anchor1.values = [0.5, 0.5]
    input_anchor1.starting_positions = np.array([[[0.5, 0.5, 0.0]]])
    input_anchor1.bound_state = False
    input_anchor1.bulk_state = False
    
    input_anchor2 = common_cv.Voronoi_cv_toy_anchor()
    input_anchor2.values = [-0.2, 0.2]
    input_anchor2.starting_positions = np.array([[[-0.2, 0.2, 0.0]]])
    input_anchor2.bound_state = False
    input_anchor2.bulk_state = False
    
    input_anchor3 = common_cv.Voronoi_cv_toy_anchor()
    input_anchor3.values = [0.0, 0.0]
    input_anchor3.starting_positions = np.array([[[0.0, 0.0, 0.0]]])
    input_anchor3.bound_state = False
    input_anchor3.bulk_state = False
    
    input_anchor4 = common_cv.Voronoi_cv_toy_anchor()
    input_anchor4.values = [0.2, -0.7]
    input_anchor4.starting_positions = np.array([[[0.2, -0.7, 0.0]]])
    input_anchor4.bound_state = False
    input_anchor4.bulk_state = False
    
    voronoi_cv.input_anchors = [input_anchor1, input_anchor2, input_anchor3, 
                                input_anchor4]
    
    stateA = common_cv.State_point()
    stateA.name = "stateA"
    stateA.location = [0.5, 0.5]
    stateB = common_cv.State_point()
    stateB.name = "stateB"
    stateB.location = [0.2, -0.7]
    voronoi_cv.state_points = [stateA, stateB]
    
    model_input.cv_inputs = [voronoi_cv]
    model_input.browndye_settings_input = None
    model_input.toy_settings_input = common_prepare.Toy_settings_input()
    model_input.toy_settings_input.potential_energy_expression \
        = "5*(x1^6+y1^6+exp(-(10*y1)^2)*(1-exp(-(10*x1)^2)))"
    model_input.toy_settings_input.num_particles = 1
    model_input.toy_settings_input.masses = np.array([10.0])
    
    return model_input

def create_ala_ala_mmvt_model_input(root_dir, ff="charmm"):
    """
    Create a generic host-guest model input object.
    """
    os.chdir(TEST_DIRECTORY)
    model_input = common_prepare.Model_input()
    model_input.calculation_type = "mmvt"
    model_input.calculation_settings = common_prepare.MMVT_input_settings()
    model_input.calculation_settings.md_output_interval = 100
    model_input.calculation_settings.md_steps_per_anchor = 1000 #1000000
    model_input.temperature = 298.15
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
    cv_input1.group1 = [6]
    cv_input1.group2 = [26]
    cv_input1.input_anchors = []
    
    radius_list = [0.7, 0.8, 0.9]
    charmm_psf_filename = os.path.abspath("data/ala_ala_ala.psf")
    charmm_param_filenames = [os.path.abspath("data/charmm22.par"),
                              os.path.abspath("data/charmm22.rtf")]
    pdb_filenames = [os.path.abspath("data/ala_ala_ala.pdb"), "", ""]
    for i, (radius, pdb_filename) in enumerate(zip(radius_list, pdb_filenames)):
        input_anchor = common_cv.Spherical_cv_anchor()
        input_anchor.radius = radius
        if ff == "charmm":
            assign_charmm_params(input_anchor, charmm_psf_filename, 
                                charmm_param_filenames, pdb_filename)
        else:
            raise Exception("ff type not supported: {}".format(ff))
        
        input_anchor.bound_state = False
        input_anchor.bulk_anchor = False
        cv_input1.input_anchors.append(input_anchor)
    
    model_input.cv_inputs = [cv_input1] 
    model_input.browndye_settings_input = None
    return model_input