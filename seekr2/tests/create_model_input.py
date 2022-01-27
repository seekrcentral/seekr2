"""
create_model_input.py

Create generic Model input objects
"""

import os

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
            assign_forcefield_params(input_anchor, forcefield_built_in_ff_list, 
                                     forcefield_custom_ff_list, pdb_filename)
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
        model_input.browndye_settings_input.num_bd_milestone_trajectories = 100
        model_input.browndye_settings_input.num_b_surface_trajectories = 10000
        model_input.browndye_settings_input.max_b_surface_trajs_to_extract = 100
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

def create_tryp_ben_mmvt_model_input(root_dir, bd=False):
    """
    Create a trypsin-benzamidine model input object.
    """
    os.chdir(TEST_DIRECTORY)
    model_input = common_prepare.Model_input()
    model_input.calculation_type = "mmvt"
    model_input.calculation_settings = common_prepare.MMVT_input_settings()
    model_input.calculation_settings.md_output_interval = 2500000
    model_input.calculation_settings.md_steps_per_anchor = 250000000
    model_input.temperature = 298.0
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
    cv_input1.group1 = [2478, 2489, 2499, 2535, 2718, 2745, 2769, 2787, 2794, 
                        2867, 2926]
    cv_input1.group2 = [3221, 3222, 3223, 3224, 3225, 3226, 3227, 3228, 3229]
    cv_input1.input_anchors = []
    
    input_anchor1 = common_cv.Spherical_cv_anchor()
    input_anchor1.radius = 0.05
    input_anchor1.lower_milestone_radius = None
    input_anchor1.upper_milestone_radius = 0.1
    input_anchor1.starting_amber_params = base.Amber_params()
    input_anchor1.starting_amber_params.prmtop_filename \
        = os.path.abspath("../data/trypsin_benzamidine_files/tryp_ben.prmtop")
    input_anchor1.starting_amber_params.pdb_coordinates_filename = \
        os.path.abspath(
            "../data/trypsin_benzamidine_files/mmvt/tryp_ben_at0.pdb")
    input_anchor1.bound_state = True
    input_anchor1.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor1)
    
    input_anchor2 = common_cv.Spherical_cv_anchor()
    input_anchor2.radius = 0.15
    input_anchor2.lower_milestone_radius = 0.1
    input_anchor2.upper_milestone_radius = 0.2
    input_anchor2.starting_amber_params = base.Amber_params()
    input_anchor2.starting_amber_params.prmtop_filename \
        = os.path.abspath("../data/trypsin_benzamidine_files/tryp_ben.prmtop")
    input_anchor2.starting_amber_params.pdb_coordinates_filename = \
        os.path.abspath(
            "../data/trypsin_benzamidine_files/mmvt/tryp_ben_at1.58.pdb")
    input_anchor2.bound_state = False
    input_anchor2.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor2)
    
    input_anchor3 = common_cv.Spherical_cv_anchor()
    input_anchor3.radius = 0.25
    input_anchor3.lower_milestone_radius = 0.2
    input_anchor3.upper_milestone_radius = 0.3
    input_anchor3.starting_amber_params = base.Amber_params()
    input_anchor3.starting_amber_params.prmtop_filename \
        = os.path.abspath("../data/trypsin_benzamidine_files/tryp_ben.prmtop")
    input_anchor3.starting_amber_params.pdb_coordinates_filename = \
        os.path.abspath(
            "../data/trypsin_benzamidine_files/mmvt/tryp_ben_at2.76.pdb")
    input_anchor3.bound_state = False
    input_anchor3.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor3)
    
    input_anchor4 = common_cv.Spherical_cv_anchor()
    input_anchor4.radius = 0.35
    input_anchor4.lower_milestone_radius = 0.3
    input_anchor4.upper_milestone_radius = 0.4
    input_anchor4.starting_amber_params = base.Amber_params()
    input_anchor4.starting_amber_params.prmtop_filename \
        = os.path.abspath("../data/trypsin_benzamidine_files/tryp_ben.prmtop")
    input_anchor4.starting_amber_params.pdb_coordinates_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/mmvt/tryp_ben_at3.69.pdb")
    input_anchor4.bound_state = False
    input_anchor4.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor4)
    
    input_anchor5 = common_cv.Spherical_cv_anchor()
    input_anchor5.radius = 0.45
    input_anchor5.lower_milestone_radius = 0.4
    input_anchor5.upper_milestone_radius = 0.6
    input_anchor5.starting_amber_params = base.Amber_params()
    input_anchor5.starting_amber_params.prmtop_filename \
        = os.path.abspath("../data/trypsin_benzamidine_files/tryp_ben.prmtop")
    input_anchor5.starting_amber_params.pdb_coordinates_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/mmvt/tryp_ben_at4.61.pdb")
    input_anchor5.bound_state = False
    input_anchor5.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor5)
    
    input_anchor6 = common_cv.Spherical_cv_anchor()
    input_anchor6.radius = 0.7
    input_anchor6.lower_milestone_radius = 0.6
    input_anchor6.upper_milestone_radius = 0.8
    input_anchor6.starting_amber_params = base.Amber_params()
    input_anchor6.starting_amber_params.prmtop_filename \
        = os.path.abspath("../data/trypsin_benzamidine_files/tryp_ben.prmtop")
    input_anchor6.starting_amber_params.pdb_coordinates_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/mmvt/tryp_ben_at7.20.pdb")
    input_anchor6.bound_state = False
    input_anchor6.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor6)
    
    input_anchor7 = common_cv.Spherical_cv_anchor()
    input_anchor7.radius = 0.9
    input_anchor7.lower_milestone_radius = 0.8
    input_anchor7.upper_milestone_radius = 1.0
    input_anchor7.starting_amber_params = base.Amber_params()
    input_anchor7.starting_amber_params.prmtop_filename \
        = os.path.abspath("../data/trypsin_benzamidine_files/tryp_ben.prmtop")
    input_anchor7.starting_amber_params.pdb_coordinates_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/mmvt/tryp_ben_at9.33.pdb")
    input_anchor7.bound_state = False
    input_anchor7.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor7)
    
    input_anchor8 = common_cv.Spherical_cv_anchor()
    input_anchor8.radius = 1.1
    input_anchor8.lower_milestone_radius = 1.0
    input_anchor8.upper_milestone_radius = 1.2
    input_anchor8.starting_amber_params = base.Amber_params()
    input_anchor8.starting_amber_params.prmtop_filename \
        = os.path.abspath("../data/trypsin_benzamidine_files/tryp_ben.prmtop")
    input_anchor8.starting_amber_params.pdb_coordinates_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/mmvt/tryp_ben_at11.06.pdb")
    input_anchor8.bound_state = False
    input_anchor8.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor8)
    
    input_anchor9 = common_cv.Spherical_cv_anchor()
    input_anchor9.radius = 1.3
    input_anchor9.lower_milestone_radius = 1.2
    input_anchor9.upper_milestone_radius = 1.4
    input_anchor9.starting_amber_params = base.Amber_params()
    input_anchor9.starting_amber_params.prmtop_filename \
        = os.path.abspath("../data/trypsin_benzamidine_files/tryp_ben.prmtop")
    input_anchor9.starting_amber_params.pdb_coordinates_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/mmvt/tryp_ben_at12.70.pdb")
    input_anchor9.bound_state = False
    input_anchor9.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor9)
    
    input_anchor10 = common_cv.Spherical_cv_anchor()
    input_anchor10.radius = 1.5
    input_anchor10.lower_milestone_radius = 1.4
    input_anchor10.upper_milestone_radius = 1.6
    input_anchor10.starting_amber_params = base.Amber_params()
    input_anchor10.starting_amber_params.prmtop_filename \
        = os.path.abspath("../data/trypsin_benzamidine_files/tryp_ben.prmtop")
    input_anchor10.starting_amber_params.pdb_coordinates_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/mmvt/tryp_ben_at14.04.pdb")
    input_anchor10.bound_state = False
    input_anchor10.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor10)
    
    input_anchor11 = common_cv.Spherical_cv_anchor()
    input_anchor11.radius = 1.7
    input_anchor11.lower_milestone_radius = 1.6
    input_anchor11.upper_milestone_radius = 1.8
    input_anchor11.starting_amber_params = base.Amber_params()
    input_anchor11.starting_amber_params.prmtop_filename \
        = os.path.abspath("../data/trypsin_benzamidine_files/tryp_ben.prmtop")
    input_anchor11.starting_amber_params.pdb_coordinates_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/mmvt/tryp_ben_at17.00.pdb")
    input_anchor11.bound_state = False
    input_anchor11.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor11)
    
    input_anchor12 = common_cv.Spherical_cv_anchor()
    input_anchor12.radius = 1.9
    input_anchor12.lower_milestone_radius = 1.8
    input_anchor12.upper_milestone_radius = None
    input_anchor12.starting_amber_params = base.Amber_params()
    input_anchor12.starting_amber_params.prmtop_filename = ""
    input_anchor12.starting_amber_params.pdb_coordinates_filename = ""
    input_anchor12.bound_state = False
    input_anchor12.bulk_anchor = True
    cv_input1.input_anchors.append(input_anchor12)
    
    model_input.cv_inputs = [cv_input1]
    
    if bd:
        model_input.browndye_settings_input \
            = common_prepare.Browndye_settings_input()
        model_input.browndye_settings_input.binary_directory = ""
        model_input.browndye_settings_input.receptor_pqr_filename \
            = os.path.abspath("../data/trypsin_benzamidine_files/trypsin.pqr")
        model_input.browndye_settings_input.ligand_pqr_filename \
            = os.path.abspath(
                "../data/trypsin_benzamidine_files/benzamidine.pqr")
        model_input.browndye_settings_input.apbs_grid_spacing = 0.5
        model_input.browndye_settings_input.receptor_indices = [
            2478, 2489, 2499, 2535, 2718, 2745, 2769, 2787, 2794, 2867, 2926]
        model_input.browndye_settings_input.ligand_indices = [
            0, 1, 2, 3, 4, 5, 6, 7, 8]
        
        ion1 = base.Ion()
        ion1.radius = 1.14
        ion1.charge = 2.0
        ion1.conc = 0.02
        ion2 = base.Ion()
        ion2.radius = 1.67
        ion2.charge = -1.0
        ion2.conc = 0.1
        ion3 = base.Ion()
        ion3.radius = 4.0
        ion3.charge = 1.0
        ion3.conc = 0.06
        model_input.browndye_settings_input.ions = [ion1, ion2, ion3]
        model_input.browndye_settings_input.num_bd_milestone_trajectories = 1000
        model_input.browndye_settings_input.num_b_surface_trajectories = 10000
        model_input.browndye_settings_input.max_b_surface_trajs_to_extract = 100
        model_input.browndye_settings_input.n_threads = 1
    else:
        model_input.browndye_settings_input  = None
    
    return model_input

def create_tryp_ben_elber_model_input(root_dir, bd=False):
    """
    Create a trypsin-benzamidine Elber model input object.
    """
    os.chdir(TEST_DIRECTORY)
    model_input = common_prepare.Model_input()
    model_input.calculation_type = "elber"
    model_input.calculation_settings = common_prepare.Elber_input_settings()
    model_input.calculation_settings.temperature_equil_progression = [
        300., 310., 320., 330., 340., 350., 340., 330., 320., 310., 300.]
    model_input.calculation_settings.num_temperature_equil_steps = 1000
    model_input.calculation_settings.num_umbrella_stage_steps = 500000000
    model_input.calculation_settings.umbrella_force_constant = 5000.0
    model_input.calculation_settings.fwd_rev_interval = 100000
    model_input.temperature = 298.0
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
    cv_input1.group1 = [2478, 2489, 2499, 2535, 2718, 2745, 2769, 2787, 2794, 
                        2867, 2926]
    cv_input1.group2 = [3221, 3222, 3223, 3224, 3225, 3226, 3227, 3228, 3229]
    cv_input1.input_anchors = []
    
    input_anchor1 = common_cv.Spherical_cv_anchor()
    input_anchor1.radius = 0.1
    input_anchor1.starting_amber_params = base.Amber_params()
    input_anchor1.starting_amber_params.prmtop_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/elber/tryp_ben.prmtop")
    input_anchor1.starting_amber_params.pdb_coordinates_filename = \
        os.path.abspath(
            "../data/trypsin_benzamidine_files/elber/tryp_ben_at1.pdb")
    input_anchor1.bound_state = True
    input_anchor1.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor1)
    
    input_anchor2 = common_cv.Spherical_cv_anchor()
    input_anchor2.radius = 0.15
    input_anchor2.starting_amber_params = base.Amber_params()
    input_anchor2.starting_amber_params.prmtop_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/elber/tryp_ben.prmtop")
    input_anchor2.starting_amber_params.pdb_coordinates_filename = \
        os.path.abspath(
            "../data/trypsin_benzamidine_files/elber/tryp_ben_at1.5.pdb")
    input_anchor2.bound_state = False
    input_anchor2.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor2)
    
    input_anchor3 = common_cv.Spherical_cv_anchor()
    input_anchor3.radius = 0.2
    input_anchor3.starting_amber_params = base.Amber_params()
    input_anchor3.starting_amber_params.prmtop_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/elber/tryp_ben.prmtop")
    input_anchor3.starting_amber_params.pdb_coordinates_filename = \
        os.path.abspath(
            "../data/trypsin_benzamidine_files/elber/tryp_ben_at2.pdb")
    input_anchor3.bound_state = False
    input_anchor3.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor3)
    
    input_anchor4 = common_cv.Spherical_cv_anchor()
    input_anchor4.radius = 0.3
    input_anchor4.starting_amber_params = base.Amber_params()
    input_anchor4.starting_amber_params.prmtop_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/elber/tryp_ben.prmtop")
    input_anchor4.starting_amber_params.pdb_coordinates_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/elber/tryp_ben_at3.pdb")
    input_anchor4.bound_state = False
    input_anchor4.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor4)
    
    input_anchor5 = common_cv.Spherical_cv_anchor()
    input_anchor5.radius = 0.4
    input_anchor5.starting_amber_params = base.Amber_params()
    input_anchor5.starting_amber_params.prmtop_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/elber/tryp_ben.prmtop")
    input_anchor5.starting_amber_params.pdb_coordinates_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/elber/tryp_ben_at4.pdb")
    input_anchor5.bound_state = False
    input_anchor5.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor5)
    
    input_anchor6 = common_cv.Spherical_cv_anchor()
    input_anchor6.radius = 0.6
    input_anchor6.starting_amber_params = base.Amber_params()
    input_anchor6.starting_amber_params.prmtop_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/elber/tryp_ben.prmtop")
    input_anchor6.starting_amber_params.pdb_coordinates_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/elber/tryp_ben_at6.pdb")
    input_anchor6.bound_state = False
    input_anchor6.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor6)
    
    input_anchor7 = common_cv.Spherical_cv_anchor()
    input_anchor7.radius = 0.8
    input_anchor7.starting_amber_params = base.Amber_params()
    input_anchor7.starting_amber_params.prmtop_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/elber/tryp_ben.prmtop")
    input_anchor7.starting_amber_params.pdb_coordinates_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/elber/tryp_ben_at8.pdb")
    input_anchor7.bound_state = False
    input_anchor7.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor7)
    
    input_anchor8 = common_cv.Spherical_cv_anchor()
    input_anchor8.radius = 1.0
    input_anchor8.starting_amber_params = base.Amber_params()
    input_anchor8.starting_amber_params.prmtop_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/elber/tryp_ben.prmtop")
    input_anchor8.starting_amber_params.pdb_coordinates_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/elber/tryp_ben_at10.pdb")
    input_anchor8.bound_state = False
    input_anchor8.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor8)
    
    input_anchor9 = common_cv.Spherical_cv_anchor()
    input_anchor9.radius = 1.2
    input_anchor9.starting_amber_params = base.Amber_params()
    input_anchor9.starting_amber_params.prmtop_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/elber/tryp_ben.prmtop")
    input_anchor9.starting_amber_params.pdb_coordinates_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/elber/tryp_ben_at12.pdb")
    input_anchor9.bound_state = False
    input_anchor9.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor9)
    
    input_anchor10 = common_cv.Spherical_cv_anchor()
    input_anchor10.radius = 1.4
    input_anchor10.starting_amber_params = base.Amber_params()
    input_anchor10.starting_amber_params.prmtop_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/elber/tryp_ben.prmtop")
    input_anchor10.starting_amber_params.pdb_coordinates_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/elber/tryp_ben_at14.pdb")
    input_anchor10.bound_state = False
    input_anchor10.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor10)
    
    input_anchor11 = common_cv.Spherical_cv_anchor()
    input_anchor11.radius = 1.6
    input_anchor11.starting_amber_params = base.Amber_params()
    input_anchor11.starting_amber_params.prmtop_filename \
        = os.path.abspath(
            "../data/trypsin_benzamidine_files/elber/tryp_ben.prmtop")
    input_anchor11.starting_amber_params.pdb_coordinates_filename \
        = ""
    input_anchor11.bound_state = False
    input_anchor11.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor11)
    
    input_anchor12 = common_cv.Spherical_cv_anchor()
    input_anchor12.radius = 1.8
    input_anchor12.starting_amber_params = base.Amber_params()
    input_anchor12.starting_amber_params.prmtop_filename = ""
    input_anchor12.starting_amber_params.pdb_coordinates_filename = ""
    input_anchor12.bound_state = False
    input_anchor12.bulk_anchor = True
    cv_input1.input_anchors.append(input_anchor12)
    
    model_input.cv_inputs = [cv_input1]
    
    if bd:
        model_input.browndye_settings_input \
            = common_prepare.Browndye_settings_input()
        model_input.browndye_settings_input.binary_directory = ""
        model_input.browndye_settings_input.receptor_pqr_filename \
            = os.path.abspath("../data/trypsin_benzamidine_files/trypsin.pqr")
        model_input.browndye_settings_input.ligand_pqr_filename \
            = os.path.abspath(
                "../data/trypsin_benzamidine_files/benzamidine.pqr")
        model_input.browndye_settings_input.apbs_grid_spacing = 0.5
        model_input.browndye_settings_input.receptor_indices = [
            2478, 2489, 2499, 2535, 2718, 2745, 2769, 2787, 2794, 2867, 2926]
        model_input.browndye_settings_input.ligand_indices = [
            0, 1, 2, 3, 4, 5, 6, 7, 8]
        
        ion1 = base.Ion()
        ion1.radius = 1.14
        ion1.charge = 2.0
        ion1.conc = 0.02
        ion2 = base.Ion()
        ion2.radius = 1.67
        ion2.charge = -1.0
        ion2.conc = 0.1
        ion3 = base.Ion()
        ion3.radius = 4.0
        ion3.charge = 1.0
        ion3.conc = 0.06
        model_input.browndye_settings_input.ions = [ion1, ion2, ion3]
        model_input.browndye_settings_input.num_bd_milestone_trajectories = 1000
        model_input.browndye_settings_input.num_b_surface_trajectories = 10000
        model_input.browndye_settings_input.max_b_surface_trajs_to_extract = 100
        model_input.browndye_settings_input.n_threads = 1
    else:
        model_input.browndye_settings_input  = None
    
    return model_input

def create_smoluchowski_mmvt_model_input(root_dir, num_input_anchors=5):
    """
    Create a Smoluchowski MMVT model input object.
    """
    os.chdir(TEST_DIRECTORY)
    model_input = common_prepare.Model_input()
    model_input.calculation_type = "mmvt"
    model_input.calculation_settings = common_prepare.MMVT_input_settings()
    model_input.calculation_settings.md_output_interval = 0
    model_input.calculation_settings.md_steps_per_anchor = 0
    model_input.temperature = 300.0
    model_input.pressure = 1.0
    model_input.ensemble = "nvt"
    model_input.root_directory = root_dir
    model_input.md_program = "toy"
    model_input.constraints = "none"
    model_input.rigidWater = False
    model_input.hydrogenMass = None
    model_input.timestep = 0.002
    model_input.nonbonded_cutoff = 0.0
    cv_input1 = common_cv.Spherical_cv_input()
    cv_input1.group1 = [0]
    cv_input1.group2 = [1]
    cv_input1.input_anchors = []
    
    for i in range(num_input_anchors):
        input_anchor = common_cv.Spherical_cv_anchor()
        input_anchor.radius = i + 0.5
        input_anchor.starting_amber_params = None
        if i == 0:
            input_anchor.bound_state = True
        else:
            input_anchor.bound_state = False
        
        if i == num_input_anchors-1:
            input_anchor.bulk_anchor = True
        else:
            input_anchor.bulk_anchor = False
        cv_input1.input_anchors.append(input_anchor)
    
    """ # TODO: remove
    input_anchor1 = common_cv.Spherical_cv_anchor()
    input_anchor1.radius = 0.5
    input_anchor1.starting_amber_params = None
    input_anchor1.bound_state = True
    input_anchor1.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor1)
    
    input_anchor2 = common_cv.Spherical_cv_anchor()
    input_anchor2.radius = 1.5
    input_anchor2.starting_amber_params = None
    input_anchor2.bound_state = False
    input_anchor2.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor2)
    
    input_anchor3 = common_cv.Spherical_cv_anchor()
    input_anchor3.radius = 2.5
    input_anchor3.starting_amber_params = None
    input_anchor3.bound_state = False
    input_anchor3.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor3)
    
    input_anchor4 = common_cv.Spherical_cv_anchor()
    input_anchor4.radius = 3.5
    input_anchor4.starting_amber_params = None
    input_anchor4.bound_state = False
    input_anchor4.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor4)
    
    input_anchor5 = common_cv.Spherical_cv_anchor()
    input_anchor5.radius = 4.5
    input_anchor5.starting_amber_params = None
    input_anchor5.bound_state = False
    input_anchor5.bulk_anchor = True
    cv_input1.input_anchors.append(input_anchor5)
    """
    
    model_input.cv_inputs = [cv_input1]
    
    model_input.browndye_settings_input  = None
    
    return model_input

def create_smoluchowski_elber_model_input(root_dir, num_input_anchors=5):
    """
    Create a Smoluchowski Elber model input object.
    """
    os.chdir(TEST_DIRECTORY)
    model_input = create_smoluchowski_mmvt_model_input(root_dir)
    model_input.calculation_type = "elber"
    model_input.calculation_settings = common_prepare.Elber_input_settings()
    model_input.calculation_settings.num_umbrella_stage_steps = 100000
    model_input.calculation_settings.umbrella_force_constant = 9000.0
    model_input.calculation_settings.fwd_rev_interval = 500
    model_input.calculation_settings.rev_output_interval = None
    model_input.calculation_settings.fwd_output_interval = None
    
    cv_input1 = common_cv.Spherical_cv_input()
    cv_input1.group1 = [0]
    cv_input1.group2 = [1]
    cv_input1.input_anchors = []
    cv_input1.input_anchors = []
    
    for i in range(num_input_anchors):
        input_anchor = common_cv.Spherical_cv_anchor()
        input_anchor.radius = i + 1.0
        input_anchor.starting_amber_params = None
        if i == 0:
            input_anchor.bound_state = True
        else:
            input_anchor.bound_state = False
        
        if i == num_input_anchors-1:
            input_anchor.bulk_anchor = True
        else:
            input_anchor.bulk_anchor = False
        cv_input1.input_anchors.append(input_anchor)
    
    """
    input_anchor1 = common_cv.Spherical_cv_anchor()
    input_anchor1.radius = 1.0
    input_anchor1.starting_amber_params = None
    input_anchor1.bound_state = True
    input_anchor1.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor1)
    
    input_anchor2 = common_cv.Spherical_cv_anchor()
    input_anchor2.radius = 2.0
    input_anchor2.starting_amber_params = None
    input_anchor2.bound_state = False
    input_anchor2.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor2)
    
    input_anchor3 = common_cv.Spherical_cv_anchor()
    input_anchor3.radius = 3.0
    input_anchor3.starting_amber_params = None
    input_anchor3.bound_state = False
    input_anchor3.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor3)
    
    input_anchor4 = common_cv.Spherical_cv_anchor()
    input_anchor4.radius = 4.0
    input_anchor4.starting_amber_params = None
    input_anchor4.bound_state = False
    input_anchor4.bulk_anchor = False
    cv_input1.input_anchors.append(input_anchor4)
    
    input_anchor5 = common_cv.Spherical_cv_anchor()
    input_anchor5.radius = 5.0
    input_anchor5.starting_amber_params = None
    input_anchor5.bound_state = False
    input_anchor5.bulk_anchor = True
    cv_input1.input_anchors.append(input_anchor5)
    """
    
    model_input.cv_inputs = [cv_input1]
    model_input.browndye_settings_input  = None
    return model_input

def create_sod_mmvt_model_input(root_dir):
    """
    Create a generic host-guest model input object.
    """
    os.chdir(TEST_DIRECTORY)
    input_xml_filename = os.path.abspath(
        "../data/superoxide_dismutase_files/input_sod.xml")
    model_input = common_prepare.Model_input()
    model_input.deserialize(input_xml_filename, user_input=True)
    model_input.root_directory = root_dir
    os.chdir(os.path.dirname(input_xml_filename))
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
    order_parameter1.group1 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 
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
    order_parameter1.group2 = [147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 
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
    
    #model_input.browndye_settings_input = None
    
    return model_input

def create_rmsd_mmvt_model_input(root_dir):
    """
    Create a RMSD TRP cage model input.
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
    model_input.nonbonded_cutoff = None
    
    cv_input1 = common_cv.RMSD_cv_input()
    cv_input1.group = [4, 16, 26, 47, 57, 74, 98, 117, 139, 151, 158, 173, 
                       179, 190, 201, 208, 240, 254, 268, 274]
    cv_input1.ref_structure = "/home/lvotapka/seekr2/seekr2/data/trp_cage_files/trp_cage.pdb"
    cv_input1.input_anchors = []
    
    values_list = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
    amber_prmtop_filename = "/home/lvotapka/seekr2/seekr2/data/trp_cage_files/trp_cage.prmtop"
    forcefield_built_in_ff_list = ["amber14/tip3pfb.xml"]
    forcefield_custom_ff_list = ["../data/hostguest_files/hostguest.xml"]
    
    pdb_filenames = ["/home/lvotapka/seekr2/seekr2/data/trp_cage_files/trp_cage.pdb", 
                     "/home/lvotapka/seekr2/seekr2/data/trp_cage_files/trp_cage_0.40.pdb", 
                     "/home/lvotapka/seekr2/seekr2/data/trp_cage_files/trp_cage_0.60.pdb", 
                     "/home/lvotapka/seekr2/seekr2/data/trp_cage_files/trp_cage_0.80.pdb", 
                     "/home/lvotapka/seekr2/seekr2/data/trp_cage_files/trp_cage_1.00.pdb", 
                     ""]
    for i, (value, pdb_filename) in enumerate(zip(values_list, pdb_filenames)):
        input_anchor = common_cv.Tiwary_cv_anchor()
        input_anchor.value = value
        assign_amber_params(input_anchor, amber_prmtop_filename, 
                            pdb_filename)
        
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
    
    return model_input