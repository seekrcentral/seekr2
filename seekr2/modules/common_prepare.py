"""
common_prepare.py

Common routines for preparing the model file for either Elber or MMVT
milestoning calculations.

The resulting model file will be used for all other stages of the 
calculation.
"""

import os
import re
import glob
import shutil
from collections import defaultdict

import numpy as np

import seekr2.modules.common_base as base
import seekr2.modules.mmvt_cvs.mmvt_cv_base as mmvt_cv_base
import seekr2.modules.mmvt_cvs.mmvt_spherical_cv as mmvt_spherical_cv
import seekr2.modules.mmvt_cvs.mmvt_tiwary_cv as mmvt_tiwary_cv
import seekr2.modules.mmvt_cvs.mmvt_planar_cv as mmvt_planar_cv
import seekr2.modules.mmvt_cvs.mmvt_rmsd_cv as mmvt_rmsd_cv
import seekr2.modules.mmvt_cvs.mmvt_closest_pair_cv as mmvt_closest_pair_cv
import seekr2.modules.mmvt_cvs.mmvt_count_contacts_cv as mmvt_count_contacts_cv
import seekr2.modules.mmvt_cvs.mmvt_external_cv as mmvt_external_cv
import seekr2.modules.mmvt_cvs.mmvt_z_distance_cv as mmvt_z_distance_cv
import seekr2.modules.mmvt_cvs.mmvt_voronoi_cv as mmvt_voronoi_cv
import seekr2.modules.elber_cvs.elber_cv_base as elber_cv_base
import seekr2.modules.elber_cvs.elber_spherical_cv as elber_spherical_cv
import seekr2.modules.elber_cvs.elber_external_cv as elber_external_cv
import seekr2.modules.common_cv as common_cv
import seekr2.modules.filetree as filetree
import seekr2.modules.common_sim_browndye2 as sim_browndye2
import seekr2.modules.runner_browndye2 as runner_browndye2
from abserdes import Serializer

KEEP_ANCHOR_RE = "hidr*|string*"

def anchor_has_files(model, anchor):
    """
    Returns True if seekr simulations have already been run in this anchor.
    """
    
    output_directory = os.path.join(
        model.anchor_rootdir, anchor.directory, anchor.production_directory)
    output_files_glob = os.path.join(output_directory, anchor.md_output_glob)
    output_restarts_list = glob.glob(output_files_glob)
    if len(output_restarts_list) > 0:
        return True
    else:
        return False
    
def anchor_has_hidr_files(model, anchor):
    """
    Returns True if HIDR has already been run in this anchor.
    """
    
    building_directory = os.path.join(
        model.anchor_rootdir, anchor.directory, anchor.building_directory)
    input_files_glob = os.path.join(building_directory, "hidr*")
    input_files_list = glob.glob(input_files_glob)
    if len(input_files_list) > 0:
        return True
    else:
        return False

def cleanse_anchor_outputs(model, anchor):
    """
    Delete all simulation outputs for an existing anchor to make way
    for new outputs.
    """
    
    if model.openmm_settings is not None:
        import seekr2.modules.runner_openmm as runner_openmm
        runner_openmm.cleanse_anchor_outputs(model, anchor)
    elif model.namd_settings is not None:
        import seekr2.modules.runner_namd as runner_namd
        runner_namd.cleanse_anchor_outputs(model, anchor)
    return

class Browndye_settings_input(Serializer):
    """
    Read and parse the inputs for the BrownDye2 program, which runs
    the BD stage of the SEEKR2 calculation
    
    Attributes:
    -----------
    binary_directory : str
        A path to the BrownDye2 programs directory. Can be left as an
        empty string if the $PATH environmental variable points to
        BrownDye's bin/ directory.
        
    receptor_pqr_filename : str
        The filename path pointing to the pqr file which will be used
        as the receptor in the BrownDye simulations.
        
    ligand_pqr_filename : str
        The filename path pointing to the pqr file which will be used
        as the ligand in the BrownDye simulations.
        
    apbs_grid_spacing : float
        The resolution (in Angstroms) of the APBS (electrostatics) 
        grid.
        
    ions : list
        A list of common_base.Ion() objects which will be passed to
        APBS.
    
    num_b_surface_trajectories : int
        The number of trajectories to run with the ligand starting at
        the b-surface.
        
    num_bd_milestone_trajectories : int
        The number of trajectories to run with the ligand starting at
        each of the BD milestones.
    
    max_b_surface_trajs_to_extract : int
        After the b-surface simulation, members of the encounter 
        complexes will be extracted to construct the FHPD. Then, these
        will be run as their own independent BD simulations. This
        parameter defines the maximum number of structures to extract
        from the FHPD and run simulations for.
    
    receptor_indices : list
        The indices of the atoms (numbered starting from zero) within 
        the receptor whose center of mass will be calculated and used
        to represent the binding site.
        
    ligand_indices : list
        The indices of the atoms (numbered starting from zero) within 
        the ligand whose center of mass will be calculated and used
        to represent the binding site.
        
    n_threads : int, Default 1
        The number of cores to use for the BrownDye calculation.
    """
    
    def __init__(self):
        self.binary_directory = ""
        self.receptor_pqr_filename = ""
        self.ligand_pqr_filename = ""
        self.apbs_grid_spacing = -1.0
        self.ions = []
        self.num_b_surface_trajectories = -1
        self.receptor_indices = []
        self.ligand_indices = []
        self.n_threads = 1

class Toy_settings_input(Serializer):
    """
    Read and parse the inputs for a toy simulation.
    
    Attributes:
    -----------
    potential_energy_expression : str
        The expression that can be used in OpenMM to direct particle 
        motions.
        
    num_particles : int
        The number of particles in the system.
    
    masses : list
        A list of particle masses
    """
    def __init__(self):
        self.potential_energy_expression = None
        self.num_particles = -1
        self.masses = []
        return

class MMVT_input_settings(Serializer):
    """
    Settings needed for running an MMVT simulation.
    
    Attributes:
    -----------
    md_output_interval : int, default 500000
        How many steps must elapse between outputs. This will effect
        the interval that trajectory frames, state outputs, and
        backup checkpoints will occur.
        
    md_steps_per_anchor : int, default 500000000
        The default number of timesteps that will be simulated at each
        anchor.
    """
    
    def __init__(self):
        self.md_output_interval = 500000
        self.md_steps_per_anchor = 500000000
        
class Elber_input_settings(Serializer):
    """
    Settings needed for running an Elber milestoning simulation.
    
    Attributes:
    -----------
    temperature_equil_progression : list
        A list of floats representing the progression of temperature
        changes (in Kelvin) to run for a temperature equilibration
        stage. If an empty list, then the temperature equilibration
        stage will be skipped.
    
    num_temperature_equil_steps : int
        The number of timesteps to run for each step of the temperature
        progression of the temperature equilibration stage.
        
    num_umbrella_stage_steps : int
        The number of steps to take per umbrella stage calculation.
        
    umbrella_force_constant : float
        The force constant (in kilojoules per mole) to use for the 
        umbrella stage.
        
    fwd_rev_interval : int
        The interval of steps in the umbrella stage at which to launch
        reversal and forward stages.
        
    rev_output_interval : int
        The interval at which to output reversal trajectory and 
        state information.
        
    fwd_output_interval : int
        The interval at which to output forward trajectory and 
        state information.
    """
    
    def __init__(self):
        #self.temperature_equil_progression = [
        #    300., 310., 320., 330., 340., 350., 340., 330., 320., 310., 300]
        self.temperature_equil_progression = []
        self.num_temperature_equil_steps = 1000
        self.num_umbrella_stage_steps = 50000
        self.umbrella_force_constant = 9000.0
        self.fwd_rev_interval = 500
        self.rev_output_interval = 500
        self.fwd_output_interval = 500

class Model_input(Serializer):
    """
    The serializable object representing parameters that would be
    input by the user.
    
    Attributes:
    -----------
    calculation_type : str
        A string to represent the calculation type.
    
    calculation_settings : MMVT_settings() or Elber_settings()
        Settings general to any MD engine, including the number of
        timesteps to run, backup intervals, etc.
    
    temperature : float, default 298.15
        The temperature (in Kelvin) at which this calculation takes 
        place. This quantity will affect the initial atomic velocities,
        the thermostat, barostat, and any quantities computed in the
        analysis stage.
        
    pressure : float, default 1.0
        If a constant pressure situation is specified, this is the 
        pressure (in bar) at which the simulation will be maintained.
        
    ensemble : str, default "nvt"
        Which ensemble to use for the MMVT simulations. Options include
        "npt" (constant T and P), "nvt" (Constant T and V), or "nve"
        (constant E and V), "npt_membrane" (Constant T and P with 
        a membrane barostat).
        
    root_directory : str
        A path to the location where the model XML file will be written, 
        along with the anchor directories and other files.
    
    md_program : str, default "openmm"
        Which MD engine to use for the MD portion of the SEEKR
        calculation. Options include "openmm" and "namd".
        
    run_minimization : bool, default False
        Whether to run minimizations of the initial structures before 
        beginning the MMVT simulations. This is not recommended, since 
        systems might deviate outside of the Voronoi cell during 
        minimizations. Instead, users should minimize or equilibrate,
        then check whether the system still remains within the Voronoi
        cell in the initial structures.
    
    hydrogenMass : float or None, Default None
        This parameter may be used for hydrogen mass repartitioning.
        If a float is provided, then mass will be shifted from every
        hydrogen atom's bonded heavy atom to the hydrogen, such that
        the hydrogen atoms will have that value (in AMU) for their
        masses. If None, then the default hydrogen mass is used.
        
    constraints : str, Default "hbonds"
        The constraints method to use for bonds and angles within the
        system. This argument is supplied to the constraints argument
        of the OpenMM System() object.
        
    rigidWater : bool, Default True
        If True, then water bonds and angles will be made rigid.
    
    integrator_type : str, Default "langevin"
        Which type of integrator to use for simulation dynamics. 
        Available options include "langevin", "langevinMiddle", ...
    
    timestep : float, Default 0.002
        The length of time taken by a simulation step (in units of 
        picoseconds).
        
    nonbonded_cutoff : float, Default 0.9
        The distance after which VDW nonbonded forces are cut off in
        units of nanometers. This argument is supplied to the 
        nonbondedCutoff argument of an OpenMM System() object.
    
    browndye_settings_input : Browndye_settings_input or None
        The Browndye_settings_input() object for this model. It 
        contains all the settings that could be used within a Browndye
        simulation.
        
    toy_settings_input : Toy_settings_input or None
        The Toy_settings_input() object for this model. It 
        contains all the settings that could be used within a toy
        simulation.
    
    cv_inputs : list
        A list containing a set of collective variable input objects.
    """
    
    def __init__(self):
        self.calculation_type = "Unknown"
        self.calculation_settings = None
        self.temperature = 298.15
        self.pressure = 1.0
        self.ensemble = "nvt"
        self.root_directory = ""
        self.md_program = "openmm"
        self.run_minimization = False
        self.hydrogenMass = None
        self.constraints = "hbonds"
        self.rigidWater = True
        self.waterModel = None
        self.integrator_type = "langevin"
        self.timestep = 0.002
        self.nonbonded_cutoff = 0.9
        self.browndye_settings_input = None
        self.toy_settings_input = None
        self.cv_inputs = []
        
    def get_cv_inputs(self):
        """
        The model_input object contains a list of CV_inputs as well as
        Combos. Flatten the CV_inputs within the model_input as well as the
        Combos, so that the model has a complete, flattened list of CVs.
        Also check for strange inputs.
        """
        flat_cv_inputs = []
        for cv_input in self.cv_inputs:
            if isinstance(cv_input, common_cv.Combo):
                for cv_input2 in cv_input.cv_inputs:
                    assert isinstance(cv_input2, common_cv.CV_input), \
                        "Combos can only contain CV_inputs in the cv_input "\
                        "field."
                    assert len(cv_input2.state_points) == 0, \
                        "The CV_inputs within a combo may not have "\
                        "state_points."
                    flat_cv_inputs.append(cv_input2)
            elif isinstance(cv_input, common_cv.CV_input):
                flat_cv_inputs.append(cv_input)
            else:
                raise Exception("Improper object detected in model_input "\
                                "cv_input field.")
                
        return flat_cv_inputs
        
def model_factory(model_input, use_absolute_directory=False):
    """
    Given the Model_input object, which contains the settings, 
    create the Model object which contains more detailed 
    information.
    """
    
    model = base.Model()
    calc_settings = model_input.calculation_settings
    if model_input.calculation_type.lower() == "mmvt":
        model.calculation_settings = mmvt_cv_base.MMVT_settings()
        model.calculation_settings.num_production_steps = \
            calc_settings.md_steps_per_anchor
        model.calculation_settings.energy_reporter_interval = \
            calc_settings.md_output_interval
        model.calculation_settings.restart_checkpoint_interval = \
            calc_settings.md_output_interval
        model.calculation_settings.trajectory_reporter_interval = \
            calc_settings.md_output_interval
    elif model_input.calculation_type.lower() == "elber":
        model.calculation_settings = elber_cv_base.Elber_settings()
        model.calculation_settings.temperature_equil_progression = \
            calc_settings.temperature_equil_progression
        model.calculation_settings.num_temperature_equil_steps = \
            calc_settings.num_temperature_equil_steps
        model.calculation_settings.num_umbrella_stage_steps = \
            calc_settings.num_umbrella_stage_steps
        model.calculation_settings.umbrella_force_constant = \
            calc_settings.umbrella_force_constant
        model.calculation_settings.fwd_rev_interval = \
            calc_settings.fwd_rev_interval
        model.calculation_settings.umbrella_energy_reporter_interval = \
            calc_settings.fwd_rev_interval
        model.calculation_settings.umbrella_trajectory_reporter_interval = \
            calc_settings.fwd_rev_interval
        model.calculation_settings.rev_energy_reporter_interval = \
            calc_settings.rev_output_interval
        model.calculation_settings.rev_trajectory_reporter_interval = \
            calc_settings.rev_output_interval
        model.calculation_settings.fwd_energy_reporter_interval = \
            calc_settings.fwd_output_interval
        model.calculation_settings.fwd_trajectory_reporter_interval = \
            calc_settings.fwd_output_interval
            
    else:
        raise Exception("Invalid calculation_type entered:", 
                        model_input.calculation_type)
    model.calculation_type = model_input.calculation_type.lower()
    temperature = model_input.temperature
    model.temperature = temperature
    if use_absolute_directory:
        model.anchor_rootdir = model_input.root_directory
    else:
        model.anchor_rootdir = "."
    if model_input.ensemble.lower() == "nvt":
        pressure = None
    elif model_input.ensemble.lower() in ["npt", "npt_membrane"]:
        pressure = model_input.pressure
    elif model_input.ensemble.lower() == "nve":
        raise Exception("NVE ensemble not yet implemented. Options must be"\
                        "NVT or NPT.")
    else:
        raise Exception("Invalid ensemble entered:", model_input.ensemble)
    
    if model_input.md_program.lower() == "openmm":
        mm_settings = base.Openmm_settings()
        mm_settings.langevin_integrator.target_temperature = temperature
        if pressure is not None:
            mm_settings.barostat = base.Barostat_settings_openmm()
            mm_settings.barostat.target_temperature = temperature
            mm_settings.barostat.target_pressure = pressure
            if model_input.ensemble.lower() == "npt_membrane":
                mm_settings.barostat.membrane = True
        mm_settings.initial_temperature = temperature
        if model_input.nonbonded_cutoff is None:
            mm_settings.nonbonded_method = "nocutoff"
        else:
            mm_settings.nonbonded_cutoff = model_input.nonbonded_cutoff
        mm_settings.run_minimization = model_input.run_minimization
        mm_settings.hydrogenMass = model_input.hydrogenMass
        mm_settings.constraints = model_input.constraints
        mm_settings.langevin_integrator.timestep = model_input.timestep
        mm_settings.langevin_integrator.integrator_type = model_input.integrator_type
        mm_settings.rigidWater = model_input.rigidWater
    
        model.openmm_settings = mm_settings
    elif model_input.md_program.lower() == "namd":
        namd_settings = base.Namd_settings()
        namd_settings.langevin_integrator.target_temperature = temperature
        if pressure is not None:
            namd_settings.barostat = base.Barostat_settings_namd()
            namd_settings.barostat.target_temperature = temperature
            namd_settings.barostat.target_pressure = pressure
        namd_settings.initial_temperature = temperature
        namd_settings.nonbonded_cutoff = model_input.nonbonded_cutoff
        namd_settings.run_minimization = model_input.run_minimization
        assert model_input.hydrogenMass is None, "hydrogen mass "\
            "repartitioning not yet implemented in NAMD SEEKR."
        #namd_settings.hydrogenMass = model_input.hydrogenMass
        namd_settings.constraints = model_input.constraints
        namd_settings.langevin_integrator.timestep = model_input.timestep
        namd_settings.rigidWater = model_input.rigidWater
        namd_settings.watermodel = model_input.waterModel
        model.namd_settings = namd_settings
    
    elif model_input.md_program.lower() == "smoluchowski":
        pass
        
    else:
        raise Exception("Invalid MD program entered:", 
                        model_input.md_program)
    
    if model_input.browndye_settings_input is not None:
        k_on_info = base.K_on_info()
        if model_input.browndye_settings_input.ions is None:
            model_input.browndye_settings_input.ions = []
        k_on_info.ions = model_input.browndye_settings_input.ions
        k_on_info.b_surface_num_trajectories \
            = model_input.browndye_settings_input.num_b_surface_trajectories
        
        model.browndye_settings = base.Browndye_settings()
        model.browndye_settings.browndye_bin_dir \
            = model_input.browndye_settings_input.binary_directory
        model.browndye_settings.receptor_pqr_filename \
            = os.path.basename(
                model_input.browndye_settings_input.receptor_pqr_filename)
        model.browndye_settings.ligand_pqr_filename \
            = os.path.basename(
                model_input.browndye_settings_input.ligand_pqr_filename)
        model.browndye_settings.apbs_grid_spacing \
            = model_input.browndye_settings_input.apbs_grid_spacing
        model.browndye_settings.n_threads \
            = model_input.browndye_settings_input.n_threads
        model.k_on_info = k_on_info
    
    if model_input.toy_settings_input is not None:
        assert model_input.md_program.lower() == "openmm", \
            "Only OpenMM can run toy systems."
        toy_settings = base.Toy_settings()
        toy_settings.potential_energy_expression \
            = model_input.toy_settings_input.potential_energy_expression
        toy_settings.num_particles \
            = model_input.toy_settings_input.num_particles
        toy_settings.masses = model_input.toy_settings_input.masses
        model.toy_settings = toy_settings
        mm_settings.cuda_platform_settings = None
        mm_settings.reference_platform = True
        
    return model

def resolve_connections(connection_flag_dict, model, associated_input_anchor,
                        root_directory, force_overwrite):
    """
    Given connection information between anchors, merge anchors
    that are connected.
    """
    anchors = model.anchors
    anchor_indices_to_remove = []
    index_reference = {}
    bulk_anchor = None
    for anchor in anchors:
        index_reference[anchor.index] = anchor.index
    
    for key in connection_flag_dict:
        anchor_connection_list = connection_flag_dict[key]
        bulk_discarded_indices = []
        if key == "bulk":
            bulk_anchor = anchor_connection_list[0]
            for anchor_discarded in anchor_connection_list:
                anchor_indices_to_remove.append(anchor_discarded.index)
                bulk_discarded_indices.append(anchor_discarded.index)
                
        if not key == "bulk":
            assert len(anchor_connection_list) > 1, \
                "Connection flag {} needs to have endpoints ".format(key) \
                +"in at least two anchors."
        anchor_kept = anchor_connection_list[0]
        for anchor_discarded in anchor_connection_list[1:]:
            input_anchor_kept = associated_input_anchor[anchor_kept.index]
            input_anchor_discarded = associated_input_anchor[anchor_discarded.index]
            
            
            assert anchor_kept.endstate == anchor_discarded.endstate, \
                "The anchors connected by connection flag {} ".format(key) \
                +"must have the same value for 'bound_state'."
            assert anchor_kept.bulkstate == anchor_discarded.bulkstate, \
                "The anchors connected by connection flag {} ".format(key) \
                +"must have the same value for 'bulk_anchor'."   
            # Ensure that both anchors don't have starting structures 
            # defined
            if anchor.__class__.__name__ \
                    in ["MMVT_toy_anchor", "Elber_toy_anchor"]:
                if input_anchor_kept.starting_positions is None:
                    input_anchor_kept.starting_positions \
                        = input_anchor_discarded.starting_positions
                
                if input_anchor_discarded.starting_positions is None:
                    input_anchor_discarded.starting_positions \
                        = input_anchor_kept.starting_positions
                print("input_anchor_kept.starting_positions:", input_anchor_kept.starting_positions)
                print("input_anchor_discarded.starting_positions:", input_anchor_discarded.starting_positions)
                assert (input_anchor_kept.starting_positions \
                      == input_anchor_discarded.starting_positions).all(), \
                    "The anchors connected by connection flag {} ".format(key) \
                    +"have inconsistent starting positions defined."
            else:
                # TODO: this will need to be handled for charmm and other inputs
                if input_anchor_kept.starting_amber_params is None \
                        or input_anchor_kept.starting_amber_params\
                        .pdb_coordinates_filename == "":
                    input_anchor_kept.starting_amber_params \
                        = input_anchor_discarded.starting_amber_params
                
                if input_anchor_discarded.starting_amber_params is None\
                        or input_anchor_discarded.starting_amber_params\
                        .pdb_coordinates_filename == "":
                    input_anchor_discarded.starting_amber_params \
                        = input_anchor_kept.starting_amber_params
                assert base.same_amber_params(
                    input_anchor_kept.starting_amber_params, 
                    input_anchor_discarded.starting_amber_params), \
                    "The anchors connected by connection flag {} ".format(key) \
                    +"have inconsistent starting positions defined."
            anchor_indices_to_remove.append(anchor_discarded.index)
            anchor_kept.milestones += anchor_discarded.milestones
            anchor_kept.variables.update(anchor_discarded.variables)
            index_reference[anchor_discarded.index] = anchor_kept.index
            
    # sort and remove all redundant anchors
    anchor_indices_to_remove = sorted(set(anchor_indices_to_remove), 
                                      reverse=True)
    for anchor_index_to_remove in anchor_indices_to_remove:
        if len(anchors) > anchor_index_to_remove+1:
            for lower_anchor in anchors[anchor_index_to_remove+1:]:
                index_reference[lower_anchor.index] -= 1
                
    for anchor_index_to_remove in anchor_indices_to_remove:
        anchors.pop(anchor_index_to_remove)
    
    if bulk_anchor is not None:
        # make sure the bulk anchor is the last anchor always
        anchors.append(bulk_anchor)
        index_reference[bulk_anchor.index] = len(anchors) - 1
        for bulk_discarded in bulk_discarded_indices:
            index_reference[bulk_discarded] = len(anchors) - 1
    
    visited_new_indices = []
    
    for old_index in index_reference:
        new_index = index_reference[old_index]
        if new_index < len(anchors):
            anchors[new_index].index = new_index
            anchors[new_index].name = "anchor_"+str(anchors[new_index].index)
            anchors[new_index].directory = anchors[new_index].name
            alias_index = 1
            milestones_to_remove = []
            for i, milestone in enumerate(anchors[new_index].milestones):
                if (model.get_type() == "mmvt") \
                        and (new_index not in visited_new_indices):
                    # don't renumber the same anchor twice
                    neighbor_id = milestone.neighbor_anchor_index
                    milestone.neighbor_anchor_index = index_reference[neighbor_id]
                    milestone.alias_index = alias_index
                    if milestone.neighbor_anchor_index \
                            == anchors[new_index].index:
                        milestones_to_remove.append(i)
                    alias_index += 1
                    
            for milestone_to_remove in sorted(milestones_to_remove, reverse=True):
                anchors[new_index].milestones.pop(milestone_to_remove)
            
        visited_new_indices.append(new_index)
            
    xml_path = os.path.join(root_directory, "model.xml")
    if os.path.exists(xml_path):
        # then a model file already exists at this location: update
        # the anchor directories.
        old_model = base.Model()
        old_model.deserialize(xml_path)
        new_anchors_with_starting_pdbs_to_keep = modify_model(
            old_model, model, root_directory, force_overwrite)
    else:
        new_anchors_with_starting_pdbs_to_keep = None
    
    for old_index in index_reference:
        new_index = index_reference[old_index]
        if new_index < len(anchors):
            input_anchor = associated_input_anchor[old_index]
            filetree.generate_filetree_by_anchor(
                anchors[new_index], root_directory)
            if (new_anchors_with_starting_pdbs_to_keep is not None) \
                    and (new_index in new_anchors_with_starting_pdbs_to_keep):
                print("Keeping starting structure for anchor:", new_index)
            else:
                filetree.copy_building_files_by_anchor(
                    anchors[new_index], input_anchor, root_directory)
                                
    return anchors

def create_cvs(model, collective_variable_inputs, root_directory):
    """
    Create the collective variable objects for the model.
    
    Parameters:
    -----------
    model : Model()
        The model object representing the entire calculation.
    
    collective_variable_inputs : list
        A list of cv input objects used to construct the anchors for
        this model.
        
    root_directory : str
        A path to the root directory of this model. Used to copy over
        RMSD CV reference structure.
    
    Returns:
    --------
    cvs : list
        A list of collective variable objects that will be placed into
        the model.
    
    """
    cvs = []
    for i, cv_input in enumerate(collective_variable_inputs):
        cv_input.index = i
        cv_input.check()
        if model.get_type() == "mmvt":
            if isinstance(cv_input, common_cv.Spherical_cv_input):
                cv = mmvt_spherical_cv.make_mmvt_spherical_cv_object(cv_input, index=i)
            elif isinstance(cv_input, common_cv.Tiwary_cv_input):
                cv = mmvt_tiwary_cv.make_mmvt_tiwary_cv_object(cv_input, index=i)
            elif isinstance(cv_input, common_cv.Planar_cv_input):
                cv = mmvt_planar_cv.make_mmvt_planar_cv_object(cv_input, index=i)
            elif isinstance(cv_input, common_cv.RMSD_cv_input):
                cv = mmvt_rmsd_cv.make_mmvt_RMSD_cv_object(
                    cv_input, index=i, root_directory=root_directory)
            elif isinstance(cv_input, common_cv.Closest_pair_cv_input):
                cv = mmvt_closest_pair_cv.make_mmvt_closest_pair_cv_object(
                    cv_input, index=i)
            elif isinstance(cv_input, common_cv.Count_contacts_cv_input):
                cv = mmvt_count_contacts_cv.make_mmvt_count_contacts_cv_object(
                    cv_input, index=i)
            elif isinstance(cv_input, common_cv.Toy_cv_input):
                cv = mmvt_external_cv.make_mmvt_external_cv_object(cv_input, index=i)
            elif isinstance(cv_input, common_cv.Z_distance_cv_input):
                cv = mmvt_z_distance_cv.make_mmvt_z_distance_cv_object(cv_input, index=i)
            elif isinstance(cv_input, common_cv.Voronoi_cv_input):
                cv = mmvt_voronoi_cv.make_mmvt_voronoi_cv_object(
                    cv_input, index=i, root_directory=root_directory)
            else:
                raise Exception("CV type not implemented in MMVT: %s" \
                                % type(cv_input))
            
        elif model.get_type() == "elber":
            if isinstance(cv_input, common_cv.Spherical_cv_input):
                cv = elber_spherical_cv.make_elber_spherical_cv_object(cv_input, index=i)
            elif isinstance(cv_input, common_cv.Toy_cv_input):
                cv = elber_external_cv.make_elber_external_cv_object(cv_input, index=i)
            else:
                raise Exception("CV type not implemented in Elber: %s" \
                                % type(cv_input))
        
        cvs.append(cv)
    return cvs

def create_anchors(model, model_input):
    """
    Create the Anchor objects for the Model from the input CVs, 
    input anchors, and Combos.
    
    Parameters:
    -----------
    model : Model()
        The model object representing the entire calculation.
    
    """
    milestone_index = 0
    anchor_index = 0
    anchors = []
    connection_flag_dict = defaultdict(list)
    associated_input_anchor = {}
    
    for cv_input in model_input.cv_inputs:
        if isinstance(cv_input, common_cv.Combo):
            this_anchors, anchor_index, milestone_index, connection_flag_dict,\
                associated_input_anchor = cv_input.make_anchors(
                model, anchor_index, milestone_index, connection_flag_dict,
                associated_input_anchor)
            anchors += this_anchors
        elif isinstance(cv_input, common_cv.Voronoi_cv_input):
            this_anchors, anchor_index, milestone_index, connection_flag_dict,\
                associated_input_anchor = cv_input.make_anchors(
                model, anchor_index, milestone_index, connection_flag_dict,
                associated_input_anchor)
            anchors += this_anchors
        else: 
            # It's a cv_input
            this_cvs_anchors = []
            for j, input_anchor in enumerate(cv_input.input_anchors):
                input_anchor.check(j, cv_input)
                anchor = common_cv.create_anchor(model, anchor_index)
                if not input_anchor.bulk_anchor:
                    anchor.md = True
                else:
                    anchor.md = False
                    anchor.bulkstate = True
                    anchor.name = "bulk"
                    
                if input_anchor.bound_state:
                    anchor.name = "bound"
                    anchor.endstate = True
                # Make the milestone objects assuming that this is a 1D CV
                variable_name = "{}_{}".format(cv_input.variable_name, 
                                               cv_input.index)
                variable_value = input_anchor.get_variable_value()
                anchor.variables[variable_name] = variable_value
                for connection_flag in input_anchor.connection_flags:
                    connection_flag_dict[connection_flag].append(anchor)
                if input_anchor.bulk_anchor:
                    connection_flag_dict["bulk"].append(anchor)
                associated_input_anchor[anchor.index] = input_anchor
                anchors.append(anchor)
                anchor_index += 1
                this_cvs_anchors.append(anchor)
                
            input_anchor_index = 0
            for i, anchor in enumerate(this_cvs_anchors):
                if model.get_type() == "mmvt":
                    this_input_anchor = cv_input.input_anchors[i]
                    if input_anchor_index < len(cv_input.input_anchors)-1:
                        upper_anchor = this_cvs_anchors[i+1]
                        upper_input_anchor = cv_input.input_anchors[i+1]
                        milestone_index \
                            = cv_input.make_mmvt_milestone_between_two_anchors(
                                anchor, upper_anchor, this_input_anchor, 
                                upper_input_anchor, milestone_index)
                    if input_anchor_index > 0:
                        lower_anchor = this_cvs_anchors[i-1]
                        lower_input_anchor = cv_input.input_anchors[i-1]
                        milestone_index \
                            = cv_input.make_mmvt_milestone_between_two_anchors(
                                lower_anchor, anchor, lower_input_anchor, 
                                this_input_anchor, milestone_index)
                        
                elif model.get_type() == "elber":
                    milestone_alias = 1
                    if isinstance(cv_input, common_cv.Spherical_cv_input):
                        milestones, milestone_alias, milestone_index = \
                            elber_spherical_cv.make_elber_milestoning_objects_spherical(
                            cv_input, milestone_alias, milestone_index, 
                            input_anchor_index, anchor.index, cv_input.input_anchors, 
                            model.calculation_settings.umbrella_force_constant)
                    elif isinstance(cv_input, common_cv.Toy_cv_input):
                        milestones, milestone_alias, milestone_index = \
                            elber_external_cv.make_elber_milestoning_objects_external(
                            cv_input, milestone_alias, milestone_index, 
                            input_anchor_index, anchor.index, cv_input.input_anchors, 
                            model.calculation_settings.umbrella_force_constant)
                            
                    anchor.milestones += milestones
                input_anchor_index += 1
                
    if model.get_type() == "elber":
        milestone_index += 1
    
    return anchors, milestone_index, connection_flag_dict, \
        associated_input_anchor

def create_bd_milestones(model, model_input):
    """
    Using all bulk state anchors, create the BD milestones.
    """
    
    model.k_on_info.bd_milestones = []
    bd_milestone_counter = 0
    for anchor in model.anchors:
        if anchor.bulkstate:
            bd_milestone = base.BD_milestone()
            bd_milestone.index = bd_milestone_counter
            bd_milestone.name = \
                "bd_milestone_%d" % bd_milestone.index
            bd_milestone.directory = bd_milestone.name
            
            if model.get_type() == "mmvt":
                bd_milestone.outer_milestone = anchor.milestones[0]
                assert "radius" in bd_milestone.outer_milestone.variables,\
                    "A BD outer milestone must be spherical."
                neighbor_anchor = model.anchors[
                    bd_milestone.outer_milestone.neighbor_anchor_index]
                for neighbor_milestone in neighbor_anchor.milestones:
                    if neighbor_milestone.index != \
                            bd_milestone.outer_milestone.index:
                        if neighbor_milestone.cv_index == \
                                bd_milestone.outer_milestone.cv_index:
                            bd_milestone.inner_milestone = \
                                neighbor_milestone
                                
                assert bd_milestone.inner_milestone is not None, "No suitable "\
                    "spherical milestone found for inner BD reaction "\
                    "criteria. Make sure that there are at least 3 spherical "\
                    "input anchors for the outermost CV."
            
            elif model.get_type() == "elber":
                bd_milestone.outer_milestone = anchor.milestones[1]
                assert "radius" in bd_milestone.outer_milestone.variables,\
                    "A BD outer milestone must be spherical."
                bd_milestone.inner_milestone = anchor.milestones[0]
            
            cv_index = bd_milestone.outer_milestone.cv_index
            cv_input = model_input.cv_inputs[cv_index]
            if len(cv_input.bd_group1)>0:
                bd_milestone.receptor_indices \
                    = base.parse_xml_list(cv_input.bd_group1)
            else:
                bd_milestone.receptor_indices  \
                    = base.parse_xml_list(
                        model_input.browndye_settings_input.receptor_indices)
            
            if len(cv_input.bd_group2)>0:
                bd_milestone.ligand_indices \
                    = base.parse_xml_list(cv_input.bd_group2)
            else:
                bd_milestone.ligand_indices \
                    = base.parse_xml_list(
                        model_input.browndye_settings_input.ligand_indices)
            
            model.k_on_info.bd_milestones.append(bd_milestone)
            bd_milestone_counter += 1
            
    assert len(model.k_on_info.bd_milestones) > 0, "No BD milestones created. "\
        "Make sure there is at least one bulk state input anchor."
                
    return

def prepare_model_cvs_and_anchors(model, model_input, force_overwrite):
    """
    Fill out the CV and anchor items within the model based on the 
    input objects.
    """
    flat_cv_inputs = model_input.get_cv_inputs()
    cvs = create_cvs(model, flat_cv_inputs, model_input.root_directory)
    model.collective_variables = cvs
    anchors, num_milestones, connection_flag_dict, associated_input_anchor\
        = create_anchors(model, model_input)
    model.anchors = anchors
    
    if model.get_type() == "mmvt":
        if model_input.md_program.lower() == "openmm":
            for anchor in model.anchors:
                anchor.md_output_glob = mmvt_cv_base.OPENMMVT_GLOB
        elif model_input.md_program.lower() == "namd":
            for anchor in model.anchors:
                anchor.md_output_glob = mmvt_cv_base.NAMDMMVT_GLOB
        model.num_milestones = num_milestones
        
    elif model.get_type() == "elber":
        if model_input.md_program.lower() == "openmm":
            for anchor in model.anchors:
                anchor.md_output_glob = elber_cv_base.OPENMM_ELBER_GLOB
        elif model_input.md_program.lower() == "namd":
            for anchor in model.anchors:
                anchor.md_output_glob = elber_cv_base.NAMD_ELBER_GLOB
        model.num_milestones = num_milestones-1
        
    if model_input.browndye_settings_input is not None:
        create_bd_milestones(model, model_input)
    
    anchors = resolve_connections(connection_flag_dict, model, 
                                  associated_input_anchor, 
                                  model_input.root_directory, force_overwrite)
    
    # check to make sure that anchors don't have repeated milestone aliases
    # or have themselves as neighbors. Also check that bulk anchors are at
    # the end of the list
    iterated_through_bulk_anchors = False
    for anchor in anchors:
        alias_indices = set()
        for milestone in anchor.milestones:
            assert milestone.alias_index not in alias_indices, \
                "Repeated alias indices detected in anchor {}".format(
                    anchor.index)
            alias_indices.add(milestone.alias_index)
            if model.get_type() == "mmvt":
                assert milestone.neighbor_anchor_index is not anchor.index, \
                    "Milestone neighbor_anchor_index cannot be its own anchor "\
                    "in MMVT. Index: {}".format(anchor.index)
        
        if anchor.bulkstate:
            iterated_through_bulk_anchors = True
        else:
            assert not iterated_through_bulk_anchors, \
                f"Non-bulk anchor {anchor.index} found after previous bulk "\
                "anchor. Make sure all bulk input anchors are listed after "\
                "non-bulk input anchors."
            
    model.num_anchors = len(anchors)
    common_cv.assign_state_points(model_input, model)
    return
    
def generate_bd_files(model, rootdir):
    """
    Create the ghost atoms in the BD files, convert PQRs to XMLs,
    make the input.xml file, as well as the reaction criteria.
    """
    
    if model.using_bd():
        b_surface_dir = os.path.join(
            rootdir, model.k_on_info.b_surface_directory)
        receptor_pqr_filename = os.path.join(
            b_surface_dir, model.browndye_settings.receptor_pqr_filename)
        ligand_pqr_filename = os.path.join(
            b_surface_dir, model.browndye_settings.ligand_pqr_filename)
        ghost_indices_rec = []
        ghost_indices_lig = []
        for bd_milestone in model.k_on_info.bd_milestones:
            #print("adding ghost atom to file:", receptor_pqr_filename)
            ghost_index_rec = \
                sim_browndye2.add_ghost_atom_to_pqr_from_atoms_center_of_mass(
                    receptor_pqr_filename, bd_milestone.receptor_indices)
            #print("adding ghost atom to file:", ligand_pqr_filename)
            ghost_index_lig = \
                sim_browndye2.add_ghost_atom_to_pqr_from_atoms_center_of_mass(
                    ligand_pqr_filename, bd_milestone.ligand_indices)
            ghost_indices_rec.append(ghost_index_rec)
            ghost_indices_lig.append(ghost_index_lig)
            
        model.browndye_settings.ghost_indices_rec = ghost_indices_rec
        model.browndye_settings.ghost_indices_lig = ghost_indices_lig
        receptor_xml_filename = sim_browndye2.make_pqrxml(
            receptor_pqr_filename, 
            browndye2_bin=model.browndye_settings.browndye_bin_dir)
        ligand_xml_filename = sim_browndye2.make_pqrxml(ligand_pqr_filename)
        debye_length, reaction_filename = \
            runner_browndye2.make_browndye_input_xml(
            model, rootdir, receptor_xml_filename, ligand_xml_filename,
            model.k_on_info.b_surface_num_trajectories)
        model.browndye_settings.debye_length = debye_length
        abs_reaction_path = os.path.join(b_surface_dir, reaction_filename)
        runner_browndye2.make_browndye_reaction_xml(model, abs_reaction_path)
        return
    
def modify_model(old_model, new_model, root_directory, force_overwrite=False):
    """
    If someone runs the prepare stage on a model with existing
    directories and simulation results, examine the changes to delete
    and rename directories properly, and only overwrite simulation files
    if forced to do so.
    """
    
    old_model.anchor_rootdir = root_directory
    anchor_pairs = []
    old_anchors_to_delete = []
    new_anchors_to_create = []
    # Look through the new and old milestones and pair them up by the value
    #  of their variables
    for alpha, anchor1 in enumerate(new_model.anchors):
        if anchor1.bulkstate:
            continue
        alpha_paired = False
        for beta, anchor2 in enumerate(old_model.anchors):
            if anchor2.bulkstate:
                continue
            if anchor1.variables == anchor2.variables:
                pair = (alpha, beta)
                anchor_pairs.append(pair)
                alpha_paired = True
                break
        if not alpha_paired:
            new_anchors_to_create.append(alpha)
        
    for beta, anchor2 in enumerate(old_model.anchors):
        if anchor2.bulkstate:
            continue
        beta_paired = False
        for pair in anchor_pairs:
            if beta == pair[1]:
                beta_paired = True
                break
        if not beta_paired:
            old_anchors_to_delete.append(beta)
        
    # Now check all the paired anchors to see if anyone's milestones
    # have changed
    old_anchors_with_changed_milestones = []
    old_anchors_with_starting_pdbs_to_keep = []
    for pair in anchor_pairs:
        (alpha, beta) = pair
        anchor1 = new_model.anchors[alpha]
        anchor2 = old_model.anchors[beta]
        milestones_changed = False
        milestone_pairs = []
        for i, milestone1 in enumerate(anchor1.milestones):
            milestone1_paired = False
            for j, milestone2 in enumerate(anchor2.milestones):
                if milestone1.variables == milestone2.variables:
                    milestone_pairs.append((i,j))
                    milestone1_paired = True
                    break
            if not milestone1_paired:
                milestones_changed = True
        for j, milestone2 in enumerate(anchor2.milestones):
            milestone2_paired = False
            for milestone_pair in milestone_pairs:
                if j == milestone_pair[1]:
                    milestone2_paired = True
                    break
            if not milestone2_paired:
                milestones_changed = True
        if milestones_changed:
            old_anchors_with_changed_milestones.append(beta)
            
        # See whether the starting structures have changed
        if not new_model.using_toy():
            old_anchor_starting = base.get_anchor_pdb_filename(anchor2)
            if re.match(KEEP_ANCHOR_RE, old_anchor_starting):
                old_anchors_with_starting_pdbs_to_keep.append(beta)
    
    warning_msg = ""
    # check whether any anchors have been given new starting structures
    if len(old_anchors_with_starting_pdbs_to_keep) > 0:
        warning_msg+="The anchor(s) {} ".format(
                old_anchors_with_starting_pdbs_to_keep)\
            +"have been given new starting structures (from HIDR or elsewhere) "\
            "and the entered command mismatches the new structure file names. "\
            "If you desire to ignore this message, then use the "\
            "--force_overwrite (-f) option, and the existing structures from HIDR "\
            "will be kept in the new model. If you wish to overwrite these, you "\
            "will need to change the model.xml file or write the root directory "\
            "to a new location.\n\n"
    
    # check whether deleting directories contain simulation files
    deleting_directories_with_files = []
    for del_anchor_index in old_anchors_to_delete:
        del_anchor = old_model.anchors[del_anchor_index]
        if anchor_has_files(old_model, del_anchor):
            deleting_directories_with_files.append(del_anchor.directory)
    
    if len(deleting_directories_with_files) > 0 :
        warning_msg+="The anchor(s) {} ".format(
                deleting_directories_with_files)\
            +"already have existing output files "\
            "and the entered command would overwrite them. If you desire "\
            "to overwrite the existing files, then use the "\
            "--force_overwrite (-f) option, and these anchors will be "\
            "deleted in directory: " + root_directory + "\n\n"
    
    # Check whether directories with changed milestones already contain
    # simulation files
    cleansing_anchors = []
    cleansing_directories = []
    for cleanse_anchor_index in old_anchors_with_changed_milestones:
        cleanse_anchor = old_model.anchors[cleanse_anchor_index]
        if anchor_has_files(old_model, cleanse_anchor):
            cleansing_anchors.append(cleanse_anchor)
            cleansing_directories.append(cleanse_anchor.directory)
    
    if len(cleansing_directories) > 0 :
        warning_msg+="The anchor(s) {} ".format(
                cleansing_directories)\
            +"already have existing output files yet the milestone locations "\
            "will be changed if this command proceeds. If you desire "\
            "to overwrite the existing files, then use the "\
            "--force_overwrite (-f) option, and these outputs will be "\
            "deleted in directory: " + root_directory + "\n\n"
    
    # Check whether deleting directories may contain HIDR starting structures
    deleting_directories_with_hidr = []
    for del_anchor_index in old_anchors_to_delete:
        del_anchor = old_model.anchors[del_anchor_index]
        if anchor_has_hidr_files(old_model, del_anchor):
            deleting_directories_with_hidr.append(del_anchor.directory)
            
    if len(deleting_directories_with_hidr) > 0 :
        warning_msg+="The anchor(s) {} ".format(
                deleting_directories_with_hidr)\
            +"already have existing HIDR files "\
            "and the entered command would delete them. If you desire "\
            "to override this warning, then use the "\
            "--force_overwrite (-f) option, and these anchors will be "\
            "deleted in directory: " + root_directory + "\n\n"
    
    if warning_msg and not force_overwrite:
        print(warning_msg)
        raise Exception("Cannot overwrite existing files.")
    
    for del_anchor_index in old_anchors_to_delete:
        del_anchor = old_model.anchors[del_anchor_index]
        full_path = os.path.join(root_directory, del_anchor.directory)
        print("removing directory: "+del_anchor.directory)
        shutil.rmtree(full_path)
    
    for cleansing_anchor in cleansing_anchors:
        print("removing output files from anchor:", cleansing_anchor.name)
        cleanse_anchor_outputs(old_model, cleansing_anchor)
    
    # Now that the equivalent anchors are known, we know which ones will
    # be deleted, and which ones are being changed, with or without
    # simulation data, 
    
    new_anchors_with_starting_pdbs_to_keep = []
    for pair in anchor_pairs:
        (alpha, beta) = pair
        anchor1 = old_model.anchors[beta]
        anchor2 = new_model.anchors[alpha]
        full_path1 = os.path.join(root_directory, anchor1.directory)
        full_path2 = os.path.join(root_directory, anchor2.directory)
        if full_path1 != full_path2:
            if os.path.exists(full_path2):
                # then we cannot rename to this without giving this one
                # a temporary directory
                for problem_anchor in old_model.anchors:
                    if problem_anchor.bulkstate:
                        continue
                    if problem_anchor.directory == anchor2.directory:
                        temporary_name = problem_anchor.directory+"temp"
                        temp_full_path1 = os.path.join(root_directory, 
                                                  problem_anchor.directory)
                        temp_full_path2 = os.path.join(root_directory, 
                                                  temporary_name)
                        print("moving directory {} to {}".format(
                            problem_anchor.directory, temporary_name))
                        os.rename(temp_full_path1, temp_full_path2)
                        problem_anchor.directory = temporary_name
                        break
            
            print("moving directory {} to {}".format(anchor1.directory, 
                                                 anchor2.directory))
            # If it gets here, it's because a new directory was made that
            # has no data in it. It can be safely deleted
            if os.path.exists(full_path2):
                shutil.rmtree(full_path2)
                
            os.rename(full_path1, full_path2)
        
        if new_model.using_toy():
            anchor2.starting_positions = anchor1.starting_positions
        else:
            anchor2.amber_params = anchor1.amber_params
            anchor2.forcefield_params = anchor1.forcefield_params
            anchor2.charmm_params = anchor1.charmm_params
            
        if anchor1.index in old_anchors_with_starting_pdbs_to_keep:
            new_anchors_with_starting_pdbs_to_keep.append(anchor2.index)
    
    for new_anchor_index_to_create in new_anchors_to_create:
        print("creating new anchor {}".format(new_anchor_index_to_create))
    #    anchor = new_model.anchors[new_anchor_index_to_create]
    #    filetree.generate_filetree_by_anchor(anchor, root_directory)
    #    copy_building_files_by_anchor(anchor, input_anchor, rootdir)
    
    # TODO: check for BD milestone changes
    if old_model.using_bd():
        bd_milestones_to_check = []
        bd_milestones_to_cleanse = []
        b_surface_directory = os.path.join(
            old_model.anchor_rootdir, 
            old_model.k_on_info.b_surface_directory)
        if not new_model.using_bd():
            warning_msg = "The new model does not have any BD settings "\
            "provided, while the old model does. This would delete all "\
            "existing BD simulation files in {}. Use the --force_overwrite "\
            "(-f) option if you wish to proceed anyway.".format(
                old_model.anchor_rootdir)
            if not force_overwrite:
                print(warning_msg)
                raise Exception("Cannot overwrite existing BD outputs.")
            
            print("Deleting folder:", b_surface_directory)
            shutil.rmtree(b_surface_directory)
            force_overwrite = False
            return new_anchors_with_starting_pdbs_to_keep
        
        if len(old_model.k_on_info.bd_milestones) \
                != len(new_model.k_on_info.bd_milestones):
            warning_msg = "A different number of bd_milestones have been defined"\
                "in this model's modifications. You should probably generate an"\
                "entirely new model for this system. Use the --force_overwrite "\
                "(-f) option if you wish to proceed anyway."
            if not force_overwrite:
                print(warning_msg)
                raise Exception("Cannot define a new number of BD milestones.")
        
            for bd_milestone in old_model.k_on_info.bd_milestones:
                bd_milestones_to_check.append(bd_milestone)
        
        else:
            # check to see if bd milestone variables changed
            warnstr = "One or more BD milestones has had variables changed. "\
                "If BD simulations have already been run, they will be "\
                "overwritten. Use the --force_overwrite (-f) option if you "\
                "wish to proceed anyway."
            for old_bd_milestone, new_bd_milestone in zip(
                    old_model.k_on_info.bd_milestones, 
                    new_model.k_on_info.bd_milestones):
                
                if old_bd_milestone.outer_milestone.variables \
                        != new_bd_milestone.outer_milestone.variables:
                    bd_milestones_to_check.append(old_bd_milestone)
                    continue
                
                if old_bd_milestone.inner_milestone.variables \
                        != new_bd_milestone.outer_milestone.variables:
                    bd_milestones_to_check.append(old_bd_milestone)
                    continue
        
        if len(bd_milestones_to_check) > 0:
            b_surface_files_present = runner_browndye2.cleanse_bd_outputs(
                b_surface_directory, check_mode=True)
            
            if b_surface_files_present:
                if not force_overwrite:
                    print(warnstr)
                    raise Exception("Cannot overwrite existing BD milestones.")
            
            b_surface_directory = os.path.join(
                old_model.anchor_rootdir, 
                old_model.k_on_info.b_surface_directory)
            runner_browndye2.cleanse_bd_outputs(b_surface_directory, 
                                                check_mode=False)
        
    return new_anchors_with_starting_pdbs_to_keep
    
