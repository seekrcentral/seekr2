"""
common_prepare.py

Common routines for preparing the model file for either Elber or MMVT
milestoning calculations.

The resulting model file will be used for all other stages of the 
calculation.
"""

import os
import glob
import shutil
from collections import defaultdict

import seekr2.modules.common_base as base
import seekr2.modules.mmvt_base as mmvt_base
import seekr2.modules.elber_base as elber_base
import seekr2.modules.common_cv as common_cv
import seekr2.modules.mmvt_cv as mmvt_cv
import seekr2.modules.elber_cv as elber_cv
import seekr2.modules.filetree as filetree
import seekr2.modules.common_sim_browndye2 as sim_browndye2
import seekr2.modules.runner_browndye2 as runner_browndye2
from abserdes import Serializer

def anchor_has_files(model, anchor):
    """
    Returns True if simulations have already been run in this anchor.
    """
    
    output_directory = os.path.join(
        model.anchor_rootdir, anchor.directory, anchor.production_directory)
    output_files_glob = os.path.join(output_directory, anchor.md_output_glob)
    output_restarts_list = glob.glob(output_files_glob)
    if len(output_restarts_list) > 0:
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
    Read and parse the outputs from the BrownDye2 program, which runs
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

class Connector(Serializer):
    """
    A Connector defines a connection between two input anchors of different
    (or possibly the same) CVs. A connection between two input anchors will
    results in a single anchor in the final model.
    """
    
    def __init__(self):
        self.index = -1
        self.flag = -1
        self.starting_amber_params = None
        self.starting_forcefield_params = None
        self.bound_state = False
        self.bulk_anchor = False

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
        (constant E and V).
        
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
    
    timestep : float, Default 0.002
        The length of time taken by a simulation step (in units of 
        picoseconds).
        
    nonbonded_cutoff : float, Default 0.9
        The distance after which VDW nonbonded forces are cut off in
        units of nanometers. This argument is supplied to the 
        nonbondedCutoff argument of an OpenMM System() object.
    
    browndye_settings_input : Browndye_settings_input
        The Browndye_settings_input() object for this model. It 
        contains all the settings that could be used within a Browndye
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
        self.timestep = 0.002
        self.nonbonded_cutoff = 0.9
        self.browndye_settings_input = None
        self.cv_inputs = []
        self.connectors = []
        
    def read_plain_input_file(self, filename):
        """
        Read a plain input file (as opposed to an XML)
        """
        
        raise Exception("Reading a plain text file is not yet implemented. "\
                        "Only an XML model input may be read at this time.")
        
def model_factory(model_input, use_absolute_directory=False):
    """
    Given the Model_input object, which contains the settings, 
    create the Model object which contains more detailed 
    information.
    """
    
    model = base.Model()
    calc_settings = model_input.calculation_settings
    if model_input.calculation_type.lower() == "mmvt":
        model.calculation_settings = mmvt_base.MMVT_settings()
        model.calculation_settings.num_production_steps = \
            calc_settings.md_steps_per_anchor
        model.calculation_settings.energy_reporter_interval = \
            calc_settings.md_output_interval
        model.calculation_settings.restart_checkpoint_interval = \
            calc_settings.md_output_interval
        model.calculation_settings.trajectory_reporter_interval = \
            calc_settings.md_output_interval
    elif model_input.calculation_type.lower() == "elber":
        model.calculation_settings = elber_base.Elber_settings()
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
    elif model_input.ensemble.lower() == "npt":
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
        mm_settings.initial_temperature = temperature
        mm_settings.nonbonded_cutoff = model_input.nonbonded_cutoff
        mm_settings.run_minimization = model_input.run_minimization
        mm_settings.hydrogenMass = model_input.hydrogenMass
        mm_settings.constraints = model_input.constraints
        mm_settings.langevin_integrator.timestep = model_input.timestep
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
        model.namd_settings = namd_settings
    
    elif model_input.md_program.lower() == "toy":
        toy_settings = base.Toy_settings()
        model.toy_settings = toy_settings
        
    else:
        raise Exception("Invalid MD program entered:", 
                        model_input.md_program)
    
    if model_input.browndye_settings_input is None:
        # Running no BD
        pass
    
    else:
        k_on_info = base.K_on_info()
        if model_input.browndye_settings_input.ions is None:
            model_input.browndye_settings_input.ions = []
        k_on_info.ions = model_input.browndye_settings_input.ions
        k_on_info.b_surface_num_trajectories = \
            model_input.browndye_settings_input.num_b_surface_trajectories
        
        model.browndye_settings = base.Browndye_settings()
        model.browndye_settings.browndye_bin_dir = \
            model_input.browndye_settings_input.binary_directory
        model.browndye_settings.receptor_pqr_filename = \
            os.path.basename(
                model_input.browndye_settings_input.receptor_pqr_filename)
        model.browndye_settings.ligand_pqr_filename = \
            os.path.basename(
                model_input.browndye_settings_input.ligand_pqr_filename)
        model.browndye_settings.apbs_grid_spacing = \
            model_input.browndye_settings_input.apbs_grid_spacing
        model.browndye_settings.n_threads = \
            model_input.browndye_settings_input.n_threads
        
        model.k_on_info = k_on_info
    
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
            assert anchor_kept.endstate == anchor_discarded.endstate, \
                "The anchors connected by connection flag {} ".format(key) \
                +"must have the same value for 'bound_state'."
            assert anchor_kept.bulkstate == anchor_discarded.bulkstate, \
                "The anchors connected by connection flag {} ".format(key) \
                +"must have the same value for 'bulk_anchor'."   
            # Ensure that both anchors don't have starting structures 
            # defined
            assert anchor_kept.amber_params is None or \
                    anchor_kept.amber_params is None, \
                "The anchors connected by connection flag {} ".format(key) \
                +"may not both have starting structures defined. One or both" \
                "must be left blank."
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
            for milestone in anchors[new_index].milestones:
                if (model.get_type() == "mmvt") \
                        and (new_index not in visited_new_indices):
                    # don't renumber the same anchor twice
                    neighbor_id = milestone.neighbor_anchor_index
                    milestone.neighbor_anchor_index = index_reference[neighbor_id]
                    milestone.alias_index = alias_index
                    alias_index += 1
            
        visited_new_indices.append(new_index)
            
    xml_path = os.path.join(root_directory, "model.xml")
    if os.path.exists(xml_path):
        # then a model file already exists at this location: update
        # the anchor directories.
        old_model = base.Model()
        old_model.deserialize(xml_path)
        modify_model(old_model, model, root_directory,
                                    force_overwrite)
    
    for old_index in index_reference:
        new_index = index_reference[old_index]
        if new_index < len(anchors):
            input_anchor = associated_input_anchor[old_index]
            filetree.generate_filetree_by_anchor(
                anchors[new_index], root_directory)
            filetree.copy_building_files_by_anchor(
                anchors[new_index], input_anchor, root_directory)
                
    return anchors

def create_cvs_and_anchors(model, collective_variable_inputs, root_directory):
    """
    Create the collective variable and Anchor objects for the Model.
    
    Parameters:
    -----------
    model : Model()
        The model object representing the entire calculation.
    
    collective_variable_inputs : list
        A list of cv input objects used to construct the anchors for
        this model.
        
    root_directory : str
        A path to the root directory of this model.
    
    Returns:
    --------
    cvs : list
        A list of collective variable objects that will be placed into
        the model.
    anchors : list
        A list of Anchor() objects. This list and its associated anchors
        may be modified before being placed into the final model.
    milestone_index : int
        A running count of the number of milestones in this model.
    connection_flag_dict : dict
        A dictionary whose keys are integers representing a connection
        index, and whose values are lists of anchors that have a connection
        together, and therefore must be combined.
    associated_input_anchor : dist
        A dictionary whose keys are anchor indices, and whose values are
        the Input_anchor() objects associated with that anchor.
    """
    
    milestone_index = 0
    anchor_index = 0
    cv_indices = []
    cvs = []
    anchors = []
    connection_flag_dict = defaultdict(list)
    associated_input_anchor = {}
    for i, cv_input in enumerate(collective_variable_inputs):
        cv_input.index = i
        cv_input.check()
        if model.get_type() == "mmvt":
            if isinstance(cv_input, common_cv.Spherical_cv_input):
                cv = mmvt_cv.make_mmvt_spherical_cv_object(cv_input, index=i)
            elif isinstance(cv_input, common_cv.Tiwary_cv_input):
                cv = mmvt_cv.make_mmvt_tiwary_cv_object(cv_input, index=i)
            elif isinstance(cv_input, common_cv.Planar_cv_input):
                cv = mmvt_cv.make_mmvt_planar_cv_object(cv_input, index=i)
            else:
                raise Exception("CV type not implemented in MMVT: %s" \
                                % type(cv_input))
            
        elif model.get_type() == "elber":
            if isinstance(cv_input, common_cv.Spherical_cv_input):
                cv = elber_cv.make_elber_spherical_cv_object(cv_input, index=i)
            else:
                raise Exception("CV type not implemented in Elber: %s" \
                                % type(cv_input))
        
        cvs.append(cv)
        cv_indices.append(i)
        input_anchor_index = 0
        for j, input_anchor in enumerate(cv_input.input_anchors):
            if model.get_type() == "mmvt":
                anchor = mmvt_base.MMVT_anchor()
            elif model.get_type() == "elber":
                anchor = elber_base.Elber_anchor()
            anchor.index = anchor_index
            anchor.name = "anchor_"+str(anchor_index)
            anchor.directory = anchor.name
            
            if not input_anchor.bulk_anchor:
                anchor.md = True
            else:
                anchor.md = False
                anchor.bulkstate = True
                
            if input_anchor.bound_state:
                anchor.endstate = True
            
            input_anchor.check(j, cv_input)
            
            milestone_alias = 1
            if model.get_type() == "mmvt":
                milestones, milestone_alias, milestone_index = \
                    cv_input.make_mmvt_milestoning_objects(
                    milestone_alias, milestone_index, input_anchor_index,
                    anchor_index)
                
            elif model.get_type() == "elber":
                milestones, milestone_alias, milestone_index = \
                    elber_cv.make_elber_milestoning_objects_spherical(
                    cv_input, milestone_alias, milestone_index, 
                    input_anchor_index, anchor_index, cv_input.input_anchors, 
                    model.calculation_settings.umbrella_force_constant)
                
            anchor.milestones += milestones
            variable_name = "{}_{}".format(cv.variable_name, i)
            variable_value = input_anchor.get_variable_value()
            anchor.variables[variable_name] = variable_value
            for connection_flag in input_anchor.connection_flags:
                connection_flag_dict[connection_flag].append(anchor)
            if input_anchor.bulk_anchor:
                connection_flag_dict["bulk"].append(anchor)
            #filetree.generate_filetree_by_anchor(anchor, root_directory)
            #filetree.copy_building_files_by_anchor(anchor, input_anchor, 
            #                                       root_directory)
            associated_input_anchor[anchor.index] = input_anchor
            anchors.append(anchor)
            anchor_index += 1
            input_anchor_index += 1
    
    if model.get_type() == "elber":
        milestone_index += 1
    
    return cvs, anchors, milestone_index, connection_flag_dict, \
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
            
            elif model.get_type() == "elber":
                bd_milestone.outer_milestone = anchor.milestones[1]
                assert "radius" in bd_milestone.outer_milestone.variables,\
                    "A BD outer milestone must be spherical."
                bd_milestone.inner_milestone = anchor.milestones[0]
            
            cv_index = bd_milestone.outer_milestone.cv_index
            cv_input = model_input.cv_inputs[cv_index]
            if len(cv_input.bd_group1)>0:
                bd_milestone.receptor_indices = cv_input.bd_group1
            else:
                bd_milestone.receptor_indices = \
                    model_input.browndye_settings_input.receptor_indices
            
            if len(cv_input.bd_group2)>0:
                bd_milestone.ligand_indices = cv_input.bd_group2
            else:
                bd_milestone.ligand_indices = \
                    model_input.browndye_settings_input.ligand_indices
            
            model.k_on_info.bd_milestones.append(bd_milestone)
            bd_milestone_counter += 1
                
    return

def prepare_model_cvs_and_anchors(model, model_input, force_overwrite):
    """
    Fill out the CV and anchor items within the model based on the 
    input objects.
    """
    
    cvs, anchors, num_milestones, connection_flag_dict, associated_input_anchor\
        = create_cvs_and_anchors(
            model, model_input.cv_inputs, model_input.root_directory)
    model.collective_variables = cvs
    model.anchors = anchors
    if model.get_type() == "mmvt":
        if model_input.md_program.lower() == "openmm":
            for anchor in model.anchors:
                anchor.md_output_glob = mmvt_base.OPENMMVT_GLOB
        elif model_input.md_program.lower() == "namd":
            for anchor in model.anchors:
                anchor.md_output_glob = mmvt_base.NAMDMMVT_GLOB
        model.num_milestones = num_milestones
        
    elif model.get_type() == "elber":
        if model_input.md_program.lower() == "openmm":
            for anchor in model.anchors:
                anchor.md_output_glob = elber_base.OPENMM_ELBER_GLOB
        elif model_input.md_program.lower() == "namd":
            for anchor in model.anchors:
                anchor.md_output_glob = elber_base.NAMD_ELBER_GLOB
        model.num_milestones = num_milestones-1
        
    if model_input.browndye_settings_input is not None:
        create_bd_milestones(model, model_input)
    
    anchors = resolve_connections(connection_flag_dict, model, 
                                  associated_input_anchor, 
                                  model_input.root_directory, force_overwrite)
    # check to make sure that anchors don't have repeated milestone aliases
    # or have themselves as neighbors.
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
        
    model.num_anchors = len(anchors)
    
    return
    
def generate_bd_files(model, rootdir):
    """
    Create the ghost atoms in the BD files, convert PQRs to XMLs,
    make the input.xml file, as well as the reaction criteria.
    """
    
    if model.k_on_info is not None and model.browndye_settings is not None:
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
    # check whether deleting directories contain simulation files
    deleting_directories_with_files = []
    for del_anchor_index in old_anchors_to_delete:
        del_anchor = old_model.anchors[del_anchor_index]
        if anchor_has_files(old_model, del_anchor):
            deleting_directories_with_files.append(del_anchor.directory)
    warning_msg = ""
    if len(deleting_directories_with_files) > 0 :
        warning_msg+="The anchor(s) {} ".format(
                deleting_directories_with_files)\
            +"already have existing output files "\
            "and the entered command would overwrite them. If you desire "\
            "to overwrite the existing files, then use the "\
            "--force_overwrite (-f) option, and these anchors will be "\
            "deleted in directory: " + root_directory + "\n"
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
            "deleted in directory: " + root_directory + "\n"
    if warning_msg and not force_overwrite:
        print(warning_msg)
        raise Exception("Cannot overwrite existing outputs.")
    
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
    
    for new_anchor_index_to_create in new_anchors_to_create:
        print("creating new anchor {}".format(new_anchor_index_to_create))
    #    anchor = new_model.anchors[new_anchor_index_to_create]
    #    filetree.generate_filetree_by_anchor(anchor, root_directory)
    #    copy_building_files_by_anchor(anchor, input_anchor, rootdir)
    
    # TODO: check for BD milestone changes
    if old_model.k_on_info is not None:
        bd_milestones_to_check = []
        bd_milestones_to_cleanse = []
        b_surface_directory = os.path.join(
            old_model.anchor_rootdir, 
            old_model.k_on_info.b_surface_directory)
        if new_model.k_on_info is None:
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
            return
        
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
        
    return
    