"""
common/prepare.py

Common routines for preparing the model file for either Elber or MMVT
milestoning calculations.

The resulting model file will be used for all other stages of the 
calculation.
"""

import os
import sys

import seekr2.common.filetree as filetree
import seekr2.common.base as base
import seekr2.mmvt.base as mmvt_base
import seekr2.elber.base as elber_base
import seekr2.common.collective_variables as common_cv
import seekr2.mmvt.collective_variables as mmvt_cv
import seekr2.elber.collective_variables as elber_cv
from seekr2.common.base import Ion
import seekr2.common.sim_browndye2 as sim_browndye2
import seekr2.common.runner_browndye2 as runner_browndye2
import seekr2.libraries.serializer.serializer as serializer

from seekr2.common.base import Forcefield_params

class Browndye_settings_input(serializer.Serializer):
    """
    Read and parse the outputs from the BrownDye program, which runs
    the BD stage of the MMVT SEEKR calculation
    
    Attributes:
    -----------
    binary_directory : str
        A path to the BrownDye programs directory. Can be left as an
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
        A list of common.base.Ion() objects which will be passed to
        APBS.
    
    num_b_surface_trajectories : int
        The number of trajectories to run with the ligand starting at
        the b-surface.
        
    num_bd_milestone_trajectories : int
        The number of trajectories to run with the ligand starting at
        each of the BD milestones.
        
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
        self.num_bd_milestone_trajectories = -1
        self.receptor_indices = []
        self.ligand_indices = []
        self.n_threads = 1

class Model_input(serializer.Serializer):
    """
    The serializable object representing parameters that would be
    input by the user.
    
    Attributes:
    -----------
    calculation_type : str
        A string to represent the calculation type.
    
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
        
    md_output_frequency : int, default 500000
        How many steps must elapse between outputs. This will effect
        the frequency that trajectory frames, state outputs, and
        backup checkpoints will occur.
        
    md_steps_per_anchor : int, default 500000000
        The default number of timesteps that will be simulated at each
        anchor.
        
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
    
    browndye_settings_input : Browndye_settings_input
        The Browndye_settings_input() object for this model. It 
        contains all the settings that could be used within a Browndye
        simulation.
    
    cv_inputs : list
        A list containing a set of collective variable input objects.
    """
    def __init__(self):
        self.calculation_type = "Unknown"
        self.temperature = 298.15
        self.pressure = 1.0
        self.ensemble = "nvt"
        self.root_directory = ""
        self.md_program = "openmm"
        self.md_output_frequency = 500000
        self.md_steps_per_anchor = 500000000
        self.run_minimization = False
        self.hydrogenMass = None
        self.constraints = "hbonds"
        self.rigidWater = True
        self.timestep = 0.002
        self.browndye_settings_input = None
        self.cv_inputs = []
        
    def read_plain_input_file(self, filename):
        """
        Read a plain input file (as opposed to an XML)
        """
        raise Exception("Reading a plain text file is not yet implemented. "\
                        "Only an XML model input may be read at this time.")
        
class Model_factory():
    """
    Create the Model object that will be used throughout the entire 
    SEEKR calculation
    """
    def __init__(self):
        return
    
    def create_model(self, model_input, use_absolute_directory=False):
        """
        Given the Model_input object, which contains the settings, 
        create the Model object which contains more detailed 
        information.
        """
        model = base.Model()
        if model_input.calculation_type.lower() == "mmvt":
            model.calculation_settings = mmvt_base.MMVT_settings()
            model.calculation_settings.num_production_steps = \
                model_input.md_steps_per_anchor
        elif model_input.calculation_type.lower() == "elber":
            model.calculation_settings = elber_base.Elber_settings()
            model.calculation_settings.num_umbrella_stage_steps = \
                model_input.md_steps_per_anchor
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
            mm_settings.energy_reporter_frequency = \
                model_input.md_output_frequency
            mm_settings.restart_checkpoint_frequency = \
                model_input.md_output_frequency
            mm_settings.trajectory_reporter_frequency = \
                model_input.md_output_frequency
            
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
            namd_settings.energy_reporter_frequency = \
                model_input.md_output_frequency
            namd_settings.restart_checkpoint_frequency = \
                model_input.md_output_frequency
            namd_settings.trajectory_reporter_frequency = \
                model_input.md_output_frequency
            namd_settings.total_simulation_length = \
                model_input.md_steps_per_anchor
            namd_settings.run_minimization = model_input.run_minimization
            assert model_input.hydrogenMass is None, "hydrogen mass "\
                "repartitioning not yet implemented in NAMD SEEKR."
            #namd_settings.hydrogenMass = model_input.hydrogenMass
            namd_settings.constraints = model_input.constraints
            namd_settings.langevin_integrator.timestep = model_input.timestep
            namd_settings.rigidWater = model_input.rigidWater
            model.namd_settings = namd_settings
            
        else:
            raise Exception("Invalid MD program entered:", 
                            model_input.md_program)
        
        if model_input.browndye_settings is None:
            # Running no BD
            pass
        
        else:
            k_on_info = base.K_on_info()
            k_on_info.ions = model_input.browndye_settings.ions
            k_on_info.b_surface_num_trajectories = \
                model_input.browndye_settings.num_b_surface_trajectories
            
            model.browndye_settings = base.Browndye_settings()
            model.browndye_settings.browndye_bin_dir = \
                model_input.browndye_settings.binary_directory
            model.browndye_settings.receptor_pqr_filename = \
                os.path.basename(
                    model_input.browndye_settings.receptor_pqr_filename)
            model.browndye_settings.ligand_pqr_filename = \
                os.path.basename(
                    model_input.browndye_settings.ligand_pqr_filename)
            model.browndye_settings.apbs_grid_spacing = \
                model_input.browndye_settings.apbs_grid_spacing
            model.browndye_settings.n_threads = \
                model_input.browndye_settings.n_threads
            
            model.k_on_info = k_on_info
        
        return model
        
def create_cvs_and_anchors(model, collective_variable_inputs):
    """
    Create the collective variable and Anchor objects for the Model.
    """
    milestone_index = 0
    anchor_index = 0
    cv_indices = []
    cvs = []
    anchors = []
    for i, cv_input in enumerate(collective_variable_inputs):
        cv_input.index = i
        cv_input.check()
        if model.get_type() == "mmvt":
            if isinstance(cv_input, common_cv.Spherical_cv_input):
                cv = mmvt_cv.make_mmvt_spherical_cv_object(cv_input, index=i)
            else:
                raise Exception("CV type not implemented: %s" % type(cv_input))
            
        elif model.get_type() == "elber":
            if isinstance(cv_input, common_cv.Spherical_cv_input):
                cv = elber_cv.make_elber_spherical_cv_object(cv_input, index=i)
            else:
                raise Exception("CV type not implemented: %s" % type(cv_input))
        
        cvs.append(cv)
        cv_indices.append(i)
        
        for input_anchor in cv_input.input_anchors:
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
            
            milestone_alias = 1
            if model.get_type() == "mmvt":
                milestones, milestone_alias, milestone_index = \
                    mmvt_cv.make_mmvt_milestoning_objects_spherical(
                    cv_input, milestone_alias, milestone_index, anchor_index, 
                    cv_input.input_anchors)
            elif model.get_type() == "elber":
                milestones, milestone_alias, milestone_index = \
                    elber_cv.make_elber_milestoning_objects_spherical(
                    cv_input, milestone_alias, milestone_index, anchor_index, 
                    cv_input.input_anchors)
            anchor.milestones += milestones
            anchors.append(anchor)
            anchor_index += 1
    
    return cvs, anchors, anchor_index, milestone_index

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
            
            bd_milestone.num_trajectories = \
                model_input.browndye_settings.num_bd_milestone_trajectories
            bd_milestone.receptor_indices = \
                model_input.browndye_settings.receptor_indices
            bd_milestone.ligand_indices = \
                model_input.browndye_settings.ligand_indices
            
            model.k_on_info.bd_milestones.append(bd_milestone)
            bd_milestone_counter += 1
                
    return

def prepare_model_cvs_and_anchors(model, model_input):
    """
    Fill out the CV and anchor items within the model based on the 
    input objects.
    """
    cvs, anchors, num_anchors, num_milestones = create_cvs_and_anchors(
        model, model_input.cv_inputs)
            
    model.collective_variables = cvs
    model.anchors = anchors
    if model.get_type() == "mmvt":
        if model_input.md_program.lower() == "openmm":
            for anchor in model.anchors:
                anchor.md_mmvt_output_glob = mmvt_base.OPENMMVT_GLOB
        elif model_input.md_program.lower() == "namd":
            for anchor in model.anchors:
                anchor.md_mmvt_output_glob = mmvt_base.NAMDMMVT_GLOB
    
    elif model.get_type() == "elber":
        if model_input.md_program.lower() == "openmm":
            for anchor in model.anchors:
                anchor.md_mmvt_output_glob = elber_base.OPENMM_ELBER_GLOB
        elif model_input.md_program.lower() == "namd":
            for anchor in model.anchors:
                anchor.md_mmvt_output_glob = elber_base.NAMD_ELBER_GLOB
    
    model.num_anchors = num_anchors
    model.num_milestones = num_milestones
    if model_input.browndye_settings is not None:
        create_bd_milestones(
            model, model_input)
    
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
            print("adding ghost atom to file:", receptor_pqr_filename)
            ghost_index_rec = \
                sim_browndye2.add_ghost_atom_to_pqr_from_atoms_center_of_mass(
                    receptor_pqr_filename, bd_milestone.receptor_indices)
            print("adding ghost atom to file:", ligand_pqr_filename)
            ghost_index_lig = \
                sim_browndye2.add_ghost_atom_to_pqr_from_atoms_center_of_mass(
                    ligand_pqr_filename, bd_milestone.ligand_indices)
            ghost_indices_rec.append(ghost_index_rec)
            ghost_indices_lig.append(ghost_index_lig)
            
        model.browndye_settings.ghost_indices_rec = ghost_indices_rec
        model.browndye_settings.ghost_indices_lig = ghost_indices_lig
        
        receptor_xml_filename = sim_browndye2.make_pqrxml(receptor_pqr_filename)
        ligand_xml_filename = sim_browndye2.make_pqrxml(ligand_pqr_filename)
        debye_length, reaction_filename = \
            runner_browndye2.make_browndye_input_xml(
            model, rootdir, receptor_xml_filename, ligand_xml_filename,
            model.k_on_info.b_surface_num_trajectories)
        model.browndye_settings.debye_length = debye_length
        
        abs_reaction_path = os.path.join(b_surface_dir, 
                                         reaction_filename)
        
        runner_browndye2.make_browndye_reaction_xml(model, abs_reaction_path)
        
        return