"""
sim_openmm.py

Create the Sim_openmm object which contains the objects needed for an
openmm simulation based on the settings provided by a user
"""

import os

from simtk import openmm
import simtk.openmm.app as openmm_app
from simtk import unit

#import openmmvt.base
import openmmvt.base as base
from mmvtplugin import MmvtLangevinIntegrator, vectori, vectord
import mmvtplugin

def make_boundary_definitions(cv, milestone):
    """
    Take a Collective_variable object and a particular milestone and
    return an OpenMM Force() object that the plugin can use to monitor
    crossings.
    
    Parameters
    ----------
    cv : Collective_variable()
        A Collective_variable object which contains all the information
        for the collective variable describine this variable. In fact,
        the boundaries are contours of the function described by cv.
        This variable contains information like the groups of atoms
        involved with the CV, and the expression which describes the
        function.
        
    milestone : Milestone()
        A Milestone object which describes the boundary between two
        Voronoi cells. This variable contains information like the 
        values of the variables which will be entered into the 
        Force() object.
        
    Returns
    -------
    myforce : openmm.Force()
        An OpenMM force object which does not affect atomic motion, but
        allows us to conveniently monitor a function of atomic 
        position.
    """
    myforce = cv.make_force_object()
    mygroup_list = []
    mygroup1 = myforce.addGroup(cv.group1)
    mygroup_list.append(mygroup1)
    mygroup2 = myforce.addGroup(cv.group2)
    mygroup_list.append(mygroup2)
        
    myforce.setForceGroup(milestone.alias_index)
    
    variable_names_list = cv.add_parameters(myforce)
    cv.add_groups_and_variables(myforce, mygroup_list, 
                                cv.get_variable_values_list(
                                    milestone))
    return myforce

class Sim_openmm:
    """
    Contain all the information necessary to run an MMVT SEEKR 
    calculation in OpenMM.
    
    system : openmm.System()
        The OpenMM System object for creating the simulation.
        
    integrator : MmvtLangevinIntegrator()
        The MMVT Langevin Integrator implemented within the plugin.
        
    simulation : openmm.Simulation()
        The OpenMM Simulation object which will be running the
        simulation itself.
        
    platform : openmm.Platform()
        The platform to use for the simulation. The CUDA platform
        is recommended, but Reference is also available.
        
    properties : dict
        A dictionary of properties passed to the Simulation() object
        that apply to the platform.
        
    traj_reporter : openmm.DCDReporter
        The OpenMM Reporter object to which the trajectory will be
        written.
        
    energy_reporter : openmm.StateDataReporter
        The OpenMM StateDataReporter to which the energies and other
        state data will be reported.
        
    timestep : float
        The timestep (in picoseconds) which was passed to the 
        integrator, but is also useful to know in other places within 
        the program.
        
    output_filename : str
        The file that the plugin will write MMVT bounces to, and which
        will be parsed by the analysis software.
    
    """
    def __init__(self):
        self.system = None
        self.integrator = None
        self.simulation = None
        self.platform = None
        self.properties = None
        self.traj_reporter = None
        self.energy_reporter = None
        self.timestep = None
        self.output_filename = ""
    
class Sim_openmm_factory:
    """
    Create the sim_openmm objects which will be used to run MMVT
    calculations in OpenMM.
    """
    def __init__(self):
        return
        
    def create_sim_openmm(self, model, anchor, output_filename, 
                          state_prefix=None):
        """
        Take all relevant model and anchor information and generate
        the necessary OpenMM objects to run the simulation.
        
        Parameters
        ----------
        model : Model()
            The Model which contains the anchor and all settings for
            the simulation that is about to be run.
            
        anchor : Anchor()
            The anchor object that this OpenMM simulation applies to.
            
        output_filename : str
            The name of the file that will be written by the plugin as
            the MMVT simulation proceeds, recording every 'bounce'.
            
        state_prefix : str or None
            The plugin can optionally save every state during a
            bounce. These can be used to seed simulations in other
            cells. This argument provides the file prefix for these
            saved states. If None, then no states will be written.
            
        Returns
        -------
        sim_openmm : Sim_openmm()
            The sim_openmm object, which contains everything needed
            to run an MMVT simulation within OpenMM.
            
        """
        sim_openmm = Sim_openmm()
        sim_openmm.output_filename = output_filename
        building_directory = os.path.join(
            model.anchor_rootdir, anchor.directory, anchor.building_directory)
        box_vectors = None
        if anchor.amber_params is not None:
            prmtop_filename = os.path.join(
                building_directory, anchor.amber_params.prmtop_filename)
            prmtop = openmm_app.AmberPrmtopFile(prmtop_filename)
            inpcrd = None
            if anchor.amber_params.inpcrd_filename is not None \
                    and anchor.amber_params.inpcrd_filename != "":
                inpcrd_filename = os.path.join(
                    building_directory, anchor.amber_params.inpcrd_filename)
                inpcrd = openmm_app.AmberInpcrdFile(inpcrd_filename)
                
            if anchor.amber_params.box_vectors is None:
                assert inpcrd is not None, "If box_vectors field is "\
                    "empty, then an inpcrd must be provided."
                assert inpcrd.boxVectors is not None, "The provided "\
                    "inpcrd file contains no box vectors: "+inpcrd_filename
                box_vectors = inpcrd.boxVectors
            else:
                box_vectors = anchor.amber_params.box_vectors
                
            if anchor.amber_params.pdb_coordinates_filename is None \
                    or anchor.amber_params.pdb_coordinates_filename == "":
                positions = inpcrd
            else:
                pdb_coordinates_filename = os.path.join(
                    building_directory, 
                    anchor.amber_params.pdb_coordinates_filename)
                positions = openmm_app.PDBFile(pdb_coordinates_filename)
            
            topology = prmtop
            
            assert box_vectors is not None, "No source of box vectors provided."
        
        elif anchor.forcefield_params is not None:
            forcefield_filenames = []
            for forcefield_filename in \
                    anchor.forcefield_params.built_in_forcefield_filenames:
                forcefield_filenames.append(forcefield_filename)
            for forcefield_filename in \
                    anchor.forcefield_params.custom_forcefield_filenames:
                forcefield_filenames.append(os.path.join(
                    building_directory, forcefield_filename))
            pdb_filename = os.path.join(building_directory, 
                                   anchor.forcefield_params.pdb_filename)
            pdb = openmm_app.PDBFile(pdb_filename)
            forcefield = openmm_app.ForceField(
                *forcefield_filenames)
            if anchor.forcefield_params.box_vectors is not None:
                box_vectors = anchor.forcefield_params.box_vectors
            
            topology = pdb
            positions = pdb
        
        elif anchor.charmm_params is not None:
            raise Exception("Charmm systems not yet implemented")
        
        else:
            raise Exception("No Amber or Charmm input settings detected.")
            
        nonbonded_method = model.openmm_settings.nonbonded_method.lower()
        if nonbonded_method == "pme":
            nonbondedMethod = openmm_app.PME
            
        elif nonbonded_method == "nocutoff":
            nonbondedMethod = openmm_app.NoCutoff
            
        elif nonbonded_method == "cutoffnonperiodic":
            nonbondedMethod = openmm_app.CutoffNonPeriodic
            
        elif nonbonded_method == "cutoffperiodic":
            nonbondedMethod = openmm_app.CutoffPeriodic
            
        elif nonbonded_method == "ewald":
            nonbondedMethod = openmm_app.Ewald
        
        else:
            raise Exception("nonbonded method not found: %s", 
                            model.openmm_settings.nonbonded_method)
        
        if model.openmm_settings.constraints is None:
            constraints_str = None
        else:
            constraints_str = model.openmm_settings.constraints
            
        if constraints_str is None:
            constraints = None
        
        elif constraints_str.lower() == "none":
            constraints = None
        
        elif constraints_str.lower() == "hbonds":
            constraints = openmm_app.HBonds
            
        elif constraints_str.lower() == "allbonds":
            constraints = openmm_app.AllBonds
            
        elif constraints_str.lower() == "hangles":
            constraints = openmm_app.HAngles
            
        else:
            raise Exception("constraints not found: %s", 
                            model.openmm_settings.constraints)
        
        hydrogenMass = model.openmm_settings.hydrogenMass
        rigidWater = model.openmm_settings.rigidWater
        
        if anchor.amber_params is not None:
            sim_openmm.system = prmtop.createSystem(
                nonbondedMethod=nonbondedMethod, 
                nonbondedCutoff=model.openmm_settings.nonbonded_cutoff, 
                constraints=constraints, hydrogenMass=hydrogenMass, 
                rigidWater=rigidWater)
        
        elif anchor.forcefield_params is not None:
            sim_openmm.system = forcefield.createSystem(
                pdb.topology, nonbondedMethod=nonbondedMethod, 
                nonbondedCutoff=model.openmm_settings.nonbonded_cutoff, 
                constraints=constraints, hydrogenMass=hydrogenMass, 
                rigidWater=rigidWater)
        
        elif anchor.charmm_params is not None:
            raise Exception("Charmm input settings not yet implemented")
            
        else:
            print("Settings for Amber or Charmm simulations not found")
        
        if model.openmm_settings.langevin_integrator is not None:
            target_temperature = \
                model.openmm_settings.langevin_integrator.target_temperature
            friction_coefficient = \
                model.openmm_settings.langevin_integrator.friction_coefficient
            random_seed = \
                model.openmm_settings.langevin_integrator.random_seed
            timestep = \
                model.openmm_settings.langevin_integrator.timestep
            rigid_constraint_tolerance = \
                model.openmm_settings.langevin_integrator\
                .rigid_tolerance
                
            sim_openmm.timestep = timestep
            
            #if random_seed is not None:
            #    raise Exception("random_seed setting not yet implemented.")
            
            sim_openmm.integrator = MmvtLangevinIntegrator(
                target_temperature*unit.kelvin, 
                friction_coefficient/unit.picoseconds, 
                timestep*unit.picoseconds, 
                sim_openmm.output_filename)
            
            if random_seed is not None:
                sim_openmm.integrator.setRandomNumberSeed(random_seed)
                
            if rigid_constraint_tolerance is not None:
                sim_openmm.integrator.setConstraintTolerance(
                    rigid_constraint_tolerance)
            
            if state_prefix is not None:
                sim_openmm.integrator.setSaveStateFileName(state_prefix)
            
        else:
            raise Exception("Settings not provided for available "\
                            "integrator type(s).")
        
        if model.openmm_settings.barostat:
            barostat = openmm.MonteCarloBarostat(
                model.openmm_settings.barostat.target_pressure, 
                model.openmm_settings.barostat.target_temperature,
                model.openmm_settings.barostat.frequency)
            sim_openmm.system.addForce(barostat)
        
        if model.openmm_settings.reference_platform:
            sim_openmm.platform = \
                openmm.Platform.getPlatformByName('Reference')
            sim_openmm.properties = {}
        elif model.openmm_settings.cuda_platform_settings is not None:
            sim_openmm.platform = openmm.Platform.getPlatformByName('CUDA')
            sim_openmm.properties = model.openmm_settings.\
                cuda_platform_settings.make_properties_dict()
        else:
            raise Exception("Platform settings not found for either CUDA "\
                            "or Reference")
        
        # create MMVT forces
        for milestone in anchor.milestones:
            cv = milestone.get_CV(model)
            myforce = make_boundary_definitions(cv, milestone)
            forcenum = sim_openmm.system.addForce(myforce)
            sim_openmm.integrator.addMilestoneGroup(milestone.alias_index)
            
        sim_openmm.simulation = openmm_app.Simulation(
            topology.topology, sim_openmm.system, 
            sim_openmm.integrator, sim_openmm.platform, sim_openmm.properties)
        
        sim_openmm.simulation.context.setPositions(positions.positions)
        if box_vectors is not None:
            sim_openmm.simulation.context.setPeriodicBoxVectors(*box_vectors)
        if model.openmm_settings.run_minimization:
            print("Warning: running minimizations. It is recommended that "\
                  "structures are minimized and verified by the user before "\
                  "running SEEKR, since minimizations might cause the system "\
                  "to drift out of the MMVT cell.")
            sim_openmm.simulation.minimizeEnergy()
        
        sim_openmm.simulation.context.setVelocitiesToTemperature(
            model.openmm_settings.initial_temperature * unit.kelvin)
        
        sim_openmm.traj_reporter = openmm_app.DCDReporter
        sim_openmm.energy_reporter = openmm_app.StateDataReporter
        
        assert sim_openmm.timestep is not None
        
        return sim_openmm