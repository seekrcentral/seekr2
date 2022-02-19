"""
mmvt_sim_openmm.py

Create the Sim_openmm object which contains the objects needed for an
openmm simulation based on the settings provided by a user. These
objects are specific to MMVT only.
"""
import os

import numpy as np
try:
    import openmm.app as openmm_app
except ImportError:
    import simtk.openmm.app as openmm_app
    
from parmed import unit

import seekr2plugin
import seekr2.modules.common_sim_openmm as common_sim_openmm

class MMVT_sim_openmm(common_sim_openmm.Common_sim_openmm):
    """
    system : The OpenMM system object for this simulation.
    
    integrator : the OpenMM integrator object for this simulation.
    
    simulation : the OpenMM simulation object.
    
    traj_reporter : openmm.DCDReporter
        The OpenMM Reporter object to which the trajectory will be
        written.
        
    energy_reporter : openmm.StateDataReporter
        The OpenMM StateDataReporter to which the energies and other
        state data will be reported.
        
    """
    def __init__(self):
        super(MMVT_sim_openmm, self).__init__()
        self.system = None
        self.integrator = None
        self.simulation = None
        self.traj_reporter = openmm_app.DCDReporter
        self.energy_reporter = openmm_app.StateDataReporter
        return
        
def add_integrator(sim_openmm, model, state_prefix=None):
    """
    Assign the proper integrator to this OpenMM simulation.
    """
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
        
        integrator_type = model.openmm_settings.langevin_integrator\
            .integrator_type
        
        if integrator_type == "langevin":
            integrator_object = seekr2plugin.MmvtLangevinIntegrator
        elif integrator_type == "langevinMiddle":
            integrator_object = seekr2plugin.MmvtLangevinMiddleIntegrator
        else:
            raise Exception("Settings not provided for available "\
                        "integrator type(s).")
        
        sim_openmm.integrator = integrator_object(
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
    return

def add_forces(sim_openmm, model, anchor):
    """
    Add the proper forces for this MMVT simulation.
    """
    curdir = os.getcwd()
    os.chdir(model.anchor_rootdir)
    for milestone in anchor.milestones:
        cv = milestone.get_CV(model)
        myforce = make_mmvt_boundary_definitions(
            cv, milestone)
        sim_openmm.integrator.addMilestoneGroup(milestone.alias_index)
        forcenum = sim_openmm.system.addForce(myforce)
        
    os.chdir(curdir)
    return

def add_simulation(sim_openmm, model, topology, positions, box_vectors):
    """
    Assign the OpenMM simulation object for MMVT.
    """
        
        
    sim_openmm.simulation = openmm_app.Simulation(
        topology, sim_openmm.system, 
        sim_openmm.integrator, sim_openmm.platform, 
        sim_openmm.properties)
    
    if positions is not None:
        sim_openmm.simulation.context.setPositions(positions)
    
    sim_openmm.simulation.context.setVelocitiesToTemperature(
        model.openmm_settings.initial_temperature * unit.kelvin)
        
    if box_vectors is not None:
        sim_openmm.simulation.context.setPeriodicBoxVectors(
            *box_vectors.to_quantity())
    if model.openmm_settings.run_minimization:
        assert positions is not None, "If states are being loaded as starting"\
            "positions, minimizations cannot be activated."
        print("Warning: running minimizations. It is recommended that "\
              "structures are minimized and verified by the user before "\
              "running SEEKR, since minimizations might cause the system "\
              "to drift out of the MMVT cell.")
        sim_openmm.simulation.minimizeEnergy()
    
    
    assert sim_openmm.timestep is not None
    return

def create_sim_openmm(model, anchor, output_filename, state_prefix=None, 
                      frame=0):
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
    
    frame : int
        Which frame of the starting positions file to retrieve.
    
    Returns
    -------
    sim_openmm : Sim_openmm()
        The sim_openmm object, which contains everything needed
        to run an MMVT simulation within OpenMM.
        
    """
    sim_openmm = MMVT_sim_openmm()
    common_sim_openmm.fill_generic_parameters(
        sim_openmm, model, anchor, output_filename)
    system, topology, positions, box_vectors, num_frames \
        = common_sim_openmm.create_openmm_system(sim_openmm, model, anchor, 
                                                 frame)
    sim_openmm.system = system
    add_integrator(sim_openmm, model, state_prefix=state_prefix)
    common_sim_openmm.add_barostat(sim_openmm, model)
    common_sim_openmm.add_platform(sim_openmm, model)
    add_forces(sim_openmm, model, anchor)
    add_simulation(sim_openmm, model, topology, positions, box_vectors)
    return sim_openmm

def make_mmvt_boundary_definitions(cv, milestone):
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
    myforce.setForceGroup(1)
    variable_names_list = cv.add_parameters(myforce)
    cv.add_groups_and_variables(myforce, cv.get_variable_values_list(
                                    milestone))
    return myforce

def get_starting_structure_num_frames(model, anchor, dummy_outfile):
    """
    For an anchor's starting structure, find and return the number of frames.
    """
    sim_openmm = MMVT_sim_openmm()
    common_sim_openmm.fill_generic_parameters(
        sim_openmm, model, anchor, dummy_outfile)
    dummy_system, dummy_topology, positions, dummy_box_vectors, num_frames \
        = common_sim_openmm.create_openmm_system(sim_openmm, model, anchor)
    
    return num_frames
    