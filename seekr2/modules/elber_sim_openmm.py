"""
elber_sim_openmm.py

Create the Sim_openmm object which contains the objects needed for an
openmm simulation based on the settings provided by a user. These
objects are specific to Elber milestoning only.
"""

# TODO: update this code and documentation once the Elber plugin is
# fixed.

import os
import re

try:
    import openmm
except ImportError:
    import simtk.openmm as openmm #
    
try:
    import openmm.app as openmm_app
except ImportError:
    import simtk.openmm.app as openmm_app
    
from parmed import unit

import seekr2plugin
import seekr2.modules.elber_cvs.elber_cv_base as elber_cv_base
import seekr2.modules.common_sim_openmm as common_sim_openmm

class Elber_sim_openmm(common_sim_openmm.Common_sim_openmm):
    """
    
    """
    def __init__(self):
        super(Elber_sim_openmm, self).__init__()
        self.umbrella_system = None
        self.umbrella_integrator = None
        self.umbrella_simulation = None
        self.umbrella_traj_reporter = openmm_app.DCDReporter
        self.umbrella_energy_reporter = openmm_app.StateDataReporter
        self.umbrella_force = None
        
        self.rev_system = None
        self.rev_integrator = None
        self.rev_simulation = None
        self.rev_traj_reporter = openmm_app.DCDReporter
        self.rev_energy_reporter = openmm_app.StateDataReporter
        self.rev_output_filename = ""
        
        self.fwd_system = None
        self.fwd_integrator = None
        self.fwd_simulation = None
        self.fwd_traj_reporter = openmm_app.DCDReporter
        self.fwd_energy_reporter = openmm_app.StateDataReporter
        self.fwd_output_filename = ""
        return
    
def add_integrators(sim_openmm, model, state_prefix=None):
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
        timestep = model.get_timestep()
        rigid_constraint_tolerance = \
            model.openmm_settings.langevin_integrator\
            .rigid_tolerance
            
        sim_openmm.timestep = timestep
        
        integrator_type = model.openmm_settings.langevin_integrator.integrator_type
        
        if random_seed == 0:
            rev_seed = 0
            fwd_seed = 0
        else:
            rev_seed = random_seed + 1
            fwd_seed = random_seed + 2
        
        if integrator_type == "langevin":
            umbrella_integrator_object = openmm.LangevinIntegrator
            fwd_rev_integrator_object = seekr2plugin.ElberLangevinIntegrator
        elif integrator_type == "langevinMiddle":
            umbrella_integrator_object = openmm.LangevinMiddleIntegrator
            fwd_rev_integrator_object \
                = seekr2plugin.ElberLangevinMiddleIntegrator
        else:
            raise Exception("Settings not provided for available "\
                        "integrator type(s).")
        
        sim_openmm.umbrella_integrator = umbrella_integrator_object(
            target_temperature*unit.kelvin, 
            friction_coefficient/unit.picoseconds, 
            timestep*unit.picoseconds)
        if random_seed is not None:
            sim_openmm.umbrella_integrator.setRandomNumberSeed(random_seed)
            
        #sim_openmm.rev_integrator = openmm.VerletIntegrator(
        #    timestep*unit.picoseconds)
        sim_openmm.rev_integrator = fwd_rev_integrator_object(
            target_temperature*unit.kelvin, 
            friction_coefficient/unit.picoseconds, 
            timestep*unit.picoseconds,
            sim_openmm.rev_output_filename)
        if random_seed is not None:
            sim_openmm.rev_integrator.setRandomNumberSeed(rev_seed)
        sim_openmm.rev_integrator.setEndOnSrcMilestone(True)
        
        #sim_openmm.fwd_integrator = openmm.VerletIntegrator(
        #    timestep*unit.picoseconds)
        sim_openmm.fwd_integrator = fwd_rev_integrator_object(
            target_temperature*unit.kelvin, 
            friction_coefficient/unit.picoseconds, 
            timestep*unit.picoseconds,
            sim_openmm.fwd_output_filename)
        if random_seed is not None:
            sim_openmm.fwd_integrator.setRandomNumberSeed(fwd_seed)
        sim_openmm.fwd_integrator.setEndOnSrcMilestone(False)
        
        if rigid_constraint_tolerance is not None:
            sim_openmm.umbrella_integrator.setConstraintTolerance(
                rigid_constraint_tolerance)
            sim_openmm.rev_integrator.setConstraintTolerance(
                rigid_constraint_tolerance)
            sim_openmm.fwd_integrator.setConstraintTolerance(
                rigid_constraint_tolerance)
        
        if state_prefix is not None:
            sim_openmm.fwd_integrator.setSaveStateFileName(state_prefix)
        
    else:
        raise Exception("Settings not provided for available "\
                        "integrator type(s).")
    return

def add_forces(sim_openmm, model, anchor):
    """
    Add the proper forces for this Elber simulation to each of their 
    respective systems. This will include the umbrella force for the
    umbrella stage, as well as the forces used to monitor for boundary
    crossings in the reversal and forward stages.
    """
    sim_openmm.umbrella_force = make_elber_umbrella_force(
        model, anchor)
    sim_openmm.umbrella_system.addForce(sim_openmm.umbrella_force)
    
    for milestone in anchor.milestones:
        cv = milestone.get_CV(model)
        my_rev_force = make_elber_boundary_definitions(
            cv, milestone)
        my_fwd_force = make_elber_boundary_definitions(
            cv, milestone)
        if milestone.is_source_milestone:
            sim_openmm.rev_integrator.addSrcMilestoneGroup(
                milestone.alias_index)
            sim_openmm.fwd_integrator.addSrcMilestoneGroup(
                milestone.alias_index)
        else:
            sim_openmm.rev_integrator.addDestMilestoneGroup(
                milestone.alias_index)
            sim_openmm.fwd_integrator.addDestMilestoneGroup(
                milestone.alias_index)
        revforcenum = sim_openmm.rev_system.addForce(my_rev_force)
        fwdforcenum = sim_openmm.fwd_system.addForce(my_fwd_force)
    return

def add_simulations(sim_openmm, model, topology, positions, box_vectors):
    """
    
    """
    sim_openmm.umbrella_simulation = openmm_app.Simulation(
        topology, sim_openmm.umbrella_system, 
        sim_openmm.umbrella_integrator, sim_openmm.platform, 
        sim_openmm.properties)
    
    #for force in self.sim_openmm.rev_system.getForces():
    #    print("force name:", force.__class__.__name__,
    #          "force.usesPeriodicBoundaryConditions():", 
    #          force.usesPeriodicBoundaryConditions())
    #exit()
    sim_openmm.rev_simulation = openmm_app.Simulation(
        topology, sim_openmm.rev_system, 
        sim_openmm.rev_integrator, sim_openmm.platform, 
        sim_openmm.properties)
    sim_openmm.fwd_simulation = openmm_app.Simulation(
        topology, sim_openmm.fwd_system, 
        sim_openmm.fwd_integrator, sim_openmm.platform, 
        sim_openmm.properties)
    
    sim_openmm.umbrella_simulation.context.setPositions(positions)
    if box_vectors is not None:
        sim_openmm.umbrella_simulation.context.setPeriodicBoxVectors(
            *box_vectors.to_quantity())
    if model.openmm_settings.run_minimization:
        sim_openmm.umbrella_simulation.minimizeEnergy()
    
    sim_openmm.umbrella_simulation.context.setVelocitiesToTemperature(
        model.openmm_settings.initial_temperature * unit.kelvin)
    assert sim_openmm.timestep is not None
    return

def create_sim_openmm(model, anchor, output_filename, state_prefix=None):
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
        The name of the directory that will be written by the plugin as
        the Elber simulation proceeds, recording every 'bounce'.
        
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
    sim_openmm = Elber_sim_openmm()
    fwd_output_filename = output_filename
    directory = os.path.dirname(output_filename)
    basename = os.path.basename(output_filename)
    suffix = re.sub(elber_cv_base.OPENMM_ELBER_BASENAME, "", basename)
    rev_basename = elber_cv_base.ELBER_REV_BASENAME+suffix
    rev_output_filename = os.path.join(directory, rev_basename)
    sim_openmm.rev_output_filename = rev_output_filename
    sim_openmm.fwd_output_filename = fwd_output_filename
    sim_openmm.output_filename = output_filename
    umbrella_system, umbrella_topology, umbrella_positions, \
        umbrella_box_vectors, umbrella_num_frames \
        = common_sim_openmm.create_openmm_system(sim_openmm, model, anchor)
    rev_system, rev_topology, rev_positions, rev_box_vectors, rev_num_frames \
        = common_sim_openmm.create_openmm_system(sim_openmm, model, anchor)
    fwd_system, fwd_topology, fwd_positions, fwd_box_vectors, fwd_num_frames \
        = common_sim_openmm.create_openmm_system(sim_openmm, model, anchor)
    sim_openmm.umbrella_system = umbrella_system
    sim_openmm.rev_system = rev_system
    sim_openmm.fwd_system = fwd_system
    add_integrators(sim_openmm, model, state_prefix=state_prefix)
    common_sim_openmm.add_barostat(umbrella_system, model)
    common_sim_openmm.add_barostat(rev_system, model)
    common_sim_openmm.add_barostat(fwd_system, model)
    common_sim_openmm.add_platform(sim_openmm, model)
    add_forces(sim_openmm, model, anchor)
    add_simulations(sim_openmm, model, umbrella_topology, umbrella_positions, 
                         umbrella_box_vectors)
    return sim_openmm

def make_elber_umbrella_force(model, anchor):
    """
    
    """
    
    for milestone in anchor.milestones:
        #if milestone.alias_index == 2:
        if milestone.is_source_milestone:
            cv = milestone.get_CV(model)
            umbrella_force = cv.make_umbrella_force_object()
            break
        
    mygroup_list = []
    #mygroup1 = umbrella_force.addGroup(cv.group1)
    #mygroup_list.append(mygroup1)
    #mygroup2 = umbrella_force.addGroup(cv.group2)
    #mygroup_list.append(mygroup2)
    for group in cv.get_atom_groups():
        mygroup = umbrella_force.addGroup(group)
        mygroup_list.append(mygroup)
    
    variable_names_list = cv.add_umbrella_parameters(umbrella_force)
    cv.add_groups_and_variables(umbrella_force, mygroup_list,
                                cv.get_variable_values_list(milestone))
    return umbrella_force

def make_elber_boundary_definitions(cv, milestone):
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
    myforce = cv.make_fwd_rev_force_object()
    myforce.setForceGroup(milestone.alias_index)
    variable_names_list, group_list = cv.add_fwd_rev_parameters(myforce)
    cv.add_groups_and_variables(myforce, group_list,
                                cv.get_variable_values_list(
                                milestone))
    return myforce