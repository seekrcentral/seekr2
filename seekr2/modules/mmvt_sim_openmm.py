"""
mmvt_sim_openmm.py

Create the Sim_openmm object which contains the objects needed for an
openmm simulation based on the settings provided by a user. These
objects are specific to MMVT only.
"""
import os
import tempfile

try:
    import openmm.app as openmm_app
except ImportError:
    import simtk.openmm.app as openmm_app
    
from parmed import unit

import seekr2plugin
import seekr2.modules.mmvt_cvs.mmvt_closest_pair_cv as mmvt_closest_pair_cv
import seekr2.modules.mmvt_cvs.mmvt_count_contacts_cv as mmvt_count_contacts_cv
import seekr2.modules.mmvt_cvs.mmvt_voronoi_cv as mmvt_voronoi_cv
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
        timestep = model.get_timestep()
        rigid_constraint_tolerance = \
            model.openmm_settings.langevin_integrator\
            .rigid_tolerance
            
        sim_openmm.timestep = timestep
        
        integrator_type = model.openmm_settings.langevin_integrator\
            .integrator_type
        
        if integrator_type == "langevin":
            # The original Langevin integrator, which is not supported
            # in OpenMM 8.2 and later.
            #integrator_object = seekr2plugin.MmvtLangevinIntegrator
            integrator_object = seekr2plugin.MmvtLangevinMiddleIntegrator
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

def assign_nonbonded_cv_info(cv, system, box_vectors):
    """
    Update the cv with the correct nonbonded information.
    """
    if system is not None:
        cv.num_system_particles = system.getNumParticles()
        forces = { force.__class__.__name__ : force for force in system.getForces() }
        reference_force = forces['NonbondedForce']
        cv.exclusion_pairs = []
        for index in range(reference_force.getNumExceptions()):
            [iatom, jatom, chargeprod, sigma, epsilon] \
                = reference_force.getExceptionParameters(index)
            cv.exclusion_pairs.append((iatom, jatom))
    
    if isinstance(cv, mmvt_closest_pair_cv.MMVT_closest_pair_CV):
        cv.cutoff_distance = 0.4 * box_vectors.get_min_length()
        
    return

def add_forces(sim_openmm, model, anchor, box_vectors):
    """
    Add the proper forces for this MMVT simulation.
    """
    curdir = os.getcwd()
    os.chdir(model.anchor_rootdir)
    for milestone in anchor.milestones:
        cv = milestone.get_CV(model)
        if isinstance(cv, mmvt_closest_pair_cv.MMVT_closest_pair_CV) \
                or isinstance(cv, mmvt_count_contacts_cv.MMVT_count_contacts_CV):
            assign_nonbonded_cv_info(cv, sim_openmm.system, box_vectors)
            
        elif isinstance(cv, mmvt_voronoi_cv.MMVT_Voronoi_CV):
            for child_cv in cv.child_cvs:
                if isinstance(child_cv, mmvt_closest_pair_cv.MMVT_closest_pair_CV) \
                        or isinstance(child_cv, mmvt_count_contacts_cv.MMVT_count_contacts_CV):
                    assign_nonbonded_cv_info(child_cv, sim_openmm.system, 
                                             box_vectors)
                #myforce = make_mmvt_boundary_definitions(
                #    child_cv, milestone)
                #forcenum = sim_openmm.system.addForce(myforce)
        
        myforce = make_mmvt_boundary_definitions(
            cv, milestone)
        forcenum = sim_openmm.system.addForce(myforce)
        sim_openmm.integrator.addMilestoneGroup(milestone.alias_index)
        
    os.chdir(curdir)
    return

def add_simulation(sim_openmm, model, topology, positions, box_vectors,
                   load_state_file=None):
    """
    Assign the OpenMM simulation object for MMVT.
    """
    sim_openmm.simulation = openmm_app.Simulation(
        topology, sim_openmm.system, 
        sim_openmm.integrator, sim_openmm.platform, 
        sim_openmm.properties)
    
    state = sim_openmm.simulation.context.getState(getPositions=True)
    old_positions = state.getPositions()
    
    if positions is not None:
        sim_openmm.simulation.context.setPositions(positions)
    
    if load_state_file is not None:
        sim_openmm.simulation.loadState(load_state_file)
        state = sim_openmm.simulation.context.getState(getPositions=True)
        positions = state.getPositions(positions)
    elif positions is not None:
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
    return positions

def create_sim_openmm(model, anchor, output_filename, state_prefix=None, 
                      frame=0, load_state_file=None, use_only_reference=False):
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
    sim_openmm.output_filename = output_filename
    system, topology, positions, box_vectors, num_frames \
        = common_sim_openmm.create_openmm_system(
            sim_openmm, model, anchor, frame, load_state_file=load_state_file)
    sim_openmm.system = system
    add_integrator(sim_openmm, model, state_prefix=state_prefix)
    common_sim_openmm.add_barostat(system, model)
    common_sim_openmm.add_platform(sim_openmm, model, use_only_reference)
    add_forces(sim_openmm, model, anchor, box_vectors)
    positions = add_simulation(
        sim_openmm, model, topology, positions, box_vectors, load_state_file)
    if anchor.__class__.__name__ == "MMVT_toy_anchor":
        out_file_name = os.path.join(model.anchor_rootdir, anchor.directory, 
                                     anchor.building_directory, "toy.pdb")
        common_sim_openmm.write_toy_pdb_file(topology, positions, out_file_name)
        
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
    myforce = cv.make_boundary_force(milestone.alias_index)
    myforce.setForceGroup(1)
    variable_names_list = cv.add_parameters(myforce)
    cv.add_groups_and_variables(myforce, cv.get_variable_values_list(
                                    milestone), milestone.alias_index)
    return myforce

def get_starting_structure_num_frames(model, anchor, dummy_outfile, 
                                      load_state_file=None):
    """
    For an anchor's starting structure, find and return the number of frames.
    """
    if load_state_file is None:
        sim_openmm = MMVT_sim_openmm()
        sim_openmm.output_filename = dummy_outfile
        dummy_system, dummy_topology, dummy_positions, dummy_box_vectors, \
            num_frames = common_sim_openmm.create_openmm_system(
                sim_openmm, model, anchor)
    else:
        num_frames = 1
    
    return num_frames

def check_if_state_in_anchor(model, anchor, state_file):
    """
    Given a state file, make sure that it is in the assigned anchor.
    """
    curdir = os.getcwd()
    os.chdir(model.anchor_rootdir)
    tmp_path = tempfile.NamedTemporaryFile()
    output_filename = tmp_path.name
    sim_openmm = create_sim_openmm(
        model, anchor, output_filename, state_prefix=None, frame=0, 
        load_state_file=state_file)
    
    for milestone in anchor.milestones:
        cv = model.collective_variables[milestone.cv_index]
        result = cv.check_openmm_context_within_boundary(
            sim_openmm.simulation.context, milestone.variables, 
            verbose=True)
        if not result:
            os.chdir(curdir)
            return False
    os.chdir(curdir)
    return True
