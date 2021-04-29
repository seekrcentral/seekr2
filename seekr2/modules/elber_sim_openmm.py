"""
elber/sim_openmm.py

TODO: description
"""

import os
import re

import openmm
import openmm.app as openmm_app
from parmed import unit

import seekr2.modules.common_base as base
import seekr2.modules.elber_base as elber_base
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
        self.rev_seekr_force = None
        
        self.fwd_system = None
        self.fwd_integrator = None
        self.fwd_simulation = None
        self.fwd_traj_reporter = openmm_app.DCDReporter
        self.fwd_energy_reporter = openmm_app.StateDataReporter
        self.fwd_output_filename = ""
        self.fwd_seekr_force = None
        return
    
class Elber_sim_openmm_factory(common_sim_openmm.Common_sim_openmm_factory):
    """
    Create the sim_openmm objects which will be used to run MMVT
    calculations in OpenMM.
    """
    def __init__(self):
        self.sim_openmm = Elber_sim_openmm()
        return
    
    def add_integrators(self, model, state_prefix=None):
        """
        
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
                
            self.sim_openmm.timestep = timestep
            
            self.sim_openmm.umbrella_integrator = openmm.LangevinIntegrator(
                target_temperature*unit.kelvin, 
                friction_coefficient/unit.picoseconds, 
                timestep*unit.picoseconds)
            if random_seed is not None:
                self.sim_openmm.umbrella_integrator.setRandomNumberSeed(random_seed)
                
            self.sim_openmm.rev_integrator = openmm.VerletIntegrator(
                timestep*unit.picoseconds)
            self.sim_openmm.fwd_integrator = openmm.VerletIntegrator(
                timestep*unit.picoseconds)
            
            if rigid_constraint_tolerance is not None:
                self.sim_openmm.umbrella_integrator.setConstraintTolerance(
                    rigid_constraint_tolerance)
                self.sim_openmm.rev_integrator.setConstraintTolerance(
                    rigid_constraint_tolerance)
                self.sim_openmm.fwd_integrator.setConstraintTolerance(
                    rigid_constraint_tolerance)
            
            #if state_prefix is not None:
            #    self.sim_openmm.integrator.setSaveStateFileName(state_prefix)
            
        else:
            raise Exception("Settings not provided for available "\
                            "integrator type(s).")
        return
    
    def add_forces(self, model, anchor, rev_data_file_name, fwd_data_file_name):
        """
        
        """
        self.sim_openmm.umbrella_force = make_elber_umbrella_force(
            model, anchor)
        self.sim_openmm.umbrella_system.addForce(self.sim_openmm.umbrella_force)
        self.sim_openmm.rev_seekr_force = make_elber_rev_force(
            model, anchor, rev_data_file_name)
        self.sim_openmm.rev_system.addForce(self.sim_openmm.rev_seekr_force)
        self.sim_openmm.fwd_seekr_force = make_elber_fwd_force(
            model, anchor, fwd_data_file_name)
        self.sim_openmm.fwd_system.addForce(self.sim_openmm.fwd_seekr_force)
        return
    
    def add_simulations(self, model, topology, positions, box_vectors):
        """
        
        """
        self.sim_openmm.umbrella_simulation = openmm_app.Simulation(
            topology.topology, self.sim_openmm.umbrella_system, 
            self.sim_openmm.umbrella_integrator, self.sim_openmm.platform, 
            self.sim_openmm.properties)
        
        #for force in self.sim_openmm.rev_system.getForces():
        #    print("force name:", force.__class__.__name__,
        #          "force.usesPeriodicBoundaryConditions():", 
        #          force.usesPeriodicBoundaryConditions())
        #exit()
        self.sim_openmm.rev_simulation = openmm_app.Simulation(
            topology.topology, self.sim_openmm.rev_system, 
            self.sim_openmm.rev_integrator, self.sim_openmm.platform, 
            self.sim_openmm.properties)
        self.sim_openmm.fwd_simulation = openmm_app.Simulation(
            topology.topology, self.sim_openmm.fwd_system, 
            self.sim_openmm.fwd_integrator, self.sim_openmm.platform, 
            self.sim_openmm.properties)
        
        self.sim_openmm.umbrella_simulation.context.setPositions(positions.positions)
        if box_vectors is not None:
            self.sim_openmm.umbrella_simulation.context.setPeriodicBoxVectors(
                *box_vectors)
        if model.openmm_settings.run_minimization:
            self.sim_openmm.umbrella_simulation.minimizeEnergy()
        
        self.sim_openmm.umbrella_simulation.context.setVelocitiesToTemperature(
            model.openmm_settings.initial_temperature * unit.kelvin)
        assert self.sim_openmm.timestep is not None
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
        fwd_output_filename = output_filename
        directory = os.path.dirname(output_filename)
        basename = os.path.basename(output_filename)
        suffix = re.sub(elber_base.OPENMM_ELBER_BASENAME, "", basename)
        rev_basename = elber_base.ELBER_REV_BASENAME+suffix
        rev_output_filename = os.path.join(directory, rev_basename)
        super(Elber_sim_openmm_factory, self).fill_generic_parameters(
            self.sim_openmm, model, anchor, output_filename)
        umbrella_system, umbrella_topology, umbrella_positions, \
            umbrella_box_vectors = super(
            Elber_sim_openmm_factory, self).create_openmm_system(
                self.sim_openmm, model, anchor)
        rev_system, rev_topology, rev_positions, rev_box_vectors = super(
            Elber_sim_openmm_factory, self).create_openmm_system(
                self.sim_openmm, model, anchor)
        fwd_system, fwd_topology, fwd_positions, fwd_box_vectors = super(
            Elber_sim_openmm_factory, self).create_openmm_system(
                self.sim_openmm, model, anchor)
        self.sim_openmm.umbrella_system = umbrella_system
        self.sim_openmm.rev_system = rev_system
        self.sim_openmm.fwd_system = fwd_system
        self.add_integrators(model, state_prefix=state_prefix)
        super(Elber_sim_openmm_factory, self).add_barostat(
            self.sim_openmm, model)
        super(Elber_sim_openmm_factory, self).add_platform(
            self.sim_openmm, model)
        self.add_forces(model, anchor, rev_output_filename, fwd_output_filename)
        self.sim_openmm.rev_output_filename = rev_output_filename
        self.sim_openmm.fwd_output_filename = fwd_output_filename
        self.add_simulations(model, umbrella_topology, umbrella_positions, 
                             umbrella_box_vectors)
        return self.sim_openmm

def make_elber_umbrella_force(model, anchor):
    """
    
    """
    for milestone in anchor.milestones:
        if milestone.alias_index == 2:
            cv = milestone.get_CV(model)
            umbrella_force = cv.make_umbrella_force_object()
            break
        
    mygroup_list = []
    mygroup1 = umbrella_force.addGroup(cv.group1)
    mygroup_list.append(mygroup1)
    mygroup2 = umbrella_force.addGroup(cv.group2)
    mygroup_list.append(mygroup2)    
    
    variable_names_list = cv.add_umbrella_parameters(umbrella_force)
    cv.add_groups_and_variables(umbrella_force, mygroup_list, 
                                cv.get_variable_values_list(
                                    milestone))
    return umbrella_force

def make_elber_rev_force(model, anchor, data_file_name):
    """
    
    """
    
    for milestone in anchor.milestones:
        if milestone.alias_index == 2:
            cv = milestone.get_CV(model)
            rev_force = cv.make_fwd_rev_force_object(anchor)
    variable_names_list = cv.add_fwd_rev_parameters(rev_force, anchor, 
                                            end_on_middle_crossing=True,
                                            data_file_name=data_file_name)
    return rev_force

def make_elber_fwd_force(model, anchor, data_file_name):
    """
    
    """
    
    for milestone in anchor.milestones:
        if milestone.alias_index == 2:
            cv = milestone.get_CV(model)
            fwd_force = cv.make_fwd_rev_force_object(anchor)
    variable_names_list = cv.add_fwd_rev_parameters(fwd_force, anchor, 
                                            end_on_middle_crossing=False,
                                            data_file_name=data_file_name)
    return fwd_force
        
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