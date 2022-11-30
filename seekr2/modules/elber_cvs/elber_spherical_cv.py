"""
Structures, methods, and functions for handling Spherical CVs.
(COM-COM distance CV). Using Elber milestoning.
"""

import numpy as np

from parmed import unit

import seekr2.modules.common_base as base
import seekr2.modules.elber_cvs.elber_cv_base as elber_cv_base
from seekr2.modules.elber_cvs.elber_cv_base import Elber_collective_variable

class Elber_spherical_CV(Elber_collective_variable):
    """
    A spherical collective variable represents the distance between two
    different groups of atoms.
    
    """+Elber_collective_variable.__doc__
    
    def __init__(self, index, groups):
        self.index = index
        self.group1 = groups[0]
        self.group2 = groups[1]
        self.name = "elber_spherical"
        self.openmm_umbrella_expression = "0.5*k*(distance(g1,g2)-radius)^2"
        self.openmm_fwd_rev_expression \
            = "step(k*(distance(g1, g2)^2 - radius^2))"
        self.num_groups = 2
        self.per_dof_variables = ["k", "radius"]
        self.global_variables = []
        self._mygroup_list = None
        self.variable_name = "r"
        return

    def __name__(self):
        return "Elber_spherical_CV"
    
    def make_umbrella_force_object(self):
        """
        Make an umbrella sampling force object, which will constrain
        the system to the milestone.
        """
        try:
            import openmm
        except ImportError:
            import simtk.openmm as openmm
            
        assert self.num_groups == 2
        return openmm.CustomCentroidBondForce(
            self.num_groups, self.openmm_umbrella_expression)
        
    def make_fwd_rev_force_object(self):
        """
        Make a list of reversal force objects, which will  be used to
        monitor milestone crossing during the reversal stage.
        """
        try:
            import openmm
        except ImportError:
            import simtk.openmm as openmm
            
        assert self.num_groups == 2
        return openmm.CustomCentroidBondForce(
            self.num_groups, self.openmm_fwd_rev_expression)
        
    def make_namd_colvar_umbrella_string(self):
        """
        This string will be put into a NAMD colvar file for applying
        an umbrella sampling force to constrain the system to the
        milestone.
        """
        serial_group1 = [str(index+1) for index in self.group1]
        serial_group2 = [str(index+1) for index in self.group2]
        serial_group1_str = " ".join(serial_group1)
        serial_group2_str = " ".join(serial_group2)
        namd_colvar_string = """
colvar {{
  name collective_variable_{0}
  outputappliedforce         off
  distance {{
    group1 {{ atomNumbers {1} }}
    group2 {{ atomNumbers {2} }}
  }}
}}
""".format(self.index, serial_group1_str, serial_group2_str)
        return namd_colvar_string
        
    def add_fwd_rev_parameters(self, force):
        """
        An OpenMM custom force object needs a list of variables
        provided to it that will occur within its expression. Both
        the per-dof and global variables are combined within the
        variable_names_list. The numerical values of these variables
        will be provided at a later step.
        """
        self._mygroup_list = []
        mygroup1 = force.addGroup(self.group1)
        self._mygroup_list.append(mygroup1)
        mygroup2 = force.addGroup(self.group2)
        self._mygroup_list.append(mygroup2)
        variable_names_list = []
        if self.per_dof_variables is not None:
            for per_dof_variable in self.per_dof_variables:
                force.addPerBondParameter(per_dof_variable)
                variable_names_list.append(per_dof_variable)
        
        if self.global_variables is not None:
            for global_variable in self.global_variables:
                force.addGlobalParameter(global_variable)
                variable_names_list.append(global_variable)
            
        return variable_names_list, self._mygroup_list
    
    def add_umbrella_parameters(self, force):
        """
        
        """
        #self._mygroup_list = []
        #for i in range(len(self.get_atom_groups())):
        #    self._mygroup_list.append(i)
        variable_names_list = []
        if self.per_dof_variables is not None:
            for per_dof_variable in self.per_dof_variables:
                force.addPerBondParameter(per_dof_variable)
                variable_names_list.append(per_dof_variable)
        
        if self.global_variables is not None:
            for global_variable in self.global_variables:
                force.addGlobalParameter(global_variable)
                variable_names_list.append(global_variable)
            
        return variable_names_list
    
    def add_groups_and_variables(self, force, group_list, variables):
        """
        Provide the custom force with additional information it needs,
        which includes a list of the groups of atoms involved with the
        CV, as well as a list of the variables' *values*.
        """
        assert len(group_list) == self.num_groups
        force.addBond(group_list, variables)
        return
    
    def get_variable_values_list(self, milestone):
        """
        Create the list of CV variables' values in the proper order
        so they can be provided to the custom force object.
        """
        assert milestone.cv_index == self.index
        values_list = []
        k = milestone.variables['k'] * unit.kilojoules_per_mole
        radius = milestone.variables['radius'] * unit.nanometers
        values_list.append(k)
        values_list.append(radius)
        
        return values_list
    
    def get_namd_fwd_rev_evaluation_string(self, milestone, cv_val_var="cv_val"):
        """
        For a given milestone, return a string that can be evaluated
        my NAMD to monitor for a crossing event. Essentially, if the 
        function defined by the string ever returns True, then a
        bounce will occur
        """
        assert milestone.cv_index == self.index
        k = milestone.variables['k']
        radius_in_nm = milestone.variables['radius'] * unit.nanometers
        radius_in_A = radius_in_nm.value_in_unit(unit.angstroms)
        eval_string = "{0} * (${1}_{2} - {3}) > 0".format(
            k, cv_val_var, self.index, radius_in_A)
        return eval_string
    
    def check_mdtraj_close_to_boundary(self, traj, milestone_variables, 
                                     verbose=False, max_avg=0.03, max_std=0.05):
        """
        
        """
        traj1 = traj.atom_slice(self.group1)
        traj2 = traj.atom_slice(self.group2)
        com1_array = elber_cv_base.traj_center_of_mass(traj1)
        com2_array = elber_cv_base.traj_center_of_mass(traj2)
        #com1_array = mdtraj.compute_center_of_mass(traj1)
        #com2_array = mdtraj.compute_center_of_mass(traj2)
        distances = []
        for frame_index in range(traj.n_frames):
            com1 = com1_array[frame_index,:]
            com2 = com2_array[frame_index,:]
            radius = np.linalg.norm(com2-com1)
            milestone_radius = milestone_variables["radius"]
            distances.append(radius - milestone_radius)
            
        avg_distance = np.mean(distances)
        std_distance = np.std(distances)
        if abs(avg_distance) > max_avg or std_distance > max_std:
            if verbose:
                warnstr = """The distance between the system and central 
    milestone were found on average to be {:.4f} nm apart.
    The standard deviation was {:.4f} nm.""".format(avg_distance, std_distance)
                print(warnstr)
            return False
            
        return True
    
    def get_atom_groups(self):
        """
        
        """
        return[self.group1, self.group2]

def make_elber_spherical_cv_object(spherical_cv_input, index):
    """
    Create a SphericalCV object to be placed into the Model.
    """
    groups = [spherical_cv_input.group1, spherical_cv_input.group2]
    cv = Elber_spherical_CV(index, groups)
    return cv
    
    
def make_elber_milestoning_objects_spherical(
        spherical_cv_input, milestone_alias, milestone_index, 
        input_anchor_index, anchor_index, input_anchors, umbrella_force_constant):
    """
    Make a set of 1-dimensional spherical Milestone objects to be put into
    an Anchor, and eventually, the Model. Each anchor needs 3 milestones.
    """
    milestones = []
    num_anchors = len(input_anchors)
    if input_anchor_index > 0:
        neighbor_index = anchor_index - 1
        neighbor_index_input = input_anchor_index - 1
        milestone1 = base.Milestone()
        milestone1.index = milestone_index-1
        milestone1.neighbor_anchor_index = neighbor_index
        #milestone1.alias_index = milestone_alias
        milestone1.alias_index = 1
        milestone1.is_source_milestone = False
        milestone1.cv_index = spherical_cv_input.index
        radius = input_anchors[neighbor_index_input].radius
        milestone1.variables = {"k": umbrella_force_constant, "radius": radius}
        #milestone_alias += 1
        #milestone_index += 1
        milestones.append(milestone1)
    
    milestone2 = base.Milestone()
    milestone2.index = milestone_index
    milestone2.neighbor_anchor_index = milestone_index
    #milestone2.alias_index = milestone_alias
    milestone2.alias_index = 2
    milestone2.is_source_milestone = True
    milestone2.cv_index = spherical_cv_input.index
    radius = input_anchors[input_anchor_index].radius
    milestone2.variables = {"k": umbrella_force_constant, "radius": radius}
    milestones.append(milestone2)
    
    if input_anchor_index < num_anchors-1:
        neighbor_index = anchor_index + 1
        neighbor_index_input = input_anchor_index + 1
        milestone3 = base.Milestone()
        milestone3.index = milestone_index+1
        milestone3.neighbor_anchor_index = neighbor_index
        #milestone3.alias_index = milestone_alias+1
        milestone3.alias_index = 3
        milestone3.is_source_milestone = False
        milestone3.cv_index = spherical_cv_input.index
        radius = input_anchors[neighbor_index_input].radius
        milestone3.variables = {"k": umbrella_force_constant, "radius": radius}
        milestones.append(milestone3)
    
    return milestones, milestone_alias, milestone_index+1
