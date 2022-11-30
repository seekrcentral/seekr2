"""
Structures, methods, and functions for handling Spherical CVs.
(COM-COM distance CV).
"""

import numpy as np

from parmed import unit

import seekr2.modules.common_base as base
import seekr2.modules.mmvt_cvs.mmvt_cv_base as mmvt_cv_base
from seekr2.modules.mmvt_cvs.mmvt_cv_base import MMVT_collective_variable

class MMVT_spherical_CV(MMVT_collective_variable):
    """
    A spherical collective variable represents the distance between two
    different groups of atoms.
    
    """+MMVT_collective_variable.__doc__
    
    def __init__(self, index, groups):
        self.index = index
        self.group1 = groups[0]
        self.group2 = groups[1]
        self.name = "mmvt_spherical"
        self.openmm_expression = "step(k*(distance(g1, g2)^2 - radius^2))"
        self.restraining_expression = "0.5*k*(distance(g1, g2) - radius)^2"
        self.cv_expression = "distance(g1, g2)"
        self.num_groups = 2
        self.per_dof_variables = ["k", "radius"]
        self.global_variables = []
        self._mygroup_list = None
        self.variable_name = "r"
        return

    def __name__(self):
        return "MMVT_spherical_CV"
    
    def make_boundary_force(self, alias_id):
        """
        Create an OpenMM force object which will be used to compute
        the value of the CV's mathematical expression as the simulation
        propagates. Each *milestone* in a cell, not each CV, will have
        a force object created.
        
        In this implementation of MMVT in OpenMM, CustomForce objects
        are used to provide the boundary definitions of the MMVT cell.
        These 'forces' are designed to never affect atomic motions,
        but are merely monitored by the plugin for bouncing event.
        So while they are technically OpenMM force objects, we don't
        refer to them as forces outside of this layer of the code,
        preferring instead the term: boundary definitions.
        """
        try:
            import openmm
        except ImportError:
            import simtk.openmm as openmm
            
        assert self.num_groups == 2
        expression_w_bitcode = "bitcode*"+self.openmm_expression
        return openmm.CustomCentroidBondForce(
            self.num_groups, expression_w_bitcode)
    
    def make_restraining_force(self, alias_id):
        """
        Create an OpenMM force object that will restrain the system to
        a given value of this CV.
        """
        try:
            import openmm
        except ImportError:
            import simtk.openmm as openmm
            
        assert self.num_groups == 2
        return openmm.CustomCentroidBondForce(
            self.num_groups, self.restraining_expression)
    
    def make_voronoi_cv_boundary_forces(self, me_val, neighbor_val, alias_id):
        """
        
        """
        try:
            import openmm
        except ImportError:
            import simtk.openmm as openmm
            
        me_expr = "(me_val - {})^2".format(self.cv_expression)
        me_force = openmm.CustomCentroidBondForce(
            self.num_groups, me_expr)
        megroup1 = me_force.addGroup(self.group1)
        megroup2 = me_force.addGroup(self.group2)
        me_force.addPerBondParameter("me_val")
        me_force.addBond([megroup1, megroup2], [me_val])
        
        neighbor_expr = "(neighbor_val - {})^2".format(self.cv_expression)
        neighbor_force = openmm.CustomCentroidBondForce(
            self.num_groups, neighbor_expr)
        neighborgroup1 = neighbor_force.addGroup(self.group1)
        neighborgroup2 = neighbor_force.addGroup(self.group2)
        neighbor_force.addPerBondParameter("neighbor_val")
        neighbor_force.addBond([neighborgroup1, neighborgroup2], [neighbor_val])
        
        return me_force, neighbor_force
    
    def update_voronoi_cv_boundary_forces(
            self, me_force, me_val, neighbor_force, neighbor_val, alias_id, 
            context):
        """
        
        """
        # the group list values of [0,1] is hacky, but probably correct
        me_force.setBondParameters(0, [0,1], [me_val])
        neighbor_force.setBondParameters(0, [0,1], [neighbor_val])
        return
    
    def make_namd_colvar_string(self):
        """
        This string will be put into a NAMD colvar file for tracking
        MMVT bounces.
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
    
    def add_groups(self, force):
        """
        
        """
        self._mygroup_list = []
        mygroup1 = force.addGroup(self.group1)
        self._mygroup_list.append(mygroup1)
        mygroup2 = force.addGroup(self.group2)
        self._mygroup_list.append(mygroup2)
        return
    
    def add_parameters(self, force):
        """
        An OpenMM custom force object needs a list of variables
        provided to it that will occur within its expression. Both
        the per-dof and global variables are combined within the
        variable_names_list. The numerical values of these variables
        will be provided at a later step.
        """
        self.add_groups(force)
        variable_names_list = []
        variable_names_list.append("bitcode")
        force.addPerBondParameter("bitcode")
        if self.per_dof_variables is not None:
            for per_dof_variable in self.per_dof_variables:
                force.addPerBondParameter(per_dof_variable)
                variable_names_list.append(per_dof_variable)
        
        if self.global_variables is not None:
            for global_variable in self.global_variables:
                force.addGlobalParameter(global_variable)
                variable_names_list.append(global_variable)
        
        return variable_names_list
    
    def add_groups_and_variables(self, force, variables, alias_id):
        """
        Provide the custom force with additional information it needs,
        which includes a list of the groups of atoms involved with the
        CV, as well as a list of the variables' *values*.
        """
        assert len(self._mygroup_list) == self.num_groups
        # TODO: add this assert to other CVs
        assert force.getNumBonds() == 0, \
            "add_groups_and_variables cannot be called twice for the same "\
            "CV. Use update_groups_and_variables."
        force.addBond(self._mygroup_list, variables)
        return
    
    def update_groups_and_variables(self, force, variables, alias_id):
        """
        Update the force's variables with a list of new values.
        
        A context must have it's reinitialize() method called after 
        this function for the changes to be applied.
        """
        force.setBondParameters(0, self._mygroup_list, variables)
        return
    
    def get_variable_values_list(self, milestone):
        """
        Create the list of CV variables' values in the proper order
        so they can be provided to the custom force object.
        """
        assert milestone.cv_index == self.index
        values_list = []
        bitcode = 2**(milestone.alias_index-1)
        k = milestone.variables["k"] * unit.kilojoules_per_mole/unit.angstrom**2
        radius = milestone.variables["radius"] * unit.nanometers
        values_list.append(bitcode)
        values_list.append(k)
        values_list.append(radius)
        return values_list
    
    def get_namd_evaluation_string(self, milestone, cv_val_var="cv_val"):
        """
        For a given milestone, return a string that can be evaluated
        my NAMD to monitor for a crossing event. Essentially, if the 
        function defined by the string ever returns True, then a
        bounce will occur
        """
        assert milestone.cv_index == self.index
        k = milestone.variables["k"]
        radius_in_nm = milestone.variables["radius"] * unit.nanometers
        radius_in_A = radius_in_nm.value_in_unit(unit.angstroms)
        eval_string = "{0} * (${1}_{2} - {3}) > 0".format(
            k, cv_val_var, self.index, radius_in_A)
        return eval_string
    
    def get_mdtraj_cv_value(self, traj, frame_index):
        """
        Determine the current CV value for an mdtraj object.
        """
        traj1 = traj.atom_slice(self.group1)
        traj2 = traj.atom_slice(self.group2)
        com1_array = mmvt_cv_base.traj_center_of_mass(traj1)
        com2_array = mmvt_cv_base.traj_center_of_mass(traj2)
        
        #if traj1.xyz.shape[1] == 1:
        #    com1_array = traj1.xyz
        #else:
        #    com1_array = mdtraj.compute_center_of_mass(traj1)
        
        #if traj2.xyz.shape[1] == 1:
        #    com2_array = traj2.xyz
        #else:
        #    com2_array = mdtraj.compute_center_of_mass(traj2)
            
        com1 = com1_array[frame_index,:]
        com2 = com2_array[frame_index,:]
        radius = np.linalg.norm(com2-com1)
        return radius
    
    def check_mdtraj_within_boundary(self, traj, milestone_variables, 
                                     verbose=False, TOL=0.001):
        """
        Check if an mdtraj Trajectory describes a system that remains
        within the expected anchor. Return True if passed, return
        False if failed.
        """
        for frame_index in range(traj.n_frames):
            radius = self.get_mdtraj_cv_value(traj, frame_index)
            result = self.check_value_within_boundary(
                radius, milestone_variables, verbose)
            if not result:
                return False
            
        return True
        
    def get_openmm_context_cv_value(self, context, positions=None, system=None):
        """
        Determine the current CV value for an openmm context.
        """
        if system is None:
            system = context.getSystem()
        if positions is None:
            state = context.getState(getPositions=True)
            positions = state.getPositions()
        com1 = base.get_openmm_center_of_mass_com(
            system, positions, self.group1)
        com2 = base.get_openmm_center_of_mass_com(
            system, positions, self.group2)
        radius = np.linalg.norm(com2-com1)
        return radius
        
    def check_openmm_context_within_boundary(
            self, context, milestone_variables, positions=None, verbose=False,
            tolerance=0.0):
        """
        Check if an OpenMM context describes a system that remains
        within the expected anchor. Return True if passed, return
        False if failed.
        """
        radius = self.get_openmm_context_cv_value(context, positions)
        result = self.check_value_within_boundary(
            radius, milestone_variables, verbose, tolerance)
        return result
        
    def check_value_within_boundary(self, radius, milestone_variables, 
                                    verbose=False, tolerance=0.0):
        """
        Check if the given CV value is within the milestone's boundaries.
        """
        milestone_k = milestone_variables["k"]
        milestone_radius = milestone_variables["radius"]
        if milestone_k*(radius - milestone_radius) > tolerance:
            if verbose:
                warnstr = """The center of masses of atom group1 and atom
group2 were found to be {:.4f} nm apart.
This distance falls outside of the milestone
boundary at {:.4f} nm.""".format(radius, milestone_radius)
                print(warnstr)
            return False
        return True
    
    def check_mdtraj_close_to_boundary(self, traj, milestone_variables, 
                                     verbose=False, max_avg=0.03, max_std=0.05):
        """
        Given an mdtraj Trajectory, check if the system lies close
        to the MMVT boundary. Return True if passed, return False if 
        failed.
        """
        traj1 = traj.atom_slice(self.group1)
        traj2 = traj.atom_slice(self.group2)
        
        com1_array = mmvt_cv_base.traj_center_of_mass(traj1)
        com2_array = mmvt_cv_base.traj_center_of_mass(traj2)
        #if traj1.xyz.shape[1] == 1:
        #    com1_array = traj1.xyz
        #else:
        #    com1_array = mdtraj.compute_center_of_mass(traj1)
        
        #if traj2.xyz.shape[1] == 1:
        #    com2_array = traj2.xyz
        #else:
        #    com2_array = mdtraj.compute_center_of_mass(traj2)
        
        distances = []
        for frame_index in range(traj.n_frames):
            com1 = com1_array[frame_index,:]
            com2 = com2_array[frame_index,:]
            radius = np.linalg.norm(com2-com1)
            milestone_radius = milestone_variables["radius"]
            distances.append(radius - milestone_radius)
            
        assert np.isfinite(distances).all(), "Non-finite numbers detected in \
            'distances'."
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
        Return a 2-list of this CV's atomic groups.
        """
        return [self.group1, self.group2]
            
    def get_variable_values(self):
        """
        This type of CV has no extra variables, so an empty list is 
        returned.
        """
        return []


def make_mmvt_spherical_cv_object(spherical_cv_input, index):
    """
    Create a SphericalCV object to be placed into the Model.
    """
    group1 = base.parse_xml_list(spherical_cv_input.group1)
    group2 = base.parse_xml_list(spherical_cv_input.group2)
    groups = [group1, group2]
    cv = MMVT_spherical_CV(index, groups)
    return cv