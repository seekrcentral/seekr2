"""
Structures, methods, and functions for handling Tiwary CVs.
(linear combinations of distances, angles, dihedrals).
"""

import numpy as np

from parmed import unit

import seekr2.modules.common_base as base
import seekr2.modules.mmvt_cvs.mmvt_cv_base as mmvt_cv_base
from seekr2.modules.mmvt_cvs.mmvt_cv_base import MMVT_collective_variable

class MMVT_tiwary_CV(MMVT_collective_variable):
    """
    A Tiwary collective variable which is a linear function of order
    parameters.
    
    """+MMVT_collective_variable.__doc__
    
    def __init__(self, index, order_parameters, order_parameter_weights):
        self.index = index
        self.order_parameters = order_parameters
        for order_parameter in self.order_parameters:
            order_parameter.initialize_groups()
        self.order_parameter_weights = order_parameter_weights
        self.name = "mmvt_tiwary"
        self.num_groups = 0
        self.openmm_expression = None
        self.restraining_expression = None
        self.cv_expression = None
        self.assign_expressions_and_variables()
        self._mygroup_list = None
        self.variable_name = "v"
        return

    def __name__(self):
        return "MMVT_tiwary_CV"
    
    def assign_expressions_and_variables(self):
        openmm_expression_list = []
        self.per_dof_variables = []
        for i, (order_parameter, order_parameter_weight) in enumerate(zip(
                self.order_parameters, self.order_parameter_weights)):
            weight_var = "c{}".format(i)
            this_expr = weight_var + "*" \
                + order_parameter.get_openmm_expression(
                    group_index=self.num_groups+1)
            
            self.per_dof_variables.append(weight_var)
            openmm_expression_list.append(this_expr)
            self.num_groups += order_parameter.get_num_groups()
        self.openmm_expression = "step(k*("+" + ".join(openmm_expression_list)\
            + " - value))"
        self.restraining_expression = "0.5*k*(" \
            + " + ".join(openmm_expression_list) + " - value)^2"
        self.cv_expression = " + ".join(openmm_expression_list)
        self.per_dof_variables.append("k")
        self.per_dof_variables.append("value")
        self.global_variables = []
        return
    
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
            
        assert self.num_groups > 1
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
        
        assert self.num_groups > 1
        return openmm.CustomCentroidBondForce(
            self.num_groups, self.restraining_expression)
    
    def make_voronoi_cv_boundary_forces(self, me_val, neighbor_val, alias_id):
        """
        
        """
        try:
            import openmm
        except ImportError:
            import simtk.openmm as openmm
            
        # Handle me force
        me_expr = "(me_val - {})^2".format(self.cv_expression)
        me_force = openmm.CustomCentroidBondForce(
            self.num_groups, me_expr)
        me_group_list = []
        me_values = [me_val]
        me_force.addPerBondParameter("me_val")
        for i, order_parameter in enumerate(self.order_parameters):
            for j in range(order_parameter.get_num_groups()):
                group = order_parameter.get_group(j)
                mygroup = me_force.addGroup(group)
                me_group_list.append(mygroup)
            me_force.addPerBondParameter("c{}".format(i))
            me_values.append(self.order_parameter_weights[i])
        
        me_force.addBond(me_group_list, me_values)
        
        # Handle neighbor force
        neighbor_expr = "(neighbor_val - {})^2".format(self.cv_expression)
        neighbor_force = openmm.CustomCentroidBondForce(
            self.num_groups, neighbor_expr)
        neighbor_group_list = []
        neighbor_values = [neighbor_val]
        neighbor_force.addPerBondParameter("neighbor_val")
        for i, order_parameter in enumerate(self.order_parameters):
            for j in range(order_parameter.get_num_groups()):
                group = order_parameter.get_group(j)
                mygroup = neighbor_force.addGroup(group)
                neighbor_group_list.append(mygroup)
            neighbor_force.addPerBondParameter("c{}".format(i))
            neighbor_values.append(self.order_parameter_weights[i])
        
        neighbor_force.addBond(neighbor_group_list, neighbor_values)
        
        return me_force, neighbor_force
    
    def update_voronoi_cv_boundary_forces(
            self, me_force, me_val, neighbor_force, neighbor_val, alias_id,
            context):
        """
        
        """
        me_values = [me_val]
        neighbor_values = [neighbor_val]
        me_group_list = []
        neighbor_group_list = []
        counter = 0
        for i, order_parameter in enumerate(self.order_parameters):
            for j in range(order_parameter.get_num_groups()):
                me_group_list.append(counter)
                neighbor_group_list.append(counter)
                counter += 1
            
            me_values.append(self.order_parameter_weights[i])
            neighbor_values.append(self.order_parameter_weights[i])
        
        me_force.setBondParameters(0, me_group_list, me_values)
        neighbor_force.setBondParameters(0, neighbor_group_list, neighbor_values)
        return
    
    def make_namd_colvar_string(self):
        """
        This string will be put into a NAMD colvar file for tracking
        MMVT bounces.
        """
        raise Exception("MMVT Tiwary CVs are not available in NAMD")
    
    def add_groups(self, force):
        """
        
        """
        self._mygroup_list = []
        for order_parameter in self.order_parameters:
            for i in range(order_parameter.get_num_groups()):
                group = order_parameter.get_group(i)
                mygroup = force.addGroup(group)
                self._mygroup_list.append(mygroup)
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
        force.addBond(self._mygroup_list, variables)
        return
    
    def update_groups_and_variables(self, force, variables, alias_id, context):
        """
        Update the force's variables with a list of new values.
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
        k = milestone.variables['k'] * unit.kilojoules_per_mole/unit.angstrom**2
        radius = milestone.variables['value'] * unit.nanometers
        values_list.append(bitcode)
        for i, order_parameter in enumerate(self.order_parameters):
            key = "c{}".format(i)
            c = milestone.variables[key]
            values_list.append(c)
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
        raise Exception("MMVT Tiwary CVs are not available in NAMD")
    
    def get_mdtraj_cv_value(self, traj, frame_index):
        """
        Determine the current CV value for an mdtraj object.
        """
        op_com_array_list = []
        for order_parameter in self.order_parameters:
            com_array_list = []
            for j in range(order_parameter.get_num_groups()):
                group = order_parameter.get_group(j)
                traj_group = traj.atom_slice(group)
                #com_array = mdtraj.compute_center_of_mass(traj_group)
                com_array = mmvt_cv_base.traj_center_of_mass(traj_group)
                com_array_list.append(com_array)
            op_com_array_list.append(com_array_list)
                
        op_value = 0.0
        for i, order_parameter in enumerate(self.order_parameters):
            com_list = []
            for j in range(order_parameter.get_num_groups()):
                com = op_com_array_list[i][j][frame_index,:]
                com_list.append(com)
            op_term = order_parameter.get_value(com_list)
            op_weight = self.order_parameter_weights[i]
            op_value += op_weight * op_term
        return op_value
    
    def check_mdtraj_within_boundary(self, traj, milestone_variables, 
                                     verbose=False, TOL=0.001):
        """
        Check if an mdtraj Trajectory describes a system that remains
        within the expected anchor. Return True if passed, return
        False if failed.
        """
        for frame_index in range(traj.n_frames):
            op_value = self.get_mdtraj_cv_value(traj, frame_index)
            result = self.check_value_within_boundary(
                op_value, milestone_variables, verbose=verbose, tolerance=TOL)
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
        op_com_list = []
        for order_parameter in self.order_parameters:
            com_list = []
            for j in range(order_parameter.get_num_groups()):
                group = order_parameter.get_group(j)
                com = base.get_openmm_center_of_mass_com(
                    system, positions, group)
                com_list.append(com)
            op_com_list.append(com_list)
                
        op_value = 0.0
        for i, order_parameter in enumerate(self.order_parameters):
            com_list = []
            for j in range(order_parameter.get_num_groups()):
                com = op_com_list[i][j]
                com_list.append(com.value_in_unit(unit.nanometer))
            op_term = order_parameter.get_value(com_list)
            op_weight = self.order_parameter_weights[i]
            op_value += op_weight * op_term
        return op_value
    
    def check_openmm_context_within_boundary(
            self, context, milestone_variables, positions=None, verbose=False,
            tolerance=0.0):
        """
        Check if an mdtraj Trajectory describes a system that remains
        within the expected anchor. Return True if passed, return
        False if failed.
        """
        op_value = self.get_openmm_context_cv_value(context, positions)
        result = self.check_value_within_boundary(op_value, milestone_variables, 
                                    verbose=verbose, tolerance=tolerance)
        return result
    
    def check_value_within_boundary(self, op_value, milestone_variables, 
                                    verbose=False, tolerance=0.0):
        """
        
        """
        milestone_k = milestone_variables["k"]
        milestone_value = milestone_variables["value"]
        if milestone_k*(op_value - milestone_value) > tolerance:
            if verbose:
                warnstr = """This system value for the Tiwary CV is {:.4f} 
but the milestone's value for the same CV is {:.4f}. The system falls outside 
the boundaries.""".format(op_value, milestone_value)
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
        print("check_mdtraj_close_to_boundary not yet implemented for "\
              "Tiwary CVs")
        return True
    
    def get_atom_groups(self):
        """
        Return a list of this CV's atomic groups.
        """
        group_list = []
        for order_parameter in self.order_parameters:
            for i in range(order_parameter.get_num_groups()):
                group = order_parameter.get_group(i)
                group_list.append(group)
            
        return group_list
    
    def get_variable_values(self):
        """
        Return the order parameter weights as the extra CV variables.
        """
        return self.order_parameter_weights

def make_mmvt_tiwary_cv_object(tiwary_cv_input, index):
    """
    Create a Tiwary CV object to be placed into the Model.
    """
    
    cv = MMVT_tiwary_CV(
        index, tiwary_cv_input.order_parameters, 
        tiwary_cv_input.order_parameter_weights)
    return cv
