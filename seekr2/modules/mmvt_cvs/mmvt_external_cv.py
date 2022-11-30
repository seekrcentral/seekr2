"""
Structures, methods, and functions for handling External CVs.
(Especially used for toy systems).
"""

import numpy as np
import math

from parmed import unit

import seekr2.modules.common_base as base
import seekr2.modules.mmvt_cvs.mmvt_cv_base as mmvt_cv_base
from seekr2.modules.mmvt_cvs.mmvt_cv_base import MMVT_collective_variable

class MMVT_external_CV(MMVT_collective_variable):
    """
    A collective variable that depends on external coordinates.
    
    """+MMVT_collective_variable.__doc__
    
    def __init__(self, index, groups):
        self.index = index
        self.groups = groups
        self.name = "mmvt_external"
        self.cv_expression = None
        self.openmm_expression = None
        self.restraining_expression = None
        self.cv_expression = None
        self.num_groups = len(groups)
        self.per_dof_variables = ["k", "value"]
        self.global_variables = []
        self._mygroup_list = None
        self.variable_name = "v"
        return

    def __name__(self):
        return "MMVT_external_CV"
    
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
        me_force = openmm.CustomCentroidBondForce(self.num_groups, me_expr)
        me_group_list = []
        for group in self.groups:
            mygroup = me_force.addGroup(group)
            me_group_list.append(mygroup)
        me_force.addPerBondParameter("me_val")
        me_force.addBond(me_group_list, [me_val])
        
        neighbor_expr = "(neighbor_val - {})^2".format(self.cv_expression)
        neighbor_force = openmm.CustomCentroidBondForce(
            self.num_groups, neighbor_expr)
        neighbor_group_list = []
        for group in self.groups:
            mygroup = neighbor_force.addGroup(group)
            neighbor_group_list.append(mygroup)
        neighbor_force.addPerBondParameter("neighbor_val")
        neighbor_force.addBond(neighbor_group_list, [neighbor_val])
        
        return me_force, neighbor_force
    
    def update_voronoi_cv_boundary_forces(
            self, me_force, me_val, neighbor_force, neighbor_val, alias_id, 
            context):
        """
        
        """
        group_list = []
        for i in range(len(self.groups)):
            group_list.append(i)
        me_force.setBondParameters(0, group_list, [me_val])
        neighbor_force.setBondParameters(0, group_list, [neighbor_val])
        return
    
    def make_namd_colvar_string(self):
        """
        This string will be put into a NAMD colvar file for tracking
        MMVT bounces.
        """
        raise Exception("MMVT External CVs are not available in NAMD")
    
    def add_groups(self, force):
        """
        
        """
        self._mygroup_list = []
        for group in self.groups:
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
    
    def update_groups_and_variables(self, force, variables, alias_id):
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
        values_list.append(bitcode)
        k_val = milestone.variables["k"]
        values_list.append(k_val)
        value_val = milestone.variables["value"]
        values_list.append(value_val)
        return values_list
    
    def get_namd_evaluation_string(self, milestone, cv_val_var="cv_val"):
        """
        For a given milestone, return a string that can be evaluated
        my NAMD to monitor for a crossing event. Essentially, if the 
        function defined by the string ever returns True, then a
        bounce will occur
        """
        raise Exception("MMVT External CVs are not available in NAMD")
    
    def get_mdtraj_cv_value(self, traj, frame_index):
        """
        Determine the current CV value for an mdtraj object.
        """
        positions = traj.xyz[frame_index, :,:] * unit.nanometers
        value = self.get_cv_value(positions)
        return value
    
    def check_mdtraj_within_boundary(self, traj, milestone_variables, 
                                     verbose=False, TOL=0.001):
        """
        For now, this will just always return True.
        """
        for frame_index in range(traj.n_frames):
            value = self.get_mdtraj_cv_value(traj, frame_index)
            result = self.check_value_within_boundary(
                value, milestone_variables, tolerance=TOL)
            if not result:
                return False
            
        return True
    
    def get_openmm_context_cv_value(self, context, positions=None, system=None):
        """
        Determine the current CV value for an openmm context.
        """
        if positions is None:
            state = context.getState(getPositions=True)
            positions = state.getPositions()
        
        value = self.get_cv_value(positions)
        return value
    
    def check_openmm_context_within_boundary(
            self, context, milestone_variables, positions=None, verbose=False, 
            tolerance=0.0):
        """
        For now, this will just always return True.
        """
        
        if positions is None:
            state = context.getState(getPositions=True)
            positions = state.getPositions()
        print("positions:", positions)
        return self.check_positions_within_boundary(
            positions.value_in_unit(unit.nanometers), milestone_variables, 
            tolerance)
    
    def get_cv_value(self, positions):
        """
        Get the value of the cv for the set of positions.
        """
        try:
            import openmm
        except ImportError:
            import simtk.openmm as openmm
            
        sqrt = math.sqrt
        exp = math.exp
        log = math.log
        sin = math.sin
        cos = math.cos
        tan = math.tan
        asin = math.asin
        acos = math.acos
        atan = math.atan
        sinh = math.sinh
        cosh = math.cosh
        tanh = math.tanh
        erf = math.erf
        erfc = math.erfc
        floor = math.floor
        ceil = math.ceil
        step = lambda x : 0 if x < 0 else 1
        delta = lambda x : 1 if x == 0 else 0
        select = lambda x, y, z : z if x == 0 else y
        expr = ""
        for i, position in enumerate(positions):
            expr_x = "x{} = {};".format(i+1, position[0].value_in_unit(
                openmm.unit.nanometer))
            expr_y = "y{} = {};".format(i+1, position[1].value_in_unit(
                openmm.unit.nanometer))
            expr_z = "z{} = {};".format(i+1, position[2].value_in_unit(
                openmm.unit.nanometer))
            expr += expr_x + expr_y + expr_z
            
        expr += base.convert_openmm_to_python_expr("result="+self.cv_expression)
        mylocals = locals()
        exec(expr, globals(), mylocals)
        result = mylocals["result"]
        return result
    
    def check_value_within_boundary(
            self, value, milestone_variables, tolerance=0.0):
        step = lambda x : 0 if x < 0 else 1
        result = step(milestone_variables["k"] \
                      * (value-milestone_variables["value"]) + tolerance)
        if result <= 0:
            return True
        else:
            return False
    
    def check_positions_within_boundary(self, positions, milestone_variables, 
                                    verbose=False, tolerance=0.0):
        """
        
        """
        try:
            import openmm
        except ImportError:
            import simtk.openmm as openmm
            
        new_positions = []
        for position in positions:
            new_positions.append(position*openmm.unit.nanometer)
        value = self.get_cv_value(new_positions)
        result = self.check_value_within_boundary(
            value, milestone_variables)
        return result
    
    def check_mdtraj_close_to_boundary(self, traj, milestone_variables, 
                                     verbose=False, max_avg=0.03, max_std=0.05):
        """
        For now, this will just always return True.
        """
        distances = []
        for frame_index in range(traj.n_frames):
            value = self.get_mdtraj_cv_value(traj, frame_index)
            milestone_value = milestone_variables["value"]
            distances.append(value - milestone_value)
            
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
        Return a list of this CV's atomic groups.
        """
        return self.groups
            
    def get_variable_values(self):
        """
        This type of CV has no extra variables, so an empty list is 
        returned.
        """
        return []

def make_mmvt_external_cv_object(external_cv_input, index):
    """
    Create a externalCV object to be placed into the Model.
    """
    
    groups = []
    for group in external_cv_input.groups:
        groups.append(base.parse_xml_list(group))
    
    cv = MMVT_external_CV(index, groups)
    cv.cv_expression = external_cv_input.cv_expression
    assert cv.cv_expression is not None, \
        "A CV expression is required to define MMVT milestones."
    if external_cv_input.openmm_expression is None:
        cv.openmm_expression = "step(k*("+cv.cv_expression+" - value))"
    else:
        cv.openmm_expression = external_cv_input.openmm_expression
    
    if external_cv_input.restraining_expression is None:
        cv.restraining_expression = "0.5*k*("+cv.cv_expression+" - value)^2"
    else:
        cv.restraining_expression = external_cv_input.restraining_expression
    return cv