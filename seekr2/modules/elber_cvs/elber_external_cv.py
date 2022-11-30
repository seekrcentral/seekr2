"""
Structures, methods, and functions for handling External CVs.
(toy CV). Using Elber milestoning.
"""

import numpy as np

from parmed import unit

import seekr2.modules.common_base as base
import seekr2.modules.elber_cvs.elber_cv_base as elber_cv_base
from seekr2.modules.elber_cvs.elber_cv_base import Elber_collective_variable

class Elber_external_CV(Elber_collective_variable):
    """
    A collective variable that depends on external coordinates.
    
    """+Elber_collective_variable.__doc__
    
    def __init__(self, index, groups):
        self.index = index
        self.groups = groups
        self.name = "elber_external"
        self.cv_expression = None
        self.openmm_umbrella_expression = None
        self.openmm_fwd_rev_expression = None
        self.num_groups = len(groups)
        self.per_dof_variables = ["k", "value"]
        self.global_variables = []
        self._mygroup_list = None
        self.variable_name = "v"
        return

    def __name__(self):
        return "Elber_external_CV"
    
    def make_umbrella_force_object(self):
        """
        Make an umbrella sampling force object, which will constrain
        the system to the milestone.
        """
        try:
            import openmm
        except ImportError:
            import simtk.openmm as openmm
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
            
        return openmm.CustomCentroidBondForce(
            self.num_groups, self.openmm_fwd_rev_expression)
        
    def make_namd_colvar_string(self):
        """
        This string will be put into a NAMD colvar file for tracking
        MMVT bounces.
        """
        raise Exception("Elber External CVs are not available in NAMD")
        
    def add_fwd_rev_parameters(self, force):
        """
        An OpenMM custom force object needs a list of variables
        provided to it that will occur within its expression. Both
        the per-dof and global variables are combined within the
        variable_names_list. The numerical values of these variables
        will be provided at a later step.
        """
        self._mygroup_list = []
        for group in self.groups:
            mygroup = force.addGroup(group)
            self._mygroup_list.append(mygroup)
                
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
        #bitcode = 2**(milestone.alias_index-1)
        #values_list.append(bitcode)
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
        raise Exception("Elber External CVs are not available in NAMD")
    
    def check_mdtraj_within_boundary(self, traj, milestone_variables, 
                                     verbose=False):
        """
        For now, this will just always return True.
        """
        return True
    
    def check_openmm_context_within_boundary(
            self, context, milestone_variables, positions=None, verbose=False):
        """
        For now, this will just always return True.
        """
        
        system = context.getSystem()
        if positions is None:
            state = context.getState(getPositions=True)
            positions = state.getPositions()
        
        return self.check_positions_within_boundary(
            positions, milestone_variables)
    
    def check_positions_within_boundary(
            self, positions, milestone_variables):
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
            
        for variable in milestone_variables:
            expr_var = "{}={};".format(variable, milestone_variables[variable])
            expr += expr_var
        
        #expr += base.convert_openmm_to_python_expr(
        #    "result="+self.openmm_fwd_rev_expression)
        expr += base.convert_openmm_to_python_expr(
            "result="+self.cv_expression)
        mylocals = locals()
        exec(expr, globals(), mylocals)
        result = mylocals["result"]
        if result <= 0:
            return True
        else:
            return False
    
    def check_value_within_boundary(self, positions, milestone_variables, 
                                    verbose=False, tolerance=0.0):
        """
        
        """
        result = self.check_positions_within_boundary(
            positions*openmm.unit.nanometer, milestone_variables)
        return result
    
    def check_mdtraj_close_to_boundary(self, traj, milestone_variables, 
                                     verbose=False, max_avg=0.03, max_std=0.05):
        """
        For now, this will just always return True.
        """
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

def make_elber_external_cv_object(external_cv_input, index):
    """
    Create a ExternalCV object to be placed into the Model.
    """
    
    groups = []
    for group in external_cv_input.groups:
        groups.append(group)
    
    cv = Elber_external_CV(index, groups)
    cv.cv_expression = external_cv_input.cv_expression
    assert cv.cv_expression is not None, \
        "A CV expression is required to define Elber milestones."
    
    cv.cv_expression = external_cv_input.cv_expression
    if external_cv_input.openmm_expression is None:
        cv.openmm_fwd_rev_expression = "step(k*("+cv.cv_expression+" - value))"
    else:
        cv.openmm_fwd_rev_expression = external_cv_input.openmm_expression
    
    if external_cv_input.restraining_expression is None:
        cv.openmm_umbrella_expression = "0.5*k*("+cv.cv_expression+" - value)^2"
    else:
        cv.openmm_umbrella_expression = external_cv_input.restraining_expression
    
    return cv
    
def make_elber_milestoning_objects_external(
        external_cv_input, milestone_alias, milestone_index, 
        input_anchor_index, anchor_index, input_anchors,
        umbrella_force_constant):
    """
    Make a set of 1-dimensional external Milestone objects to be put into
    an Anchor, and eventually, the Model.
    """
    
    milestones = []
    num_input_anchors = len(input_anchors)
    if input_anchor_index > 0:
        neighbor_index = anchor_index - 1
        neighbor_index_input = input_anchor_index - 1
        milestone1 = base.Milestone()
        milestone1.index = milestone_index - 1
        milestone1.neighbor_anchor_index = neighbor_index
        milestone1.alias_index = 1
        milestone1.is_source_milestone = False
        milestone1.cv_index = external_cv_input.index
        if input_anchors[input_anchor_index].lower_milestone_value is None:
            value = 0.5 * (input_anchors[input_anchor_index].value \
                            + input_anchors[neighbor_index_input].value)
        else:
            value = input_anchors[input_anchor_index].lower_milestone_value
        milestone1.variables = {"k": umbrella_force_constant, "value": value}
        milestones.append(milestone1)
    
    milestone2 = base.Milestone()
    milestone2.index = milestone_index
    milestone2.neighbor_anchor_index = milestone_index
    #milestone2.alias_index = milestone_alias
    milestone2.alias_index = 2
    milestone2.is_source_milestone = True
    milestone2.cv_index = external_cv_input.index
    value = input_anchors[input_anchor_index].value
    milestone2.variables = {"k": umbrella_force_constant, "value": value}
    milestones.append(milestone2)
    
    if input_anchor_index < num_input_anchors-1:
        neighbor_index = anchor_index + 1
        neighbor_index_input = input_anchor_index + 1
        milestone3 = base.Milestone()
        milestone3.index = milestone_index + 1
        milestone3.neighbor_anchor_index = neighbor_index
        milestone3.alias_index = 3
        milestone3.is_source_milestone = False
        milestone3.cv_index = external_cv_input.index
        if input_anchors[input_anchor_index].upper_milestone_value is None:
            value = 0.5 * (input_anchors[input_anchor_index].value \
                            + input_anchors[neighbor_index_input].value)
        else:
            value = input_anchors[input_anchor_index].upper_milestone_value
        milestone3.variables = {"k": umbrella_force_constant, "value": value}
        milestones.append(milestone3)
            
    return milestones, milestone_alias, milestone_index+1