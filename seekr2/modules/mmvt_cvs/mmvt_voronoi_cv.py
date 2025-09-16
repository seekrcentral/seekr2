"""
Structures, methods, and functions for handling Voronoi CVs.
"""

import time

import numpy as np

from parmed import unit

import seekr2.modules.common_base as base
import seekr2.modules.common_cv as common_cv
import seekr2.modules.mmvt_cvs.mmvt_cv_base as mmvt_cv_base
from seekr2.modules.mmvt_cvs.mmvt_cv_base import MMVT_collective_variable
import seekr2.modules.mmvt_cvs.mmvt_spherical_cv as mmvt_spherical_cv
import seekr2.modules.mmvt_cvs.mmvt_tiwary_cv as mmvt_tiwary_cv
import seekr2.modules.mmvt_cvs.mmvt_planar_cv as mmvt_planar_cv
import seekr2.modules.mmvt_cvs.mmvt_rmsd_cv as mmvt_rmsd_cv
import seekr2.modules.mmvt_cvs.mmvt_closest_pair_cv as mmvt_closest_pair_cv
import seekr2.modules.mmvt_cvs.mmvt_count_contacts_cv as mmvt_count_contacts_cv
import seekr2.modules.mmvt_cvs.mmvt_external_cv as mmvt_external_cv

class MMVT_Voronoi_CV(MMVT_collective_variable):
    """
    A Voronoi collective variable which contains multidimensional functions
    to produce a Voronoi Tessellation.
    
    """+MMVT_collective_variable.__doc__
    
    def __init__(self, index):
        self.index = index
        self.child_cvs = []
        self.name = "mmvt_voronoi"
        self.openmm_expression = None
        self.restraining_expression = None
        self.cv_expression = None
        self._mygroup_list = None
        self.num_groups = 0
        return

    def __name__(self):
        return "MMVT_Voronoi_CV"
    
    def make_cv_expr(self, alias_id):
        """
        
        """
        me_expr_list = []
        neighbor_expr_list = []
        num_groups = 0
        for i, child_cv in enumerate(self.child_cvs):
            me_cv_expr = "me_cv_{}_alias_{}".format(i, alias_id)
            me_expr_list.append(me_cv_expr)
            neighbor_cv_expr = "neighbor_cv_{}_alias_{}".format(i, alias_id)
            neighbor_expr_list.append(neighbor_cv_expr)
            num_groups += child_cv.num_groups
        
        me_expr = " + ".join(me_expr_list)
        neighbor_expr = " + ".join(neighbor_expr_list)
        cv_expr = "(" + me_expr + ") - (" + neighbor_expr + ")"
        self.cv_expression = cv_expr
        return cv_expr
    
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
        
        cv_expression = self.make_cv_expr(alias_id)
        self.openmm_expression = "step(" + cv_expression + ")"
        expression_w_bitcode = "bitcode_{}*".format(alias_id)+self.openmm_expression
        return openmm.CustomCVForce(expression_w_bitcode)
    
    def make_restraining_force(self, alias_id):
        """
        Create an OpenMM force object that will restrain the system to
        a given value of this CV.
        """
        try:
            import openmm
        except ImportError:
            import simtk.openmm as openmm
        
        cv_expression = self.make_cv_expr(alias_id)
        self.restraining_expression = "0.5*k_{}*(".format(alias_id) + cv_expression + ")^2"
        
        return openmm.CustomCVForce(self.restraining_expression)
    
    def make_cv_force(self, alias_id):
        raise Exception("Not yet implemented.")
    
    def make_namd_colvar_string(self):
        """
        This string will be put into a NAMD colvar file for tracking
        MMVT bounces.
        """
        raise Exception("MMVT Voronoi CVs are not available in NAMD")
    
    def add_groups(self, force):
        """
        
        """
        return
    
    def add_parameters(self, force):
        """
        An OpenMM custom force object needs a list of variables
        provided to it that will occur within its expression. Both
        the per-dof and global variables are combined within the
        variable_names_list. The numerical values of these variables
        will be provided at a later step.
        """
        return
    
    def add_groups_and_variables(self, force, variables, alias_id):
        """
        Provide the custom force with additional information it needs,
        which includes a list of the groups of atoms involved with the
        CV, as well as a list of the variables' *values*.
        """
        force.addGlobalParameter("bitcode_{}".format(alias_id), variables[0])
        force.addGlobalParameter("k_{}".format(alias_id), variables[1])
        for i, child_cv in enumerate(self.child_cvs):
            me_force, neighbor_force = child_cv.make_voronoi_cv_boundary_forces(
                variables[2*(i+1)], variables[2*(i+1)+1], alias_id)
            force.addCollectiveVariable("me_cv_{}_alias_{}".format(i, alias_id), me_force)
            force.addCollectiveVariable("neighbor_cv_{}_alias_{}".format(i, alias_id), neighbor_force)
            
        return
    
    def update_groups_and_variables(self, force, variables, alias_id, context):
        """
        Update the force's variables with a list of new values.
        """
        raise Exception(
            "update_groups_and_variables not yet implemented for Voronoi CVs.")
            
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
        values_list.append(bitcode)
        values_list.append(k)
        for i, child_cv in enumerate(self.child_cvs):
            me_i_value = milestone.variables["me_{}".format(i)]
            neighbor_i_value = milestone.variables["neighbor_{}".format(i)]
            values_list.append(me_i_value)
            values_list.append(neighbor_i_value)
            
        return values_list
    
    def get_namd_evaluation_string(self, milestone, cv_val_var="cv_val"):
        """
        For a given milestone, return a string that can be evaluated
        my NAMD to monitor for a crossing event. Essentially, if the 
        function defined by the string ever returns True, then a
        bounce will occur
        """
        raise Exception("MMVT RMSD CVs are not available in NAMD")
    
    def get_mdtraj_cv_value(self, traj, frame_index):
        """
        Determine the current CV value for an mdtraj object.
        """
        values = []
        for i, child_cv in enumerate(self.child_cvs):
            cv_value = child_cv.get_mdtraj_cv_value(traj, frame_index)
            values.append(cv_value)
        
        return values
    
    def check_mdtraj_within_boundary(self, traj, milestone_variables, 
                                     verbose=False, TOL=0.001):
        """
        Check if an mdtraj Trajectory describes a system that remains
        within the expected anchor. Return True if passed, return
        False if failed.
        """
        for frame_index in range(traj.n_frames):
            value = self.get_mdtraj_cv_value(traj, frame_index)
            result = self.check_value_within_boundary(
                value, milestone_variables, tolerance=TOL)
            if not result:
                return False
            
        return True
    
    def get_openmm_context_cv_value(
            self, context, positions=None, system=None):
        """
        
        """
        values = []
        for i, child_cv in enumerate(self.child_cvs):
            cv_value = child_cv.get_openmm_context_cv_value(
                context, positions=positions, system=system)
            values.append(cv_value)
            
        return values
    
    def get_cv_value(self, positions):
        # Only works if the children are all External_CVs used in
        # toy systems
        values = []
        for i, child_cv in enumerate(self.child_cvs):
            cv_value = child_cv.get_cv_value(positions=positions)
            values.append(cv_value)
            assert np.isfinite(cv_value), "Non-finite value detected."
            
        return values
    
    def check_openmm_context_within_boundary(
            self, context, milestone_variables, positions=None, 
            verbose=False, tolerance=0.0):
        """
        Check if an OpenMM context describes a system that remains
        within the expected anchor. Return True if passed, return
        False if failed.
        """
        values = self.get_openmm_context_cv_value(context, positions)
        result = self.check_value_within_boundary(
            values, milestone_variables, tolerance=tolerance)
        return result
    
    
    def check_value_within_boundary(self, value, milestone_variables, 
                                    tolerance=0.0):
        me_total_value = 0.0
        neighbor_total_value = 0.0
        for i, child_cv in enumerate(self.child_cvs):
            cv_value = value[i]
            me_i_value = milestone_variables["me_{}".format(i)]
            neighbor_i_value = milestone_variables["neighbor_{}".format(i)]
            me_this_value = (me_i_value - cv_value)**2
            neighbor_this_value = (neighbor_i_value - cv_value)**2
            me_total_value += me_this_value
            neighbor_total_value += neighbor_this_value
                
        distance = me_total_value - neighbor_total_value
        if (distance-tolerance) > 0.0: 
            return False
        
        return True
    
    def check_positions_within_boundary(self, positions, milestone_variables, 
                                        tolerance=0.0):
        me_total_value = 0.0
        neighbor_total_value = 0.0
        for i, child_cv in enumerate(self.child_cvs):
            cv_value = child_cv.get_cv_value(positions)
            me_i_value = milestone_variables["me_{}".format(i)]
            neighbor_i_value = milestone_variables["neighbor_{}".format(i)]
            me_this_value = (me_i_value - cv_value)**2
            neighbor_this_value = (neighbor_i_value - cv_value)**2
            me_total_value += me_this_value
            neighbor_total_value += neighbor_this_value
            
        distance = me_total_value - neighbor_total_value
        if (distance-tolerance) > 0.0: 
            return False
        
        return True
        
    def check_mdtraj_close_to_boundary(self, traj, milestone_variables, 
                                     verbose=False, max_avg=0.03, max_std=0.05):
        """
        Given an mdtraj Trajectory, check if the system lies close
        to the MMVT boundary. Return True if passed, return False if 
        failed.
        """
        TOL = 0.001
        diffs = []
        for frame_index in range(traj.n_frames):
            values = self.get_mdtraj_cv_value(traj, frame_index)
            me_dist = 0.0
            neighbor_dist = 0.0
            for i, child_cv in enumerate(self.child_cvs):
                me_val = milestone_variables["me_{}".format(i)]
                neighbor_val = milestone_variables["neighbor_{}".format(i)]
                me_dist += (me_val - values[i])**2
                neighbor_dist += (neighbor_val - values[i])**2
            
            diff = me_dist - neighbor_dist
            diffs.append(diff)
            
        avg_diff = np.mean(diffs)
        std_diff = np.std(diffs)
        if abs(avg_diff) > max_avg or std_diff > max_std:
            if verbose:
                warnstr = """The distance between the system and central 
    milestone were found on average to be {:.4f} nm apart.
    The standard deviation was {:.4f} nm.""".format(avg_diff, std_diff)
                print(warnstr)
            return False
            
        return True
    
    def get_atom_groups(self):
        """
        Return a 1-list of this CV's atomic group.
        """
        groups = []
        for i, child_cv in enumerate(self.child_cvs):
            groups_cv = child_cv.get_atom_groups()
            groups += groups_cv
        return groups
            
    def get_variable_values(self):
        """
        Return the child variable values.
        """
        return []


def make_mmvt_voronoi_cv_object(voronoi_cv_input, index, root_directory):
    """
    Create a RMSD CV object to be placed into the Model.
    """
    # TODO: these need to be imported and called from different modules
    cv = MMVT_Voronoi_CV(index)
    for i, cv_input in enumerate(voronoi_cv_input.cv_inputs):
        cv_input.check()
        if isinstance(cv_input, common_cv.Spherical_cv_input):
            child_cv = mmvt_spherical_cv.make_mmvt_spherical_cv_object(cv_input, index=i)
        elif isinstance(cv_input, common_cv.Tiwary_cv_input):
            child_cv = mmvt_tiwary_cv.make_mmvt_tiwary_cv_object(cv_input, index=i)
        elif isinstance(cv_input, common_cv.Planar_cv_input):
            child_cv = mmvt_planar_cv.make_mmvt_planar_cv_object(cv_input, index=i)
        elif isinstance(cv_input, common_cv.RMSD_cv_input):
            child_cv = mmvt_rmsd_cv.make_mmvt_RMSD_cv_object(
                cv_input, index=i, root_directory=root_directory)
        elif isinstance(cv_input, common_cv.Closest_pair_cv_input):
            child_cv = mmvt_closest_pair_cv.make_mmvt_closest_pair_cv_object(cv_input, index=i)
        elif isinstance(cv_input, common_cv.Toy_cv_input):
            child_cv = mmvt_external_cv.make_mmvt_external_cv_object(cv_input, index=i)
        else:
            raise Exception("CV type not available for Voronoi CV: %s" \
                            % type(cv_input))
        
        cv.child_cvs.append(child_cv)
    
    return cv
