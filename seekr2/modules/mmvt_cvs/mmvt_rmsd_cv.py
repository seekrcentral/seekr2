"""
Structures, methods, and functions for handling RMSD CVs.
"""

import os
import shutil

import numpy as np
from scipy.spatial import transform
import mdtraj

from parmed import unit

import seekr2.modules.common_base as base
import seekr2.modules.mmvt_cvs.mmvt_cv_base as mmvt_cv_base
from seekr2.modules.mmvt_cvs.mmvt_cv_base import MMVT_collective_variable

class MMVT_RMSD_CV(MMVT_collective_variable):
    """
    A RMSD collective variable which is the root mean squared deviation
    of the system from a reference structure.
    
    """+MMVT_collective_variable.__doc__
    
    def __init__(self, index, group, ref_structure, align_group=None):
        self.index = index
        self.group = group
        self.align_group = align_group
        self.ref_structure = ref_structure
        self.name = "mmvt_rmsd"
        self.openmm_expression = None
        self.restraining_expression = None
        self.cv_expression = "RMSD"
        self.num_groups = 1
        self.per_dof_variables = []
        self.global_variables = [] #["k", "value"]
        self._mygroup_list = None
        self.variable_name = "v"
        return

    def __name__(self):
        return "MMVT_RMSD_CV"
    
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
            
        assert self.num_groups == 1
        self.openmm_expression = "step(k_{}*(RMSD - value_{}))".format(alias_id, alias_id)
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
        
        assert self.num_groups == 1
        self.restraining_expression = "0.5*k_{}*(RMSD - value_{})^2".format(alias_id, alias_id)
        return openmm.CustomCVForce(self.restraining_expression)
    
    def make_cv_force(self, alias_id):
        try:
            import openmm
        except ImportError:
            import simtk.openmm as openmm
        return openmm.CustomCentroidBondForce(
            self.num_groups, self.cv_expression)
    
    def make_voronoi_cv_boundary_forces(self, me_val, neighbor_val, alias_id):
        """
        
        """
        try:
            import openmm
        except ImportError:
            import simtk.openmm as openmm
        
        try:
            import openmm.app as openmm_app
        except ImportError:
            import simtk.openmm.app as openmm_app
            
        pdb_file = openmm_app.PDBFile(self.ref_structure)
        
        if self.align_group is None:
            rmsd_me_force = openmm.RMSDForce(pdb_file.positions, self.group)
            rmsd_neighbor_force = openmm.RMSDForce(pdb_file.positions, self.group)
        else:
            try:
                import rmsdplusforceplugin
            except ImportError:
                print("Unable to load RMSDPlusForcePlugin. Please install "\
                      "from: https://github.com/seekrcentral/"\
                      "rmsdplusforceplugin.git")
                exit()
            
            rmsd_me_force = rmsdplusforceplugin.RMSDPlusForce(
                pdb_file.positions, self.group, self.align_group)
            rmsd_neighbor_force = rmsdplusforceplugin.RMSDPlusForce(
                pdb_file.positions, self.group, self.align_group)
        
        
        me_expr = "(me_val_{}_alias_{} - {})^2".format(self.index, alias_id, self.cv_expression)
        me_force = openmm.CustomCVForce(me_expr)
        me_force.addGlobalParameter("me_val_{}_alias_{}".format(self.index, alias_id), me_val)
        me_force.addCollectiveVariable("RMSD", rmsd_me_force)
        
        neighbor_expr = "(neighbor_val_{}_alias_{} - {})^2".format(self.index, alias_id, self.cv_expression)
        neighbor_force = openmm.CustomCVForce(neighbor_expr)
        neighbor_force.addGlobalParameter("neighbor_val_{}_alias_{}".format(self.index, alias_id), neighbor_val)
        neighbor_force.addCollectiveVariable("RMSD", rmsd_neighbor_force)
        
        return me_force, neighbor_force
    
    def update_voronoi_cv_boundary_forces(
            self, me_force, me_val, neighbor_force, neighbor_val, alias_id, 
            context):
        """
        
        """
        #me_force.setGlobalParameterDefaultValue(0, me_val)
        #neighbor_force.setGlobalParameterDefaultValue(0, neighbor_val)
        context.setParameter("me_val_{}_alias_{}".format(self.index, alias_id), me_val)
        context.setParameter("neighbor_val_{}_alias_{}".format(self.index, alias_id), neighbor_val)
        return
    
    def make_namd_colvar_string(self):
        """
        This string will be put into a NAMD colvar file for tracking
        MMVT bounces.
        """
        raise Exception("MMVT RMSD CVs are not available in NAMD")
    
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
        try:
            import openmm
        except ImportError:
            import simtk.openmm as openmm
            
        try:
            import openmm.app as openmm_app
        except ImportError:
            import simtk.openmm.app as openmm_app
            
        pdb_file = openmm_app.PDBFile(self.ref_structure)
        if self.align_group is None:
            rmsd_force = openmm.RMSDForce(pdb_file.positions, self.group)
        else:
            try:
                import rmsdplusforceplugin
            except ImportError:
                print("Unable to load RMSDPlusForcePlugin. Please install "\
                      "from: https://github.com/seekrcentral/"\
                      "rmsdplusforceplugin.git")
                exit()
            
            rmsd_force = rmsdplusforceplugin.RMSDPlusForce(
                pdb_file.positions, self.group, self.align_group)
            
        force.addCollectiveVariable("RMSD", rmsd_force)
        
        return
    
    def add_groups_and_variables(self, force, variables, alias_id):
        """
        Provide the custom force with additional information it needs,
        which includes a list of the groups of atoms involved with the
        CV, as well as a list of the variables' *values*.
        """
        force.addGlobalParameter("bitcode_{}".format(alias_id), variables[0])
        force.addGlobalParameter("k_{}".format(alias_id), variables[1])
        force.addGlobalParameter("value_{}".format(alias_id), variables[2])
        return
    
    def update_groups_and_variables(self, force, variables, alias_id, context):
        """
        Update the force's variables with a list of new values.
        """
        #force.setGlobalParameterDefaultValue(0, variables[0])
        #force.setGlobalParameterDefaultValue(1, variables[1])
        #force.setGlobalParameterDefaultValue(2, variables[2])
        context.setParameter("bitcode_{}".format(alias_id), variables[0])
        context.setParameter("k_{}".format(alias_id), variables[1])
        context.setParameter("value_{}".format(alias_id), variables[2])
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
        value = milestone.variables["value"] * unit.nanometers
        values_list.append(bitcode)
        values_list.append(k)
        values_list.append(value)
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
        if self.align_group is None:
            align_group = self.group
        else:
            align_group = self.align_group
        assert os.path.exists(self.ref_structure), \
            "File {} does not exist. Make sure ".format(self.ref_structure) \
            +"that any programs using the get_mdtraj_cv_value() method " \
            "within an API is performed in the model directory."
        ref_traj = mdtraj.load(self.ref_structure)
        traj.superpose(ref_traj, atom_indices=align_group)
        ref_xyz = ref_traj.xyz[0]
        traj_xyz = traj.xyz[frame_index]
        value = 0.0
        for i in self.group:
            increment = (traj_xyz[i,0] - ref_xyz[i,0])**2 \
                + (traj_xyz[i,1] - ref_xyz[i,1])**2 \
                + (traj_xyz[i,2] - ref_xyz[i,2])**2
            value += increment
        
        rmsd = np.sqrt(value/len(self.group))
        return rmsd
    
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
                value, milestone_variables, verbose=verbose, tolerance=TOL)
            if not result:
                return False
            
        return True
    
    def get_openmm_context_cv_value(
            self, context, positions=None, ref_positions=None, verbose=False, 
            system=None, tolerance=0.0):
        """
        
        """
        try:
            import openmm.app as openmm_app
        except ImportError:
            import simtk.openmm.app as openmm_app
            
        if system is None:
            system = context.getSystem()
        if positions is None:
            state = context.getState(getPositions=True)
            positions = state.getPositions()
            
        if ref_positions is None:
            pdb_filename = self.ref_structure
            pdb_file = openmm_app.PDBFile(pdb_filename)
            ref_positions = pdb_file.positions
        
        pos_subset = []
        pos_subset_rmsd = []
        ref_subset = []
        ref_subset_rmsd = []
        if self.align_group is None:
            align_group = self.group
        else:
            align_group = self.align_group
        for atom_index in align_group:
            pos_subset.append(positions[atom_index].value_in_unit(unit.nanometers))
            ref_subset.append(ref_positions[atom_index].value_in_unit(unit.nanometers))
        
        for atom_index in self.group:
            pos_subset_rmsd.append(positions[atom_index].value_in_unit(unit.nanometers))
            ref_subset_rmsd.append(ref_positions[atom_index].value_in_unit(unit.nanometers))
        
        pos_subset = np.array(pos_subset)
        ref_subset = np.array(ref_subset)
        pos_subset_rmsd = np.array(pos_subset_rmsd)
        ref_subset_rmsd = np.array(ref_subset_rmsd)
        
        avg_pos = np.array([float(np.mean(pos_subset[:,0])), 
                            float(np.mean(pos_subset[:,1])), 
                            float(np.mean(pos_subset[:,2]))])
        avg_ref = np.array([np.mean(ref_subset[:,0]), 
                            np.mean(ref_subset[:,1]), 
                            np.mean(ref_subset[:,2])])        
        pos_subset_centered = pos_subset - avg_pos
        ref_subset_centered = ref_subset - avg_ref
        pos_subset_rmsd_centered = pos_subset_rmsd - avg_pos
        ref_subset_rmsd_centered = ref_subset_rmsd - avg_ref
        rotation, rssd_align = transform.Rotation.align_vectors(
            pos_subset_centered, ref_subset_centered)
        new_pos_subset_rmsd = rotation.apply(pos_subset_rmsd_centered)
        
        value = 0.0
        for i in range(len(self.group)):
            increment = (new_pos_subset_rmsd[i,0] - ref_subset_rmsd_centered[i,0])**2 \
                + (new_pos_subset_rmsd[i,1] - ref_subset_rmsd_centered[i,1])**2 \
                + (new_pos_subset_rmsd[i,2] - ref_subset_rmsd_centered[i,2])**2
            value += increment
            
        rmsd = np.sqrt(value / len(self.group))
        assert np.isfinite(rmsd), "Non-finite value detected."
        return rmsd
    
    def check_openmm_context_within_boundary(
            self, context, milestone_variables, positions=None, 
            ref_positions=None, verbose=False, tolerance=0.0):
        """
        Check if an OpenMM context describes a system that remains
        within the expected anchor. Return True if passed, return
        False if failed.
        """
        
        value = self.get_openmm_context_cv_value(
            context, positions=positions, 
            ref_positions=ref_positions, verbose=verbose, tolerance=tolerance)
        result = self.check_value_within_boundary(
            value, milestone_variables, verbose=verbose, tolerance=tolerance)
        return result
    
    
    def check_value_within_boundary(self, value, milestone_variables, 
                                    verbose=False, tolerance=0.0):
        milestone_k = milestone_variables["k"]
        milestone_value = milestone_variables["value"]
        if milestone_k*(value - milestone_value) > tolerance:
            if verbose:
                warnstr = """The RMSD value of atom group 
was found to be {:.4f} nm apart.
This distance falls outside of the milestone
boundary at {:.4f} nm.""".format(value, milestone_value)
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
        
        TOL = 0.001
        diffs = []
        for frame_index in range(traj.n_frames):
            value = self.get_mdtraj_cv_value(traj, frame_index)
            milestone_value = milestone_variables["value"]
            diff = value - milestone_value
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
        return [self.group]
            
    def get_variable_values(self):
        """
        This type of CV has no extra variables, so an empty list is 
        returned.
        """
        return []

def make_mmvt_RMSD_cv_object(RMSD_cv_input, index, root_directory):
    """
    Create a RMSD CV object to be placed into the Model.
    """
    RMSD_ref_pdb = "rmsd_reference_cv_{}.pdb".format(index)
    group = base.parse_xml_list(RMSD_cv_input.group)
    align_group = base.parse_xml_list(RMSD_cv_input.align_group)
    absolute_RMSD_ref_pdb = os.path.join(root_directory, RMSD_ref_pdb)
    shutil.copyfile(RMSD_cv_input.ref_structure, absolute_RMSD_ref_pdb)
    cv = MMVT_RMSD_CV(index, group, RMSD_ref_pdb, align_group)
    return cv
