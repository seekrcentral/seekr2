"""
Structures, methods, and functions for handling closest-pair CVs.
"""

import numpy as np

from parmed import unit

import seekr2.modules.common_base as base
import seekr2.modules.mmvt_cvs.mmvt_cv_base as mmvt_cv_base
from seekr2.modules.mmvt_cvs.mmvt_cv_base import MMVT_collective_variable

class MMVT_closest_pair_CV(MMVT_collective_variable):
    """
    A collective variable represents the closest distance between a pair
    of atoms out of two sets of atoms.
    
    """+MMVT_collective_variable.__doc__
    
    def __init__(self, index, groups):
        self.index = index
        self.group1 = groups[0]
        self.group2 = groups[1]
        self.name = "mmvt_closest_pair"
        self.openmm_expression = None
        self.restraining_expression = None
        self.cv_expression = "closest"
        self.num_groups = 2
        self.per_dof_variables = []
        self.global_variables = []
        self._mygroup_list = None
        self.variable_name = "v"
        self.exponent = 50.0
        self.num_system_particles = 0
        self.exclusion_pairs = []
        self.cutoff_distance = 10.0
        return

    def __name__(self):
        return "MMVT_closest_pair_CV"
    
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
        self.openmm_expression \
            = "step(k_{}*((closest + (1/{})^exponent)^(-1.0/exponent) - value_{}))".format(
                alias_id, self.cutoff_distance, alias_id)
        expression_w_bitcode = "bitcode_{}*".format(alias_id)\
            +self.openmm_expression
        force = openmm.CustomCVForce(expression_w_bitcode)
        return force
    
    def make_restraining_force(self, alias_id):
        """
        Create an OpenMM force object that will restrain the system to
        a given value of this CV.
        """
        raise Exception("Restraining force not available for Closest pair CV.")
        """
        try:
            import openmm
        except ImportError:
            import simtk.openmm as openmm
            
        assert self.num_groups == 2
        self.restraining_expression \
            = "0.5*k_{}*((closest + (1/{})^exponent)^(-1.0/exponent) - value_{})^2"\
                .format(alias_id, self.cutoff_distance, alias_id)
        force = openmm.CustomCVForce(self.restraining_expression)
        return force
        """
    
    def make_voronoi_cv_boundary_forces(self, me_val, neighbor_val, alias_id):
        """
        
        """
        try:
            import openmm
        except ImportError:
            import simtk.openmm as openmm
            
        #pdb_file = openmm_app.PDBFile(self.ref_structure)
        #rmsd_me_force = openmm.RMSDForce(pdb_file.positions, self.group)
        #rmsd_neighbor_force = openmm.RMSDForce(pdb_file.positions, self.group)
        closest_force_me = openmm.CustomNonbondedForce("(1/r)^exponent")
        closest_force_me.addGlobalParameter("exponent", self.exponent)
        closest_force_me.addInteractionGroup(self.group1, self.group2)
        closest_force_me.setUseSwitchingFunction(False)
        closest_force_me.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
        closest_force_me.setCutoffDistance(self.cutoff_distance*unit.nanometers)
        closest_force_me.setUseLongRangeCorrection(False)
        for i in range(self.num_system_particles):
            closest_force_me.addParticle([])
        for (iatom, jatom) in self.exclusion_pairs:
            closest_force_me.addExclusion(iatom, jatom)
            
        closest_force_neighbor = openmm.CustomNonbondedForce("(1/r)^exponent")
        closest_force_neighbor.addGlobalParameter("exponent", self.exponent)
        closest_force_neighbor.addInteractionGroup(self.group1, self.group2)
        closest_force_neighbor.setUseSwitchingFunction(False)
        closest_force_neighbor.setNonbondedMethod(
            openmm.CustomNonbondedForce.CutoffPeriodic)
        closest_force_neighbor.setCutoffDistance(self.cutoff_distance*unit.nanometers)
        closest_force_neighbor.setUseLongRangeCorrection(False)
        for i in range(self.num_system_particles):
            closest_force_neighbor.addParticle([])
        for (iatom, jatom) in self.exclusion_pairs:
            closest_force_neighbor.addExclusion(iatom, jatom)
        
        me_expr = "(({} + (1/{})^exponent)^(-1.0/exponent) - me_val_{}_alias_{})^2".format(
            self.cv_expression, self.cutoff_distance, self.index, alias_id)
        me_force = openmm.CustomCVForce(me_expr)
        me_force.addGlobalParameter("me_val_{}_alias_{}".format(
            self.index, alias_id), me_val)
        me_force.addGlobalParameter("exponent", self.exponent)
        me_force.addCollectiveVariable(self.cv_expression, closest_force_me)
        
        neighbor_expr \
            = "(({} + (1/{})^exponent)^(-1.0/exponent) - neighbor_val_{}_alias_{})^2".format(
            self.cv_expression, self.cutoff_distance, self.index, alias_id)
        neighbor_force = openmm.CustomCVForce(neighbor_expr)
        neighbor_force.addGlobalParameter("neighbor_val_{}_alias_{}".format(
            self.index, alias_id), neighbor_val)
        neighbor_force.addGlobalParameter("exponent", self.exponent)
        neighbor_force.addCollectiveVariable(
            self.cv_expression, closest_force_neighbor)
        
        return me_force, neighbor_force
    
    def update_voronoi_cv_boundary_forces(
            self, me_force, me_val, neighbor_force, neighbor_val, alias_id, 
            context):
        """
        
        """
        context.setParameter("me_val_{}_alias_{}".format(self.index, alias_id), me_val)
        context.setParameter("neighbor_val_{}_alias_{}".format(self.index, alias_id), neighbor_val)
        return
    
    def make_namd_colvar_string(self):
        """
        This string will be put into a NAMD colvar file for tracking
        MMVT bounces.
        """
        raise Exception("MMVT closest pair CVs are not available in NAMD")
    
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
            
        closest_force = openmm.CustomNonbondedForce("(1/r)^exponent")
        closest_force.addGlobalParameter("exponent", self.exponent)
        closest_force.addInteractionGroup(self.group1, self.group2)
        closest_force.setUseSwitchingFunction(False)
        closest_force.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
        closest_force.setCutoffDistance(self.cutoff_distance*unit.nanometers)
        closest_force.setUseLongRangeCorrection(False)
        for i in range(self.num_system_particles):
            closest_force.addParticle([])
        for (iatom, jatom) in self.exclusion_pairs:
            closest_force.addExclusion(iatom, jatom)
        
        force.addCollectiveVariable("closest", closest_force)
        self.add_groups(force)
        variable_names_list = []
        return variable_names_list
    
    def add_groups_and_variables(self, force, variables, alias_id):
        """
        Provide the custom force with additional information it needs,
        which includes a list of the groups of atoms involved with the
        CV, as well as a list of the variables' *values*.
        """
        force.addGlobalParameter("bitcode_{}".format(alias_id), variables[0])
        force.addGlobalParameter("k_{}".format(alias_id), variables[1])
        force.addGlobalParameter("value_{}".format(alias_id), variables[2])
        force.addGlobalParameter("exponent", self.exponent)
        return
    
    def update_groups_and_variables(self, force, variables, alias_id, context):
        """
        Update the force's variables with a list of new values.
        """
        #force.setGlobalParameter(0, "bitcode_{}".format(alias_id), variables[0])
        #force.setGlobalParameter(1, "k_{}".format(alias_id), variables[1])
        #force.setGlobalParameter(2, "value_{}".format(alias_id), variables[2])
        #force.setGlobalParameter(3, "exponent".format(alias_id), variables[3])
        context.setParameter("bitcode_{}".format(alias_id), variables[0])
        context.setParameter("k_{}".format(alias_id), variables[1])
        context.setParameter("value_{}".format(alias_id), variables[2])
        context.setParameter("exponent", variables[3])
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
        raise Exception("MMVT Closest pair CVs are not available in NAMD")
    
    def get_mdtraj_cv_value(self, traj, frame_index):
        """
        Determine the current CV value for an mdtraj object.
        """
        atom_set_list = []
        for index in self.group1:
            atom_set_list.append(set([traj.topology.atom(index)]))
        #for index in self.group2:
        #    atom_set_list.append(set([traj.topology.atom(index)]))
        traj.image_molecules(inplace=True, anchor_molecules=atom_set_list)
        sum = (1.0 / self.cutoff_distance) ** self.exponent
        for atom_index1 in self.group1:
            for atom_index2 in self.group2:
                com1 = traj.xyz[frame_index, atom_index1,:]
                com2 = traj.xyz[frame_index, atom_index2,:]
                dist = np.linalg.norm(com2-com1)
                if dist < self.cutoff_distance:
                    sum += (1.0/dist) ** self.exponent
                
        min_value = sum ** (-1.0/self.exponent)
        return min_value
    
    def check_mdtraj_within_boundary(self, traj, milestone_variables, 
                                     verbose=False, TOL=0.001):
        """
        Check if an mdtraj Trajectory describes a system that remains
        within the expected anchor. Return True if passed, return
        False if failed.
        """
        atom_set_list = []
        for index in self.group1:
            atom_set_list.append(set([traj.topology.atom(index)]))
        #for index in self.group2:
        #    atom_set_list.append(set([traj.topology.atom(index)]))
        traj.image_molecules(inplace=True, anchor_molecules=atom_set_list)
        #traj.image_molecules(inplace=True, anchor_molecules=[self.group1])
        for frame_index in range(traj.n_frames):
            value = self.get_mdtraj_cv_value(traj, frame_index)
            result = self.check_value_within_boundary(
                value, milestone_variables, verbose, tolerance=TOL)
            if not result:
                return False
            
        return True
        
    def get_openmm_context_cv_value(self, context, positions=None, system=None):
        """
        Determine the current CV value for an openmm context.
        """
        # TODO: this system will not be properly imaged
        if system is None:
            system = context.getSystem()
        if positions is None:
            state = context.getState(getPositions=True)
            positions = state.getPositions()
        sum = (1.0 / self.cutoff_distance) ** self.exponent
        for atom_index1 in self.group1:
            for atom_index2 in self.group2:
                com1 = positions[atom_index1]
                com2 = positions[atom_index2]
                dist = np.linalg.norm(com2-com1)
                if dist < self.cutoff_distance:
                    sum += (1.0/dist) ** self.exponent
                
        min_value = sum ** (-1.0/self.exponent)
        return min_value
        
    def check_openmm_context_within_boundary(
            self, context, milestone_variables, positions=None, verbose=False,
            tolerance=0.0):
        """
        Check if an OpenMM context describes a system that remains
        within the expected anchor. Return True if passed, return
        False if failed.
        """
        value = self.get_openmm_context_cv_value(context, positions)
        result = self.check_value_within_boundary(
            value, milestone_variables, verbose, tolerance)
        return result
        
    def check_value_within_boundary(self, value, milestone_variables, 
                                    verbose=False, tolerance=0.0):
        """
        Check if the given CV value is within the milestone's boundaries.
        """
        milestone_k = milestone_variables["k"]
        milestone_value = milestone_variables["value"]
        if milestone_k*((1.0/milestone_value)**self.exponent - (1.0/value)**self.exponent) > tolerance:
            if verbose:
                warnstr = """The pair distance of atom group1 and atom
group2 were found to be {:.4f} nm apart.
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
        atom_set_list = []
        for index in self.group1:
            atom_set_list.append(set([traj.atoms[index]]))
        traj.image_molecules(inplace=True, anchor_molecules=atom_set_list)
        #traj.image_molecules(inplace=True, anchor_molecules=[self.group1])
        distances = []
        for frame_index in range(traj.n_frames):
            sum = (1.0 / self.cutoff_distance) ** self.exponent
            for atom_index1 in self.group1:
                for atom_index2 in self.group2:
                    com1 = traj.xyz[frame_index, atom_index1,:]
                    com2 = traj.xyz[frame_index, atom_index2,:]
                    value = np.linalg.norm(com2-com1)
                    if value < self.cutoff_distance:
                        sum += (1.0/value) ** self.exponent
                    
            milestone_value = milestone_variables["value"]
            min_value = sum ** (-1.0/self.exponent)
            distances.append(min_value - milestone_value)
            
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

def make_mmvt_closest_pair_cv_object(closest_pair_cv_input, index):
    """
    Create a closest pair CV object to be placed into the Model.
    """
    group1 = base.parse_xml_list(closest_pair_cv_input.group1)
    group2 = base.parse_xml_list(closest_pair_cv_input.group2)
    groups = [group1, group2]
    cv = MMVT_closest_pair_CV(index, groups)
    cv.exponent = closest_pair_cv_input.exponent
    return cv