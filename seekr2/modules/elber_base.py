"""
elber_base.py

Base classes, objects, and constants used in multiple stages of the
Elber milestoning calculations.
"""

# TODO: update this code and documentation once the Elber plugin is
# fixed.

import numpy as np
from parmed import unit
import mdtraj

#import seekr2.libraries.serializer.serializer as serializer
from abserdes import Serializer

OPENMM_ELBER_BASENAME = "forward"
OPENMM_ELBER_EXTENSION = "out"
OPENMM_ELBER_GLOB = "%s*.%s" % (OPENMM_ELBER_BASENAME, OPENMM_ELBER_EXTENSION)
NAMD_ELBER_BASENAME = "forward"
NAMD_ELBER_EXTENSION = "out" # TODO: consolidate OpenMM vs. NAMD output: not necessary
NAMD_ELBER_GLOB = "%s*.%s*" % (NAMD_ELBER_BASENAME, NAMD_ELBER_EXTENSION)

ELBER_UMBRELLA_BASENAME = "umbrella"
ELBER_FWD_BASENAME = "reverse"
ELBER_FWD_EXTENSION = "out"
ELBER_FWD_GLOB = "%s*.%s" % (ELBER_FWD_BASENAME, ELBER_FWD_EXTENSION)
ELBER_REV_BASENAME = "reverse"
ELBER_REV_EXTENSION = "out"
ELBER_REV_GLOB = "%s*.%s" % (ELBER_REV_BASENAME, ELBER_REV_EXTENSION)

class Elber_settings(Serializer):
    """
    Settings that are specific to an Elber milestoning calculation.
    
    Attributes:
    -----------
    temperature_equil_progression : list
        A list of temperatures (in Kelvin) to warm the simulation to
        during the temperature equilibration stage.
    num_temperature_equil_steps : int
        The number of steps to do per temperature during the
        temperature equilibration stage.
    temperature_equil_trajectory_interval : int or None
        The interval to write trajectory frames during the temperature
        equilibration stage. If None, then the trajectory won't be 
        written
    num_umbrella_stage_steps : int
        The number of steps to take within a given MMVT production
        run for a Voronoi cell.
    umbrella_stage_trajectory_interval : int
        The interval to write trajectory frames during the umbrella
        stage.
    """
    #num_equilibration_steps : int
    #    The number of steps to take during an equilibration run, where
    #    no statistics will be reported
    def __init__(self):
        #self.temperature_equil_progression = [
        #    300., 310., 320., 330., 340., 350., 340., 330., 320., 310., 300]
        self.temperature_equil_progression = []
        self.num_temperature_equil_steps = 1000
        self.num_umbrella_stage_steps = 50000
        self.umbrella_force_constant = 9000.0
        self.fwd_rev_interval = 500
        self.num_rev_launches = 1
        self.umbrella_energy_reporter_interval = None
        self.umbrella_trajectory_reporter_interval = None
        self.rev_energy_reporter_interval = None
        self.rev_trajectory_reporter_interval = None
        self.fwd_energy_reporter_interval = None
        self.fwd_trajectory_reporter_interval = None

class Elber_collective_variable(Serializer):
    """
    Collective variables represent the function of system positions
    and velocities so that Elber milestones can be defined
    
    Attributes:
    -----------
    index : int
        Every collective variable needs an index so that it may be
        quickly and easily referenced by one of the many milestones
        in the model.
        
    name : str
        Each type of collective variable has a shorthand 'name' for
        quick reference and identification. Example: 'elber_spherical'.
        
    openmm_umbrella_expression : str
        In order to restrain a system along a milestone, an umbrella 
        sampling potential energy expression must be applied.
        
    num_groups : int
        The number of atomic groups that are needed for the function
        describing this collective variable. Example: 2 for spherical
        CVs because a distance requires two points.
        
    groups : list
        A list of lists of integers. The length of the outer list is
        equal to self.num_groups. The inner lists contain integer
        values representing the indices of atoms in that group.
        
    per_dof_variables : list
        A list of strings of the names of variables used in 
        self.expression that apply to individual degrees of 
        freedom.
        
    global_variables : list
        A list of strings of the names of variables used in
        self.expression that apply globally, regardless of the degrees
        of freedom.
        
    """
    def __init__(self, index, groups):
        self.index = index
        self.groups = groups
        return

    def __name__(self):
        return "elber_baseCV"
    
    def make_force_object(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
    
    def make_namd_colvar_umbrella_string(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
        
    def add_parameters(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
        
    def add_groups_and_variables(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
        
    def get_variable_values_list(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
        
    def get_namd_evaluation_string(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
        
    def check_mdtraj_within_boundary(self, parmed_structure, 
                                               milestone_variables):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
    
    def get_atom_groups(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")

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
        com1_array = mdtraj.compute_center_of_mass(traj1)
        com2_array = mdtraj.compute_center_of_mass(traj2)
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

class Elber_anchor(Serializer):
    """
    An anchor object for representing a Voronoi cell in an Elber 
    milestoning calculation.
    
    Attributes
    ----------
    index : int
        The index of this anchor (cell) within the model.
    
    directory : str
        The directory (within the model's root directory) that contains
        the information and calculations for this Voronoi cell.
        
    amber_params : Amber_params
        Settings if this anchor starts the simulation using the
        AMBER forcefield and files.
    
    charmm_params : Charmm_params
        Settings if this anchor starts the simulation using the
        CHARMM forcefield and files.
        
    forcefield_params : Forcefield_params
        Settings if this anchor starts the simulation using an XML
        forcefield file and a PDB.
    
    md_directory : str or None
        The directory within the 'directory' argument above which 
        contains the MD simulation information. If None, then no MD
        is performed for this anchor.
        
    bd_directory : str or None
        The directory within the 'directory' argument above which
        contains the BD simulation information. If None, then no BD
        is performed for this anchor.
        
    production_directory : str
        The directory within the MD or BD directory above in which the
        simulations will be performed.
        
    md_output_glob : str
        A glob to select all the MD output files within the production
        directory above.
        
    name : str
        A unique name for this anchor.
        
    md : bool
        A boolean of whether MD is performed in this Voronoi cell.
        
    bd : bool
        A boolean of whether BD is performed in this Voronoi cell.
        
    endstate : bool
        A boolean of whether this is an end state or not - does it
        act as the bulk or a bound state or another state of interest?
        All end states will have kinetics calculated to all other
        end states.
        
    bulkstate : bool
        A boolean of whether this state acts as the bulk state (That
        is, the state represents a large separation distance between
        ligand and receptor.
        
    milestones : list
        A list of Milestone() objects, which are the boundaries 
        bordering this cell.
    """
    def __init__(self):
        self.index = 0
        self.directory = ""
        self.amber_params = None
        self.charmm_params = None
        self.forcefield_params = None
        self.building_directory = "building"
        self.production_directory = "prod"
        self.md_output_glob = OPENMM_ELBER_GLOB
        self.name = ""
        self.md = False
        self.endstate = False
        self.bulkstate = False
        self.milestones = []
        self.variables = {}
        return
    
    def _make_milestone_collection(self):
        """
        Make the dictionaries that allow for easy access of milestone
        indices, aliases, and neighboring indices.
        """
        id_key_alias_value_dict = {}
        alias_key_id_value_dict = {}
        neighbor_id_key_alias_value_dict = {}
        
        for milestone in self.milestones:
                index = milestone.index
                neighbor_index = milestone.neighbor_anchor_index
                alias_index = milestone.alias_index
                id_key_alias_value_dict[index] = alias_index
                neighbor_id_key_alias_value_dict[neighbor_index] = alias_index
                alias_key_id_value_dict[alias_index] = index
        
        return id_key_alias_value_dict, alias_key_id_value_dict, \
            neighbor_id_key_alias_value_dict
    
    def id_from_alias(self, alias_id):
        """
        Accept the alias index of a milestone and return the model-wide
        index.
        """
        id_key_alias_value_dict, alias_key_id_value_dict, \
            neighbor_id_key_alias_value_dict = self._make_milestone_collection()
        if alias_id in alias_key_id_value_dict:
            return alias_key_id_value_dict[alias_id]
        else:
            return None
    
    def alias_from_id(self, my_id):
        """
        Accept the model-wide index and return the milestone's alias
        index.
        """
        id_key_alias_value_dict, alias_key_id_value_dict, \
            neighbor_id_key_alias_value_dict = self._make_milestone_collection()
        if my_id in id_key_alias_value_dict:
            return id_key_alias_value_dict[my_id]
        else:
            return None
        
    def alias_from_neighbor_id(self, neighbor_id):
        """
        Take the index of the neighbor anchor's index and provide the
        milestone's alias index.
        """
        id_key_alias_value_dict, alias_key_id_value_dict, \
            neighbor_id_key_alias_value_dict = self._make_milestone_collection()
        if neighbor_id in neighbor_id_key_alias_value_dict:
            return neighbor_id_key_alias_value_dict[neighbor_id]
        else:
            return None
        
    def get_ids(self):
        """
        Return a list of model-wide incides.
        """
        id_key_alias_value_dict, alias_key_id_value_dict, \
            neighbor_id_key_alias_value_dict = self._make_milestone_collection()
        return id_key_alias_value_dict.keys()
    
class Elber_toy_anchor(Elber_anchor):
    """
    An anchor object for representing a Voronoi cell in an Elber 
    milestoning within a toy system.
    
    Attributes
    ----------
    index : int
        The index of this anchor (cell) within the model.
    
    directory : str
        The directory (within the model's root directory) that contains
        the information and calculations for this Voronoi cell.
        
    
    
    md_directory : str or None
        The directory within the 'directory' argument above which 
        contains the MD simulation information. If None, then no MD
        is performed for this anchor.
        
    bd_directory : str or None
        The directory within the 'directory' argument above which
        contains the BD simulation information. If None, then no BD
        is performed for this anchor.
        
    production_directory : str
        The directory within the MD or BD directory above in which the
        simulations will be performed.
        
    md_output_glob : str
        A glob to select all the MD output files within the production
        directory above.
        
    name : str
        A unique name for this anchor.
        
    md : bool
        A boolean of whether MD is performed in this Voronoi cell.
        
    bd : bool
        A boolean of whether BD is performed in this Voronoi cell.
        
    endstate : bool
        A boolean of whether this is an end state or not - does it
        act as the bulk or a bound state or another state of interest?
        All end states will have kinetics calculated to all other
        end states.
        
    bulkstate : bool
        A boolean of whether this state acts as the bulk state (That
        is, the state represents a large separation distance between
        ligand and receptor.
        
    milestones : list
        A list of Milestone() objects, which are the boundaries 
        bordering this cell.
    """
    
    def __init__(self):
        self.index = 0
        self.directory = ""
        self.starting_positions = []
        self.building_directory = "building"
        self.production_directory = "prod"
        self.md_output_glob = OPENMM_ELBER_GLOB
        self.name = ""
        self.md = False
        self.endstate = False
        self.bulkstate = False
        self.milestones = []
        self.variables = {}
        return