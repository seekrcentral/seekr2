"""
mmvt_cv_base.py

Base classes, objects, and constants used when handling CVs of the
MMVT calculations.
"""

import mdtraj

from abserdes import Serializer

OPENMMVT_BASENAME = "mmvt"
OPENMMVT_EXTENSION = "out"
OPENMMVT_GLOB = "%s*.%s" % (OPENMMVT_BASENAME, OPENMMVT_EXTENSION)
NAMDMMVT_BASENAME = "namdmmvt"
NAMDMMVT_EXTENSION = "out"
NAMDMMVT_GLOB = "%s*.%s*" % (NAMDMMVT_BASENAME, NAMDMMVT_EXTENSION)

def traj_center_of_mass(traj):
    """
    Returns a center of mass array by frames for a traj. Avoids NaN
    errors that come when using mdtraj's .
    """
    if traj.xyz.shape[1] == 1:
        com_array = traj.xyz[:,0,:]
    else:
        com_array = mdtraj.compute_center_of_mass(traj)
    return com_array

class MMVT_settings(Serializer):
    """
    Settings that are specific to an MMVT calculation.
    
    Attributes:
    -----------
    num_production_steps : int
        The number of steps to take within a given MMVT production
        run for a Voronoi cell.
        
    energy_reporter_interval : int or None, Default None
        The interval to report system state information, including
        energy. If set to None, then this information is never printed.
        
    restart_checkpoint_interval : int or None, Default None
        The interval to write a restart checkpoint for the system to
        be easily restarted. None means do not write backups.
        
    trajectory_reporter_interval : int or None, Default None
        The interval to write frames of the trajectory to file.
        None means do not write to a trajectory file.
    """
    
    def __init__(self):
        self.num_production_steps = 30000
        self.energy_reporter_interval = None
        self.trajectory_reporter_interval = None
        self.restart_checkpoint_interval = None
        
class MMVT_collective_variable(Serializer):
    """
    Collective variables represent the function of system positions
    and velocities along with the MMVT cells and boundaries can be
    defined.
    
    Attributes:
    -----------
    index : int
        Every collective variable needs an index so that it may be
        quickly and easily referenced by one of the many milestones
        in the model.
        
    groups : list
        A list of lists of integers. The length of the outer list is
        equal to self.num_groups. The inner lists contain integer
        values representing the indices of atoms in that group.
        
    """
    def __init__(self, index, groups):
        self.index = index
        self.groups = groups
        return

    def __name__(self):
        return "mmvt_base_CV"
    
    def make_boundary_force(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
    
    def make_restraining_force(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable restraining force.")
    
    def make_cv_force(self, alias_id):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable force.")
    
    def make_voronoi_cv_boundary_force(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
    
    def update_voronoi_cv_boundary_force(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
    
    def make_namd_colvar_string(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
    
    def add_groups(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
    
    def add_parameters(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
        
    def add_groups_and_variables(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
    
    def update_groups_and_variables(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
        
    def get_variable_values_list(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
        
    def get_namd_evaluation_string(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
    
    def get_mdtraj_cv_value(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
    
    def check_mdtraj_within_boundary(self, parmed_structure, 
                                               milestone_variables):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
    
    def get_openmm_context_cv_value(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
    
    def check_openmm_context_within_boundary(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
    
    def check_value_within_boundary(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
    
    def check_mdtraj_close_to_boundary(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
    
    def get_atom_groups(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
        
    def get_variable_values(self):
        raise Exception("This base class cannot be used for creating a "\
                        "collective variable boundary definition.")
        
class MMVT_anchor(Serializer):
    """
    An anchor object for representing a Voronoi cell in an MMVT 
    calculation.
    
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
    
    production_directory : str
        The directory within the MD or BD directory above in which the
        simulations will be performed.
        
    building_directory : str
        The directory in which the prepared parameter and starting 
        structure files are stored.
        
    md_output_glob : str
        A glob to select all the MD output files within the production
        directory above.
        
    name : str
        A unique name for this anchor.
        
    md : bool
        A boolean of whether MD is performed in this Voronoi cell.
        
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
        self.production_directory = "prod"
        self.building_directory = "building"
        self.md_output_glob = OPENMMVT_GLOB
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
        index of the milestone.
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
        return list(id_key_alias_value_dict.keys())
    
class MMVT_toy_anchor(MMVT_anchor):
    """
    An anchor object for representing a Voronoi cell in an MMVT 
    calculation within a toy system.
    
    Attributes
    ----------
    index : int
        The index of this anchor (cell) within the model.
    
    directory : str
        The directory (within the model's root directory) that contains
        the information and calculations for this Voronoi cell.
        
    starting_positions : list
        A list of lists for each particle's starting positions
    
    production_directory : str
        The directory within the MD or BD directory above in which the
        simulations will be performed.
        
    building_directory : str
        The directory in which the prepared parameter and starting 
        structure files are stored.
        
    md_output_glob : str
        A glob to select all the MD output files within the production
        directory above.
        
    name : str
        A unique name for this anchor.
        
    md : bool
        A boolean of whether MD is performed in this Voronoi cell.
        
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
        self.production_directory = "prod"
        self.building_directory = "building"
        self.md_output_glob = OPENMMVT_GLOB
        self.name = ""
        self.md = False
        self.endstate = False
        self.bulkstate = False
        self.milestones = []
        self.variables = {}
        return