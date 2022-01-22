"""
common_cv.py

Define collective variable superclasses (also known as milestone 
shapes) that might be used in SEEKR2 calculations.
"""

import numpy as np
from abserdes import Serializer

import seekr2.modules.mmvt_cv as mmvt_cv

class Spherical_cv_anchor(Serializer):
    """
    This object represents an anchor within the concentric spherical
    CV. Used for input purposes only.
    
    Attributes:
    -----------
    radius : float
        The radius of this spherical anchor in units of nanometers.
    
    lower_milestone_radius : float
        Optionally define the locations of the milestones for each
        anchor. This is the radius of the lower milestone.
        
    upper_milestone_radius : float
        Optionally define the locations of the milestones for each
        anchor. This is the radius of the lower milestone.
    
    starting_amber_params : Amber_params or None
        If Amber inputs are used for this anchor, this object contains
        the necessary inputs to start a new simulation.
        
    starting_forcefield_params : Forcefield_params or None
        If Forcefield XML inputs are used for this anchor, this object
        contains the necessary inputs to start a new simulation.
        
    bound_state : bool
        Whether this anchor represents the bound state of a ligand-
        receptor system.
        
    bulk_anchor : bool
        Whether this anchor acts as a bulk state of a ligand-receptor
        system.
    """
    
    def __init__(self):
        self.radius = 0.0
        self.lower_milestone_radius = None
        self.upper_milestone_radius = None
        self.starting_amber_params = None
        self.starting_forcefield_params = None
        self.bound_state = False
        self.bulk_anchor = False
        self.connection_flags = []
        
    def check(self, j, cv_input):
        if self.lower_milestone_radius is not None:
            assert j > 0, "lower_milestone_radius must be None for lowest "\
                "anchor in cv."
            assert self.lower_milestone_radius \
                == cv_input.input_anchors[j-1].upper_milestone_radius,\
                "If lower_milestone_radius is defined for anchor "\
                "{}, the anchor below (number {}).".format(j, j-1)\
                +" must have a corresponding upper_milestone_radius."
                
        if self.upper_milestone_radius is not None:
            assert j < len(cv_input.input_anchors), \
                "upper_milestone_radius must be None for highest anchor "\
                "in cv."
            assert self.upper_milestone_radius \
                == cv_input.input_anchors[j+1].lower_milestone_radius,\
                "If upper_milestone_radius is defined for anchor "\
                "{} at value {:.3f}, the anchor above ".format(j, 
                self.upper_milestone_radius)\
                +"(number {}).".format(j+1)\
                +" must have a corresponding lower_milestone_radius, "\
                "current value: {:.3f}.".format(
                    cv_input.input_anchors[j+1].lower_milestone_radius)
        return
    
    def get_variable_value(self):
        return self.radius
        
class Spherical_cv_input(Serializer):
    """
    Inputs by the user resulting in concentric spherical anchors
    with milestones and the collective variable (CV).
    
    Attributes:
    -----------
    index : int
        The index of this CV input object in the Model_input object.
        
    group1 : list
        A list of ints representing atom indices whose center of mass
        is one end of the CV distance vector.
        
    group2 : list
        A list of ints representing atom indices whose center of mass
        is the other end of the CV distance vector.
        
    input_anchors : list
        A list of Spherical_cv_anchor objects which specify inputs for
        the spherical anchors.
    """
    
    def __init__(self):
        self.index = 0
        self.group1 = []
        self.group2 = []
        self.bd_group1 = []
        self.bd_group2 = []
        self.input_anchors = []
        return
        
    def read_plain_input(self, inputs):
        """
        Read a plain input file (as opposed to an XML)
        TODO: REMOVE?
        """
        
        raise Exception("Reading a plain text file is not yet implemented. "\
                        "Only an XML CV input may be read at this time.")
        return
    
    def check(self):
        """
        Check user inputs to ensure they have been entered properly.
        """
        
        last_radius = -1e9
        found_bulk_anchor = False
        for i, input_anchor in enumerate(self.input_anchors):
            radius = input_anchor.radius
            assert radius >= 0.0, "A radius must be greater than "\
                "or equal to zero."
            assert radius > last_radius, "Each subsequent radius "\
                "argument must be greater than the last (sorted)."
            
            if input_anchor.bound_state is None:
                input_anchor.bound_state = False
            
            assert input_anchor.bound_state in [True, False], \
                "bound_state must be a boolean"
                
            if input_anchor.bulk_anchor is None:
                input_anchor.bulk_anchor = False
                
            assert input_anchor.bulk_anchor in [True, False], \
                "bulk_anchor must be a boolean"
            
            if input_anchor.bulk_anchor:
                assert not found_bulk_anchor, "Only one bulk anchor allowed "\
                    "per set of anchors in a CV."
                found_bulk_anchor = False
            else:
                assert not found_bulk_anchor, "Only the outermost anchor "\
                    "should be the bulk anchor."
            
            if i > 0:
                assert not input_anchor.bound_state, "Only the lowest"\
                    "anchor can be the bound state."
                    
            assert len(self.input_anchors) > 1, "A CV must contain "\
                "more than one anchor."
        return
    
    def make_mmvt_milestoning_objects(self, milestone_alias, milestone_index, 
                                      input_anchor_index, anchor_index):
        milestones, milestone_alias, milestone_index = \
            mmvt_cv.make_mmvt_milestoning_objects_spherical(
            self, milestone_alias, milestone_index, input_anchor_index, 
            anchor_index, self.input_anchors)
        return milestones, milestone_alias, milestone_index
    
class Tiwary_cv_anchor(Serializer):
    """
    This object represents an anchor within a Tiwary-style CV,
    which is composed of a linear sum of order parameters. Used 
    for input purposes only.
    
    Attributes:
    -----------
    value : float
        The value of the anchor (center of Voronoi cell)
    
    lower_milestone_value : float
        Optionally define the locations of the milestones for each
        anchor. This is the value of the lower milestone.
        
    upper_milestone_value : float
        Optionally define the locations of the milestones for each
        anchor. This is the value of the lower milestone.
    
    starting_amber_params : Amber_params or None
        If Amber inputs are used for this anchor, this object contains
        the necessary inputs to start a new simulation.
        
    starting_forcefield_params : Forcefield_params or None
        If Forcefield XML inputs are used for this anchor, this object
        contains the necessary inputs to start a new simulation.
        
    bound_state : bool
        Whether this anchor represents the bound state of a ligand-
        receptor system.
        
    bulk_anchor : bool
        Whether this anchor acts as a bulk state of a ligand-receptor
        system.
    """
    
    def __init__(self):
        self.value = 0.0
        self.lower_milestone_value = None
        self.upper_milestone_value = None
        self.starting_amber_params = None
        self.starting_forcefield_params = None
        self.bound_state = False
        self.bulk_anchor = False
        self.connection_flags = []
        return
        
    def check(self, j, cv_input):
        if self.lower_milestone_value is not None:
            assert j > 0, "lower_milestone_value must be None for lowest "\
                "anchor in cv."
            assert self.lower_milestone_value \
                == cv_input.input_anchors[j-1].upper_milestone_value,\
                "If lower_milestone_value is defined for anchor "\
                "{}, the anchor below (number {}).".format(j, j-1)\
                +" must have a corresponding upper_milestone_value."
                
        if self.upper_milestone_value is not None:
            assert j < len(cv_input.input_anchors), \
                "upper_milestone_value must be None for highest anchor "\
                "in cv."
            assert self.upper_milestone_value \
                == cv_input.input_anchors[j+1].lower_milestone_value,\
                "If upper_milestone_value is defined for anchor "\
                "{} at value {:.3f}, the anchor above ".format(j, 
                self.upper_milestone_value)\
                +"(number {}).".format(j+1)\
                +" must have a corresponding lower_milestone_value, "\
                "current value: {:.3f}.".format(
                    cv_input.input_anchors[j+1].lower_milestone_value)
        return
    
    def get_variable_value(self):
        return self.value

class Tiwary_cv_distance_order_parameter(Serializer):
    """
    An order parameter object for a Tiwary CV the represents the
    distance between two groups of atoms.
    """
    
    def __init__(self):
        self.index = 0
        self.group1 = []
        self.group2 = []
        return
    
    def get_openmm_expression(self, group_index):
        """
        
        """
        return "distance(g{0}, g{1})".format(group_index, group_index+1)
    
    def get_num_groups(self):
        """
        
        """
        return 2
    
    def get_group(self, index):
        """
        
        """
        if index == 0:
            return self.group1
        elif index == 1:
            return self.group2
        else:
            raise Exception("Invalid index: {}".format(index))
        return
    
    def get_value(self, coms):
        """
        
        """
        assert len(coms) == 2
        distance = np.linalg.norm(coms[1]-coms[0])
        return distance
    
class Tiwary_cv_angle_order_parameter(Serializer):
    """
    An order parameter object for a Tiwary CV the represents the
    distance between two groups of atoms.
    """
    
    def __init__(self):
        self.index = 0
        self.group1 = []
        self.group2 = []
        self.group3 = []
        return
    
    def get_openmm_expression(self, group_index):
        """
        
        """
        return "angle(g{0}, g{1}, g{2})".format(group_index, group_index+1,
                                                group_index+2)
    
    def get_num_groups(self):
        """
        
        """
        return 3
    
    def get_group(self, index):
        """
        
        """
        if index == 0:
            return self.group1
        elif index == 1:
            return self.group2
        elif index == 2:
            return self.group3
        else:
            raise Exception("Invalid index: {}".format(index))
        
    def get_value(self, coms):
        """
        
        """
        assert len(coms) == 3
        vec1 = coms[0] - coms[1]
        vec2 = coms[2] - coms[1]
        angle = np.arccos(np.dot(vec1/np.linalg.norm(vec1), 
                                 vec2/np.linalg.norm(vec2)))
        return angle
    
class Tiwary_cv_torsion_order_parameter(Serializer):
    """
    An order parameter object for a Tiwary CV the represents the
    distance between two groups of atoms.
    """
    
    def __init__(self):
        self.index = 0
        self.group1 = []
        self.group2 = []
        self.group3 = []
        self.group4 = []
        return
    
    def get_openmm_expression(self, group_index):
        return "dihedral(g{0}, g{1}, g{2}, g{3})".format(
            group_index, group_index+1, group_index+2, group_index+3)
    
    def get_num_groups(self):
        return 4
    
    def get_group(self, index):
        if index == 0:
            return self.group1
        elif index == 1:
            return self.group2
        elif index == 2:
            return self.group3
        elif index == 3:
            return self.group4
        else:
            raise Exception("Invalid index: {}".format(index))
        
    def get_value(self, coms):
        """
        
        """
        assert len(coms) == 4
        x1 = coms[0][0]; y1 = coms[0][1]; z1 = coms[0][2]
        x2 = coms[1][0]; y2 = coms[1][1]; z2 = coms[1][2]
        x3 = coms[2][0]; y3 = coms[2][1]; z3 = coms[2][2]
        x4 = coms[3][0]; y4 = coms[3][1]; z4 = coms[3][2]
        vec1 = np.array([x2-x1,y2-y1,z2-z1])
        axis = np.array([x3-x2,y3-y2,z3-z2])
        vec2 = np.array([x4-x3,y4-y3,z4-z3])
        cross1 = np.cross(vec1,axis)
        cross2 = np.cross(axis,vec2)
        cross1 = cross1/np.linalg.norm(cross1)
        cross2 = cross2/np.linalg.norm(cross2)
        x = np.dot(cross1, cross2)
        y = np.dot(np.cross(cross1, axis/np.linalg.norm(axis)), cross2)
        phi = -np.arctan2(y,x) % (2.0*np.pi)
        return phi

class Tiwary_cv_input(Serializer):
    """
    Inputs by the user resulting in Tiwary-style anchors
    with milestones and the collective variable (CV).
    
    TODO: update docstring
    
    Attributes:
    -----------
    index : int
        The index of this CV input object in the Model_input object.
        
    group1 : list
        A list of ints representing atom indices whose center of mass
        is one end of the CV distance vector.
        
    group2 : list
        A list of ints representing atom indices whose center of mass
        is the other end of the CV distance vector.
        
    input_anchors : list
        A list of Spherical_cv_anchor objects which specify inputs for
        the spherical anchors.
    """
    
    def __init__(self):
        self.index = 0
        self.order_parameters = []
        self.order_parameter_weights = []
        self.input_anchors = []
        return
        
    def read_plain_input(self, inputs):
        """
        Read a plain input file (as opposed to an XML)
        """
        
        raise Exception("Reading a plain text file is not yet implemented. "\
                        "Only an XML CV input may be read at this time.")
        return
    
    def check(self):
        """
        Check user inputs to ensure they have been entered properly.
        TODO: REMOVE?
        """
        
        last_value = -1e9
        num_bulk_anchors = 0
        for i, input_anchor in enumerate(self.input_anchors):
            value = input_anchor.value
            
            assert value > last_value, "Each subsequent value "\
                "argument must be greater than the last (sorted)."
            
            if input_anchor.bound_state is None:
                input_anchor.bound_state = False
            
            assert input_anchor.bound_state in [True, False], \
                "bound_state must be a boolean"
                
            if input_anchor.bulk_anchor is None:
                input_anchor.bulk_anchor = False
                
            if input_anchor.bulk_anchor:
                num_bulk_anchors += 1
            
            if i > 0:
                assert not input_anchor.bound_state, "Only the lowest"\
                    "anchor can be the bound state."
                    
            assert len(self.input_anchors) > 1, "A CV must contain "\
                "more than one anchor."
        
        assert num_bulk_anchors <= 1, "There cannot be more than one bulk "\
            "anchor."
        return
    
    def make_mmvt_milestoning_objects(self, milestone_alias, milestone_index, 
                                      input_anchor_index, anchor_index):
        milestones, milestone_alias, milestone_index = \
            mmvt_cv.make_mmvt_milestoning_objects_tiwary(
            self, milestone_alias, milestone_index, input_anchor_index, 
            anchor_index, self.input_anchors)
        return milestones, milestone_alias, milestone_index
    
class Planar_cv_anchor(Serializer):
    """
    This object represents an anchor within the planar incremental
    CV. Used for input purposes only.
    
    TODO: update docstring
    
    Attributes:
    -----------
    radius : float
        The radius of this spherical anchor in units of nanometers.
    
    lower_milestone_radius : float
        Optionally define the locations of the milestones for each
        anchor. This is the radius of the lower milestone.
        
    upper_milestone_radius : float
        Optionally define the locations of the milestones for each
        anchor. This is the radius of the lower milestone.
    
    starting_amber_params : Amber_params or None
        If Amber inputs are used for this anchor, this object contains
        the necessary inputs to start a new simulation.
        
    starting_forcefield_params : Forcefield_params or None
        If Forcefield XML inputs are used for this anchor, this object
        contains the necessary inputs to start a new simulation.
        
    bound_state : bool
        Whether this anchor represents the bound state of a ligand-
        receptor system.
        
    bulk_anchor : bool
        Whether this anchor acts as a bulk state of a ligand-receptor
        system.
    """
    
    def __init__(self):
        self.value = 0.0
        self.lower_milestone_value = None
        self.upper_milestone_value = None
        self.starting_amber_params = None
        self.starting_forcefield_params = None
        self.bound_state = False
        self.bulk_anchor = False
        self.connection_flags = []
        return
        
    def check(self, j, cv_input):
        if self.lower_milestone_value is not None:
            assert j > 0, "lower_milestone_value must be None for lowest "\
                "anchor in cv."
            assert self.lower_milestone_value \
                == cv_input.input_anchors[j-1].upper_milestone_value,\
                "If lower_milestone_value is defined for anchor "\
                "{}, the anchor below (number {}).".format(j, j-1)\
                +" must have a corresponding upper_milestone_value."
                
        if self.upper_milestone_value is not None:
            assert j < len(cv_input.input_anchors), \
                "upper_milestone_value must be None for highest anchor "\
                "in cv."
            assert self.upper_milestone_value \
                == cv_input.input_anchors[j+1].lower_milestone_value,\
                "If upper_milestone_value is defined for anchor "\
                "{} at value {:.3f}, the anchor above ".format(j, 
                self.upper_milestone_value)\
                +"(number {}).".format(j+1)\
                +" must have a corresponding lower_milestone_value, "\
                "current value: {:.3f}.".format(
                    cv_input.input_anchors[j+1].lower_milestone_value)
        return
    
    def get_variable_value(self):
        return self.value
        
class Planar_cv_input(Serializer):
    """
    Inputs by the user resulting in planar incremental anchors
    with milestones and the collective variable (CV).
    
    Attributes:
    -----------
    index : int
        The index of this CV input object in the Model_input object.
        
    group1 : list
        A list of ints representing atom indices whose center of mass
        is one end of the CV distance vector.
        
    group2 : list
        A list of ints representing atom indices whose center of mass
        is the other end of the CV distance vector.
        
    input_anchors : list
        A list of Spherical_cv_anchor objects which specify inputs for
        the spherical anchors.
    """
    
    def __init__(self):
        self.index = 0
        self.start_group = []
        self.end_group = []
        self.mobile_group = []
        self.input_anchors = []
        return
        
    def read_plain_input(self, inputs):
        """
        Read a plain input file (as opposed to an XML)
        TODO: REMOVE?
        """
        
        raise Exception("Reading a plain text file is not yet implemented. "\
                        "Only an XML CV input may be read at this time.")
        return
    
    def check(self):
        """
        Check user inputs to ensure they have been entered properly.
        """
        
        last_value = -1e9
        num_bulk_anchors = 0
        for i, input_anchor in enumerate(self.input_anchors):
            value = input_anchor.value
            assert value > last_value, "Each subsequent value "\
                "argument must be greater than the last (sorted)."
            
            if input_anchor.bound_state is None:
                input_anchor.bound_state = False
            
            assert input_anchor.bound_state in [True, False], \
                "bound_state must be a boolean"
                
            if input_anchor.bulk_anchor is None:
                input_anchor.bulk_anchor = False
                
            if input_anchor.bulk_anchor:
                num_bulk_anchors += 1
            
            if i > 0:
                assert not input_anchor.bound_state, "Only the lowest"\
                    "anchor can be the bound state."
                    
            assert len(self.input_anchors) > 1, "A CV must contain "\
                "more than one anchor."
        assert num_bulk_anchors <= 1, "There cannot be more than one bulk "\
            "anchor."
            
        return
    
    def make_mmvt_milestoning_objects(self, milestone_alias, milestone_index, 
                                      input_anchor_index, anchor_index):
        milestones, milestone_alias, milestone_index = \
            mmvt_cv.make_mmvt_milestoning_objects_planar(
            self, milestone_alias, milestone_index, input_anchor_index,
            anchor_index, self.input_anchors)
        return milestones, milestone_alias, milestone_index

class RMSD_cv_anchor(Serializer):
    """
    This object represents an anchor within an RMSD CV. Used for input 
    purposes only.
    
    Attributes:
    -----------
    value : float
        The value of this RMSD anchor in units of nanometers.
    
    lower_milestone_value : float
        Optionally define the locations of the milestones for each
        anchor. This is the radius of the lower milestone.
        
    upper_milestone_value : float
        Optionally define the locations of the milestones for each
        anchor. This is the radius of the lower milestone.
    
    starting_amber_params : Amber_params or None
        If Amber inputs are used for this anchor, this object contains
        the necessary inputs to start a new simulation.
        
    starting_forcefield_params : Forcefield_params or None
        If Forcefield XML inputs are used for this anchor, this object
        contains the necessary inputs to start a new simulation.
        
    bound_state : bool
        Whether this anchor represents the bound state of a ligand-
        receptor system.
        
    bulk_anchor : bool
        Whether this anchor acts as a bulk state of a ligand-receptor
        system.
    """
    
    def __init__(self):
        self.value = 0.0
        self.lower_milestone_value = None
        self.upper_milestone_value = None
        self.starting_amber_params = None
        self.starting_forcefield_params = None
        self.bound_state = False
        self.bulk_anchor = False
        self.connection_flags = []
        
    def check(self, j, cv_input):
        if self.lower_milestone_value is not None:
            assert j > 0, "lower_milestone_value must be None for lowest "\
                "anchor in cv."
            assert self.lower_milestone_value \
                == cv_input.input_anchors[j-1].upper_milestone_value,\
                "If lower_milestone_value is defined for anchor "\
                "{}, the anchor below (number {}).".format(j, j-1)\
                +" must have a corresponding upper_milestone_value."
                
        if self.upper_milestone_value is not None:
            assert j < len(cv_input.input_anchors), \
                "upper_milestone_value must be None for highest anchor "\
                "in cv."
            assert self.upper_milestone_value \
                == cv_input.input_anchors[j+1].lower_milestone_value,\
                "If upper_milestone_value is defined for anchor "\
                "{} at value {:.3f}, the anchor above ".format(j, 
                self.upper_milestone_value)\
                +"(number {}).".format(j+1)\
                +" must have a corresponding lower_milestone_value, "\
                "current value: {:.3f}.".format(
                    cv_input.input_anchors[j+1].lower_milestone_value)
        return
    
    def get_variable_value(self):
        return self.value
        
class RMSD_cv_input(Serializer):
    """
    Inputs by the user resulting in RMSD anchors with milestones and the 
    collective variable (CV).
    
    Attributes:
    -----------
    index : int
        The index of this CV input object in the Model_input object.
        
    group : list
        A list of ints representing atom indices whose center of mass
        is one end of the CV distance vector.
        
    input_anchors : list
        A list of Spherical_cv_anchor objects which specify inputs for
        the spherical anchors.
    """
    
    def __init__(self):
        self.index = 0
        self.group = []
        self.input_anchors = []
        return
        
    def read_plain_input(self, inputs):
        """
        Read a plain input file (as opposed to an XML)
        TODO: REMOVE?
        """
        
        raise Exception("Reading a plain text file is not yet implemented. "\
                        "Only an XML CV input may be read at this time.")
        return
    
    def check(self):
        """
        Check user inputs to ensure they have been entered properly.
        """
        
        last_value = -1e9
        found_bulk_anchor = False
        for i, input_anchor in enumerate(self.input_anchors):
            value = input_anchor.value
            assert value >= 0.0, "A value must be greater than "\
                "or equal to zero."
            assert value > last_value, "Each subsequent value "\
                "argument must be greater than the last (sorted)."
            
            if input_anchor.bound_state is None:
                input_anchor.bound_state = False
            
            assert input_anchor.bound_state in [True, False], \
                "bound_state must be a boolean"
                
            if input_anchor.bulk_anchor is None:
                input_anchor.bulk_anchor = False
                
            assert input_anchor.bulk_anchor in [True, False], \
                "bulk_anchor must be a boolean"
            
            assert not input_anchor.bulk_anchor, "An RMSD CV must not have "\
                "a bulk anchor."
                    
            assert len(self.input_anchors) > 1, "A CV must contain "\
                "more than one anchor."
        return
    
    def make_mmvt_milestoning_objects(self, milestone_alias, milestone_index, 
                                      input_anchor_index, anchor_index):
        milestones, milestone_alias, milestone_index = \
            mmvt_cv.make_mmvt_milestoning_objects_RMSD(
            self, milestone_alias, milestone_index, input_anchor_index, 
            anchor_index, self.input_anchors)
        return milestones, milestone_alias, milestone_index
    

class Multidimensional_cv(Serializer):
    """
    An object for the input to define when input CVs should be combined to
    make a multidimensional anchor - an anchor with more than one dimension
    of milestones.
    """
    def __init__(self):
        self.cv_inputs = []
        
        return