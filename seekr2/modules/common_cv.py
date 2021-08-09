"""
common_cv.py

Define collective variable superclasses (also known as milestone 
shapes) that might be used in SEEKR2 calculations.
"""

#import seekr2.libraries.serializer.serializer as serializer
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
                "{} at value {}, the anchor above ".format(j, 
                self.upper_milestone_radius)\
                +"(number {}).".format(j+1)\
                +" must have a corresponding lower_milestone_radius, "\
                "current value: {}.".format(
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
                                      anchor_index):
        milestones, milestone_alias, milestone_index = \
            mmvt_cv.make_mmvt_milestoning_objects_spherical(
            self, milestone_alias, milestone_index, anchor_index, 
            self.input_anchors)
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
                "{} at value {}, the anchor above ".format(j, 
                self.upper_milestone_value)\
                +"(number {}).".format(j+1)\
                +" must have a corresponding lower_milestone_value, "\
                "current value: {}.".format(
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
        return "distance(g{0}, g{1})".format(group_index, group_index+1)
    
    def get_num_groups(self):
        return 2
    
    def get_group(self, index):
        if index == 0:
            return self.group1
        elif index == 1:
            return self.group2
        else:
            raise Exception("Invalid index: {}".format(index))
    
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
        return "angle(g{0}, g{1}, g{2})".format(group_index, group_index+1,
                                                group_index+2)
    
    def get_num_groups(self):
        return 3
    
    def get_group(self, index):
        if index == 0:
            return self.group1
        elif index == 1:
            return self.group2
        elif index == 2:
            return self.group3
        else:
            raise Exception("Invalid index: {}".format(index))
    
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
        return "angle(g{0}, g{1}, g{2}, g{3})".format(
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
                
            assert input_anchor.bulk_anchor == False, \
                "Tiwary anchors may not have their 'bulk_anchor' property "\
                "set to True. This shape is currently not supported in "\
                "Brownian dynamics software."
            
            if i > 0:
                assert not input_anchor.bound_state, "Only the lowest"\
                    "anchor can be the bound state."
                    
            assert len(self.input_anchors) > 1, "A CV must contain "\
                "more than one anchor."
        return
    
    def make_mmvt_milestoning_objects(self, milestone_alias, milestone_index, 
                                      anchor_index):
        milestones, milestone_alias, milestone_index = \
            mmvt_cv.make_mmvt_milestoning_objects_tiwary(
            self, milestone_alias, milestone_index, anchor_index, 
            self.input_anchors)
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
                "{} at value {}, the anchor above ".format(j, 
                self.upper_milestone_value)\
                +"(number {}).".format(j+1)\
                +" must have a corresponding lower_milestone_value, "\
                "current value: {}.".format(
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
            assert value > last_value, "Each subsequent value "\
                "argument must be greater than the last (sorted)."
            
            if input_anchor.bound_state is None:
                input_anchor.bound_state = False
            
            assert input_anchor.bound_state in [True, False], \
                "bound_state must be a boolean"
                
            if input_anchor.bulk_anchor is None:
                input_anchor.bulk_anchor = False
                
            assert input_anchor.bulk_anchor == False, \
                "Planar anchors may not have their 'bulk_anchor' property "\
                "set to True. This shape is currently not supported in "\
                "Brownian dynamics software."
            
            if i > 0:
                assert not input_anchor.bound_state, "Only the lowest"\
                    "anchor can be the bound state."
                    
            assert len(self.input_anchors) > 1, "A CV must contain "\
                "more than one anchor."
        return
    
    def make_mmvt_milestoning_objects(self, milestone_alias, milestone_index, 
                                      anchor_index):
        milestones, milestone_alias, milestone_index = \
            mmvt_cv.make_mmvt_milestoning_objects_planar(
            self, milestone_alias, milestone_index, anchor_index, 
            self.input_anchors)
        return milestones, milestone_alias, milestone_index
    