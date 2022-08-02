"""
common_cv.py

Define collective variable superclasses (also known as milestone 
shapes) that might be used in SEEKR2 calculations.
"""

from abc import abstractmethod

import numpy as np
from abserdes import Serializer

import seekr2.modules.common_base as base
import seekr2.modules.mmvt_base as mmvt_base
import seekr2.modules.elber_base as elber_base

def create_anchor(model, anchor_index):
    """
    Make the anchor object and fill out its generic attributes.
    """
    if model.get_type() == "mmvt":
        if not model.using_toy():
            anchor = mmvt_base.MMVT_anchor()
        else:
            anchor = mmvt_base.MMVT_toy_anchor()
    elif model.get_type() == "elber":
        if not model.using_toy():
            anchor = elber_base.Elber_anchor()
        else:
            anchor = elber_base.Elber_toy_anchor()
    
    anchor.index = anchor_index
    anchor.name = "anchor_"+str(anchor_index)
    anchor.directory = anchor.name
    return anchor

def assign_state_point(state_point, model):
    """
    Given a State_point object, find which anchor it belongs in, and
    assign that anchor's name and end_state accordingly.
    """
    already_found_anchor = False
    for anchor in model.anchors:
        in_all_milestones = True
        for milestone in anchor.milestones:
            cv = model.collective_variables[milestone.cv_index]
            result = cv.check_value_within_boundary(
                state_point.location[cv.index], milestone.variables)
            if not result:
                in_all_milestones = False
                break
        
        if in_all_milestones:
            assert not already_found_anchor, \
                "State point named {} found in multiple anchors.".format(
                    state_point.name)
            
            anchor.name = state_point.name
            if state_point.name.lower() == "bulk":
                anchor.bulkstate = True
            else:
                anchor.endstate = True
                
            already_found_anchor = True
    
    if not already_found_anchor:
        raise Exception(
            "Cannot assign state point named {}: no suitable anchors found. "\
            "The point might be located on a milestone.".format(
                state_point.name))
    return

def assign_state_points(model_input, model):
    """
    
    """
    for cv_input in model_input.cv_inputs:
        for state_point in cv_input.state_points:
            state_point.expand_state_point(cv_input)
            assign_state_point(state_point, model)
    return

class State_point(Serializer):
    """
    An object that can be put into the input file to indicate the location
    of a state within a CV or Combo. The name of the anchor containing this
    point will be changed, and end_state will be set to True.
    """
    def __init__(self):
        self.name = ""
        self.location = []
        return
        
    def expand_state_point(self, cv_input):
        if isinstance(cv_input, Toy_cv_input):
            self.location = [self.location]
        else:
            self.location = [self.location]
        if isinstance(cv_input, Combo):
            old_location = self.location[0]
            self.location = []
            for inner_cv_input in cv_input.cv_inputs:
                if isinstance(inner_cv_input, Toy_cv_input):
                    self.location.append(old_location)
                else:
                    self.location = old_location
        
        return

class CV_input(Serializer):
    """
    A base class for CV_inputs.
    """
    def make_mmvt_milestone_between_two_anchors(
            self, anchor1, anchor2, input_anchor1, input_anchor2, 
            milestone_index):
        neighbor_index1 = anchor2.index
        neighbor_index2 = anchor1.index
        assert anchor1.index != anchor2.index
        milestone1 = base.Milestone()
        milestone1.index = milestone_index
        milestone1.neighbor_anchor_index = neighbor_index1
        milestone1.alias_index = len(anchor1.milestones)+1
        milestone1.cv_index = self.index
        if input_anchor1.upper_milestone_value is None:
            value = 0.5 * (input_anchor1.value + input_anchor2.value)
        else:
            value = input_anchor2.upper_milestone_value
        milestone2 = base.Milestone()
        milestone2.index = milestone_index
        milestone2.neighbor_anchor_index = neighbor_index2
        milestone2.alias_index = len(anchor2.milestones)+1
        milestone2.cv_index = self.index
        if input_anchor2.lower_milestone_value is None:
            value = 0.5 * (input_anchor1.value + input_anchor2.value)
        else:
            value = input_anchor2.lower_milestone_value
        # assign k
        milestone1.variables = {"k": 1.0, "value": value}
        milestone2.variables = {"k": -1.0, "value": value}
        for milestone in anchor1.milestones:
            if milestone.neighbor_anchor_index == neighbor_index1:
                # Then anchor1 already has this milestone added to it
                return milestone_index
            
        anchor1.milestones.append(milestone1)
        anchor2.milestones.append(milestone2)
        milestone_index += 1
        return milestone_index
    
class CV_anchor(Serializer):
    """
    A base class for CV anchors.
    """
    pass

class Spherical_cv_anchor(CV_anchor):
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
        self.name = None
        self.radius = 0.0
        self.lower_milestone_radius = None
        self.upper_milestone_radius = None
        self.starting_amber_params = None
        self.starting_forcefield_params = None
        self.bound_state = False
        self.bulk_anchor = False
        self.connection_flags = []
        return
        
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
        
class Spherical_cv_input(CV_input):
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
        self.variable_name = "r"
        self.state_points = []
        return
        
    def check(self):
        """
        Check user inputs to ensure they have been entered properly.
        """
        
        last_radius = -1e9
        found_bulk_anchor = False
        for i, input_anchor in enumerate(self.input_anchors):
            assert input_anchor.__class__.__name__ == "Spherical_cv_anchor"
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
        
        assert len(self.group1) > 0, "Any input CV groups must contain atoms."
        assert len(self.group2) > 0, "Any input CV groups must contain atoms."
        return
    
    def make_mmvt_milestone_between_two_anchors(
            self, anchor1, anchor2, input_anchor1, input_anchor2, 
            milestone_index):
        neighbor_index1 = anchor2.index
        neighbor_index2 = anchor1.index        
        milestone1 = base.Milestone()
        milestone1.index = milestone_index
        milestone1.neighbor_anchor_index = neighbor_index1
        milestone1.alias_index = len(anchor1.milestones)+1
        milestone1.cv_index = self.index
        if input_anchor1.upper_milestone_radius is None:
            value = 0.5 * (input_anchor1.radius + input_anchor2.radius)
        else:
            value = input_anchor2.upper_milestone_radius
        milestone2 = base.Milestone()
        milestone2.index = milestone_index
        milestone2.neighbor_anchor_index = neighbor_index2
        milestone2.alias_index = len(anchor2.milestones)+1
        milestone2.cv_index = self.index
        if input_anchor2.lower_milestone_radius is None:
            value = 0.5 * (input_anchor1.radius + input_anchor2.radius)
        else:
            value = input_anchor2.lower_milestone_radius
        # assign k
        milestone1.variables = {"k": 1.0, "radius": value}
        milestone2.variables = {"k": -1.0, "radius": value}
        for milestone in anchor1.milestones:
            if milestone.neighbor_anchor_index == neighbor_index1:
                # Then anchor1 already has this milestone added to it
                return milestone_index
        anchor1.milestones.append(milestone1)
        anchor2.milestones.append(milestone2)
        milestone_index += 1
        return milestone_index
    
class Tiwary_cv_anchor(CV_anchor):
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
        self.name = None
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

class Tiwary_cv_input(CV_input):
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
        self.variable_name = "v"
        self.state_points = []
        return
        
    def check(self):
        """
        Check user inputs to ensure they have been entered properly.
        """
        
        last_value = -1e9
        num_bulk_anchors = 0
        for i, input_anchor in enumerate(self.input_anchors):
            assert input_anchor.__class__.__name__ == "Tiwary_cv_anchor"
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
    
    def make_mmvt_milestone_between_two_anchors(
            self, anchor1, anchor2, input_anchor1, input_anchor2, 
            milestone_index):
        milestone_index = super(Tiwary_cv_input, self)\
            .make_mmvt_milestone_between_two_anchors(
                anchor1, anchor2, input_anchor1, input_anchor2, milestone_index)
        
        for i, order_parameter_weight in enumerate(
                self.order_parameter_weights):
            key = "c{}".format(i)
            anchor1.milestones[-1].variables[key] = order_parameter_weight
            anchor2.milestones[-1].variables[key] = order_parameter_weight
        
        return milestone_index
    
class Planar_cv_anchor(CV_anchor):
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
        self.name = None
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
        
class Planar_cv_input(CV_input):
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
        A list of Planar_cv_anchor objects which specify inputs for
        the planar anchors.
    """
    
    def __init__(self):
        self.index = 0
        self.start_group = []
        self.end_group = []
        self.mobile_group = []
        self.input_anchors = []
        self.variable_name = "v"
        self.state_points = []
        return
        
    def check(self):
        """
        Check user inputs to ensure they have been entered properly.
        """
        
        last_value = -1e9
        num_bulk_anchors = 0
        for i, input_anchor in enumerate(self.input_anchors):
            assert input_anchor.__class__.__name__ == "Planar_cv_anchor"
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
        
        assert len(self.start_group) > 0, "Any input CV groups must contain atoms."
        assert len(self.end_group) > 0, "Any input CV groups must contain atoms."
        assert len(self.mobile_group) > 0, "Any input CV groups must contain atoms."
        return
    
    def make_mmvt_milestone_between_two_anchors(
            self, anchor1, anchor2, input_anchor1, input_anchor2, 
            milestone_index):
        milestone_index = super(Planar_cv_input, self)\
            .make_mmvt_milestone_between_two_anchors(
                anchor1, anchor2, input_anchor1, input_anchor2, milestone_index)
        return milestone_index
    
class RMSD_cv_anchor(CV_anchor):
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
        self.name = None
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
        
class RMSD_cv_input(CV_input):
    """
    Inputs by the user resulting in RMSD anchors with milestones and the 
    collective variable (CV).
    
    Attributes:
    -----------
    index : int
        The index of this CV input object in the Model_input object.
        
    group : list
        A list of ints representing atom indices whose RMSD will be 
        computed.
        
    ref_structure : str
        A PDB file name for the structure that will be the reference 
        for the RMSD calculation.
        
    input_anchors : list
        A list of Spherical_cv_anchor objects which specify inputs for
        the spherical anchors.
    """
    
    def __init__(self):
        self.index = 0
        self.group = []
        self.ref_structure = ""
        self.input_anchors = []
        self.variable_name = "v"
        self.state_points = []
        return
        
    def check(self):
        """
        Check user inputs to ensure they have been entered properly.
        """
        
        last_value = -1e9
        found_bulk_anchor = False
        for i, input_anchor in enumerate(self.input_anchors):
            assert input_anchor.__class__.__name__ == "RMSD_cv_anchor"
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
            
            #assert not input_anchor.bulk_anchor, "An RMSD CV must not have "\
            #    "a bulk anchor."
                    
            assert len(self.input_anchors) > 1, "A CV must contain "\
                "more than one anchor."
        
        assert len(self.group) > 0, "Any input CV groups must contain atoms."
        return
    
    def make_mmvt_milestone_between_two_anchors(
            self, anchor1, anchor2, input_anchor1, input_anchor2, 
            milestone_index):
        milestone_index = super(RMSD_cv_input, self)\
            .make_mmvt_milestone_between_two_anchors(
                anchor1, anchor2, input_anchor1, input_anchor2, milestone_index)
        return milestone_index
    
class Toy_cv_anchor(CV_anchor):
    """
    This object represents an anchor within a Toy
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
        self.name = None
        self.value = 0.0
        self.lower_milestone_value = None
        self.upper_milestone_value = None
        self.starting_positions = None
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
        
class Toy_cv_input(CV_input):
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
        self.groups = []
        self.cv_expression = None
        self.openmm_expression = None
        self.restraining_expression = None
        self.input_anchors = []
        self.variable_name = "v"
        self.state_points = []
        return
        
    def check(self):
        """
        Check user inputs to ensure they have been entered properly.
        """
        
        last_value = -1e9
        found_bulk_anchor = False
        for i, input_anchor in enumerate(self.input_anchors):
            assert input_anchor.__class__.__name__ == "Toy_cv_anchor"
            value = input_anchor.value
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
        
        for group in self.groups:
            assert len(group) > 0, \
                "Any input CV groups must contain atoms."
        return
    
    def make_mmvt_milestone_between_two_anchors(
            self, anchor1, anchor2, input_anchor1, input_anchor2, 
            milestone_index):
        milestone_index = super(Toy_cv_input, self)\
            .make_mmvt_milestone_between_two_anchors(
                anchor1, anchor2, input_anchor1, input_anchor2, milestone_index)
        return milestone_index
    
class Combo(Serializer):
    """
    A combo superclass - a way to combine CVs in a multidimensional
    arrangement
    """
    def __init__(self):
        self.cv_inputs = []
        self.state_points = []
        return
    
    @abstractmethod
    def make_anchors(self):
        raise NotImplementedError(
            "Combo superclass does not define make_anchors().")

class Grid_combo(Combo):
    """
    An object for the input to define when input CVs should be combined to
    make a multidimensional anchor - an anchor with more than one dimension
    of milestones.
    """
    def __init__(self):
        super(Grid_combo, self).__init__()
        return
    
    def make_anchors(self, model, anchor_index, milestone_index, 
                     connection_flag_dict, associated_input_anchor):
        assert model.get_type() == "mmvt", \
            "Only MMVT models may use combos."
        anchors = []
        num_anchors = 1
        cv_widths = []
        for cv_input in self.cv_inputs:
            cv_width = len(cv_input.input_anchors)
            cv_widths.append(cv_width)
            num_anchors *= cv_width
            for j, input_anchor in enumerate(cv_input.input_anchors):
                input_anchor.check(j, cv_input)
        
        for i in range(num_anchors):
            anchor = create_anchor(model, anchor_index)
            input_anchor_index_list = []
            input_anchor_list = []
            
            divider = 1
            bulkstate = False
            endstate = False
            for j, cv_input in enumerate(self.cv_inputs):
                input_anchor_index = (i // divider) % cv_widths[j]
                input_anchor_index_list.append(input_anchor_index)
                divider *= cv_widths[j]
                input_anchor = cv_input.input_anchors[input_anchor_index]
                input_anchor_list.append(input_anchor)
                if input_anchor.bulk_anchor:
                    bulkstate = True
                if input_anchor.bound_state:
                    endstate = True                
                # TODO: is there a better way to assign the 
                #  associated_input_anchor object? For now, we are assuming
                #  that the assignment is arbitrary
                associated_input_anchor[anchor.index] = input_anchor
                # Grid objects with assigned starting positions are meaningless
                if anchor.__class__.__name__ == "MMVT_toy_anchor":
                    starting_positions = input_anchor.starting_positions
                    assert starting_positions is None, \
                        "Grid combo cv's with anchors cannot have starting "\
                        "positions. Use HIDR to assign starting positions."
                else:
                    starting_positions = input_anchor.starting_amber_params
                    if anchor.amber_params is None:
                        anchor.amber_params = starting_positions
                    else:
                        assert input_anchor.starting_amber_params is not None, \
                            "All input_anchors that apply to the same anchor "\
                            "must have matching amber_params."
                        assert (anchor.amber_params.prmtop_filename \
                            == input_anchor.starting_amber_params\
                            .prmtop_filename) and (anchor.amber_params.\
                            box_vectors == input_anchor.starting_amber_params.\
                            box_vectors) and (anchor.amber_params.\
                            pdb_coordinates_filename == input_anchor.\
                            starting_amber_params.pdb_coordinates_filename), \
                            "Multiple input anchors must have matching "\
                            "starting_amber_params if they apply to the same "\
                            "anchor."
                            
                
                variable_name = "{}_{}".format(cv_input.variable_name, 
                                           cv_input.index)
                variable_value = input_anchor.get_variable_value()
                anchor.variables[variable_name] = variable_value
                
            # Assign md, bulkstate, and endstate
            if not bulkstate:
                anchor.md = True
            else:
                anchor.md = False
                anchor.bulkstate = True
                
            if endstate:
                anchor.endstate = True
            # Make the variable names
            
            # Handle the connection flags
            assert len(input_anchor.connection_flags) == 0, \
                "Connection flags within combos not currently supported."
            if input_anchor.bulk_anchor:
                connection_flag_dict["bulk"].append(anchor)
            anchors.append(anchor)
            anchor_index += 1
        
        for i, anchor in enumerate(anchors):
            input_anchor_index_list = []
            input_anchor_list = []
            divider = 1
            for j, cv_input in enumerate(self.cv_inputs):
                input_anchor_index = (i // divider) % cv_widths[j]
                input_anchor_index_list.append(input_anchor_index)
                input_anchor = cv_input.input_anchors[input_anchor_index]
                this_input_anchor = cv_input.input_anchors[input_anchor_index]
                
                if input_anchor_index < len(cv_input.input_anchors)-1:
                    upper_anchor_index = i + divider
                    upper_anchor = anchors[upper_anchor_index]
                    upper_input_anchor = cv_input.input_anchors[
                        input_anchor_index+1]
                    milestone_index \
                        = cv_input.make_mmvt_milestone_between_two_anchors(
                            anchor, upper_anchor, this_input_anchor, 
                            upper_input_anchor, milestone_index)
                        
                if input_anchor_index > 0:
                    lower_anchor_index = i - divider
                    lower_anchor = anchors[lower_anchor_index]
                    lower_input_anchor = cv_input.input_anchors[
                        input_anchor_index-1]
                    milestone_index \
                        = cv_input.make_mmvt_milestone_between_two_anchors(
                            lower_anchor, anchor, lower_input_anchor, 
                            this_input_anchor, milestone_index)
                
                divider *= cv_widths[j]
            
        return anchors, anchor_index, milestone_index, connection_flag_dict,\
            associated_input_anchor
        
        