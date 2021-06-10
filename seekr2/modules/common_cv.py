"""
common_cv.py

Define collective variable superclasses (also known as milestone 
shapes) that might be used in SEEKR2 calculations.
"""

#import seekr2.libraries.serializer.serializer as serializer
from abserdes import Serializer

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
    
