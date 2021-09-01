"""
elber/collective_variables.py

Define any type of collective variable (or milestone shape) that might
be used in an Elber milestoning calculation.
"""

# TODO: update this code and documentation once the Elber plugin is
# fixed.

import seekr2.modules.common_base as base
import seekr2.modules.elber_base as elber_base
#import seekr2.libraries.serializer.serializer as serializer
#from abserdes import Serializer

def make_elber_spherical_cv_object(spherical_cv_input, index):
    """
    Create a SphericalCV object to be placed into the Model.
    """
    groups = [spherical_cv_input.group1, spherical_cv_input.group2]
    cv = elber_base.Elber_spherical_CV(index, groups)
    return cv
    
    
def make_elber_milestoning_objects_spherical(
        spherical_cv_input, milestone_alias, milestone_index, 
        input_anchor_index, anchor_index, input_anchors, umbrella_force_constant):
    """
    Make a set of 1-dimensional spherical Milestone objects to be put into
    an Anchor, and eventually, the Model. Each anchor needs 3 milestones.
    """
    milestones = []
    num_anchors = len(input_anchors)
    if input_anchor_index > 0:
        neighbor_index = anchor_index - 1
        neighbor_index_input = input_anchor_index - 1
        milestone1 = base.Milestone()
        milestone1.index = milestone_index-1
        milestone1.neighbor_anchor_index = neighbor_index
        #milestone1.alias_index = milestone_alias
        milestone1.alias_index = 1
        milestone1.is_source_milestone = False
        milestone1.cv_index = spherical_cv_input.index
        radius = input_anchors[neighbor_index_input].radius
        milestone1.variables = {"k": umbrella_force_constant, "radius": radius}
        #milestone_alias += 1
        #milestone_index += 1
        milestones.append(milestone1)
    
    milestone2 = base.Milestone()
    milestone2.index = milestone_index
    milestone2.neighbor_anchor_index = None
    #milestone2.alias_index = milestone_alias
    milestone2.alias_index = 2
    milestone2.is_source_milestone = True
    milestone2.cv_index = spherical_cv_input.index
    radius = input_anchors[input_anchor_index].radius
    milestone2.variables = {"k": umbrella_force_constant, "radius": radius}
    milestones.append(milestone2)
    
    if input_anchor_index < num_anchors-1:
        neighbor_index = anchor_index + 1
        neighbor_index_input = input_anchor_index + 1
        milestone3 = base.Milestone()
        milestone3.index = milestone_index+1
        milestone3.neighbor_anchor_index = neighbor_index
        #milestone3.alias_index = milestone_alias+1
        milestone3.alias_index = 3
        milestone3.is_source_milestone = False
        milestone3.cv_index = spherical_cv_input.index
        radius = input_anchors[neighbor_index_input].radius
        milestone3.variables = {"k": umbrella_force_constant, "radius": radius}
        milestones.append(milestone3)
    
    return milestones, milestone_alias, milestone_index+1
