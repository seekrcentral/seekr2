"""
elber/collective_variables.py

Define any type of collective variable (or milestone shape) that might
be used in an Elber milestoning calculation.
"""
import seekr2.common.base as base
import seekr2.elber.base as elber_base
import seekr2.libraries.serializer.serializer as serializer

def make_elber_spherical_cv_object(spherical_cv_input, index):
    """
    Create a SphericalCV object to be placed into the Model.
    """
    groups = [spherical_cv_input.group1, spherical_cv_input.group2]
    cv = elber_base.Elber_spherical_CV(index, groups)
    return cv
    
    
def make_elber_milestoning_objects_spherical(
        spherical_cv_input, milestone_alias, milestone_index, 
        index, input_anchors):
    """
    Make a set of 1-dimensional spherical Milestone objects to be put into
    an Anchor, and eventually, the Model. Each anchor needs 3 milestones.
    """
    milestones = []
    num_anchors = len(input_anchors)
    if index > 0:
        neighbor_index = index - 1
        milestone1 = base.Milestone()
        milestone1.index = milestone_index-1
        milestone1.neighbor_anchor_index = neighbor_index
        milestone1.alias_index = milestone_alias
        milestone1.cv_index = spherical_cv_input.index
        radius = input_anchors[neighbor_index].radius
        milestone1.variables = {"k": -1.0, "radius": radius}
        milestone_alias += 1
        #milestone_index += 1
        milestones.append(milestone1)
    
    milestone2 = base.Milestone()
    milestone2.index = milestone_index
    milestone2.neighbor_anchor_index = None
    milestone2.alias_index = milestone_alias
    milestone2.cv_index = spherical_cv_input.index
    radius = input_anchors[index].radius
    milestone2.variables = {"k": 1.0, "radius": radius}
    milestones.append(milestone2)
    
    if index < num_anchors-1:
        neighbor_index = index + 1
        milestone3 = base.Milestone()
        milestone3.index = milestone_index+1
        milestone3.neighbor_anchor_index = neighbor_index
        milestone3.alias_index = milestone_alias+1
        milestone3.cv_index = spherical_cv_input.index
        radius = input_anchors[neighbor_index].radius
        milestone3.variables = {"k": 1.0, "radius": radius}
        milestones.append(milestone3)
    
    return milestones, milestone_alias, milestone_index+1
