"""
mmvt/collective_variables.py

Define any type of collective variable (or milestone shape) that might
be used in an MMVT calculation.
"""
import seekr2.modules.common_base as base
import seekr2.modules.mmvt_base as mmvt_base
import seekr2.libraries.serializer.serializer as serializer

def make_mmvt_spherical_cv_object(spherical_cv_input, index):
    """
    Create a SphericalCV object to be placed into the Model.
    """
    groups = [spherical_cv_input.group1, spherical_cv_input.group2]
    cv = mmvt_base.MMVT_spherical_CV(index, groups)
    return cv
    
    
def make_mmvt_milestoning_objects_spherical(
        spherical_cv_input, milestone_alias, milestone_index, 
        index, input_anchors):
    """
    Make a set of 1-dimensional spherical Milestone objects to be put into
    an Anchor, and eventually, the Model.
    """
    milestones = []
    num_anchors = len(input_anchors)
    if index > 0:
        neighbor_index = index - 1
        milestone1 = base.Milestone()
        milestone1.index = milestone_index
        milestone1.neighbor_anchor_index = neighbor_index
        milestone1.alias_index = milestone_alias
        milestone1.cv_index = spherical_cv_input.index
        radius = 0.5 * (input_anchors[index].radius \
                        + input_anchors[neighbor_index].radius)
        milestone1.variables = {"k": -1.0, "radius": radius}
        milestone_alias += 1
        milestone_index += 1
        milestones.append(milestone1)
        
    if index < num_anchors-1:
        neighbor_index = index + 1
        milestone2 = base.Milestone()
        milestone2.index = milestone_index
        milestone2.neighbor_anchor_index = neighbor_index
        milestone2.alias_index = milestone_alias
        milestone2.cv_index = spherical_cv_input.index
        radius = 0.5 * (input_anchors[index].radius \
                        + input_anchors[neighbor_index].radius)
        milestone2.variables = {"k": 1.0, "radius": radius}
        milestones.append(milestone2)
    
    return milestones, milestone_alias, milestone_index
