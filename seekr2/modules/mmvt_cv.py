"""
mmvt_cv.py

Define any type of collective variable (or milestone shape) that might
be used in an MMVT calculation.
"""
import seekr2.modules.common_base as base
import seekr2.modules.mmvt_base as mmvt_base

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
        if input_anchors[index].lower_milestone_radius is None:
            radius = 0.5 * (input_anchors[index].radius \
                            + input_anchors[neighbor_index].radius)
        else:
            radius = input_anchors[index].lower_milestone_radius
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
        if input_anchors[index].upper_milestone_radius is None:
            radius = 0.5 * (input_anchors[index].radius \
                            + input_anchors[neighbor_index].radius)
        else:
            radius = input_anchors[index].upper_milestone_radius
        milestone2.variables = {"k": 1.0, "radius": radius}
        milestones.append(milestone2)
    
    return milestones, milestone_alias, milestone_index

def make_mmvt_tiwary_cv_object(tiwary_cv_input, index):
    """
    Create a Tiwary CV object to be placed into the Model.
    """
    
    #groups = [spherical_cv_input.group1, spherical_cv_input.group2]
    #cv = mmvt_base.MMVT_spherical_CV(index, groups)
    cv = mmvt_base.MMVT_tiwary_CV(
        index, tiwary_cv_input.order_parameters, 
        tiwary_cv_input.order_parameter_weights)
    return cv

def make_mmvt_milestoning_objects_tiwary(
        tiwary_cv_input, milestone_alias, milestone_index, 
        index, input_anchors):
    """
    Make a set of 1-dimensional Tiwary Milestone objects to be put into
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
        milestone1.cv_index = tiwary_cv_input.index
        if input_anchors[index].lower_milestone_value is None:
            value = 0.5 * (input_anchors[index].value \
                            + input_anchors[neighbor_index].value)
        else:
            value = input_anchors[index].lower_milestone_value
        milestone1.variables = {"k": -1.0, "value": value}
        # Loop thru weights
        for i, order_parameter_weight in enumerate(
                tiwary_cv_input.order_parameter_weights):
            key = "c{}".format(i)
            milestone1.variables[key] = order_parameter_weight
            
        milestone_alias += 1
        milestone_index += 1
        milestones.append(milestone1)
        
    if index < num_anchors-1:
        neighbor_index = index + 1
        milestone2 = base.Milestone()
        milestone2.index = milestone_index
        milestone2.neighbor_anchor_index = neighbor_index
        milestone2.alias_index = milestone_alias
        milestone2.cv_index = tiwary_cv_input.index
        if input_anchors[index].upper_milestone_value is None:
            value = 0.5 * (input_anchors[index].value \
                            + input_anchors[neighbor_index].value)
        else:
            value = input_anchors[index].upper_milestone_value
        milestone2.variables = {"k": 1.0, "value": value}
        for i, order_parameter_weight in enumerate(
                tiwary_cv_input.order_parameter_weights):
            key = "c{}".format(i)
            milestone2.variables[key] = order_parameter_weight
            
        milestones.append(milestone2)
    
    return milestones, milestone_alias, milestone_index

def make_mmvt_planar_cv_object(planar_cv_input, index):
    """
    Create a PlanarCV object to be placed into the Model.
    """
    cv = mmvt_base.MMVT_planar_CV(
        index, planar_cv_input.start_group, planar_cv_input.end_group, 
              planar_cv_input.mobile_group)
    return cv
    
def make_mmvt_milestoning_objects_planar(
        planar_cv_input, milestone_alias, milestone_index, 
        index, input_anchors):
    """
    Make a set of 1-dimensional planar Milestone objects to be put into
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
        milestone1.cv_index = planar_cv_input.index
        if input_anchors[index].lower_milestone_value is None:
            value = 0.5 * (input_anchors[index].value \
                            + input_anchors[neighbor_index].value)
        else:
            value = input_anchors[index].lower_milestone_value
        milestone1.variables = {"k": -1.0, "value": value}
        milestone_alias += 1
        milestone_index += 1
        milestones.append(milestone1)
        
    if index < num_anchors-1:
        neighbor_index = index + 1
        milestone2 = base.Milestone()
        milestone2.index = milestone_index
        milestone2.neighbor_anchor_index = neighbor_index
        milestone2.alias_index = milestone_alias
        milestone2.cv_index = planar_cv_input.index
        if input_anchors[index].upper_milestone_value is None:
            value = 0.5 * (input_anchors[index].value \
                            + input_anchors[neighbor_index].value)
        else:
            value = input_anchors[index].upper_milestone_value
        milestone2.variables = {"k": 1.0, "value": value}
        milestones.append(milestone2)
    
    return milestones, milestone_alias, milestone_index
