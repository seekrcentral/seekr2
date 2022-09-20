"""
mmvt_cv.py

Define any type of collective variable (or milestone shape) that might
be used in an MMVT calculation.
"""
import os
import shutil

import seekr2.modules.common_base as base
import seekr2.modules.mmvt_base as mmvt_base
import seekr2.modules.common_cv as common_cv

def make_mmvt_spherical_cv_object(spherical_cv_input, index):
    """
    Create a SphericalCV object to be placed into the Model.
    """
    group1 = base.parse_xml_list(spherical_cv_input.group1)
    group2 = base.parse_xml_list(spherical_cv_input.group2)
    groups = [group1, group2]
    cv = mmvt_base.MMVT_spherical_CV(index, groups)
    return cv
    
def make_mmvt_tiwary_cv_object(tiwary_cv_input, index):
    """
    Create a Tiwary CV object to be placed into the Model.
    """
    
    cv = mmvt_base.MMVT_tiwary_CV(
        index, tiwary_cv_input.order_parameters, 
        tiwary_cv_input.order_parameter_weights)
    return cv

def make_mmvt_planar_cv_object(planar_cv_input, index):
    """
    Create a PlanarCV object to be placed into the Model.
    """
    start_group = base.parse_xml_list(planar_cv_input.start_group)
    end_group = base.parse_xml_list(planar_cv_input.end_group)
    mobile_group = base.parse_xml_list(planar_cv_input.mobile_group)
    cv = mmvt_base.MMVT_planar_CV(index, start_group, end_group, mobile_group)
    return cv

def make_mmvt_RMSD_cv_object(RMSD_cv_input, index, root_directory):
    """
    Create a RMSD CV object to be placed into the Model.
    """
    RMSD_ref_pdb = "rmsd_reference_cv_{}.pdb".format(index)
    group = base.parse_xml_list(RMSD_cv_input.group)
    absolute_RMSD_ref_pdb = os.path.join(root_directory, RMSD_ref_pdb)
    shutil.copyfile(RMSD_cv_input.ref_structure, absolute_RMSD_ref_pdb)
    cv = mmvt_base.MMVT_RMSD_CV(index, group, RMSD_ref_pdb)
    return cv

def make_mmvt_closest_pair_cv_object(closest_pair_cv_input, index, root_directory):
    """
    Create a closest pair CV object to be placed into the Model.
    """
    group1 = base.parse_xml_list(closest_pair_cv_input.group1)
    group2 = base.parse_xml_list(closest_pair_cv_input.group2)
    groups = [group1, group2]
    cv = mmvt_base.MMVT_closest_pair_CV(index, groups)
    return cv

def make_mmvt_external_cv_object(external_cv_input, index):
    """
    Create a SphericalCV object to be placed into the Model.
    """
    
    groups = []
    for group in external_cv_input.groups:
        groups.append(base.parse_xml_list(group))
    
    cv = mmvt_base.MMVT_external_CV(index, groups)
    cv.cv_expression = external_cv_input.cv_expression
    assert cv.cv_expression is not None, \
        "A CV expression is required to define MMVT milestones."
    if external_cv_input.openmm_expression is None:
        cv.openmm_expression = "step(k*("+cv.cv_expression+" - value))"
    else:
        cv.openmm_expression = external_cv_input.openmm_expression
    
    if external_cv_input.restraining_expression is None:
        cv.restraining_expression = "0.5*k*("+cv.cv_expression+" - value)^2"
    else:
        cv.restraining_expression = external_cv_input.restraining_expression
    return cv
    
def make_mmvt_voronoi_cv_object(voronoi_cv_input, index, root_directory):
    """
    Create a RMSD CV object to be placed into the Model.
    """
    cv = mmvt_base.MMVT_Voronoi_CV(index)
    for i, cv_input in enumerate(voronoi_cv_input.cv_inputs):
        cv_input.check()
        if isinstance(cv_input, common_cv.Spherical_cv_input):
            child_cv = make_mmvt_spherical_cv_object(cv_input, index=i)
        elif isinstance(cv_input, common_cv.Tiwary_cv_input):
            child_cv = make_mmvt_tiwary_cv_object(cv_input, index=i)
        elif isinstance(cv_input, common_cv.Planar_cv_input):
            child_cv = make_mmvt_planar_cv_object(cv_input, index=i)
        elif isinstance(cv_input, common_cv.RMSD_cv_input):
            child_cv = make_mmvt_RMSD_cv_object(
                cv_input, index=i, root_directory=root_directory)
        if isinstance(cv_input, common_cv.Closest_pair_cv_input):
            child_cv = make_mmvt_closest_pair_cv_object(cv_input, index=i)
        elif isinstance(cv_input, common_cv.Toy_cv_input):
            child_cv = make_mmvt_external_cv_object(cv_input, index=i)
        else:
            raise Exception("CV type not available for Voronoi CV: %s" \
                            % type(cv_input))
        
        cv.child_cvs.append(child_cv)
    
    return cv
