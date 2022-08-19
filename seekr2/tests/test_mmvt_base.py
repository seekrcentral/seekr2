"""
test_mmvt_base.py
"""

import os
import shutil

import numpy as np

from seekr2.modules import mmvt_base
from seekr2.modules import mmvt_sim_openmm

def test_MMVT_anchor_id_from_alias(host_guest_mmvt_model):
    """
    Test the anchor's ability to return information about itself and its
    milestones - id_from_alias.
    """
    anchor0 = host_guest_mmvt_model.anchors[0]
    anchor1 = host_guest_mmvt_model.anchors[1]
    # test id_from_alias
    assert anchor0.id_from_alias(1) == 0
    assert anchor1.id_from_alias(1) == 0
    assert anchor1.id_from_alias(2) == 1
    return

def test_MMVT_anchor_alias_from_id(host_guest_mmvt_model):
    """
    Test the anchor's ability to return information about itself and its
    milestones - alias_from_id.
    """
    anchor0 = host_guest_mmvt_model.anchors[0]
    anchor1 = host_guest_mmvt_model.anchors[1]
    # test id_from_alias
    assert anchor0.alias_from_id(0) == 1
    assert anchor1.alias_from_id(0) == 1
    assert anchor1.alias_from_id(1) == 2
    return

def test_MMVT_anchor_alias_from_neighbor_id(host_guest_mmvt_model):
    """
    Test the anchor's ability to return information about itself and its
    milestones - alias_from_neighbor_id.
    """
    anchor0 = host_guest_mmvt_model.anchors[0]
    anchor1 = host_guest_mmvt_model.anchors[1]
    # test id_from_alias
    assert anchor0.alias_from_neighbor_id(1) == 1
    assert anchor1.alias_from_neighbor_id(0) == 1
    assert anchor1.alias_from_neighbor_id(2) == 2
    return

def test_MMVT_anchor_get_ids(host_guest_mmvt_model):
    """
    Test the anchor's ability to return information about itself and its
    milestones - id_from_alias.
    """
    anchor0 = host_guest_mmvt_model.anchors[0]
    anchor1 = host_guest_mmvt_model.anchors[1]
    # test id_from_alias
    assert anchor0.get_ids() == [0]
    assert anchor1.get_ids() == [0, 1]
    return

def test_check_openmm_context_within_boundary(host_guest_mmvt_model, tmp_path):
    """
    Test whether the check can find systems that exist outside the proper
    Voronoi cells.
    """
    # Test check success: system starting within Voronoi cells
    anchor = host_guest_mmvt_model.anchors[0]
    output_file = os.path.join(tmp_path, "output.txt")
    host_guest_mmvt_model.openmm_settings.cuda_platform_settings = None
    host_guest_mmvt_model.openmm_settings.reference_platform = True
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(
        host_guest_mmvt_model, anchor, output_file)
    context = my_sim_openmm.simulation.context
    
    final_result = True
    for milestone in anchor.milestones:
        cv = host_guest_mmvt_model.collective_variables[milestone.cv_index]
        result = cv.check_openmm_context_within_boundary(
            context, milestone.variables, verbose=True)
        final_result = final_result and result
        
    assert final_result
    
    # Test check failure(s): system outside Voronoi cells
    anchor_pdb1 \
        = host_guest_mmvt_model.anchors[0].amber_params.pdb_coordinates_filename
    anchor_pdb2 \
        = host_guest_mmvt_model.anchors[1].amber_params.pdb_coordinates_filename
    anchor_pdb_src_path1 = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, 
        host_guest_mmvt_model.anchors[0].directory,
        host_guest_mmvt_model.anchors[0].building_directory, 
        anchor_pdb1)
    anchor_pdb_src_path2 = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, 
        host_guest_mmvt_model.anchors[1].directory,
        host_guest_mmvt_model.anchors[1].building_directory, 
        anchor_pdb2)
    anchor_pdb_dest_path1 = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, 
        host_guest_mmvt_model.anchors[0].directory,
        host_guest_mmvt_model.anchors[0].building_directory, 
        anchor_pdb2)
    anchor_pdb_dest_path2 = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, 
        host_guest_mmvt_model.anchors[1].directory,
        host_guest_mmvt_model.anchors[1].building_directory, 
        anchor_pdb1)
    shutil.copyfile(anchor_pdb_src_path1, anchor_pdb_dest_path2)
    shutil.copyfile(anchor_pdb_src_path2, anchor_pdb_dest_path1)
    host_guest_mmvt_model.anchors[0].amber_params.pdb_coordinates_filename \
        = anchor_pdb2
    host_guest_mmvt_model.anchors[1].amber_params.pdb_coordinates_filename \
        = anchor_pdb1
        
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(
        host_guest_mmvt_model, anchor, output_file)
    context = my_sim_openmm.simulation.context
    
    final_result = True
    for milestone in anchor.milestones:
        cv = host_guest_mmvt_model.collective_variables[milestone.cv_index]
        result = cv.check_openmm_context_within_boundary(
            context, milestone.variables, verbose=True)
        final_result = final_result and result
        
    assert not final_result
    
    return

def test_check_toy_openmm_context_within_boundary(toy_mmvt_model, tmp_path):
    """
    Test whether the check can find systems that exist outside the proper
    Voronoi cells.
    """
    # Test check success: system starting within Voronoi cells
    anchor = toy_mmvt_model.anchors[0]
    output_file = os.path.join(tmp_path, "output.txt")
    toy_mmvt_model.openmm_settings.cuda_platform_settings = None
    toy_mmvt_model.openmm_settings.reference_platform = True
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(
        toy_mmvt_model, anchor, output_file)
    context = my_sim_openmm.simulation.context
    
    final_result = True
    for milestone in anchor.milestones:
        cv = toy_mmvt_model.collective_variables[milestone.cv_index]
        result = cv.check_openmm_context_within_boundary(
            context, milestone.variables, verbose=True)
        final_result = final_result and result
        
    assert final_result
    
    # Test check failure(s): system outside Voronoi cells
    toy_mmvt_model.anchors[0].starting_positions = np.array(
        [[[0.0, -0.7, 0.0]], [[0.2, -0.5, 0.0]]])
    toy_mmvt_model.anchors[1].starting_positions = np.array(
        [[[0.0, -0.3, 0.0]]])
    
    anchor = toy_mmvt_model.anchors[0]
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(
        toy_mmvt_model, anchor, output_file, frame=1)
    context = my_sim_openmm.simulation.context
    
    final_result = True
    for milestone in anchor.milestones:
        cv = toy_mmvt_model.collective_variables[milestone.cv_index]
        result = cv.check_openmm_context_within_boundary(
            context, milestone.variables, verbose=True)
        final_result = final_result and result
        
    assert not final_result
    
    anchor = toy_mmvt_model.anchors[1]
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(
        toy_mmvt_model, anchor, output_file)
    context = my_sim_openmm.simulation.context
    
    final_result = True
    for milestone in anchor.milestones:
        cv = toy_mmvt_model.collective_variables[milestone.cv_index]
        result = cv.check_openmm_context_within_boundary(
            context, milestone.variables, verbose=True)
        final_result = final_result and result
        
    assert not final_result
    
    return