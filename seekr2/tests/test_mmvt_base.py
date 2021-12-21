"""
test_mmvt_base.py
"""

import os
import shutil

from seekr2.modules import mmvt_sim_openmm

def test_MMVT_anchor_id_from_alias(tryp_ben_mmvt_model):
    """
    Test the anchor's ability to return information about itself and its
    milestones - id_from_alias.
    """
    anchor0 = tryp_ben_mmvt_model.anchors[0]
    anchor1 = tryp_ben_mmvt_model.anchors[1]
    # test id_from_alias
    assert anchor0.id_from_alias(1) == 0
    assert anchor1.id_from_alias(1) == 0
    assert anchor1.id_from_alias(2) == 1
    return

def test_MMVT_anchor_alias_from_id(tryp_ben_mmvt_model):
    """
    Test the anchor's ability to return information about itself and its
    milestones - alias_from_id.
    """
    anchor0 = tryp_ben_mmvt_model.anchors[0]
    anchor1 = tryp_ben_mmvt_model.anchors[1]
    # test id_from_alias
    assert anchor0.alias_from_id(0) == 1
    assert anchor1.alias_from_id(0) == 1
    assert anchor1.alias_from_id(1) == 2
    return

def test_MMVT_anchor_alias_from_neighbor_id(tryp_ben_mmvt_model):
    """
    Test the anchor's ability to return information about itself and its
    milestones - alias_from_neighbor_id.
    """
    anchor0 = tryp_ben_mmvt_model.anchors[0]
    anchor1 = tryp_ben_mmvt_model.anchors[1]
    # test id_from_alias
    assert anchor0.alias_from_neighbor_id(1) == 1
    assert anchor1.alias_from_neighbor_id(0) == 1
    assert anchor1.alias_from_neighbor_id(2) == 2
    return

def test_MMVT_anchor_get_ids(tryp_ben_mmvt_model):
    """
    Test the anchor's ability to return information about itself and its
    milestones - id_from_alias.
    """
    anchor0 = tryp_ben_mmvt_model.anchors[0]
    anchor1 = tryp_ben_mmvt_model.anchors[1]
    # test id_from_alias
    assert anchor0.get_ids() == [0]
    assert anchor1.get_ids() == [0, 1]
    return

def test_check_openmm_context_within_boundary(tryp_ben_mmvt_model, tmp_path):
    """
    Test whether the check can find systems that exist outside the proper
    Voronoi cells.
    """
    # Test check success: system starting within Voronoi cells
    anchor = tryp_ben_mmvt_model.anchors[0]
    output_file = os.path.join(tmp_path, "output.txt")
    tryp_ben_mmvt_model.openmm_settings.cuda_platform_settings = None
    tryp_ben_mmvt_model.openmm_settings.reference_platform = True
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(
        tryp_ben_mmvt_model, anchor, output_file)
    context = my_sim_openmm.simulation.context
    
    for milestone in anchor.milestones:
        cv = tryp_ben_mmvt_model.collective_variables[milestone.cv_index]
        result = cv.check_openmm_context_within_boundary(
            context, milestone.variables, verbose=True)
        assert result
    
    # Test check failure(s): system outside Voronoi cells
    anchor_pdb1 \
        = tryp_ben_mmvt_model.anchors[0].amber_params.pdb_coordinates_filename
    anchor_pdb2 \
        = tryp_ben_mmvt_model.anchors[1].amber_params.pdb_coordinates_filename
    anchor_pdb_src_path1 = os.path.join(
        tryp_ben_mmvt_model.anchor_rootdir, 
        tryp_ben_mmvt_model.anchors[0].directory,
        tryp_ben_mmvt_model.anchors[0].building_directory, 
        anchor_pdb1)
    anchor_pdb_src_path2 = os.path.join(
        tryp_ben_mmvt_model.anchor_rootdir, 
        tryp_ben_mmvt_model.anchors[1].directory,
        tryp_ben_mmvt_model.anchors[1].building_directory, 
        anchor_pdb2)
    anchor_pdb_dest_path1 = os.path.join(
        tryp_ben_mmvt_model.anchor_rootdir, 
        tryp_ben_mmvt_model.anchors[0].directory,
        tryp_ben_mmvt_model.anchors[0].building_directory, 
        anchor_pdb2)
    anchor_pdb_dest_path2 = os.path.join(
        tryp_ben_mmvt_model.anchor_rootdir, 
        tryp_ben_mmvt_model.anchors[1].directory,
        tryp_ben_mmvt_model.anchors[1].building_directory, 
        anchor_pdb1)
    shutil.copyfile(anchor_pdb_src_path1, anchor_pdb_dest_path2)
    shutil.copyfile(anchor_pdb_src_path2, anchor_pdb_dest_path1)
    tryp_ben_mmvt_model.anchors[0].amber_params.pdb_coordinates_filename \
        = anchor_pdb2
    tryp_ben_mmvt_model.anchors[1].amber_params.pdb_coordinates_filename \
        = anchor_pdb1
        
    my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(
        tryp_ben_mmvt_model, anchor, output_file)
    context = my_sim_openmm.simulation.context
    
    for milestone in anchor.milestones:
        cv = tryp_ben_mmvt_model.collective_variables[milestone.cv_index]
        result = cv.check_openmm_context_within_boundary(
            context, milestone.variables, verbose=True)
        assert not result
    
    
    return