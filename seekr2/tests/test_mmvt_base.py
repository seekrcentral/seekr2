"""
test_mmvt_base.py
"""

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