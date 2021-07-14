"""
test_common_prepare.py
"""

"""
test_prepare_1d_spherical.py
"""

import os
import sys
import glob

import pytest

import seekr2.modules.common_base as base
import seekr2.modules.mmvt_base as mmvt_base
import seekr2.prepare as prepare
from seekr2.modules import common_prepare
import seekr2.modules.common_cv as common_cv 
from seekr2.modules.common_cv import Spherical_cv_input

TEST_DIRECTORY = os.path.dirname(__file__)

def test_namd_model_input(host_guest_mmvt_model_input):
    """
    Test the ability to create a NAMD model input system.
    """
    model_input = host_guest_mmvt_model_input
    host_guest_mmvt_model_input.browndye_settings_input = None
    model_input.md_program = "namd"
    os.chdir(TEST_DIRECTORY)
    model, xml_path = prepare.generate_seekr2_model_and_filetree(
        model_input, True)
    return

def test_move_add_delete_input_anchors(tmp_path, tryp_ben_mmvt_model_input):
    """
    Test common_prepare.py's ability to move, add, and delete anchors
    in the model.
    """
    tmp_bd_settings = tryp_ben_mmvt_model_input.browndye_settings_input
    tryp_ben_mmvt_model_input.browndye_settings_input = None
    tryp_ben_mmvt_model_input.root_directory = os.path.join(tmp_path, "tryp_ben_mmvt_2")
    tryp_ben_mmvt_model_input.cv_inputs[0].input_anchors[0].radius = 0.06
    tryp_ben_mmvt_model_input.cv_inputs[0].input_anchors[0].upper_milestone_radius = 0.11
    tryp_ben_mmvt_model_input.cv_inputs[0].input_anchors[1].lower_milestone_radius = 0.11
    os.chdir(TEST_DIRECTORY)
    model2, model2_xml_path \
        = prepare.generate_seekr2_model_and_filetree(
            tryp_ben_mmvt_model_input, force_overwrite=True)
    model_dir = os.path.dirname(model2_xml_path)
    model2.anchor_rootdir = os.path.abspath(model_dir)
    assert model2.anchors[0].variables["r_0"] == 0.06
    assert model2.anchors[0].milestones[0].variables['radius'] == 0.11
    assert model2.anchors[1].milestones[0].variables['radius'] == 0.11
    assert len(model2.anchors) == 12
    assert model2.num_anchors == 12
    assert len(glob.glob(os.path.join(model2.anchor_rootdir, "anchor_*"))) == 11
    
    # Add anchor
    new_input_anchor = common_cv.Spherical_cv_anchor()
    new_input_anchor.radius = 0.1
    new_input_anchor.starting_amber_params = base.Amber_params()
    tryp_ben_mmvt_model_input.cv_inputs[0].input_anchors.insert(1, new_input_anchor)
    tryp_ben_mmvt_model_input.cv_inputs[0].input_anchors[0].upper_milestone_radius = None
    tryp_ben_mmvt_model_input.cv_inputs[0].input_anchors[2].lower_milestone_radius = None
    tryp_ben_mmvt_model_input.root_directory = os.path.join(tmp_path, "tryp_ben_mmvt_3")
    os.chdir(TEST_DIRECTORY)
    model3, model3_xml_path \
        = prepare.generate_seekr2_model_and_filetree(
            tryp_ben_mmvt_model_input, force_overwrite=False)
    model_dir = os.path.dirname(model3_xml_path)
    model3.anchor_rootdir = os.path.abspath(model_dir)
    assert model3.anchors[1].variables["r_0"] == 0.1
    assert model3.anchors[1].milestones[0].variables['radius'] == 0.08
    assert model3.anchors[1].milestones[1].variables['radius'] == 0.125
    assert len(model3.anchors) == 13
    assert model3.num_anchors == 13
    assert len(glob.glob(os.path.join(model3.anchor_rootdir, "anchor_*"))) == 12
    
    # Delete anchor
    tryp_ben_mmvt_model_input.cv_inputs[0].input_anchors.pop(1)
    tryp_ben_mmvt_model_input.browndye_settings_input = tmp_bd_settings
    tryp_ben_mmvt_model_input.root_directory = os.path.join(tmp_path, "tryp_ben_mmvt_4")
    os.chdir(TEST_DIRECTORY)
    model4, model4_xml_path \
        = prepare.generate_seekr2_model_and_filetree(
            tryp_ben_mmvt_model_input, force_overwrite=False)
    model_dir = os.path.dirname(model4_xml_path)
    model4.anchor_rootdir = os.path.abspath(model_dir)
    assert model4.anchors[0].milestones[0].variables['radius'] == 0.105
    assert model4.anchors[1].milestones[0].variables['radius'] == 0.105
    assert model4.anchors[1].milestones[1].variables['radius'] == 0.2
    assert len(model4.anchors) == 12
    assert model4.num_anchors == 12
    assert len(glob.glob(os.path.join(model4.anchor_rootdir, "anchor_*"))) == 11
    
    # adjust BD milestone location
    tryp_ben_mmvt_model_input.cv_inputs[0].input_anchors[-3].upper_milestone_radius = 1.65
    tryp_ben_mmvt_model_input.cv_inputs[0].input_anchors[-2].lower_milestone_radius = 1.65
    tryp_ben_mmvt_model_input.cv_inputs[0].input_anchors[-2].upper_milestone_radius = 1.85
    tryp_ben_mmvt_model_input.cv_inputs[0].input_anchors[-1].lower_milestone_radius = 1.85
    tryp_ben_mmvt_model_input.root_directory = os.path.join(tmp_path, "tryp_ben_mmvt_5")
    os.chdir(TEST_DIRECTORY)
    model5, model5_xml_path \
        = prepare.generate_seekr2_model_and_filetree(
            tryp_ben_mmvt_model_input, force_overwrite=False)
    model_dir = os.path.dirname(model5_xml_path)
    model5.anchor_rootdir = os.path.abspath(model_dir)
    assert model5.k_on_info.bd_milestones[0].outer_milestone.variables["radius"] == 1.85
    assert model5.k_on_info.bd_milestones[0].inner_milestone.variables["radius"] == 1.65
    
    # delete BD calculation
    tryp_ben_mmvt_model_input.browndye_settings_input = None
    tryp_ben_mmvt_model_input.root_directory = os.path.join(tmp_path, "tryp_ben_mmvt_6")
    os.chdir(TEST_DIRECTORY)
    model6, model6_xml_path \
        = prepare.generate_seekr2_model_and_filetree(
            tryp_ben_mmvt_model_input, force_overwrite=True)
    model_dir = os.path.dirname(model6_xml_path)
    model6.anchor_rootdir = os.path.abspath(model_dir)
    
    return