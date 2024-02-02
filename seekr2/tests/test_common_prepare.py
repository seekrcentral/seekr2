"""
test_common_prepare.py
"""

import os
import sys
import glob
import pytest

import seekr2.modules.common_base as base
import seekr2.prepare as prepare
from seekr2.modules import common_prepare
import seekr2.modules.common_cv as common_cv 
from seekr2.modules.common_cv import Spherical_cv_input
import seekr2.modules.check as check 

TEST_DIRECTORY = os.path.dirname(__file__)

def test_namd_model_input(tmp_path, host_guest_mmvt_model_input):
    """
    Test the ability to create a NAMD model input system.
    """
    model_input = host_guest_mmvt_model_input
    host_guest_mmvt_model_input.browndye_settings_input = None
    host_guest_mmvt_model_input.root_directory = os.path.join(
        tmp_path, "host_guest_mmvt_namd")
    model_input.md_program = "namd"
    os.chdir(TEST_DIRECTORY)
    model, xml_path = prepare.prepare(model_input, True)
    return

def test_move_add_delete_input_anchors(tmp_path, host_guest_mmvt_model_input):
    """
    Test common_prepare.py's ability to move, add, and delete anchors
    in the model.
    """
    tmp_bd_settings = host_guest_mmvt_model_input.browndye_settings_input
    host_guest_mmvt_model_input.browndye_settings_input = None
    host_guest_mmvt_model_input.root_directory = os.path.join(
        tmp_path, "host_guest_mmvt_2")
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors[0].radius = 0.06
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors[0]\
        .upper_milestone_radius = 0.11
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors[1]\
        .lower_milestone_radius = 0.11
    os.chdir(TEST_DIRECTORY)
    model2, model2_xml_path \
        = prepare.prepare(host_guest_mmvt_model_input, force_overwrite=True)
    model_dir = os.path.dirname(model2_xml_path)
    model2.anchor_rootdir = os.path.abspath(model_dir)
    assert model2.anchors[0].variables["r_0"] == 0.06
    assert model2.anchors[0].milestones[0].variables['radius'] == 0.11
    assert model2.anchors[1].milestones[0].variables['radius'] == 0.11
    assert len(model2.anchors) == 14
    assert model2.num_anchors == 14
    assert len(glob.glob(os.path.join(model2.anchor_rootdir, "anchor_*"))) == 13
    
    # Add anchor
    new_input_anchor = common_cv.Spherical_cv_anchor()
    new_input_anchor.radius = 0.1
    new_input_anchor.starting_amber_params = base.Amber_params()
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors.insert(
        1, new_input_anchor)
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors[0]\
        .upper_milestone_radius = None
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors[2]\
        .lower_milestone_radius = None
    host_guest_mmvt_model_input.root_directory = os.path.join(
        tmp_path, "host_guest_mmvt_3")
    os.chdir(TEST_DIRECTORY)
    model3, model3_xml_path \
        = prepare.prepare(host_guest_mmvt_model_input, force_overwrite=False)
    model_dir = os.path.dirname(model3_xml_path)
    model3.anchor_rootdir = os.path.abspath(model_dir)
    assert model3.anchors[1].variables["r_0"] == 0.1
    assert model3.anchors[1].milestones[0].variables['radius'] == 0.08
    assert model3.anchors[1].milestones[1].variables['radius'] == 0.125
    assert len(model3.anchors) == 15
    assert model3.num_anchors == 15
    assert len(glob.glob(os.path.join(model3.anchor_rootdir, "anchor_*"))) == 14
    
    # Delete anchor
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors.pop(1)
    host_guest_mmvt_model_input.browndye_settings_input = tmp_bd_settings
    host_guest_mmvt_model_input.root_directory = os.path.join(
        tmp_path, "host_guest_mmvt_4")
    os.chdir(TEST_DIRECTORY)
    model4, model4_xml_path \
        = prepare.prepare(host_guest_mmvt_model_input, force_overwrite=False)
    model_dir = os.path.dirname(model4_xml_path)
    model4.anchor_rootdir = os.path.abspath(model_dir)
    assert model4.anchors[0].milestones[0].variables['radius'] == 0.105
    assert model4.anchors[1].milestones[0].variables['radius'] == 0.105
    assert model4.anchors[1].milestones[1].variables['radius'] == 0.2
    assert len(model4.anchors) == 14
    assert model4.num_anchors == 14
    assert len(glob.glob(os.path.join(model4.anchor_rootdir, "anchor_*"))) == 13
    
    # adjust BD milestone location
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors[-3]\
        .upper_milestone_radius = 1.65
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors[-2]\
        .lower_milestone_radius = 1.65
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors[-2]\
        .radius = 1.75
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors[-2]\
        .upper_milestone_radius = 1.85
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors[-1]\
        .lower_milestone_radius = 1.85
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors[-1]\
        .radius = 1.95
    host_guest_mmvt_model_input.root_directory = os.path.join(
        tmp_path, "host_guest_mmvt_5")
    os.chdir(TEST_DIRECTORY)
    model5, model5_xml_path \
        = prepare.prepare(host_guest_mmvt_model_input, force_overwrite=False)
    model_dir = os.path.dirname(model5_xml_path)
    model5.anchor_rootdir = os.path.abspath(model_dir)
    assert model5.k_on_info.bd_milestones[0].outer_milestone\
        .variables["radius"] == 1.85
    assert model5.k_on_info.bd_milestones[0].inner_milestone\
        .variables["radius"] == 1.65
    
    # delete BD calculation
    host_guest_mmvt_model_input.browndye_settings_input = None
    host_guest_mmvt_model_input.root_directory = os.path.join(
        tmp_path, "host_guest_mmvt_6")
    os.chdir(TEST_DIRECTORY)
    model6, model6_xml_path \
        = prepare.prepare(host_guest_mmvt_model_input, force_overwrite=True)
    model_dir = os.path.dirname(model6_xml_path)
    model6.anchor_rootdir = os.path.abspath(model_dir)
    
    return

def test_modify_model(tmp_path, host_guest_mmvt_model_input):
    """
    Test the ability to modify models from the input and save the
    resulting anchors and files.
    """
    host_guest_mmvt_model_input.root_directory = os.path.join(
        tmp_path, "host_guest_mmvt_2")
    os.chdir(TEST_DIRECTORY)
    model2, model2_xml_path \
        = prepare.prepare(host_guest_mmvt_model_input, force_overwrite=False)
    model_dir = os.path.dirname(model2_xml_path)
    model2.anchor_rootdir = os.path.abspath(model_dir)
    check.check_systems_within_Voronoi_cells(model2)
    
    # First, move an input milestone
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors[0].radius = 0.06
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors[0]\
        .upper_milestone_radius = 0.11
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors[1]\
        .lower_milestone_radius = 0.11
    os.chdir(TEST_DIRECTORY)
    model2, model2_xml_path \
        = prepare.prepare(host_guest_mmvt_model_input, force_overwrite=False)
    model_dir = os.path.dirname(model2_xml_path)
    model2.anchor_rootdir = os.path.abspath(model_dir)
    check.check_systems_within_Voronoi_cells(model2)
        
    # Next, add an input milestone
    new_input_anchor = common_cv.Spherical_cv_anchor()
    new_input_anchor.radius = 0.1
    new_input_anchor.starting_amber_params = base.Amber_params()
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors.insert(
        1, new_input_anchor)
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors[0]\
        .upper_milestone_radius = None
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors[2]\
        .lower_milestone_radius = None
    os.chdir(TEST_DIRECTORY)
    model2, model2_xml_path \
        = prepare.prepare(host_guest_mmvt_model_input, force_overwrite=False)
    model_dir = os.path.dirname(model2_xml_path)
    model2.anchor_rootdir = os.path.abspath(model_dir)
    check.check_systems_within_Voronoi_cells(model2)
        
    # Add an output file to a relevant anchor
    new_anchor_directory = os.path.join(model_dir, "anchor_1")
    output_file_name = os.path.join(new_anchor_directory, "prod/mmvt1.out")
    with open(output_file_name, "w") as f:
        f.write("MMVT output")
    
    # Test that attempting to delete that anchor will fail
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors.pop(1)
    os.chdir(TEST_DIRECTORY)
    with pytest.raises(Exception) as my_exception:
        model2, model2_xml_path \
            = prepare.prepare(host_guest_mmvt_model_input, 
                              force_overwrite=False)
    
    # Test that attempting to move that anchor will fail
    new_input_anchor = common_cv.Spherical_cv_anchor()
    new_input_anchor.radius = 0.105
    new_input_anchor.starting_amber_params = base.Amber_params()
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors.insert(
        1, new_input_anchor)
    os.chdir(TEST_DIRECTORY)
    with pytest.raises(Exception) as my_exception:
        model2, model2_xml_path \
            = prepare.prepare(host_guest_mmvt_model_input, 
                              force_overwrite=False)
    
    # Force delete the anchor
    os.chdir(TEST_DIRECTORY)
    model2, model2_xml_path \
        = prepare.prepare(host_guest_mmvt_model_input, force_overwrite=True)
    model_dir = os.path.dirname(model2_xml_path)
    model2.anchor_rootdir = os.path.abspath(model_dir)
    check.check_systems_within_Voronoi_cells(model2)
    
    # Add an output file to the BD directory
    b_surface_directory = os.path.join(model_dir, "b_surface")
    output_file_name = os.path.join(b_surface_directory, "results1.xml")
    with open(output_file_name, "w") as f:
        f.write("BD output")
    
    # Modify the BD anchor
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors[-1].radius = 1.37
    os.chdir(TEST_DIRECTORY)
    with pytest.raises(Exception) as my_exception:
        model2, model2_xml_path \
            = prepare.prepare(host_guest_mmvt_model_input, 
                              force_overwrite=False)
    model_dir = os.path.dirname(model2_xml_path)
    model2.anchor_rootdir = os.path.abspath(model_dir)
    check.check_systems_within_Voronoi_cells(model2)
    
    # Force modification of BD milestone
    host_guest_mmvt_model_input.cv_inputs[0].input_anchors[-1].radius = 1.38
    os.chdir(TEST_DIRECTORY)
    model2, model2_xml_path \
        = prepare.prepare(host_guest_mmvt_model_input, force_overwrite=True)
    model_dir = os.path.dirname(model2_xml_path)
    model2.anchor_rootdir = os.path.abspath(model_dir)
    check.check_systems_within_Voronoi_cells(model2)