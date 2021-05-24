"""
test_common_prepare.py
"""

"""
test_prepare_1d_spherical.py
"""

import os
import sys

import seekr2.modules.common_base as base
import seekr2.modules.mmvt_base as mmvt_base
import seekr2.prepare as prepare
from seekr2.modules import common_prepare
import seekr2.modules.common_cv as common_cv 
from seekr2.modules.common_cv import Spherical_cv_input

def create_model_input(root_dir):
    """
    TEMPORARY
    """
    this_dir = os.path.dirname(os.path.abspath(__file__))
    model_input = common_prepare.Model_input()
    model_input.calculation_settings = mmvt_base.MMVT_settings()
    model_input.root_directory = root_dir
    model_input.temperature = 282.55
    model_input.calculation_type = "mmvt"
    model_input.calculation_settings.md_output_interval = 12345
    model_input.calculation_settings.md_steps_per_anchor = 76543
    cv_input1 = common_cv.Spherical_cv_input()
    cv_input1.group1 = list(range(147))
    cv_input1.group2 = list(range(147, 162))
    cv_input1.input_anchors = []
    input_anchor1 = common_cv.Spherical_cv_anchor()
    input_anchor1.radius = 0.05
    input_anchor1.starting_amber_params = base.Amber_params()
    input_anchor1.starting_amber_params.prmtop_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest.parm7")
    #input_anchor1.starting_amber_params.inpcrd_filename = \
    #    os.path.join(this_dir, "../data/hostguest_files/hostguest.rst7")
    input_anchor1.starting_amber_params.pdb_coordinates_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest_at0.5.pdb")
    input_anchor1.bound_state = True
    cv_input1.input_anchors.append(input_anchor1)
    
    input_anchor2 = common_cv.Spherical_cv_anchor()
    input_anchor2.radius = 0.15
    input_anchor2.starting_amber_params = base.Amber_params()
    input_anchor2.starting_amber_params.prmtop_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest.parm7")
    #input_anchor2.starting_amber_params.inpcrd_filename = \
    #    os.path.join(this_dir, "../data/hostguest_files/hostguest.rst7")
    input_anchor2.starting_amber_params.pdb_coordinates_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest_at1.5.pdb")
    input_anchor2.bound_state = False
    cv_input1.input_anchors.append(input_anchor2)
    
    input_anchor3 = common_cv.Spherical_cv_anchor()
    input_anchor3.radius = 0.25
    input_anchor3.starting_amber_params = base.Amber_params()
    input_anchor3.starting_amber_params.prmtop_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest.parm7")
    #input_anchor3.starting_amber_params.inpcrd_filename = \
    #    os.path.join(this_dir, "../data/hostguest_files/hostguest.rst7")
    input_anchor3.starting_amber_params.pdb_coordinates_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest_at2.5.pdb")
    input_anchor3.bound_state = False
    cv_input1.input_anchors.append(input_anchor3)
    
    input_anchor4 = common_cv.Spherical_cv_anchor()
    input_anchor4.radius = 0.35
    input_anchor4.starting_amber_params = base.Amber_params()
    input_anchor4.starting_amber_params.prmtop_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest.parm7")
    #input_anchor4.starting_amber_params.inpcrd_filename = \
    #    os.path.join(this_dir, "../data/hostguest_files/hostguest.rst7")
    input_anchor4.starting_amber_params.pdb_coordinates_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest_at3.5.pdb")
    input_anchor4.bound_state = False
    cv_input1.input_anchors.append(input_anchor4)
    
    input_anchor5 = common_cv.Spherical_cv_anchor()
    input_anchor5.radius = 0.45
    input_anchor5.starting_amber_params = base.Amber_params()
    input_anchor5.starting_amber_params.prmtop_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest.parm7")
    #input_anchor5.starting_amber_params.inpcrd_filename = \
    #    os.path.join(this_dir, "../data/hostguest_files/hostguest.rst7")
    input_anchor5.starting_amber_params.pdb_coordinates_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest_at4.5.pdb")
    input_anchor5.bound_state = False
    cv_input1.input_anchors.append(input_anchor5)
    
    input_anchor6 = common_cv.Spherical_cv_anchor()
    input_anchor6.radius = 0.55
    input_anchor6.starting_amber_params = base.Amber_params()
    input_anchor6.starting_amber_params.prmtop_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest.parm7")
    #input_anchor6.starting_amber_params.inpcrd_filename = \
    #    os.path.join(this_dir, "../data/hostguest_files/hostguest.rst7")
    input_anchor6.starting_amber_params.pdb_coordinates_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest_at5.5.pdb")
    input_anchor6.bound_state = False
    cv_input1.input_anchors.append(input_anchor6)
    
    input_anchor7 = common_cv.Spherical_cv_anchor()
    input_anchor7.radius = 0.65
    input_anchor7.starting_amber_params = base.Amber_params()
    input_anchor7.starting_amber_params.prmtop_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest.parm7")
    #input_anchor7.starting_amber_params.inpcrd_filename = \
    #    os.path.join(this_dir, "../data/hostguest_files/hostguest.rst7")
    input_anchor7.starting_amber_params.pdb_coordinates_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest_at6.5.pdb")
    input_anchor7.bound_state = False
    cv_input1.input_anchors.append(input_anchor7)
    
    input_anchor8 = common_cv.Spherical_cv_anchor()
    input_anchor8.radius = 0.75
    input_anchor8.starting_amber_params = base.Amber_params()
    input_anchor8.starting_amber_params.prmtop_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest.parm7")
    #input_anchor8.starting_amber_params.inpcrd_filename = \
    #    os.path.join(this_dir, "../data/hostguest_files/hostguest.rst7")
    input_anchor8.starting_amber_params.pdb_coordinates_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest_at7.5.pdb")
    input_anchor8.bound_state = False
    cv_input1.input_anchors.append(input_anchor8)
    
    input_anchor9 = common_cv.Spherical_cv_anchor()
    input_anchor9.radius = 0.85
    input_anchor9.starting_amber_params = base.Amber_params()
    input_anchor9.starting_amber_params.prmtop_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest.parm7")
    #input_anchor9.starting_amber_params.inpcrd_filename = \
    #    os.path.join(this_dir, "../data/hostguest_files/hostguest.rst7")
    input_anchor9.starting_amber_params.pdb_coordinates_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest_at8.5.pdb")
    input_anchor9.bound_state = False
    cv_input1.input_anchors.append(input_anchor9)
    
    input_anchor10 = common_cv.Spherical_cv_anchor()
    input_anchor10.radius = 0.95
    input_anchor10.starting_amber_params = base.Amber_params()
    input_anchor10.starting_amber_params.prmtop_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest.parm7")
    #input_anchor10.starting_amber_params.inpcrd_filename = \
    #    os.path.join(this_dir, "../data/hostguest_files/hostguest.rst7")
    input_anchor10.starting_amber_params.pdb_coordinates_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest_at9.5.pdb")
    input_anchor10.bound_state = False
    cv_input1.input_anchors.append(input_anchor10)
    
    input_anchor11 = common_cv.Spherical_cv_anchor()
    input_anchor11.radius = 1.05
    input_anchor11.starting_amber_params = base.Amber_params()
    input_anchor11.starting_amber_params.prmtop_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest.parm7")
    #input_anchor11.starting_amber_params.inpcrd_filename = \
    #    os.path.join(this_dir, "../data/hostguest_files/hostguest.rst7")
    input_anchor11.starting_amber_params.pdb_coordinates_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest_at10.5.pdb")
    input_anchor11.bound_state = False
    cv_input1.input_anchors.append(input_anchor11)
    
    input_anchor12 = common_cv.Spherical_cv_anchor()
    input_anchor12.radius = 1.15
    input_anchor12.starting_amber_params = base.Amber_params()
    input_anchor12.starting_amber_params.prmtop_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest.parm7")
    #input_anchor12.starting_amber_params.inpcrd_filename = \
    #    os.path.join(this_dir, "../data/hostguest_files/hostguest.rst7")
    input_anchor12.starting_amber_params.pdb_coordinates_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest_at11.5.pdb")
    input_anchor12.bound_state = False
    cv_input1.input_anchors.append(input_anchor12)
    
    input_anchor13 = common_cv.Spherical_cv_anchor()
    input_anchor13.radius = 1.25
    input_anchor13.starting_amber_params = base.Amber_params()
    input_anchor13.starting_amber_params.prmtop_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest.parm7")
    #input_anchor13.starting_amber_params.inpcrd_filename = \
    #    os.path.join(this_dir, "../data/hostguest_files/hostguest.rst7")
    input_anchor13.starting_amber_params.pdb_coordinates_filename = \
        os.path.join(this_dir, "../data/hostguest_files/hostguest_at12.5.pdb")
    input_anchor13.bound_state = False
    cv_input1.input_anchors.append(input_anchor13)
    
    input_anchor14 = common_cv.Spherical_cv_anchor()
    input_anchor14.radius = 1.35
    input_anchor14.starting_amber_params = base.Amber_params()
    input_anchor14.starting_amber_params.prmtop_filename = ""
    #input_anchor14.starting_amber_params.inpcrd_filename = ""
    input_anchor14.starting_amber_params.pdb_coordinates_filename = ""
    input_anchor14.bound_state = False
    input_anchor14.bulk_anchor = True
    cv_input1.input_anchors.append(input_anchor14)
    
    model_input.cv_inputs = [cv_input1]
    return model_input

def test_Input_model(tmp_path):
    model_input = create_model_input(tmp_path)
    assert model_input.temperature == 282.55
    assert model_input.calculation_settings.md_output_interval == 12345
    assert model_input.calculation_settings.md_steps_per_anchor == 76543
    
    assert len(model_input.cv_inputs) == 1
    assert len(model_input.cv_inputs[0].input_anchors) == 14
    print("model_input.root_directory:", model_input.root_directory)
    
    model, xml_path = prepare.generate_seekr2_model_and_filetree(
        model_input, True)
    model.anchor_rootdir = tmp_path
    assert os.path.exists(os.path.join(model.anchor_rootdir, "anchor_1"))
    assert os.path.exists(os.path.join(model.anchor_rootdir, "anchor_1", 
                                       "prod"))
    assert os.path.exists(os.path.join(model.anchor_rootdir, "anchor_1", 
                                       "building"))
    assert os.path.exists(os.path.join(model.anchor_rootdir, "anchor_1", 
                                       "building", "hostguest.parm7"))
    #assert os.path.exists(os.path.join(model.anchor_rootdir, "anchor_1", 
    #                                   "building", "hostguest.rst7"))
    assert os.path.exists(os.path.join(model.anchor_rootdir, "anchor_1", 
                                       "building", "hostguest_at1.5.pdb"))
    assert model.temperature == 282.55
    assert model.num_anchors == 14
    assert model.num_milestones == 13
    assert len(model.anchors) == 14
    assert len(model.collective_variables) == 1
    assert model.openmm_settings.initial_temperature == 282.55
    assert model.calculation_settings.energy_reporter_interval == 12345
    assert model.calculation_settings.restart_checkpoint_interval == 12345
    assert model.calculation_settings.trajectory_reporter_interval == 12345
    assert model.calculation_settings.num_production_steps == 76543

""" # this test needs the serializer to be able to find the Spherical_cv_input 
    #class
def test_read_simplified_input(tmp_path):
    this_dir = os.path.dirname(os.path.abspath(__file__))
    model_input_filename = os.path.join(this_dir, "../data/sample_input.xml")
    model_input = prepare_1d_spherical.Model_input()
    model_input.deserialize(model_input_filename, user_input = True)
    model_input.root_directory = os.path.join(tmp_path, "test_mmvt")
    prepare_1d_spherical.generate_openmmvt_model_and_filetree(model_input)
"""