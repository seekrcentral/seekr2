"""
test_common_base.py
"""
import os
import pytest
import random

import numpy as np
from parmed import unit

import seekr2.modules.common_base as base

TEST_DIRECTORY = os.path.dirname(__file__)

def test_strBool():
    assert base.strBool('True') == True
    assert base.strBool('true') == True
    assert base.strBool('TRUE') == True
    assert base.strBool('False') == False
    assert base.strBool('false') == False
    assert base.strBool('FALSE') == False
    with pytest.raises(Exception):
        base.strBool('balderdash')
    return

def test_order_files_numerically():
    string_list = ["/path/to/anchor0/output0_0", "/path/to/anchor0/output0_1",
                   "/path/to/anchor0/output0_2", "/path/to/anchor0/output1_0",
                   "/path/to/anchor0/output1_1", "/path/to/anchor0/output1_2",
                   "/path/to/anchor1/output0_0", "/path/to/anchor1/output0_1",
                   "/path/to/anchor1/output2_0", "/path/to/anchor1/output10_0"]
    desired_list = string_list[:]
    random.shuffle(string_list)
    ordered_list = base.order_files_numerically(string_list)
    
    for item1, item2 in zip(ordered_list, desired_list):
        assert item1==item2
    
    scrambled_list = ["me2.345.pdb", "me1.234.pdb", "me-3.456.pdb"]
    correct_list = ["me-3.456.pdb", "me1.234.pdb", "me2.345.pdb"]
    ordered_list = base.order_files_numerically(scrambled_list, func=float)
    for item1, item2 in zip(ordered_list, correct_list):
        assert item1==item2
    
    string_list = ["mmvt.restart1.out", "mmvt.restart2.out"]
    desired_list = string_list[:]
    random.shuffle(string_list)
    ordered_list = base.order_files_numerically(string_list)
    
    for item1, item2 in zip(ordered_list, desired_list):
        assert item1==item2
    
    string_list = ["/path/to/anchor_0.1/output0_0", 
                   "/path/to/anchor_0.2/output0_1"]
    
    desired_list = string_list[:]
    random.shuffle(string_list)
    ordered_list = base.order_files_numerically(string_list)
    
    return

def test_box_vectors():
    box_vector_q = unit.Quantity(
        [[64.9127105, 0.0, 0.0], 
         [-21.6375684, 61.2002909, 0.0], 
         [-21.6375684, -30.6001417, 53.0010088]], 
                      unit=unit.angstrom)
    box_vector = base.Box_vectors()
    box_vector.from_quantity(box_vector_q)
    assert np.isclose(box_vector.ax, 6.49127105)
    assert np.isclose(box_vector.ay, 0.0)
    assert np.isclose(box_vector.az, 0.0)
    assert np.isclose(box_vector.bx, -2.16375684)
    assert np.isclose(box_vector.by, 6.12002909)
    assert np.isclose(box_vector.bz, 0.0)
    assert np.isclose(box_vector.cx, -2.16375684)
    assert np.isclose(box_vector.cy, -3.06001417)
    assert np.isclose(box_vector.cz, 5.30010088)
    box_vector_q2 = box_vector.to_quantity()
    assert np.isclose(box_vector_q2.value_in_unit(unit.nanometers)[0][0], 
                      6.49127105)
    assert np.isclose(box_vector_q2.value_in_unit(unit.nanometers)[0][1], 0.0)
    assert np.isclose(box_vector_q2.value_in_unit(unit.nanometers)[0][2], 0.0)
    assert np.isclose(box_vector_q2.value_in_unit(unit.nanometers)[1][0], 
                      -2.16375684)
    assert np.isclose(box_vector_q2.value_in_unit(unit.nanometers)[1][1], 
                      6.12002909)
    assert np.isclose(box_vector_q2.value_in_unit(unit.nanometers)[1][2], 0.0)
    assert np.isclose(box_vector_q2.value_in_unit(unit.nanometers)[2][0], 
                      -2.16375684)
    assert np.isclose(box_vector_q2.value_in_unit(unit.nanometers)[2][1], 
                      -3.06001417)
    assert np.isclose(box_vector_q2.value_in_unit(unit.nanometers)[2][2], 
                      5.30010088)
    box_6_vector = [61.239410, 61.239410, 61.239410, 109.471, 109.471, 109.471]
    box_vector2 = base.Box_vectors()
    box_vector2.from_6_vector(box_6_vector)
    assert np.isclose(box_vector2.ax, 6.1239410, atol=0.001)
    assert np.isclose(box_vector2.ay, 0.0)
    assert np.isclose(box_vector2.az, 0.0)
    assert np.isclose(box_vector2.bx, -2.041313667, atol=0.001)
    assert np.isclose(box_vector2.by, 5.773707003, atol=0.001)
    assert np.isclose(box_vector2.bz, 0.0)
    assert np.isclose(box_vector2.cx, -2.041313667, atol=0.001)
    assert np.isclose(box_vector2.cy, -2.886853152, atol=0.001)
    assert np.isclose(box_vector2.cz, 5.00017714, atol=0.001)
    expected_box_6_vector = [6.1239410, 6.1239410, 6.1239410, 109.47122, 
                             109.47122, 109.47122]
    resulting_box_6_vector = box_vector2.to_6_vector()
    assert np.isclose(np.array(expected_box_6_vector), 
                      np.array(resulting_box_6_vector)).all()
    return
    
def test_Barostat_settings():
    """
    Initialize the Barostat_settings_openmm and Barostat_settings_namd
    object to get more coverage.
    """
    dummy1 = base.Barostat_settings_openmm()
    dummy2 = base.Barostat_settings_namd()
    return
    
def test_Cuda_platform_settings():
    """
    Initialize the Cuda_platforms_settings() object.
    """
    cuda_platform_settings = base.Cuda_platform_settings()
    properties = cuda_platform_settings.make_properties_dict()
    assert properties["CudaDeviceIndex"] == "0"
    assert properties["CudaPrecision"] == "mixed"
    return

def test_model_get_type(host_guest_mmvt_model, host_guest_elber_model):
    """
    Assure that the model.get_type() function is working properly
    """
    assert host_guest_mmvt_model.get_type() == "mmvt"
    assert host_guest_elber_model.get_type() == "elber"
    host_guest_mmvt_model.calculation_type = "wrong"
    with pytest.raises(Exception) as e_info:
        host_guest_mmvt_model.get_type()
    return

def test_model_get_timestep(host_guest_mmvt_model):
    """
    Assure that the model.get_type() function is working properly
    """
    assert host_guest_mmvt_model.get_timestep() == 0.002
    return

def test_model_get_bulk_index(host_guest_mmvt_model, toy_multi_model):
    assert host_guest_mmvt_model.get_bulk_index() == 13
    assert toy_multi_model.get_bulk_index() == None
    return

def test_get_box_vectors_from_pdb():
    expected_box_vectors = unit.Quantity(
        [[40.142, 0.0, 0.0], 
         [0.0, 40.329, 0.0], 
         [0.0, 0.0, 32.472]], 
                      unit=unit.angstrom)
    test_pdb_filename = os.path.join(
        TEST_DIRECTORY, 
        "../data/hostguest_files/hostguest_at0.5.pdb")
    result = base.get_box_vectors_from_pdb(test_pdb_filename)
    assert np.isclose(expected_box_vectors.value_in_unit(unit.angstroms), 
                      result.value_in_unit(unit.angstroms)).all()
                      
def test_parse_xml_list():
    input_list1 = [3,4,5,6]
    assert base.parse_xml_list(input_list1) == input_list1
    input_range1 = "range(5)"
    assert base.parse_xml_list(input_range1) == [0,1,2,3,4]
    input_range1 = "range(2,7)"
    assert base.parse_xml_list(input_range1) == [2,3,4,5,6]
    input_range1 = "range(3,9,2)"
    assert base.parse_xml_list(input_range1) == [3,5,7]
    with pytest.raises(Exception):
        base.parse_xml_list(2)
    with pytest.raises(Exception):
        base.parse_xml_list("balderdash")
    return

def test_save_new_model(host_guest_mmvt_model):
    """
    
    """
    MODEL_GLOB = "model_pre_test_*.xml"
    MODEL_BASE = "model_pre_test_{}.xml"
    root = host_guest_mmvt_model.anchor_rootdir
    base.save_new_model(host_guest_mmvt_model, MODEL_GLOB, MODEL_BASE)
    old_model_path = os.path.join(root, MODEL_BASE.format(0))
    print("old_model_path:", old_model_path)
    assert(os.path.exists(old_model_path))
    new_model_path = os.path.join(root, "model.xml")
    assert(os.path.exists(new_model_path))
    
    old_model_path2 = os.path.join(root, MODEL_BASE.format(1))
    base.save_new_model(host_guest_mmvt_model, MODEL_GLOB, MODEL_BASE)
    assert(os.path.exists(old_model_path2))
    assert(os.path.exists(new_model_path))
    return