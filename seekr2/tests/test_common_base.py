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
        
    return

def test_box_vectors():
    box_vector_q = unit.Quantity(
        [[64.0, 0.0, 0.0], 
         [-21.0, 61.0, 0.0], 
         [-21.0, -30.0, 53.0]], 
                      unit=unit.angstrom)
    box_vector = base.Box_vectors()
    box_vector.from_quantity(box_vector_q)
    assert np.isclose(box_vector.ax, 6.4)
    assert np.isclose(box_vector.ay, 0.0)
    assert np.isclose(box_vector.az, 0.0)
    assert np.isclose(box_vector.bx, -2.1)
    assert np.isclose(box_vector.by, 6.1)
    assert np.isclose(box_vector.bz, 0.0)
    assert np.isclose(box_vector.cx, -2.1)
    assert np.isclose(box_vector.cy, -3.0)
    assert np.isclose(box_vector.cz, 5.3)
    box_vector_q2 = box_vector.to_quantity()
    assert np.isclose(box_vector_q2.value_in_unit(unit.nanometers)[0][0], 6.4)
    assert np.isclose(box_vector_q2.value_in_unit(unit.nanometers)[0][1], 0.0)
    assert np.isclose(box_vector_q2.value_in_unit(unit.nanometers)[0][2], 0.0)
    assert np.isclose(box_vector_q2.value_in_unit(unit.nanometers)[1][0], -2.1)
    assert np.isclose(box_vector_q2.value_in_unit(unit.nanometers)[1][1], 6.1)
    assert np.isclose(box_vector_q2.value_in_unit(unit.nanometers)[1][2], 0.0)
    assert np.isclose(box_vector_q2.value_in_unit(unit.nanometers)[2][0], -2.1)
    assert np.isclose(box_vector_q2.value_in_unit(unit.nanometers)[2][1], -3.0)
    assert np.isclose(box_vector_q2.value_in_unit(unit.nanometers)[2][2], 5.3)
    
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

def test_model_get_type(tryp_ben_mmvt_model, tryp_ben_elber_model):
    """
    Assure that the model.get_type() function is working properly
    """
    assert tryp_ben_mmvt_model.get_type() == "mmvt"
    assert tryp_ben_elber_model.get_type() == "elber"
    tryp_ben_mmvt_model.calculation_type = "wrong"
    with pytest.raises(Exception) as e_info:
        tryp_ben_mmvt_model.get_type()
    return
    
def test_get_box_vectors_from_pdb():
    expected_box_vectors = unit.Quantity(
        [[61.359, 0.0, 0.0], 
         [-20.451767557539473, 57.850255701875476, 0.0], 
         [-20.451767557539473, -28.92251350474921, 50.101300355778946]], 
                      unit=unit.angstrom)
    test_pdb_filename = os.path.join(
        TEST_DIRECTORY, 
        "../data/trypsin_benzamidine_files/mmvt/tryp_ben_at0.pdb")
    result = base.get_box_vectors_from_pdb(test_pdb_filename)
    assert np.isclose(expected_box_vectors.value_in_unit(unit.angstroms), 
                      result.value_in_unit(unit.angstroms)).all()
    