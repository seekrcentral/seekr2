"""
test_common_base.py
"""

import pytest
import random

import numpy as np
from parmed import unit

import seekr2.modules.common_base as base

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
    