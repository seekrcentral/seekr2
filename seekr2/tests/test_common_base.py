"""
test_common_base.py
"""

import pytest
import random

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
