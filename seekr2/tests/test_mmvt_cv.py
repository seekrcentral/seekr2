"""
Test the mmvt_cv.py module
"""

import numpy as np

import seekr2.modules.common_cv as common_cv
import seekr2.modules.check as check

def test_tiwary_input_and_model(tiwary_mmvt_model_input, tiwary_mmvt_model):
    """
    
    """
    assert len(tiwary_mmvt_model_input.cv_inputs[0].order_parameters) == 3
    check.check_pre_simulation_all(tiwary_mmvt_model)
    return

def test_distance_order_parameter_get_value():
    """
    
    """
    op = common_cv.Tiwary_cv_distance_order_parameter()
    com1 = np.array([1.0, 1.0, 1.0])
    com2 = np.array([4.0, 5.0, 1.0])
    expected_dist = 5.0
    dist = op.get_value([com1, com2])
    assert np.isclose(dist, expected_dist)
    return

def test_angle_order_parameter_get_value():
    """
    
    """
    op = common_cv.Tiwary_cv_angle_order_parameter()
    com1 = np.array([1.0, 0.0, 0.0])
    com2 = np.array([0.0, 0.0, 0.0])
    com3 = np.array([0.0, 0.0, 1.0])
    expected_angle = 0.5 * np.pi
    angle = op.get_value([com1, com2, com3])
    assert np.isclose(angle, expected_angle)
    return

def test_torsion_order_parameter_get_value():
    """
    
    """
    op = common_cv.Tiwary_cv_torsion_order_parameter()
    com1 = np.array([0.0, 0.0, 1.0])
    com2 = np.array([0.0, 0.0, 0.0])
    com3 = np.array([1.0, 0.0, 0.0])
    com4 = np.array([1.0, 1.0, 0.0])
    expected_phi = 1.5 * np.pi
    phi = op.get_value([com1, com2, com3, com4])
    assert np.isclose(phi, expected_phi)
    return