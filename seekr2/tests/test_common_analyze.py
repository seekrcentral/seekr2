"""
test_common_analyze.py

Testing common_analyze.py
"""

import os

import numpy as np
import pytest

import seekr2.modules.common_analyze as common_analyze
from seekr2.tests.conftest import compare_dicts

TEST_DIRECTORY = os.path.dirname(__file__)

def test_Q_to_K():
    """
    Test the conversion of a rate matrix to a probability transition 
    matrix.
    """
    Q = np.array(
        [[-0.5, 0.5, 0.0, 0.0],
         [0.1, -0.3, 0.2, 0.0],
         [0.0, 0.15, -0.3, 0.15],
         [0.0, 0.0, 0.3, -0.4]])
    expected_K = np.array(
        [[0.0, 1.0, 0.0, 0.0],
         [1.0/3.0, 0.0, 2.0/3.0, 0.0],
         [0.0, 0.5, 0.0, 0.5],
         [0.0, 0.0, 0.75, 0.0]])
    K = common_analyze.Q_to_K(Q)
    error = np.linalg.norm(expected_K - K)
    assert error < 1e-8
    return
    
@pytest.mark.parametrize("a, b, c",
                         [(2.0, 0.0, 2.0),
                          (0.0, 3.5, 3.5),
                          (3.0, 4.0, 5.0),
                          (0.125, 0.0036, 0.1250518293)])
def test_quadriture(a, b, c):
    c_output = common_analyze.quadriture(a, b)
    assert np.isclose(c, c_output)
    return
    
@pytest.mark.parametrize("before, index, after",
                         [(np.array([1,2,3,4]), 0, np.array([2,3,4])),
                          (np.array([1,2,3,4]), 1, np.array([1,3,4])),
                          (np.array([1,2,3,4]), 3, np.array([1,2,3])),
                          ])
def test_minor1d(before, index, after):
    """
    Test the function which removes an element from a 1d array.
    """
    output = common_analyze.minor1d(before, index)
    assert np.all(output == after)
    return
    
@pytest.mark.parametrize("before, index_i, index_j, after",
                         [(np.array([[1,2,3],[4,5,6],[7,8,9]]), 0, 0, 
                           np.array([[5,6],[8,9]])),
                          (np.array([[1,2,3],[4,5,6],[7,8,9]]), 1, 0, 
                           np.array([[2,3],[8,9]])),
                          (np.array([[1,2,3],[4,5,6],[7,8,9]]), 1, 2, 
                           np.array([[1,2],[7,8]])),
                          ])
def test_minor2d(before, index_i, index_j, after):
    """
    Test the function which removes both a single row and a single 
    column from a 2x2 array.
    """
    output = common_analyze.minor2d(before, index_i, index_j)
    assert np.all(output == after)
    return
    
def test_pretty_string_value_error():
    """
    Test the function to print 'pretty' scientific notation to the
    terminal.
    """
    mystr = common_analyze.pretty_string_value_error(
        5.6e-2, 2.0e-3, error_digits=1, use_unicode=False)
    expectedstr = "(5.6 +/- 0.2) * 10^-02"
    assert(mystr == expectedstr)
    mystr = common_analyze.pretty_string_value_error(
        5.6e-2, 2.0e-1, error_digits=1, use_unicode=False)
    expectedstr = "(5.6 +/- 20.0) * 10^-02"
    assert(mystr == expectedstr)
    mystr = common_analyze.pretty_string_value_error(
        1.23456789e8, 4.5678e5, error_digits=2, use_unicode=False)
    expectedstr = "(1.2346 +/- 0.0046) * 10^+08"
    assert(mystr == expectedstr)
    return

def test_browndye_run_compute_rate_constant():
    """
    Test the function which extracts Browndye results from the 
    b-surface simulations.
    """
    test_bd_results_filename = os.path.join(TEST_DIRECTORY, 
                                            "data/sample_bd_results_file.xml")
    k_ons, k_on_errors, reaction_probabilities, \
        reaction_probability_errors, transition_counts \
        = common_analyze.browndye_run_compute_rate_constant(
            "compute_rate_constant", [test_bd_results_filename], 
            sample_error_from_normal=False)
    
    expected_k_ons = {10: 1.18338e+09, "escaped": 2.480783458e+10, "stuck": 0}
    compare_dicts(k_ons, expected_k_ons)
    expected_k_on_errors = {10: 1.7539737e+07, "escaped": 8.0298973e+07, 
                            "stuck": 1e+99}
    compare_dicts(k_on_errors, expected_k_on_errors)
    expected_reaction_probabilities = {10: 0.04553, "escaped": 0.95447, 
                                        "stuck": 0.0}
    compare_dicts(reaction_probabilities, expected_reaction_probabilities)
    expected_reaction_probability_errors = {10: 0.0006748333290290023, 
                                            "escaped": 0.0030894659734406714, 
                                            "stuck": 1e+99}
    compare_dicts(reaction_probability_errors, 
                  expected_reaction_probability_errors)
    expected_transition_counts = {10: 4553, "escaped": 95447, "stuck": 0}
    compare_dicts(transition_counts, expected_transition_counts)
    return

def test_init_Data_sample():
    """
    Initialize Data_sample() to get more code coverage, though this
    initialize is never likely to be directly called because the
    Data_sample() class will only be inherited.
    """
    dummy = common_analyze.Data_sample(None)
    return

# CANNOT MAKE FULL UNIT TESTS OF DATA_SAMPLE OBJECT UNTIL TOY ENGINE 
# IMPLEMENTED
# TODO: use toy engines to make full suite of data sample tests.