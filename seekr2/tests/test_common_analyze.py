"""
test_common_analyze.py

Testing common_analyze.py
"""

import os

import numpy as np
import pytest
import shutil

import seekr2.modules.common_analyze as common_analyze
from seekr2.tests.conftest import compare_dicts

TEST_DIRECTORY = os.path.dirname(__file__)

def test_solve_rate_matrix():
    """
    Test the alternative, more stable way of solving a rate matrix 
    used in SEEKR2.
    """
    Q = np.array(
        [[-0.5, 0.5, 0.0, 0.0],
         [0.1, -0.3, 0.2, 0.0],
         [0.0, 0.15, -0.3, 0.15],
         [0.0, 0.0, 0.3, -0.4]])
    
    K = np.zeros(Q.shape, dtype=np.longdouble)
    for i in range(Q.shape[0]):
        for j in range(Q.shape[0]):
            if i == j:
                K[i,j] = 0.0
            else:
                K[i,j] = -Q[i,j] / Q[i,i]
         
    for i in range(K.shape[0]-1):
        my_sum = sum(K[i,:])
        for j in range(K.shape[0]):
            K[i,j] = K[i,j] / my_sum
            
    test_times_1 = common_analyze.solve_rate_matrix(Q)
    
    one_vector = np.ones((Q.shape[0]))
    test_times_2 = np.linalg.solve(Q, -one_vector)
    
    error = np.linalg.norm(test_times_2 - test_times_1)
    assert error < 1e-8
    return

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

# TODO: remove test since it's for old BD milestoning way
@pytest.mark.skip()
def test_browndye_parse_bd_milestone_results():
    """
    Test the function which extracts Browndye results from the BD
    milestone.
    """
    test_bd_results_filename = os.path.join(
        TEST_DIRECTORY, "data/sample_bd_milestone_results.xml")
    transition_probabilities, transition_counts \
        = common_analyze.browndye_parse_bd_milestone_results(
            [test_bd_results_filename], sample_error_from_normal=False)
    expected_transition_probabilities = {9: 0.6, 
                                         "escaped": 0.4, 
                                         "stuck": 0}
    compare_dicts(transition_probabilities, 
                  expected_transition_probabilities)
    expected_transition_counts = {9: 600, "escaped": 400, "stuck": 0}
    compare_dicts(transition_counts, expected_transition_counts)

#TODO: modify for writing FHPD structures without old BD milestoning way
@pytest.mark.skip()
def test_combine_fhpd_results(host_guest_mmvt_model):
    test_bd_results_filename = os.path.join(
        TEST_DIRECTORY, "data/sample_bd_milestone_results.xml")
    lig_filenames = ["lig1_0_0", "lig1_0_1", "lig1_0_2"]
    bd_milestone_dir = os.path.join(
        host_guest_mmvt_model.anchor_rootdir,
        host_guest_mmvt_model.k_on_info.bd_milestones[0].directory)
    fhpd_dir = os.path.join(
        bd_milestone_dir,
        tryp_ben_mmvt_model.k_on_info.bd_milestones[0].fhpd_directory)
    if not os.path.exists(fhpd_dir):
        os.mkdir(fhpd_dir)
    lig_filenames_full_paths = []
    for lig_filename in lig_filenames:
        lig_filename_full_path = os.path.join(fhpd_dir, lig_filename)
        lig_filenames_full_paths.append(lig_filename_full_path)
        if not os.path.exists(lig_filename_full_path):
            os.mkdir(lig_filename_full_path)
        dest_filename = os.path.join(lig_filename_full_path, "results1.xml")
        shutil.copyfile(test_bd_results_filename, dest_filename)
    
    combined_results_filename = os.path.join(bd_milestone_dir, "results.xml")
    common_analyze.combine_fhpd_results(
        tryp_ben_mmvt_model.k_on_info.bd_milestones[0], 
        lig_filenames_full_paths, 
        combined_results_filename)
    
    transition_probabilities, transition_counts \
        = common_analyze.browndye_parse_bd_milestone_results(
            [combined_results_filename], sample_error_from_normal=False)
    expected_transition_probabilities = {9: 0.6, 
                                         "escaped": 0.4, 
                                         "stuck": 0}
    compare_dicts(transition_probabilities, 
                  expected_transition_probabilities)
    expected_transition_counts = {9: 1800, "escaped": 1200, "stuck": 0}
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