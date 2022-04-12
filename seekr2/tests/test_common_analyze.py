"""
test_common_analyze.py

Testing common_analyze.py
"""

import os

import numpy as np
import pytest
from shutil import copyfile

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
    
    expected_k_ons = {10: 1.18338e+09, "escaped": 2.480783458e+10, "stuck": 0, 
                      "total": 2.599121458e+10}
    compare_dicts(k_ons, expected_k_ons)
    expected_k_on_errors = {10: 1.7539737e+07, "escaped": 8.0298973e+07, 
                            "stuck": 1e+99, "total":8.21918481e+7}
    compare_dicts(k_on_errors, expected_k_on_errors)
    expected_reaction_probabilities = {10: (4553/100000), 
                                       "escaped": (95447/100000), 
                                       "stuck": 0.0, "total": 1.0}
    compare_dicts(reaction_probabilities, expected_reaction_probabilities)
    expected_reaction_probability_errors = {10: 0.0006748333290290023, 
                                            "escaped": 0.0030894659734406714, 
                                            "stuck": 1e+99,
                                            "total": 0.0031622934716752666}
    compare_dicts(reaction_probability_errors, 
                  expected_reaction_probability_errors)
    expected_transition_counts = {10: 4553, "escaped": 95447, "stuck": 0,
                                  "total": 100000}
    compare_dicts(transition_counts, expected_transition_counts)
    return

def test_Data_sample_parse_browndye_results(host_guest_mmvt_model):
    test_bd_results_filename = os.path.join(
        TEST_DIRECTORY, "data/sample_bd_results_file2.xml")
    b_surface_directory = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, "b_surface")
    b_surface_results_file = os.path.join(b_surface_directory, "results1.xml")
    data_sample = common_analyze.Data_sample(host_guest_mmvt_model)
    copyfile(test_bd_results_filename, b_surface_results_file)
    data_sample.parse_browndye_results()
    assert len(data_sample.bd_transition_counts) == 2
    expected_counts_b_surface = {"total": 100000, "escaped": 54070, 12: 65239,
                                 11: 45930, "stuck": 0}
    counts_b_surface = data_sample.bd_transition_counts["b_surface"]
    compare_dicts(counts_b_surface, expected_counts_b_surface)
    expected_counts_bd_milestone0 = {"total": 65239, "escaped": 65239-45930,
                                     11: 45930}
    counts_bd_milestone0 = data_sample.bd_transition_counts[0]
    compare_dicts(counts_bd_milestone0, expected_counts_bd_milestone0)
    expected_b_surface_probabilities = {"total": 1.0, "escaped": 54070/100000,
                                        12: 65239/100000, 11: 45930/100000, 
                                        "stuck": 0.0}
    b_surface_probabilities \
        = data_sample.bd_transition_probabilities["b_surface"]
    compare_dicts(b_surface_probabilities, expected_b_surface_probabilities)
    k_on = data_sample.b_surface_k_ons_src["total"]
    expected_k_ons = {"total": k_on, "escaped": k_on * 54070/100000,
                      "stuck": 0.0, 12: k_on * 65239/100000, 
                      11: k_on * 45930/100000}
    k_ons = data_sample.b_surface_k_ons_src
    compare_dicts(expected_k_ons, k_ons)
    k_on_err = data_sample.b_surface_b_surface_k_on_errors_src
    expected_k_ons_err = {"total": k_on / np.sqrt(100000-1), 
                          "escaped": k_on * (54070/100000) / np.sqrt(54070-1),
                          "stuck": 1e99, 
                          12: k_on * (65239/100000) / np.sqrt(65239-1), 
                          11: k_on * (45930/100000) / np.sqrt(45930-1)}
    compare_dicts(k_on_err, expected_k_ons_err)

def test_Data_sample_compute_rate_matrix(toy_mmvt_model):
    data_sample = common_analyze.Data_sample(toy_mmvt_model)
    n = toy_mmvt_model.num_milestones
    data_sample.N_ij = np.zeros((n,n))
    data_sample.R_i = np.zeros(n)
    for i in range(n):
        data_sample.R_i[i] = 0.5
        if i > 0:
            data_sample.N_ij[i, i-1] = 10
        if i < n-1:
            data_sample.N_ij[i, i+1] = 10
    
    data_sample.compute_rate_matrix()
    for i in range(n):
        if i > 0:
            assert data_sample.Q[i, i-1] == 10 / 0.5
        else:
            assert data_sample.K[i, i+1] == 1.0
            
        if i < n-1:
            assert data_sample.Q[i, i+1] == 10 / 0.5
        else:
            assert data_sample.K[i, i-1] == 1.0
        
        if (i > 0) and (i < n-1):
            assert data_sample.Q[i,i] == -20 / 0.5
            assert data_sample.K[i, i+1] == 0.5
            assert data_sample.K[i, i-1] == 0.5

        else:
            assert data_sample.Q[i,i] == -10 / 0.5
        
        assert data_sample.K[i,i] == 0.0
        for j in range(n):
            if (j == i) or (j == i-1) or (j == i+1):
                continue
            else:
                assert data_sample.Q[i,j] == 0.0
                assert data_sample.K[i,j] == 0.0
    return

def fill_out_simple_statistics(model, data_sample):
    n = model.num_milestones
    data_sample.Q = np.zeros((n,n))
    data_sample.N_ij = np.zeros((n,n))
    data_sample.R_i = np.zeros(n)
    for i in range(n):
        if i == 0:
            data_sample.Q[i, i+1] = 20.
            data_sample.N_ij[i, i+1] = 10
            data_sample.R_i[i] = 0.5
            data_sample.Q[i, i] = -20.
        elif i == n-1:
            data_sample.Q[i, i-1] = 20.
            data_sample.N_ij[i, i-1] = 10
            data_sample.R_i[i] = 0.5
            data_sample.Q[i, i] = -20.
        else:
            data_sample.Q[i, i-1] = 20.
            data_sample.Q[i, i+1] = 20.
            data_sample.N_ij[i, i+1] = 10
            data_sample.N_ij[i, i-1] = 10
            data_sample.R_i[i] = 0.5
            data_sample.Q[i, i] = -40.
    
    data_sample.K = common_analyze.Q_to_K(data_sample.Q)
    return

def test_Data_sample_calculate_thermodynamics(toy_mmvt_model):
    data_sample = common_analyze.Data_sample(toy_mmvt_model)
    fill_out_simple_statistics(toy_mmvt_model, data_sample)
    data_sample.calculate_thermodynamics()
    # assert all values of p_i are equal
    assert np.all(np.isclose(data_sample.p_i, data_sample.p_i[0]))
    assert np.all(np.isclose(data_sample.free_energy_profile, 0.0))
    return
    
def test_Data_sample_calculate_kinetics(toy_mmvt_model):
    data_sample = common_analyze.Data_sample(toy_mmvt_model)
    fill_out_simple_statistics(toy_mmvt_model, data_sample)
    data_sample.calculate_thermodynamics()
    data_sample.calculate_kinetics()
    assert data_sample.MFPTs is not None
    assert data_sample.k_off is not None
    return