"""
test_mmvt_analyze.py
"""
# TODO: fill out these tests when the toy engine is implemented.

import os

import numpy as np
import pytest

import seekr2.modules.mmvt_analyze as mmvt_analyze

TEST_DIRECTORY = os.path.dirname(os.path.realpath(__file__))

test_output_filename = os.path.join(TEST_DIRECTORY, 
                                    "data/test_analyze_outputfile.txt")
test_statistics_filename = os.path.join(TEST_DIRECTORY, 
                                        "data/test_analyze_statistics.txt")

def test_read_output_file():
    N_i_j_alpha, R_i_alpha_total, N_alpha_beta, \
        T_alpha_total, existing_lines \
            = mmvt_analyze.openmm_read_output_file_list(
                [test_output_filename])
    
    N_i_j_alpha_dict1 = N_i_j_alpha
    R_i_alpha_dict1 = R_i_alpha_total
    N_alpha_beta_dict1 = N_alpha_beta
    T_alpha1 = T_alpha_total
    
    N_i_j_alpha_dict2 = {(1, 2): 52, (2, 1): 52}
    R_i_alpha_dict2 = {1: 1658.696, 2: 198.912}
    N_alpha_beta_dict2 = {1: 2423, 2: 98}
    T_alpha2 = 1954.760
    
    for key in N_i_j_alpha_dict1:
        assert key in N_i_j_alpha_dict2
        assert np.isclose(N_i_j_alpha_dict1[key], N_i_j_alpha_dict2[key])
        
    for key in R_i_alpha_dict1:
        assert key in R_i_alpha_dict2
        assert np.isclose(R_i_alpha_dict1[key], R_i_alpha_dict2[key])
        
    for key in N_alpha_beta_dict1:
        assert key in N_alpha_beta_dict2
        assert np.isclose(N_alpha_beta_dict1[key], N_alpha_beta_dict2[key])
    
    assert np.isclose(T_alpha1, T_alpha2)
    return

@pytest.mark.skip()
def test_openmm_read_statistics_file():
    N_i_j_alpha_dict1, R_i_alpha_dict1, N_alpha_beta_dict1, T_alpha1 = \
        mmvt_analyze.openmm_read_statistics_file(test_statistics_filename)
    
    N_i_j_alpha_dict2 = {(1, 2): 52, (2, 1): 52}
    R_i_alpha_dict2 = {1: 1658.696, 2: 198.912}
    N_alpha_beta_dict2 = {1: 2423, 2: 98}
    T_alpha2 = 1954.760
    
    for key in N_i_j_alpha_dict1:
        assert key in N_i_j_alpha_dict2
        assert np.isclose(N_i_j_alpha_dict1[key], N_i_j_alpha_dict2[key])
        
    for key in R_i_alpha_dict1:
        assert key in R_i_alpha_dict2
        assert np.isclose(R_i_alpha_dict1[key], R_i_alpha_dict2[key])
        
    for key in N_alpha_beta_dict1:
        assert key in N_alpha_beta_dict2
        assert np.isclose(N_alpha_beta_dict1[key], N_alpha_beta_dict2[key])
    
    assert np.isclose(T_alpha1, T_alpha2)
    
    return

def test_namd_read_output_file_list():
    pass

def test_init_objects():
    """
    Initialize the MMVT_anchor_statistics() object just to cover it.
    Additional testing of it will be performed with the toy engine.
    """
    stats = mmvt_analyze.MMVT_anchor_statistics(0)
    stats.print_stats()
    return

def test_find_nonzero_matrix_entries():
    """
    Test whether non-diagonal non-zero entries in a matrix will be 
    found.
    """
    M0 = np.zeros((3,3))
    result0 = mmvt_analyze.find_nonzero_matrix_entries(M0)
    assert len(result0) == 0
    
    M1 = np.zeros((3,3))
    M1[0,0] = 2.0
    M1[0,1] = 1.0
    M1[2,1] = 3.0
    result1 = mmvt_analyze.find_nonzero_matrix_entries(M1)
    assert len(result1) == 2
    assert not((0,0) in result1)
    assert (0,1) in result1
    assert (2,1) in result1
    return