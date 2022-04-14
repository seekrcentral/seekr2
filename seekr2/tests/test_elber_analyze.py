"""
test_elber_analyze.py

Testing elber_analyze.py
"""

import os

import numpy as np
import pytest
from shutil import copyfile

import seekr2.modules.common_base as base
import seekr2.modules.elber_base as elber_base
import seekr2.modules.common_analyze as common_analyze
import seekr2.modules.elber_analyze as elber_analyze
from seekr2.tests.conftest import compare_dicts
import seekr2.modules.markov_chain_monte_carlo as markov_chain_monte_carlo

TEST_DIRECTORY = os.path.dirname(os.path.realpath(__file__))

test_output_filename = os.path.join(TEST_DIRECTORY, 
                                    "data/test_analyze_output_elber.txt")

def test_openmm_read_output_file_list():
    N_i_j, R_i_total, lines = elber_analyze.openmm_read_output_file_list(
        [test_output_filename])
    
    N_i_j_exp = {1: 3, 3: 5}
    R_i_exp = 1.324
    for key in N_i_j:
        assert key in N_i_j_exp
        assert np.isclose(N_i_j[key], N_i_j_exp[key])
        
    assert np.isclose(R_i_total, R_i_exp)
    
    N_i_j, R_i_list, lines = elber_analyze.openmm_read_output_file_list(
        [test_output_filename], min_time=1)
    assert len(lines) == 10
    N_i_j, R_i_list, lines = elber_analyze.openmm_read_output_file_list(
        [test_output_filename], max_time=9)
    assert len(lines) == 10
    return

def test_Elber_anchor_statistics_read_output_file_list(toy_elber_model):
    """
    Test the read_output_file_list method of the MMVT_anchor_statistics
    object.
    """
    stats = elber_analyze.Elber_anchor_statistics(1)
    engine = "openmm"
    output_file_list = [test_output_filename]
    min_time = None
    max_time = None
    stats.read_output_file_list(engine, output_file_list, min_time, 
                              max_time, None, None)
    
    N_i_j_exp = {1: 3, 3: 5}
    R_i_exp = 1.324
    for key in stats.N_i_j:
        assert key in N_i_j_exp
        assert np.isclose(stats.N_i_j[key], N_i_j_exp[key])
        
    assert np.isclose(stats.R_i_total, R_i_exp)
    stats.print_stats()
    return

def make_simple_model():
    n_anchors = 3
    n_milestones = 3
    
    # generate data to feed directly into MMVT_data_sample()
    model = base.Model()
    model.calculation_type = "elber"
    model.num_milestones = n_milestones
    model.num_anchors = n_anchors
    anchor0 = elber_base.Elber_toy_anchor()
    anchor0.index = 0
    anchor0.endstate = True
    milestone_0_0 = base.Milestone()
    milestone_0_0.index = 0
    milestone_0_0.neighbor_anchor_index = 0
    milestone_0_0.alias_index = 2
    milestone_0_1 = base.Milestone()
    milestone_0_1.index = 1
    milestone_0_1.neighbor_anchor_index = 1
    milestone_0_1.alias_index = 3
    anchor0.milestones = [milestone_0_0, milestone_0_1]
    
    anchor1 = elber_base.Elber_toy_anchor()
    anchor1.index = 1
    milestone_1_0 = base.Milestone()
    milestone_1_0.index = 0
    milestone_1_0.neighbor_anchor_index = 0
    milestone_1_0.alias_index = 1
    milestone_1_1 = base.Milestone()
    milestone_1_1.index = 1
    milestone_1_1.neighbor_anchor_index = 1
    milestone_1_1.alias_index = 2
    milestone_1_2 = base.Milestone()
    milestone_1_2.index = 2
    milestone_1_2.neighbor_anchor_index = 2
    milestone_1_2.alias_index = 3
    anchor1.milestones = [milestone_1_0, milestone_1_1, milestone_1_2]
    
    anchor2 = elber_base.Elber_toy_anchor()
    anchor2.index = 2
    milestone_2_0 = base.Milestone()
    milestone_2_0.index = 1
    milestone_2_0.neighbor_anchor_index = 1
    milestone_2_0.alias_index = 1
    milestone_2_1 = base.Milestone()
    milestone_2_1.index = 2
    milestone_2_1.neighbor_anchor_index = 2
    milestone_2_1.alias_index = 2
    
    anchor2.milestones = [milestone_2_0, milestone_2_1]
    anchor2.bulkstate = True
    model.anchors = [anchor0, anchor1, anchor2]
    return model

def test_Elber_data_sample_fill_out_data_quantities():
    model = make_simple_model()
    N_i_j_list = [{(0,1): 4}, {(1,0): 4, (1,2): 2}, {(2,1): 2}]
    R_i_list = [2.4, 2.4]
    N_i_j_exp = {(0,1): 4, (1,0): 4, (1,2): 2}
    data_sample = elber_analyze.Elber_data_sample(model, N_i_j_list, R_i_list)
    data_sample.fill_out_data_quantities()
    compare_dicts(N_i_j_exp, data_sample.N_ij)
    assert len(R_i_list) == len(data_sample.R_i)
    for item1, item2 in zip(R_i_list, data_sample.R_i):
        assert item1 == data_sample.R_i[item2]
    return

def test_monte_carlo_milestoning_error():
    model = make_simple_model()
    N_i_j_list = [{(0,1): 4}, {(1,0): 4, (1,2): 2}, {(2,1): 2, (2,3): 2}]
    R_i_list = [2.4, 2.4, 1.2]
    N_i_j_exp = {(0,1): 4, (1,0): 4, (1,2): 2, (2,1): 2}
    data_sample = elber_analyze.Elber_data_sample(model, N_i_j_list, R_i_list)
    data_sample.fill_out_data_quantities()
    data_sample.compute_rate_matrix()
    data_sample.calculate_thermodynamics()
    data_sample.calculate_kinetics()
    
    # Elber matrix sampler
    num = 1000
    stride = 10
    skip = 10
    n_milestones = 3
    n_anchors = 3
    data_sample_list, p_i_error, free_energy_profile_err, MFPTs_error, \
        k_off_error, k_ons_error = elber_analyze.monte_carlo_milestoning_error(
        data_sample, num=num, stride=stride, skip=skip, verbose=True)
    
    result_q1_distribution = []
    result_q2_distribution = []
    for data_sample in data_sample_list:
        result_q1_distribution.append(data_sample.Q[0,1])
        result_q2_distribution.append(data_sample.Q[1,0])
    
    N_i_j = {(0,1): 4, (1,0): 4, (1,2): 2, (2,1): 2}
    R_i = {0: 2.4, 1: 2.4, 2: 1.2}
    elber_N = np.array([[0, 4, 0],
                        [4, 0, 2],
                        [0, 2, 0]])
    elber_R = np.array([[2.4],
                        [2.4],
                        [1.2]])
    
    elber_Q = np.zeros((n_milestones, n_milestones))
    for i in range(n_milestones):
        for j in range(n_milestones):
            key = (i,j)
            if key in N_i_j:
                elber_Q[i,j] = N_i_j[key] / R_i[i]
    
    for i in range(n_milestones):
        elber_Q[i,i] = -np.sum(elber_Q[i,:])
    
    # Elber matrix sampler
    elber_q1_distribution = []
    elber_q2_distribution = []
    for counter in range(num * (stride) + skip):
        #if verbose: print("MCMC stepnum: ", counter)
        elber_Qnew = markov_chain_monte_carlo\
            .irreversible_stochastic_matrix_algorithm_sample(
                elber_Q, elber_N, elber_R)
        if counter > skip and counter % stride == 0:
            elber_q1_distribution.append(elber_Q[0,1])
            elber_q2_distribution.append(elber_Q[1,0])
        
        elber_Q = elber_Qnew
    
    assert np.isclose(np.mean(elber_q1_distribution), 
                      np.mean(result_q1_distribution), rtol=0.5, atol=0.01)
    assert np.isclose(np.std(elber_q1_distribution), 
                      np.std(result_q1_distribution), rtol=0.5, atol=0.01)
    assert np.isclose(np.mean(elber_q2_distribution), 
                      np.mean(result_q2_distribution), rtol=0.5, atol=0.01)
    assert np.isclose(np.std(elber_q2_distribution), 
                      np.std(result_q2_distribution), rtol=0.5, atol=0.01)
    return