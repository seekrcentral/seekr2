"""
test_mmvt_analyze.py
"""
# TODO: fill out these tests when the toy engine is implemented.

import os

import numpy as np
import pytest

import seekr2.modules.common_base as base
import seekr2.modules.mmvt_base as mmvt_base
import seekr2.modules.runner_openmm as runner_openmm
import seekr2.modules.mmvt_analyze as mmvt_analyze
from seekr2.tests.conftest import compare_dicts
import seekr2.run as run
import seekr2.analyze as analyze

TEST_DIRECTORY = os.path.dirname(os.path.realpath(__file__))

test_output_filename = os.path.join(TEST_DIRECTORY, 
                                    "data/test_analyze_outputfile.txt")
test_output_filename_restart1 = os.path.join(TEST_DIRECTORY, 
                                    "data/test_analyze_outputfile_restart1.txt")
test_output_filename_restart2 = os.path.join(TEST_DIRECTORY, 
                                    "data/test_analyze_outputfile_restart2.txt")
test_output_filename_checkpoint1 = os.path.join(TEST_DIRECTORY, 
                                "data/test_analyze_outputfile_checkpoint1.txt")
test_output_filename_checkpoint2 = os.path.join(TEST_DIRECTORY, 
                                "data/test_analyze_outputfile_checkpoint2.txt")
test_statistics_filename = os.path.join(TEST_DIRECTORY, 
                                        "data/test_analyze_statistics.txt")

test_output_filename_namd = os.path.join(TEST_DIRECTORY, 
                                    "data/test_analyze_outputfile_namd.txt")
test_output_filename_namd_restart1 = os.path.join(TEST_DIRECTORY, 
                            "data/test_analyze_outputfile_namd_restart1.txt")
test_output_filename_namd_restart2 = os.path.join(TEST_DIRECTORY, 
                            "data/test_analyze_outputfile_namd_restart2.txt")
test_output_filename_namd_checkpoint1 = os.path.join(TEST_DIRECTORY, 
                            "data/test_analyze_outputfile_namd_checkpoint1.txt")
test_output_filename_namd_checkpoint2 = os.path.join(TEST_DIRECTORY, 
                            "data/test_analyze_outputfile_namd_checkpoint2.txt")

def test_flux_matrix_to_K():
    """
    Test the conversion of a rate matrix to a probability transition 
    matrix.
    """
    M = np.array(
        [[-0.5, 0.5, 0.0, 0.0],
         [0.1, -0.3, 0.2, 0.0],
         [0.0, 0.15, -0.3, 0.15],
         [0.0, 0.0, 0.3, -0.4]])
    expected_K = np.array(
        [[0.0, 1.0, 0.0],
         [1.0/3.0, 0.0, 2.0/3.0],
         [0.0, 0.5, 0.0]])
    K = mmvt_analyze.flux_matrix_to_K(M)
    error = np.linalg.norm(expected_K - K)
    assert error < 1e-8
    return

def test_make_new_Nij_alpha():
    Q = np.array(
        [[-0.5, 0.5, 0.0, 0.0],
         [0.1, -0.3, 0.2, 0.0],
         [0.0, 0.15, -0.3, 0.15],
         [0.0, 0.0, 0.3, -0.4]])
    R = np.array(
        [[0.1],
         [0.2],
         [0.3],
         [0.5]])
    Nij = mmvt_analyze.make_new_Nij_alpha(Q, R)
    expected_Nij = np.array(
        [[0.0, 0.05, 0.0, 0.0],
         [0.02, 0.0, 0.04, 0.0],
         [0.0, 0.045, 0.0, 0.045],
         [0.0, 0.0, 0.15, 0.0]])
    assert np.all(np.isclose(Nij, expected_Nij))

def test_openmm_read_output_file_list_normal():
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

def test_openmm_read_output_file_list_restart():
    N_i_j_alpha, R_i_alpha_total, N_alpha_beta, \
        T_alpha_total, existing_lines \
            = mmvt_analyze.openmm_read_output_file_list(
                [test_output_filename_restart1, test_output_filename_restart2])
    assert len(existing_lines) == 10
    return

def test_openmm_read_output_file_list_checkpoint():
    N_i_j_alpha, R_i_alpha_total, N_alpha_beta, \
        T_alpha_total, existing_lines \
            = mmvt_analyze.openmm_read_output_file_list(
                [test_output_filename_checkpoint1, 
                 test_output_filename_checkpoint2])
    assert len(existing_lines) == 10
    return

def test_openmm_read_output_file_list_min_time_max_time():
    min_time = 45.0
    N_i_j_alpha, R_i_alpha_total, N_alpha_beta, \
        T_alpha_total, existing_lines \
            = mmvt_analyze.openmm_read_output_file_list(
                [test_output_filename], min_time=min_time)
    N_alpha_beta_dict1 = N_alpha_beta
    N_alpha_beta_dict2 = {1: 2418, 2: 97}
    for key in N_alpha_beta_dict1:
        assert key in N_alpha_beta_dict2
        assert np.isclose(N_alpha_beta_dict1[key], N_alpha_beta_dict2[key])
        
    max_time = 1993.0
    N_i_j_alpha, R_i_alpha_total, N_alpha_beta, \
        T_alpha_total, existing_lines \
            = mmvt_analyze.openmm_read_output_file_list(
                [test_output_filename], max_time=max_time)
    N_alpha_beta_dict1 = N_alpha_beta
    N_alpha_beta_dict2 = {1: 2420, 2: 98}
    for key in N_alpha_beta_dict1:
        assert key in N_alpha_beta_dict2
        assert np.isclose(N_alpha_beta_dict1[key], N_alpha_beta_dict2[key])

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

def test_namd_read_output_file_list(toy_mmvt_model):
    N_i_j_alpha, R_i_alpha_total, N_alpha_beta, T_alpha_total, \
        existing_lines = mmvt_analyze.namd_read_output_file_list(
            [test_output_filename_namd], toy_mmvt_model.anchors[1],
            0.002)
    
    N_i_j_alpha_dict1 = N_i_j_alpha
    R_i_alpha_dict1 = R_i_alpha_total
    N_alpha_beta_dict1 = N_alpha_beta
    T_alpha1 = T_alpha_total
    
    N_i_j_alpha_dict2 = {(1, 2): 1, (2, 1): 1}
    R_i_alpha_dict2 = {1: 0.01, 2: 0.01}
    N_alpha_beta_dict2 = {1: 8, 2: 1}
    T_alpha2 = 0.09
    
    assert len(existing_lines) == 10
    
    for key in N_i_j_alpha_dict1:
        assert key in N_i_j_alpha_dict2
        assert np.isclose(N_i_j_alpha_dict1[key], N_i_j_alpha_dict2[key])
    
    for key in N_alpha_beta_dict1:
        assert key in N_alpha_beta_dict2
        assert np.isclose(N_alpha_beta_dict1[key], N_alpha_beta_dict2[key])
    
    for key in R_i_alpha_dict1:
        assert key in R_i_alpha_dict2
        assert np.isclose(R_i_alpha_dict1[key], R_i_alpha_dict2[key])
        
    assert np.isclose(T_alpha1, T_alpha2)
    return

def test_namd_read_output_file_list_restart(toy_mmvt_model):
    N_i_j_alpha, R_i_alpha_total, N_alpha_beta, T_alpha_total, \
        existing_lines = mmvt_analyze.namd_read_output_file_list(
            [test_output_filename_namd_restart1, 
             test_output_filename_namd_restart2], 
            toy_mmvt_model.anchors[1], 0.002)
    
    assert len(existing_lines) == 10
    return

def test_namd_read_output_file_list_checkpoint(toy_mmvt_model):
    N_i_j_alpha, R_i_alpha_total, N_alpha_beta, T_alpha_total, \
        existing_lines = mmvt_analyze.namd_read_output_file_list(
            [test_output_filename_namd_checkpoint1, 
             test_output_filename_namd_checkpoint2], 
            toy_mmvt_model.anchors[1], 0.002)
    assert len(existing_lines) == 10
    return

def test_namd_read_output_file_list_min_time_max_time(toy_mmvt_model):
    min_time = 120 * 0.002
    N_i_j_alpha, R_i_alpha_total, N_alpha_beta, T_alpha_total, \
        existing_lines = mmvt_analyze.namd_read_output_file_list(
            [test_output_filename_namd], toy_mmvt_model.anchors[1],
            0.002, min_time=min_time)
    N_alpha_beta_dict1 = N_alpha_beta
    N_alpha_beta_dict2 = {1: 7, 2: 0}
    for key in N_alpha_beta_dict1:
        assert key in N_alpha_beta_dict2
        assert np.isclose(N_alpha_beta_dict1[key], N_alpha_beta_dict2[key])
        
    max_time = 150 * 0.002
    N_i_j_alpha, R_i_alpha_total, N_alpha_beta, T_alpha_total, \
        existing_lines = mmvt_analyze.namd_read_output_file_list(
            [test_output_filename_namd], toy_mmvt_model.anchors[1],
            0.002, max_time=max_time)
    N_alpha_beta_dict1 = N_alpha_beta
    N_alpha_beta_dict2 = {1: 7, 2: 1}
    for key in N_alpha_beta_dict1:
        assert key in N_alpha_beta_dict2
        assert np.isclose(N_alpha_beta_dict1[key], N_alpha_beta_dict2[key])
    return

def test_MMVT_anchor_statistics_read_output_file_list(toy_mmvt_model):
    """
    Test the read_output_file_list method of the MMVT_anchor_statistics
    object.
    """
    stats = mmvt_analyze.MMVT_anchor_statistics(1)
    engine = "openmm"
    output_file_list = [test_output_filename]
    min_time = None
    max_time = None
    anchor = toy_mmvt_model.anchors[1]
    timestep = toy_mmvt_model.get_timestep()
    stats.read_output_file_list(engine, output_file_list, min_time, 
                              max_time, anchor, timestep)
    
    N_i_j_alpha_dict2 = {(1, 2): 52, (2, 1): 52}
    R_i_alpha_dict2 = {1: 1658.696, 2: 198.912}
    N_alpha_beta_dict2 = {1: 2423, 2: 98}
    T_alpha2 = 1954.760
    k_alpha_beta_dict2 = {1: 2423/1954.760, 2: 98/1954.760}
    compare_dicts(N_i_j_alpha_dict2, stats.N_i_j_alpha)
    compare_dicts(R_i_alpha_dict2, stats.R_i_alpha_total)
    compare_dicts(N_alpha_beta_dict2, stats.N_alpha_beta)
    assert np.isclose(T_alpha2, stats.T_alpha_total)
    compare_dicts(k_alpha_beta_dict2, stats.k_alpha_beta)
    stats.print_stats()
    return

def make_mmvt_data_sample(model):
    num_steps = 100000
    model.openmm_settings.cuda_platform_settings = None
    model.openmm_settings.reference_platform = True
    model.openmm_settings.langevin_integrator.friction_coefficient \
        = 100.0
    model.calculation_settings.num_production_steps = num_steps
    model.calculation_settings.energy_reporter_interval = num_steps
    for anchor in model.anchors:
        runner_openmm.cleanse_anchor_outputs(
            model, anchor, skip_umbrella_files=False)
    run.run(model, "1", force_overwrite = True)
    analysis = analyze.Analysis(model)
    analysis.extract_data()
    analysis.fill_out_data_samples()
    return analysis.main_data_sample

def test_MMVT_data_sample_init(toy_mmvt_model):
    data_sample = make_mmvt_data_sample(toy_mmvt_model)
    for alpha, anchor in enumerate(toy_mmvt_model.anchors):
        if anchor.bulkstate:
            continue
        if alpha == 1:
            assert data_sample.N_alpha[alpha] > 0
            # Because it's low-flow - it's been sampled
            assert data_sample.k_alpha[alpha] < 1000.0
        else:
            assert data_sample.N_alpha[alpha] == 0
            # high flow - it hasn't been sampled
            assert data_sample.k_alpha[alpha] >= 1000.0
    return

def make_simple_model():
    n_anchors = 4
    n_milestones = 3
    
    # generate data to feed directly into MMVT_data_sample()
    model = base.Model()
    model.num_milestones = n_milestones
    model.num_anchors = n_anchors
    anchor0 = mmvt_base.MMVT_toy_anchor()
    anchor0.index = 0
    milestone_0_0 = base.Milestone()
    milestone_0_0.index = 0
    milestone_0_0.neighbor_anchor_index = 1
    milestone_0_0.alias_index = 1
    anchor0.milestones = [milestone_0_0]
    
    anchor1 = mmvt_base.MMVT_toy_anchor()
    anchor1.index = 1
    milestone_1_0 = base.Milestone()
    milestone_1_0.index = 0
    milestone_1_0.neighbor_anchor_index = 0
    milestone_1_0.alias_index = 1
    milestone_1_1 = base.Milestone()
    milestone_1_1.index = 1
    milestone_1_1.neighbor_anchor_index = 2
    milestone_1_1.alias_index = 2
    anchor1.milestones = [milestone_1_0, milestone_1_1]
    
    anchor2 = mmvt_base.MMVT_toy_anchor()
    anchor2.index = 2
    milestone_2_0 = base.Milestone()
    milestone_2_0.index = 1
    milestone_2_0.neighbor_anchor_index = 1
    milestone_2_0.alias_index = 1
    milestone_2_1 = base.Milestone()
    milestone_2_1.index = 2
    milestone_2_1.neighbor_anchor_index = 3
    milestone_2_1.alias_index = 2
    anchor2.milestones = [milestone_2_0, milestone_2_1]
    
    anchor3 = mmvt_base.MMVT_toy_anchor()
    anchor3.index = 3
    milestone_3_0 = base.Milestone()
    milestone_3_0.index = 2
    milestone_3_0.neighbor_anchor_index = 2
    milestone_3_0.alias_index = 1
    anchor3.milestones = [milestone_3_0]
    
    model.anchors = [anchor0, anchor1, anchor2, anchor3]
    return model

def test_MMVT_data_sample_calculate_pi_alpha():
    model = make_simple_model()
    N_alpha_beta = {(0,1):10, (1,0):10,
                    (1,2):10, (2,1):10,
                    (2,3):10,  (3,2):10}
    k_alpha_beta = {(0,1):4.0, (1,0):4.0,
                    (1,2):4.0,  (2,1):4.0,
                    (2,3):4.0, (3,2):4.0}
    N_i_j_alpha = [{},
                   {(0,1):3, (1,0):3}, 
                   {(1,2):3, (2,1):3},
                   {}]
    R_i_alpha_total = [{0: 8.0},
                       {0: 8.0, 1:8.0},
                       {1: 8.0, 2:8.0},
                       {2: 8.0}]
    T_alpha_total = [2.5,
                     2.5,
                     2.5,
                     2.5]
    
    main_data_sample = mmvt_analyze.MMVT_data_sample(
            model, N_alpha_beta, k_alpha_beta, N_i_j_alpha, 
            R_i_alpha_total, T_alpha_total)
    main_data_sample.calculate_pi_alpha()
    assert np.all(np.isclose(main_data_sample.pi_alpha[:-1], 
                             main_data_sample.pi_alpha[0]))
    assert np.isclose(main_data_sample.pi_alpha[-1], 0.0)
    return

def test_MMVT_data_sample_fill_out_data_quantities():
    model = make_simple_model()
    N_alpha_beta = {(0,1):12, (1,0):12,
                    (1,2):12, (2,1):12,
                    (2,3):6,  (3,2):6}
    k_alpha_beta = {(0,1):20.0, (1,0):10.0,
                    (1,2):10.0,  (2,1):(40.0/3.0),
                    (2,3):(20.0/3.0), (3,2):20.0}
    N_i_j_alpha = [{},
                   {(0,1):4, (1,0):4}, 
                   {(1,2):2, (2,1):2},
                   {}]
    R_i_alpha_total = [{0: 1.2},
                       {0: 1.2, 1:1.2},
                       {1: 1.2, 2:0.6},
                       {2: 0.6}]
    T_alpha_total = [1.2,
                     2.4,
                     1.8,
                     0.6]
    main_data_sample = mmvt_analyze.MMVT_data_sample(
            model, N_alpha_beta, k_alpha_beta, N_i_j_alpha, 
            R_i_alpha_total, T_alpha_total)
    main_data_sample.calculate_pi_alpha()
    main_data_sample.fill_out_data_quantities()
    N_ij = {(0,1): 1, (1,0): 1, (1,2): 0.5, (2,1): 0.5}
    R_i = {0: 0.6, 1: 0.6, 2: 0.3}
    compare_dicts(N_ij, main_data_sample.N_ij)
    compare_dicts(R_i, main_data_sample.R_i)
    return

def test_MMVT_data_sample_make_mcmc_quantities():
    model = make_simple_model()
    N_alpha_beta = {(0,1):12, (1,0):12,
                    (1,2):12, (2,1):12,
                    (2,3):6,  (3,2):6}
    k_alpha_beta = {(0,1):20.0, (1,0):10.0,
                    (1,2):10.0,  (2,1):(40.0/3.0),
                    (2,3):(20.0/3.0), (3,2):20.0}
    N_i_j_alpha = [{},
                   {(0,1):4, (1,0):4}, 
                   {(1,2):2, (2,1):2},
                   {}]
    R_i_alpha_total = [{0: 1.2},
                       {0: 1.2, 1:1.2},
                       {1: 1.2, 2:0.6},
                       {2: 0.6}]
    T_alpha_total = [1.2,
                     2.4,
                     1.8,
                     0.6]
    main_data_sample = mmvt_analyze.MMVT_data_sample(
            model, N_alpha_beta, k_alpha_beta, N_i_j_alpha, 
            R_i_alpha_total, T_alpha_total)
    main_data_sample.calculate_pi_alpha()
    k_alpha_beta_matrix, N_alpha_beta_matrix, T_alpha_matrix,\
        mmvt_Nij_alpha, mmvt_Ri_alpha, mmvt_Qij_alpha, T \
        = main_data_sample.make_mcmc_quantities()
    k_alpha_beta_matrix_exp = np.array([
        [-20.0, 20.0, 0.0, 0.0],
        [10.0, -20.0, 10.0, 0.0],
        [0.0, (40.0/3.0), -(40.0/3.0)-(20.0/3.0), (20.0/3.0)],
        [0.0, 0.0, 20.0, -20.0]])
    N_alpha_beta_matrix_exp = np.array([
        [0, 12, 0, 0],
        [12, 0, 12, 0],
        [0, 12, 0, 6],
        [0, 0, 6, 0]])
    T_alpha_matrix_exp = np.array([T_alpha_total]).T
    mmvt_Nij_alpha_exp = [
        np.array([
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]]),
        np.array([
        [0, 4, 0],
        [4, 0, 0],
        [0, 0, 0]]),
        np.array([
        [0, 0, 0],
        [0, 0, 2],
        [0, 2, 0]]),
        np.array([
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]])]
    mmvt_Ri_alpha_exp = [
        np.array([
        [1.2],
        [0],
        [0]]),
        np.array([
        [1.2],
        [1.2],
        [0]]),
        np.array([
        [0],
        [1.2],
        [0.6]]),
        np.array([
        [0],
        [0],
        [0.6]])]
    mmvt_Qij_alpha_exp = [
        np.array([
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]]),
        np.array([
        [-4/1.2, 4/1.2, 0],
        [4/1.2, -4/1.2, 0],
        [0, 0, 0]]),
        np.array([
        [0, 0, 0],
        [0, -2/1.2, 2/1.2],
        [0, 2/0.6, -2/0.6]]),
        np.array([
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]])]   
    
    assert np.all(np.isclose(k_alpha_beta_matrix, k_alpha_beta_matrix_exp))
    assert np.all(np.isclose(N_alpha_beta_matrix, N_alpha_beta_matrix_exp))
    assert np.all(np.isclose(T_alpha_matrix, T_alpha_matrix_exp))
    assert len(mmvt_Nij_alpha_exp) == len(mmvt_Nij_alpha)
    for mmvt_Nij, mmvt_Nij_exp in zip(mmvt_Nij_alpha, mmvt_Nij_alpha_exp):
        assert np.all(np.isclose(mmvt_Nij, mmvt_Nij_exp))
    assert len(mmvt_Ri_alpha_exp) == len(mmvt_Ri_alpha)
    for mmvt_Ri, mmvt_Ri_exp in zip(mmvt_Ri_alpha, mmvt_Ri_alpha_exp):
        assert np.all(np.isclose(mmvt_Ri, mmvt_Ri_exp))
    assert len(mmvt_Qij_alpha_exp) == len(mmvt_Qij_alpha)
    for mmvt_Qij, mmvt_Qij_exp in zip(mmvt_Qij_alpha, mmvt_Qij_alpha_exp):
        assert np.all(np.isclose(mmvt_Qij, mmvt_Qij_exp))
    return

def test_MMVT_data_sample_fill_k_from_matrices():
    model = make_simple_model()
    k_alpha_beta_exp = {(0,1):20.0, (1,0):10.0,
                    (1,2):10.0,  (2,1):(40.0/3.0),
                    (2,3):(20.0/3.0), (3,2):20.0}
    k_alpha_beta_matrix = np.array([
        [-20.0, 20.0, 0.0, 0.0],
        [10.0, -20.0, 10.0, 0.0],
        [0.0, (40.0/3.0), -(40.0/3.0)-(20.0/3.0), (20.0/3.0)],
        [0.0, 0.0, 20.0, -20.0]])
    main_data_sample = mmvt_analyze.MMVT_data_sample(
            model, None, None, None, 
            None, None)
    main_data_sample.fill_k_from_matrices(k_alpha_beta_matrix)
    compare_dicts(k_alpha_beta_exp, main_data_sample.k_alpha_beta)
    return

def test_MMVT_data_sample_fill_N_R_alpha_from_matrices():
    model = make_simple_model()
    N_i_j_alpha_exp = [{},
                   {(0,1):4, (1,0):4}, 
                   {(1,2):2, (2,1):2},
                   {}]
    R_i_alpha_total_exp = [{0: 1.2},
                       {0: 1.2, 1:1.2},
                       {1: 1.2, 2:0.6},
                       {2: 0.6}]
    mmvt_Nij_alpha = [
        np.array([
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]]),
        np.array([
        [0, 4, 0],
        [4, 0, 0],
        [0, 0, 0]]),
        np.array([
        [0, 0, 0],
        [0, 0, 2],
        [0, 2, 0]]),
        np.array([
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]])]
    mmvt_Ri_alpha = [
        np.array([
        [1.2],
        [0],
        [0]]),
        np.array([
        [1.2],
        [1.2],
        [0]]),
        np.array([
        [0],
        [1.2],
        [0.6]]),
        np.array([
        [0],
        [0],
        [0.6]])]
    N_i_j_alpha_start = [{},
                   {(0,1):3, (1,0):3}, 
                   {(1,2):3, (2,1):3},
                   {}]
    R_i_alpha_total_start = [{0: 1.0},
                       {0: 1.0, 1:1.0},
                       {1: 1.0, 2:0.4},
                       {2: 0.4}]
    T_alpha_total_start = [1.2,
                     2.4,
                     1.8,
                     0.6]
    main_data_sample = mmvt_analyze.MMVT_data_sample(
            model, None, None, N_i_j_alpha_start, 
            R_i_alpha_total_start, T_alpha_total_start)
    main_data_sample.fill_N_R_alpha_from_matrices(mmvt_Nij_alpha, mmvt_Ri_alpha)
    assert len(N_i_j_alpha_exp) == len(main_data_sample.N_i_j_alpha)
    for alpha in range(len(N_i_j_alpha_exp)):
        mmvt_Nij_exp = N_i_j_alpha_exp[alpha]
        mmvt_Nij = main_data_sample.N_i_j_alpha[alpha]
        compare_dicts(mmvt_Nij_exp, mmvt_Nij)
    assert len(R_i_alpha_total_exp) == len(main_data_sample.R_i_alpha)
    for alpha in range(len(R_i_alpha_total_exp)):
        mmvt_Ri_exp = R_i_alpha_total_exp[alpha]
        mmvt_Ri = main_data_sample.R_i_alpha[alpha]
        compare_dicts(mmvt_Ri_exp, mmvt_Ri)
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

def test_monte_carlo_milestoning_error():
    # This test is already implemented in test_markov_chain_monte_carlo.py
    pass