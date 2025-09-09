"""
test_analyze.py

Test the analyze.py program.
"""

import numpy as np
import pytest

import seekr2.modules.common_analyze as common_analyze
import seekr2.modules.mmvt_analyze as mmvt_analyze
import seekr2.modules.elber_analyze as elber_analyze
import seekr2.run as run
import seekr2.analyze as analyze
import seekr2.tests.test_mmvt_analyze as test_mmvt_analyze
import seekr2.tests.test_elber_analyze as test_elber_analyze
from seekr2.tests.conftest import compare_dicts

def test_check_graph_connected():
    graph1 = {0:[1,2], 1:[0,3], 2:[0], 3:[1]}
    index1 = 0
    index2 = 3
    assert analyze.check_graph_connected(graph1, index1, index2)
    graph1 = {0:[1], 1:[0], 2:[3], 3:[2]}
    index1 = 0
    index2 = 3
    assert not analyze.check_graph_connected(graph1, index1, index2)
    return

@pytest.mark.skip(reason="This test requires too much sampling to pass.")
def test_Analysis_elber_check_anchor_stats(toy_elber_model):
    num_steps = 100000
    fwd_rev_interval = 10
    toy_elber_model.openmm_settings.cuda_platform_settings = None
    toy_elber_model.openmm_settings.reference_platform = True
    toy_elber_model.openmm_settings.langevin_integrator.friction_coefficient \
        = 100.0
    toy_elber_model.calculation_settings.num_umbrella_stage_steps = num_steps
    toy_elber_model.calculation_settings.fwd_rev_interval = fwd_rev_interval
    run.run(toy_elber_model, "0")
    run.run(toy_elber_model, "1")
    run.run(toy_elber_model, "3")
    run.run(toy_elber_model, "4")
    analysis = analyze.Analysis(toy_elber_model, num_error_samples=100)
    analysis.extract_data()
    with pytest.raises(common_analyze.MissingStatisticsError):
        analysis.check_extraction()
    
    run.run(toy_elber_model, "any")
    analysis = analyze.Analysis(toy_elber_model, num_error_samples=100)
    analysis.extract_data()
    analysis.check_extraction()
    return

def test_Analysis_mmvt_check_anchor_stats(toy_mmvt_model):
    num_steps = 1000000
    long_timescale_residence_time_in_ps = 1263.06
    toy_mmvt_model.openmm_settings.cuda_platform_settings = None
    toy_mmvt_model.openmm_settings.reference_platform = True
    toy_mmvt_model.openmm_settings.langevin_integrator.friction_coefficient \
        = 100.0
    toy_mmvt_model.calculation_settings.num_production_steps = num_steps
    toy_mmvt_model.calculation_settings.energy_reporter_interval = num_steps
    run.run(toy_mmvt_model, "0")
    run.run(toy_mmvt_model, "1")
    run.run(toy_mmvt_model, "3")
    run.run(toy_mmvt_model, "4")
    analysis = analyze.Analysis(toy_mmvt_model, num_error_samples=100)
    analysis.extract_data()
    with pytest.raises(common_analyze.MissingStatisticsError):
        analysis.check_extraction()
    
    run.run(toy_mmvt_model, "any")
    analysis = analyze.Analysis(toy_mmvt_model, num_error_samples=100)
    analysis.extract_data()
    analysis.check_extraction()
    return

def test_Analysis_fill_out_data_samples_mmvt():
    model = test_mmvt_analyze.make_simple_model()
    N_alpha_beta = [{1:12}, 
                    {1:12, 2:12}, 
                    {1:12, 2:6},
                    {1:6}]
    k_alpha_beta = [{1:20.0}, 
                    {1:10.0, 2:10.0}, 
                    {1:(40.0/3.0), 2:(20.0/3.0)}, 
                    {1:20.0}]
    N_i_j_alpha_list = [{},
                   {(1,2):4, (2,1):4}, 
                   {(1,2):2, (2,1):2},
                   {}]
    R_i_alpha_total_list = [{1: 1.2},
                       {1: 1.2, 2:1.2},
                       {1: 1.2, 2:0.6},
                       {1: 0.6}]
    T_alpha_total = [1.2,
                     2.4,
                     1.8,
                     0.6]
    analysis = analyze.Analysis(model)
    for alpha in range(model.num_anchors):
        anchor_stats = mmvt_analyze.MMVT_anchor_statistics(alpha)
        anchor_stats.N_i_j_alpha = N_i_j_alpha_list[alpha]
        anchor_stats.R_i_alpha_total = R_i_alpha_total_list[alpha]
        anchor_stats.T_alpha_total = T_alpha_total[alpha]
        anchor_stats.N_alpha_beta = N_alpha_beta[alpha]
        anchor_stats.k_alpha_beta = k_alpha_beta[alpha]
        analysis.anchor_stats_list.append(anchor_stats)
        
    analysis.fill_out_data_samples_mmvt()
    N_alpha_beta_exp = {(0,1):12, (1,0):12,
                    (1,2):12, (2,1):12,
                    (2,3):6,  (3,2):6}
    k_alpha_beta_exp = {(0,1):20.0, (1,0):10.0,
                    (1,2):10.0,  (2,1):(40.0/3.0),
                    (2,3):(20.0/3.0), (3,2):20.0}
    N_i_j_alpha_exp = [{},
                   {(0,1):4, (1,0):4}, 
                   {(1,2):2, (2,1):2},
                   {}]
    R_i_alpha_total_exp = [{0: 1.2},
                       {0: 1.2, 1:1.2},
                       {1: 1.2, 2:0.6},
                       {2: 0.6}]
    T_alpha_total_exp = [1.2,
                     2.4,
                     1.8,
                     0.6]
    compare_dicts(analysis.main_data_sample.N_alpha_beta, N_alpha_beta_exp)
    compare_dicts(analysis.main_data_sample.k_alpha_beta, k_alpha_beta_exp)
    for dict1, dict2 in zip(analysis.main_data_sample.N_i_j_alpha, 
                            N_i_j_alpha_exp):
        compare_dicts(dict1, dict2)
    for dict1, dict2 in zip(analysis.main_data_sample.R_i_alpha, 
                            R_i_alpha_total_exp):
        compare_dicts(dict1, dict2)
    for val1, val2 in zip(analysis.main_data_sample.T_alpha, 
                          T_alpha_total_exp):
        assert np.isclose(val1, val2)
    
    return

def test_Analysis_fill_out_data_samples_elber():
    model = test_elber_analyze.make_simple_model()
    N_i_j_list = [{3: 4}, {1: 4, 3: 2}, {1: 2}]
    R_i_list = [2.4, 2.4, 1.2]
    N_i_j_list_exp = [{(0,1): 4}, {(1,0): 4, (1,2): 2}]
    R_i_list_exp = [2.4, 2.4]
    analysis = analyze.Analysis(model)
    for alpha in range(model.num_anchors):
        anchor_stats = elber_analyze.Elber_anchor_statistics(alpha)
        anchor_stats.N_i_j = N_i_j_list[alpha]
        anchor_stats.R_i_total = R_i_list[alpha]
        analysis.anchor_stats_list.append(anchor_stats)
        
    analysis.fill_out_data_samples_elber()
    for alpha in range(model.num_anchors):
        if model.anchors[alpha].bulkstate:
            continue
        compare_dicts(N_i_j_list_exp[alpha], 
                      analysis.main_data_sample.N_i_j_list[alpha])
        assert analysis.main_data_sample.R_i_list[alpha] == R_i_list_exp[alpha]
    return

def test_entropy_barrier_timescale_mmvt(toy_mmvt_model):
    """
    Test entropy barrier system for if it recreates reasonable
    kinetics timescales.
    """
    num_steps = 1000000
    long_timescale_residence_time_in_ps = 1263.06
    toy_mmvt_model.openmm_settings.cuda_platform_settings = None
    toy_mmvt_model.openmm_settings.reference_platform = True
    toy_mmvt_model.openmm_settings.langevin_integrator.friction_coefficient \
        = 100.0
    toy_mmvt_model.calculation_settings.num_production_steps = num_steps
    toy_mmvt_model.calculation_settings.energy_reporter_interval = num_steps
    run.run(toy_mmvt_model, "any")
    analysis = analyze.analyze(toy_mmvt_model, num_error_samples=100, 
                               skip_checks=False)
    #assert np.isclose(analysis.MFPTs[('anchor_0', 'bulk')], 
    #                  long_timescale_residence_time_in_ps,
    #                  rtol=0.5)
    assert analysis.k_off is not None
    analysis.print_results()
    image_directory = common_analyze.make_image_directory(toy_mmvt_model, None)
    analysis.save_plots(image_directory)
    return

@pytest.mark.skip(reason="This test requires too much sampling to pass.")
def test_entropy_barrier_timescale_elber(toy_elber_model):
    """
    Test entropy barrier system for if it recreates reasonable
    kinetics timescales.
    """
    num_steps = 1000000
    fwd_rev_interval = 100
    long_timescale_residence_time_in_ps = 1263.06
    toy_elber_model.openmm_settings.cuda_platform_settings = None
    toy_elber_model.openmm_settings.reference_platform = True
    toy_elber_model.openmm_settings.langevin_integrator.friction_coefficient \
        = 100.0
    toy_elber_model.calculation_settings.num_umbrella_stage_steps = num_steps
    toy_elber_model.calculation_settings.fwd_rev_interval = fwd_rev_interval
    run.run(toy_elber_model, "any")
    analysis = analyze.analyze(toy_elber_model, num_error_samples=100, 
                               skip_checks=False)
    #assert np.isclose(analysis.MFPTs[('anchor_0', 'bulk')], 
    #                  long_timescale_residence_time_in_ps,
    #                  rtol=0.5)
    assert analysis.MFPTs is not None
    assert analysis.k_off is not None
    analysis.print_results()
    image_directory = common_analyze.make_image_directory(toy_elber_model, None)
    analysis.save_plots(image_directory)
    return