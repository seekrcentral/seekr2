"""
test_common_converge.py

Testing common_converge.py
"""

import os
from shutil import copyfile

import numpy as np

import seekr2.modules.runner_openmm as runner_openmm
import seekr2.modules.common_analyze as common_analyze
import seekr2.modules.common_converge as common_converge
import seekr2.converge as converge
from seekr2.tests.conftest import compare_dicts

TEST_DIRECTORY = os.path.dirname(os.path.realpath(__file__))

test_output_filename = os.path.join(TEST_DIRECTORY, 
                                    "data/test_analyze_outputfile.txt")
test_output_filename_elber = os.path.join(TEST_DIRECTORY, 
                                    "data/test_analyze_output_elber.txt")

def test_analyze_bd_only(host_guest_mmvt_model):
    """
    Test function analyze_bd_only(). This one also tests 
    get_bd_transition_counts().
    """
    test_bd_results_filename = os.path.join(
        TEST_DIRECTORY, "data/sample_bd_results_file2.xml")
    b_surface_directory = os.path.join(
        host_guest_mmvt_model.anchor_rootdir, "b_surface")
    b_surface_results_file = os.path.join(b_surface_directory, "results1.xml")
    data_sample = common_analyze.Data_sample(host_guest_mmvt_model)
    copyfile(test_bd_results_filename, b_surface_results_file)
    common_converge.analyze_bd_only(host_guest_mmvt_model, data_sample)
    assert len(data_sample.bd_transition_counts) == 1
    expected_counts_b_surface = {"total": 100000, "escaped": 54070, 12: 65239,
                                 11: 45930, "stuck": 0}
    counts_b_surface = data_sample.bd_transition_counts["b_surface"]
    compare_dicts(counts_b_surface, expected_counts_b_surface)
    return

def test_mmvt_max_steps(toy_mmvt_model):
    #max_steps = 998000
    max_time = 997991 * 0.002
    for anchor in toy_mmvt_model.anchors:
        runner_openmm.cleanse_anchor_outputs(toy_mmvt_model, anchor)
    anchor1 = toy_mmvt_model.anchors[1]
    anchor1_output_filename = os.path.join(
        toy_mmvt_model.anchor_rootdir, anchor1.directory, 
        anchor1.production_directory, "mmvt.restart1.out")
    copyfile(test_output_filename, anchor1_output_filename)
    anchor_max_times = common_converge.get_mmvt_max_steps(toy_mmvt_model)[-1]
    assert anchor_max_times[1] == max_time
    return

def test_elber_max_steps(toy_elber_model):
    #max_steps = 100
    max_time = 0.036
    for anchor in toy_elber_model.anchors:
        runner_openmm.cleanse_anchor_outputs(toy_elber_model, anchor)
    anchor1 = toy_elber_model.anchors[1]
    anchor1_output_filename = os.path.join(
        toy_elber_model.anchor_rootdir, anchor1.directory, 
        anchor1.production_directory, "forward.restart1.out")
    copyfile(test_output_filename_elber, anchor1_output_filename)
    anchor_max_times =  common_converge.get_elber_max_steps(toy_elber_model)[-1]
    assert np.isclose(anchor_max_times[1], max_time)
    return

def test_calc_transition_steps_mmvt(toy_mmvt_model):
    for anchor in toy_mmvt_model.anchors:
        runner_openmm.cleanse_anchor_outputs(toy_mmvt_model, anchor)
    anchor1 = toy_mmvt_model.anchors[1]
    anchor1_output_filename = os.path.join(
        toy_mmvt_model.anchor_rootdir, anchor1.directory, 
        anchor1.production_directory, "mmvt.restart1.out")
    copyfile(test_output_filename, anchor1_output_filename)
    image_directory = common_analyze.make_image_directory(
        toy_mmvt_model, None)
    data_sample_list, times_dict = converge.converge(
        toy_mmvt_model, 0, image_directory=image_directory,
        verbose=True)
    print("data_sample_list[-1].N_i_j_alpha[1]:", data_sample_list[-1].N_i_j_alpha[1])
    transition_minima, transition_prob_results, transition_time_results \
        = common_converge.calc_transition_steps(
        toy_mmvt_model, data_sample_list[-1])
    for alpha, anchor in enumerate(toy_mmvt_model.anchors):
        if anchor.bulkstate:
            continue
        if anchor.index == 1:
            assert transition_minima[alpha] == 52
            assert transition_prob_results[alpha][(0,1)] == 52
            assert transition_prob_results[alpha][(1,0)] == 52
            assert np.isclose(transition_time_results[alpha][0], (1658.696/52))
            assert np.isclose(transition_time_results[alpha][1], (198.912/52))
        else:
            assert transition_minima[alpha] == 0
            for key in transition_prob_results[alpha]:
                assert transition_prob_results[alpha][key] == 0
            for key in transition_time_results[alpha]:
                assert transition_time_results[alpha][key] == 0

    return

def test_calc_transition_steps_elber(toy_elber_model):
    for anchor in toy_elber_model.anchors:
        runner_openmm.cleanse_anchor_outputs(toy_elber_model, anchor)
    anchor1 = toy_elber_model.anchors[1]
    anchor1_output_filename = os.path.join(
        toy_elber_model.anchor_rootdir, anchor1.directory, 
        anchor1.production_directory, "forward.restart1.out")
    copyfile(test_output_filename_elber, anchor1_output_filename)
    image_directory = common_analyze.make_image_directory(
        toy_elber_model, None)
    data_sample_list, times_dict = converge.converge(
        toy_elber_model, 0, image_directory=image_directory,
        verbose=True)
    transition_minima, transition_prob_results, transition_time_results \
        = common_converge.calc_transition_steps(
        toy_elber_model, data_sample_list[-1])
        
    for alpha, anchor in enumerate(toy_elber_model.anchors):
        if anchor.bulkstate:
            continue
        if anchor.index == 1:
            assert transition_minima[alpha] == 3
            assert transition_prob_results[alpha][(1,2)] == 5
            assert transition_prob_results[alpha][(1,0)] == 3
            assert np.isclose(transition_time_results[alpha][1], (1.324/5))
        else:
            assert transition_minima[alpha] == 0
            for key in transition_prob_results[alpha]:
                assert transition_prob_results[alpha][key] == 0
            for key in transition_time_results[alpha]:
                assert transition_time_results[alpha][key] == 0
                
    return

def test_calc_window_rmsd():
    conv_values = [1.0, 2.0, 3.0]
    rmsd, average = common_converge.calc_window_rmsd(conv_values)
    expected_rmsd = np.sqrt(2.0/3.0)
    expected_average = 2.0
    assert np.isclose(rmsd, expected_rmsd)
    assert np.isclose(average, expected_average)
    
    conv_values = [5.0, 9.0]
    rmsd, average = common_converge.calc_window_rmsd(conv_values)
    expected_rmsd = 2.0
    expected_average = 7.0
    assert np.isclose(rmsd, expected_rmsd)
    assert np.isclose(average, expected_average)
    return

def test_calc_RMSD_conv_amount_mmvt(toy_mmvt_model):
    for anchor in toy_mmvt_model.anchors:
        runner_openmm.cleanse_anchor_outputs(toy_mmvt_model, anchor)
    anchor1 = toy_mmvt_model.anchors[1]
    anchor1_output_filename = os.path.join(
        toy_mmvt_model.anchor_rootdir, anchor1.directory, 
        anchor1.production_directory, "mmvt.restart1.out")
    copyfile(test_output_filename, anchor1_output_filename)
    image_directory = common_analyze.make_image_directory(
        toy_mmvt_model, None)
    data_sample_list, times_dict = converge.converge(
        toy_mmvt_model, 0, image_directory=image_directory,
        verbose=True)
    convergence_results = common_converge.calc_RMSD_conv_amount(
        toy_mmvt_model, data_sample_list)
    for alpha, anchor in enumerate(toy_mmvt_model.anchors):
        if anchor.bulkstate:
            continue
        if anchor.index == 1:
            assert convergence_results[alpha] < 1.0
        else:
            assert convergence_results[alpha] == 1e99
            
    return

def test_calc_RMSD_conv_amount_elber(toy_elber_model):
    for anchor in toy_elber_model.anchors:
        runner_openmm.cleanse_anchor_outputs(toy_elber_model, anchor)
    anchor1 = toy_elber_model.anchors[1]
    anchor1_output_filename = os.path.join(
        toy_elber_model.anchor_rootdir, anchor1.directory, 
        anchor1.production_directory, "forward.restart1.out")
    copyfile(test_output_filename_elber, anchor1_output_filename)
    image_directory = common_analyze.make_image_directory(
        toy_elber_model, None)
    data_sample_list, times_dict = converge.converge(
        toy_elber_model, 0, image_directory=image_directory,
        verbose=True)
    convergence_results = common_converge.calc_RMSD_conv_amount(
        toy_elber_model, data_sample_list)
    
    for alpha, anchor in enumerate(toy_elber_model.anchors):
        if anchor.bulkstate:
            continue
        if anchor.index == 1:
            assert convergence_results[alpha] < 1.0
        else:
            assert convergence_results[alpha] == 1e99
    
    return

def test_array_to_dict():
    a1 = np.array([1.0, 2.0, 3.0])
    d1 = common_converge.array_to_dict(a1)
    assert d1[(0,)] == 1.0
    assert d1[(1,)] == 2.0
    assert d1[(2,)] == 3.0
    
    a2 = np.array([[5.0, 6.0],
                   [7.0, 8.0]])
    d2 = common_converge.array_to_dict(a2)
    assert d2[(0,0)] == 5.0
    assert d2[(0,1)] == 6.0
    assert d2[(1,0)] == 7.0
    assert d2[(1,1)] == 8.0
    return

def test_collapse_list_of_dicts():
    L1 = [{(1,2):1.0}, 
          {(2,3):2.0},
          {(3,4):3.0}]
    d1 = common_converge.collapse_list_of_dicts(L1)
    assert d1[(1,2,0)] ==  1.0
    assert d1[(2,3,1)] ==  2.0
    assert d1[(3,4,2)] ==  3.0
    return

