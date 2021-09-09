"""
test_converge.py

Testing converge.py
"""
import os
import glob

import pytest

import seekr2.converge as converge
import seekr2.tests.test_common_converge as test_common_converge
import seekr2.modules.common_analyze as common_analyze
import seekr2.modules.common_converge as common_converge

def test_converge_default(smoluchowski_mmvt_model):
    steps = 1000
    cutoff = 0.1
    minimum_anchor_transitions = 100
    image_directory = common_analyze.make_image_directory(
        smoluchowski_mmvt_model, None)
    k_on_state = 0
    smoluchowski_mmvt_model, standard_time, standard_k_on\
        = test_common_converge.make_model_completed(smoluchowski_mmvt_model, 
                                                    steps)
    data_sample_list = converge.converge(
        smoluchowski_mmvt_model, k_on_state, image_directory=image_directory,
        verbose=True)
    rmsd_convergence_results = common_converge.calc_RMSD_conv_amount(
        smoluchowski_mmvt_model, data_sample_list)
    transition_minima, transition_prob_results, transition_time_results \
        = common_converge.calc_transition_steps(
        smoluchowski_mmvt_model, data_sample_list[-1])

    bd_transition_counts = data_sample_list[-1].bd_transition_counts
        
    converge.print_convergence_results(
        smoluchowski_mmvt_model, rmsd_convergence_results, cutoff, 
        transition_prob_results, transition_time_results,
        minimum_anchor_transitions, bd_transition_counts)
    return

def test_converge_no_bd(smoluchowski_mmvt_model):
    steps = 1000
    cutoff = 0.1
    minimum_anchor_transitions = 100
    k_on_state = 0
    smoluchowski_mmvt_model, standard_time, standard_k_on\
        = test_common_converge.make_model_completed(smoluchowski_mmvt_model, 
                                                    steps)
    smoluchowski_mmvt_model.browndye_settings = None
    smoluchowski_mmvt_model.k_on_info = None
    data_sample_list = converge.converge(
        smoluchowski_mmvt_model, k_on_state, image_directory=None,
        verbose=False)
    rmsd_convergence_results = common_converge.calc_RMSD_conv_amount(
        smoluchowski_mmvt_model, data_sample_list)
    transition_minima, transition_prob_results, transition_time_results \
        = common_converge.calc_transition_steps(
        smoluchowski_mmvt_model, data_sample_list[-1])

    bd_transition_counts = data_sample_list[-1].bd_transition_counts
    converge.print_convergence_results(
        smoluchowski_mmvt_model, rmsd_convergence_results, cutoff, 
        transition_prob_results, transition_time_results,
        minimum_anchor_transitions, bd_transition_counts)
    return

def test_converge_missing_bd(smoluchowski_mmvt_model):
    steps = 1000
    cutoff = 0.1
    minimum_anchor_transitions = 100
    k_on_state = 0
    smoluchowski_mmvt_model, standard_time, standard_k_on\
        = test_common_converge.make_model_completed(smoluchowski_mmvt_model, 
                                                    steps)
    b_surface_dir_path = os.path.join(
        smoluchowski_mmvt_model.anchor_rootdir, 
        smoluchowski_mmvt_model.k_on_info.b_surface_directory)
    results_glob = os.path.join(b_surface_dir_path, "results*.xml")
    for results_file in glob.glob(results_glob):
        print("removing file:", results_file)
        os.remove(results_file)
    data_sample_list = converge.converge(
        smoluchowski_mmvt_model, k_on_state, image_directory=None,
        verbose=False)
    rmsd_convergence_results = common_converge.calc_RMSD_conv_amount(
        smoluchowski_mmvt_model, data_sample_list)
    transition_minima, transition_prob_results, transition_time_results \
        = common_converge.calc_transition_steps(
        smoluchowski_mmvt_model, data_sample_list[-1])

    bd_transition_counts = data_sample_list[-1].bd_transition_counts
    converge.print_convergence_results(
        smoluchowski_mmvt_model, rmsd_convergence_results, cutoff, 
        transition_prob_results, transition_time_results,
        minimum_anchor_transitions, bd_transition_counts)
    return

def test_converge_missing_anchor_stats(smoluchowski_mmvt_model):
    steps = 1000
    cutoff = 0.1
    minimum_anchor_transitions = 100
    k_on_state = 0
    smoluchowski_mmvt_model, standard_time, standard_k_on\
        = test_common_converge.make_model_completed(smoluchowski_mmvt_model, 
                                                    steps)
    anchor_dir_path = os.path.join(
        smoluchowski_mmvt_model.anchor_rootdir, 
        smoluchowski_mmvt_model.anchors[1].directory,
        smoluchowski_mmvt_model.anchors[1].production_directory)
    results_glob = os.path.join(anchor_dir_path, "mmvt*.out")
    for results_file in glob.glob(results_glob):
        os.remove(results_file)
    data_sample_list = converge.converge(
        smoluchowski_mmvt_model, k_on_state, image_directory=None,
        verbose=False)
    rmsd_convergence_results = common_converge.calc_RMSD_conv_amount(
        smoluchowski_mmvt_model, data_sample_list)
    transition_minima, transition_prob_results, transition_time_results \
        = common_converge.calc_transition_steps(
        smoluchowski_mmvt_model, data_sample_list[-1])

    bd_transition_counts = data_sample_list[-1].bd_transition_counts
    converge.print_convergence_results(
        smoluchowski_mmvt_model, rmsd_convergence_results, cutoff, 
        transition_prob_results, transition_time_results,
        minimum_anchor_transitions, bd_transition_counts)
    return

def test_converge_sparse_anchor_stats(smoluchowski_mmvt_model):
    steps = 1000
    cutoff = 0.1
    minimum_anchor_transitions = 100
    k_on_state = 0
    smoluchowski_mmvt_model, standard_time, standard_k_on\
        = test_common_converge.make_model_completed(smoluchowski_mmvt_model, 
                                                    steps)
    anchor_dir_path = os.path.join(
        smoluchowski_mmvt_model.anchor_rootdir, 
        smoluchowski_mmvt_model.anchors[1].directory,
        smoluchowski_mmvt_model.anchors[1].production_directory)
    results_glob = os.path.join(anchor_dir_path, "mmvt*.out")
    for results_file in glob.glob(results_glob):
        with open(results_file, "w") as f:
            f.write('#"Bounced boundary ID","bounce index","total time (ps)"\n')
            f.write('1,0,0.442\n')
        
    data_sample_list = converge.converge(
        smoluchowski_mmvt_model, k_on_state, image_directory=None,
        verbose=False)
    rmsd_convergence_results = common_converge.calc_RMSD_conv_amount(
        smoluchowski_mmvt_model, data_sample_list)
    transition_minima, transition_prob_results, transition_time_results \
        = common_converge.calc_transition_steps(
        smoluchowski_mmvt_model, data_sample_list[-1])

    bd_transition_counts = data_sample_list[-1].bd_transition_counts
    converge.print_convergence_results(
        smoluchowski_mmvt_model, rmsd_convergence_results, cutoff, 
        transition_prob_results, transition_time_results,
        minimum_anchor_transitions, bd_transition_counts)
    return