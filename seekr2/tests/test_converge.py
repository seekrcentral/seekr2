"""
test_converge.py

Testing converge.py
"""
import os
import glob

import pytest

import seekr2.run as run
import seekr2.converge as converge
#import seekr2.tests.test_common_converge as test_common_converge
import seekr2.modules.common_analyze as common_analyze
import seekr2.modules.common_converge as common_converge
import seekr2.modules.runner_browndye2 as runner_browndye2
import seekr2.modules.runner_openmm as runner_openmm
# remove a warning about too many open figures.
common_converge.plt.rcParams.update({'figure.max_open_warning': 0})

def test_converge_default_md_only_mmvt(toy_mmvt_model):
    num_steps = 100000
    toy_mmvt_model.openmm_settings.cuda_platform_settings = None
    toy_mmvt_model.openmm_settings.reference_platform = True
    toy_mmvt_model.openmm_settings.langevin_integrator.friction_coefficient \
        = 100.0
    toy_mmvt_model.calculation_settings.num_production_steps = num_steps
    toy_mmvt_model.calculation_settings.energy_reporter_interval = num_steps
    run.run(toy_mmvt_model, "any", force_overwrite=True)
    
    cutoff = 0.1
    minimum_anchor_transitions = 100
    image_directory = common_analyze.make_image_directory(
        toy_mmvt_model, None)
    k_on_state = None
    data_sample_list, times_dict = converge.converge(
        toy_mmvt_model, k_on_state, image_directory=image_directory,
        verbose=True)
    rmsd_convergence_results = common_converge.calc_RMSD_conv_amount(
        toy_mmvt_model, data_sample_list)
    transition_minima, transition_prob_results, transition_time_results \
        = common_converge.calc_transition_steps(
            toy_mmvt_model, data_sample_list[-1])
    bd_transition_counts = data_sample_list[-1].bd_transition_counts
    converge.print_convergence_results(
        toy_mmvt_model, rmsd_convergence_results, cutoff, 
        transition_prob_results, transition_time_results,
        minimum_anchor_transitions, bd_transition_counts, times_dict)
    return

def test_converge_default_bd_only_mmvt(host_guest_mmvt_model):
    model = host_guest_mmvt_model
    bd_directory = os.path.join(host_guest_mmvt_model.anchor_rootdir, 
                                "b_surface")
    runner_browndye2.run_bd_top(model.browndye_settings.browndye_bin_dir, 
               bd_directory, force_overwrite=True)
    runner_browndye2.run_nam_simulation(
        model.browndye_settings.browndye_bin_dir, bd_directory, 
        model.k_on_info.bd_output_glob)
    cutoff = 0.1
    minimum_anchor_transitions = 100
    image_directory = common_analyze.make_image_directory(
        host_guest_mmvt_model, None)
    k_on_state = 0
    data_sample_list, times_dict = converge.converge(
        host_guest_mmvt_model, k_on_state, image_directory=image_directory,
        verbose=True)
    rmsd_convergence_results = common_converge.calc_RMSD_conv_amount(
        host_guest_mmvt_model, data_sample_list)
    transition_minima, transition_prob_results, transition_time_results \
        = common_converge.calc_transition_steps(
            host_guest_mmvt_model, data_sample_list[-1])
    bd_transition_counts = data_sample_list[-1].bd_transition_counts
    converge.print_convergence_results(
        host_guest_mmvt_model, rmsd_convergence_results, cutoff, 
        transition_prob_results, transition_time_results,
        minimum_anchor_transitions, bd_transition_counts, times_dict)
    return

def test_converge_incomplete_md_only_mmvt(toy_mmvt_model):
    num_steps = 100000
    toy_mmvt_model.openmm_settings.cuda_platform_settings = None
    toy_mmvt_model.openmm_settings.reference_platform = True
    toy_mmvt_model.openmm_settings.langevin_integrator.friction_coefficient \
        = 100.0
    toy_mmvt_model.calculation_settings.num_production_steps = num_steps
    toy_mmvt_model.calculation_settings.energy_reporter_interval = num_steps
    run.run(toy_mmvt_model, "any", force_overwrite=True)
    runner_openmm.cleanse_anchor_outputs(
        toy_mmvt_model, toy_mmvt_model.anchors[3])
    
    cutoff = 0.1
    minimum_anchor_transitions = 100
    image_directory = common_analyze.make_image_directory(
        toy_mmvt_model, None)
    k_on_state = None
    data_sample_list, times_dict = converge.converge(
        toy_mmvt_model, k_on_state, image_directory=image_directory,
        verbose=True)
    rmsd_convergence_results = common_converge.calc_RMSD_conv_amount(
        toy_mmvt_model, data_sample_list)
    transition_minima, transition_prob_results, transition_time_results \
        = common_converge.calc_transition_steps(
            toy_mmvt_model, data_sample_list[-1])
    bd_transition_counts = data_sample_list[-1].bd_transition_counts
    converge.print_convergence_results(
        toy_mmvt_model, rmsd_convergence_results, cutoff, 
        transition_prob_results, transition_time_results,
        minimum_anchor_transitions, bd_transition_counts, times_dict)
    return

def test_converge_default_md_only_elber(toy_elber_model):
    num_steps = 100000
    fwd_rev_interval = 100
    toy_elber_model.openmm_settings.cuda_platform_settings = None
    toy_elber_model.openmm_settings.reference_platform = True
    toy_elber_model.openmm_settings.langevin_integrator.friction_coefficient \
        = 100.0
    toy_elber_model.calculation_settings.num_umbrella_stage_steps = num_steps
    toy_elber_model.calculation_settings.fwd_rev_interval = fwd_rev_interval
    run.run(toy_elber_model, "any")
    
    cutoff = 0.1
    minimum_anchor_transitions = 100
    image_directory = common_analyze.make_image_directory(
        toy_elber_model, None)
    k_on_state = None
    data_sample_list, times_dict = converge.converge(
        toy_elber_model, k_on_state, image_directory=image_directory,
        verbose=True)
    rmsd_convergence_results = common_converge.calc_RMSD_conv_amount(
        toy_elber_model, data_sample_list)
    transition_minima, transition_prob_results, transition_time_results \
        = common_converge.calc_transition_steps(
            toy_elber_model, data_sample_list[-1])
    bd_transition_counts = data_sample_list[-1].bd_transition_counts
    converge.print_convergence_results(
        toy_elber_model, rmsd_convergence_results, cutoff, 
        transition_prob_results, transition_time_results,
        minimum_anchor_transitions, bd_transition_counts, times_dict)
    return

def test_converge_incomplete_md_only_elber(toy_elber_model):
    num_steps = 100000
    fwd_rev_interval = 100
    toy_elber_model.openmm_settings.cuda_platform_settings = None
    toy_elber_model.openmm_settings.reference_platform = True
    toy_elber_model.openmm_settings.langevin_integrator.friction_coefficient \
        = 100.0
    toy_elber_model.calculation_settings.num_umbrella_stage_steps = num_steps
    toy_elber_model.calculation_settings.fwd_rev_interval = fwd_rev_interval
    run.run(toy_elber_model, "any")
    runner_openmm.cleanse_anchor_outputs(
        toy_elber_model, toy_elber_model.anchors[3])
    
    cutoff = 0.1
    minimum_anchor_transitions = 100
    image_directory = common_analyze.make_image_directory(
        toy_elber_model, None)
    k_on_state = None
    data_sample_list, times_dict = converge.converge(
        toy_elber_model, k_on_state, image_directory=image_directory,
        verbose=True)
    rmsd_convergence_results = common_converge.calc_RMSD_conv_amount(
        toy_elber_model, data_sample_list)
    transition_minima, transition_prob_results, transition_time_results \
        = common_converge.calc_transition_steps(
            toy_elber_model, data_sample_list[-1])
    bd_transition_counts = data_sample_list[-1].bd_transition_counts
    converge.print_convergence_results(
        toy_elber_model, rmsd_convergence_results, cutoff, 
        transition_prob_results, transition_time_results,
        minimum_anchor_transitions, bd_transition_counts, times_dict)
    return

def test_converge_default_bd_only_elber(host_guest_elber_model):
    model = host_guest_elber_model
    bd_directory = os.path.join(host_guest_elber_model.anchor_rootdir, 
                                "b_surface")
    runner_browndye2.run_bd_top(model.browndye_settings.browndye_bin_dir, 
               bd_directory, force_overwrite=True)
    runner_browndye2.run_nam_simulation(
        model.browndye_settings.browndye_bin_dir, bd_directory, 
        model.k_on_info.bd_output_glob)
    cutoff = 0.1
    minimum_anchor_transitions = 100
    image_directory = common_analyze.make_image_directory(
        host_guest_elber_model, None)
    k_on_state = 0
    data_sample_list, times_dict = converge.converge(
        host_guest_elber_model, k_on_state, image_directory=image_directory,
        verbose=True)
    rmsd_convergence_results = common_converge.calc_RMSD_conv_amount(
        host_guest_elber_model, data_sample_list)
    transition_minima, transition_prob_results, transition_time_results \
        = common_converge.calc_transition_steps(
            host_guest_elber_model, data_sample_list[-1])
    bd_transition_counts = data_sample_list[-1].bd_transition_counts
    converge.print_convergence_results(
        host_guest_elber_model, rmsd_convergence_results, cutoff, 
        transition_prob_results, transition_time_results,
        minimum_anchor_transitions, bd_transition_counts, times_dict)
    return