"""
test_common_converge.py

Testing common_converge.py
"""

import os

import numpy as np
import pytest

import seekr2.modules.common_base as base
import seekr2.modules.common_converge as common_converge
import seekr2.converge as converge
import seekr2.toy.smoluchowski as smol
import seekr2.tests.test_analyze as test_analyze

def make_model_completed(smoluchowski_model, steps):
    """
    Using the Smoluchowski engine, create a model that has all its
    statistics filled out.
    """
    smoluchowski_model.openmm_settings = base.Openmm_settings()
    smoluchowski_model.openmm_settings.langevin_integrator \
        = base.Langevin_integrator_settings()
    D = 1.0
    b_surface_dir_path = os.path.join(
        smoluchowski_model.anchor_rootdir, 
        smoluchowski_model.k_on_info.b_surface_directory)
    smoluchowski_model.k_on_info.bd_milestones[0].directory \
        = "bd_milestone_{}".format(
            smoluchowski_model.k_on_info.bd_milestones[0].index)
    if not os.path.exists(b_surface_dir_path):
        os.mkdir(b_surface_dir_path)
    k_on_src, transition_probs = test_analyze.make_bd_statistics_flat(
        smoluchowski_model, D)
    potential_energy_function = smol.FlatPotentialEnergyFunction()
    if smoluchowski_model.get_type() == "mmvt":
        mmvt_time, k_on = test_analyze.make_mmvt_calculation_based_on_output_files(
            smoluchowski_model, potential_energy_function, k_on_src, transition_probs, steps=steps)
    elif smoluchowski_model.get_type() == "elber":
        mmvt_time, k_on = test_analyze.make_elber_calculation_based_on_output_files(
            smoluchowski_model, potential_energy_function, k_on_src, transition_probs, steps=steps)
    bound_state_radius = smoluchowski_model.anchors[0].milestones[0].variables["radius"]
    standard_k_on = 4.0 * np.pi * bound_state_radius * D
    standard_time = test_analyze.make_smoluchowski_standard_for_k_off(
        smoluchowski_model, potential_energy_function)
    return smoluchowski_model, standard_time, standard_k_on

def test_analyze_bd_only(smoluchowski_mmvt_model):
    steps = 10000
    cutoff = 0.1
    minimum_anchor_transitions = 100
    image_directory = None
    k_on_state = 0
    smoluchowski_mmvt_model, standard_time, standard_k_on\
        = make_model_completed(smoluchowski_mmvt_model, steps)
    data_sample_list = converge.converge(
        smoluchowski_mmvt_model, k_on_state, image_directory=image_directory,
        verbose=True)
    common_converge.analyze_bd_only(smoluchowski_mmvt_model, 
                                    data_sample_list[-1])
    
    return

def test_mmvt_max_steps(smoluchowski_mmvt_model):
    steps = 10000
    cutoff = 0.1
    minimum_anchor_transitions = 100
    image_directory = None
    k_on_state = 0
    dt = 0.002
    smoluchowski_mmvt_model, standard_time, standard_k_on\
        = make_model_completed(smoluchowski_mmvt_model, steps)
    max_step_list = common_converge.get_mmvt_max_steps(
        smoluchowski_mmvt_model, dt)
    return

def test_elber_max_steps(smoluchowski_elber_model):
    steps = 10000
    cutoff = 0.1
    minimum_anchor_transitions = 100
    image_directory = None
    k_on_state = 0
    dt = 0.002
    smoluchowski_mmvt_model, standard_time, standard_k_on\
        = make_model_completed(smoluchowski_elber_model, steps)
    max_step_list = common_converge.get_mmvt_max_steps(
        smoluchowski_mmvt_model, dt)
    return

def test_check_milestone_convergence(smoluchowski_mmvt_model):
    steps = 10000
    k_on_state = 0
    smoluchowski_mmvt_model, standard_time, standard_k_on\
        = make_model_completed(smoluchowski_mmvt_model, steps)
    k_on_conv, k_off_conv, N_ij_conv, R_i_conv, max_step_list, \
        timestep_in_ns, data_sample_list \
        = common_converge.check_milestone_convergence(
            smoluchowski_mmvt_model, k_on_state=k_on_state, 
            pre_equilibrium_approx=False, verbose=True)
    return
    
def test_check_milestone_convergence_pre_equilibrium_approx(
        smoluchowski_mmvt_model):
    steps = 10000
    k_on_state = 0
    smoluchowski_mmvt_model, standard_time, standard_k_on\
        = make_model_completed(smoluchowski_mmvt_model, steps)
    k_on_conv, k_off_conv, N_ij_conv, R_i_conv, max_step_list, \
        timestep_in_ns, data_sample_list \
        = common_converge.check_milestone_convergence(
            smoluchowski_mmvt_model, k_on_state=k_on_state, 
            pre_equilibrium_approx=True, verbose=True)
    return

def test_plot_scalar_conv(smoluchowski_mmvt_model):
    steps = 10000
    k_on_state = 0
    smoluchowski_mmvt_model, standard_time, standard_k_on\
        = make_model_completed(smoluchowski_mmvt_model, steps)
    k_on_conv, k_off_conv, N_ij_conv, R_i_conv, max_step_list, \
        timestep_in_ns, data_sample_list \
        = common_converge.check_milestone_convergence(
            smoluchowski_mmvt_model, k_on_state=k_on_state, 
            pre_equilibrium_approx=False, verbose=True)
    
    k_off_fig, ax = common_converge.plot_scalar_conv(
        k_off_conv, max_step_list[0,:], title="$k_{off}$ Convergence", 
        label="k_{off} (s^{-1})", timestep_in_ns=timestep_in_ns)
    return

def test_plot_dict_conv(smoluchowski_mmvt_model):
    steps = 10000
    k_on_state = 0
    smoluchowski_mmvt_model, standard_time, standard_k_on\
        = make_model_completed(smoluchowski_mmvt_model, steps)
    k_on_conv, k_off_conv, N_ij_conv, R_i_conv, max_step_list, \
        timestep_in_ns, data_sample_list \
        = common_converge.check_milestone_convergence(
            smoluchowski_mmvt_model, k_on_state=k_on_state, 
            pre_equilibrium_approx=False, verbose=True)
    
    N_ij_fig_list, ax, N_ij_title_list, N_ij_name_list \
        = common_converge.plot_dict_conv(
        N_ij_conv, max_step_list[0,:], label_base="N", 
        timestep_in_ns=timestep_in_ns)
        
    R_i_fig_list, ax, R_i_title_list, R_i_name_list \
        = common_converge.plot_dict_conv(
        R_i_conv, max_step_list[0,:], label_base="R", 
        timestep_in_ns=timestep_in_ns)
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
    
def test_calc_transition_steps(smoluchowski_mmvt_model):
    steps = 10000
    cutoff = 0.1
    minimum_anchor_transitions = 100
    image_directory = None
    k_on_state = 0
    smoluchowski_mmvt_model, standard_time, standard_k_on\
        = make_model_completed(smoluchowski_mmvt_model, steps)
    data_sample_list = converge.converge(
        smoluchowski_mmvt_model, k_on_state, image_directory=image_directory,
        verbose=True)
    transition_minima, transition_prob_results, transition_time_results \
        = common_converge.calc_transition_steps(
        smoluchowski_mmvt_model, data_sample_list[-1])
    return

def test_calc_RMSD_conv_amount(smoluchowski_mmvt_model):
    steps = 10000
    cutoff = 0.1
    minimum_anchor_transitions = 100
    image_directory = None
    k_on_state = 0
    smoluchowski_mmvt_model, standard_time, standard_k_on\
        = make_model_completed(smoluchowski_mmvt_model, steps)
    data_sample_list = converge.converge(
        smoluchowski_mmvt_model, k_on_state, image_directory=image_directory,
        verbose=True)
    convergence_results = common_converge.calc_RMSD_conv_amount(
        smoluchowski_mmvt_model, data_sample_list, window_size=30,
        number_of_convergence_windows=20)
    return