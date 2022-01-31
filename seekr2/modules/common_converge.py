"""
common_converge.py

Contain functions used by converge.py script to determine the
convergence of SEEKR2 calculations.
"""

import os
import glob
import math
from collections import defaultdict
import functools
from itertools import islice

import numpy as np
import matplotlib.pyplot as plt
from parmed import unit

import seekr2.analyze as analyze
import seekr2.modules.common_base as base
import seekr2.modules.common_analyze as common_analyze

# The default number of points to include in convergence plots
DEFAULT_NUM_POINTS = 100

# How many steps to skip before computing the convergence values
DEFAULT_SKIP = 0

# The threshold beneath which to skip plotting the convergence
MIN_PLOT_NORM = 1e-8

# The interval between which to update the user on convergence script progress
PROGRESS_UPDATE_INTERVAL = DEFAULT_NUM_POINTS // 10

def analyze_bd_only(model, data_sample):
    """
    If there are missing MD statistics, then perhaps only a BD analysis
    should be performed. This function only performs a BD analysis on
    a particular data sample.
    """
    bulk_milestones = []
    MFPTs = {}
    for alpha, anchor in enumerate(model.anchors):
        if anchor.bulkstate:
            for milestone_id in anchor.get_ids():
                bulk_milestones.append(milestone_id)
                
    output_file_glob = os.path.join(
        model.anchor_rootdir, model.k_on_info.b_surface_directory, 
        model.k_on_info.bd_output_glob)
    output_file_list = glob.glob(output_file_glob)
    output_file_list = base.order_files_numerically(output_file_list)
    if len(output_file_list) > 0:
        if model.browndye_settings is not None:
            k_ons_src, k_on_errors_src, reaction_probabilities, \
                reaction_probability_errors, transition_counts = \
                common_analyze.browndye_run_compute_rate_constant(os.path.join(
                    model.browndye_settings.browndye_bin_dir,
                    "compute_rate_constant"), output_file_list, 
                    sample_error_from_normal=False)
            data_sample.bd_transition_counts["b_surface"] = transition_counts
        else:
            raise Exception("No valid BD program settings provided.")
        
        if len(bulk_milestones) > 0:
            bulk_milestone = bulk_milestones[0]
            for bd_milestone in model.k_on_info.bd_milestones:
                bd_results_file = os.path.join(
                    model.anchor_rootdir, bd_milestone.directory, 
                    "results.xml")
                
                if not os.path.exists(bd_results_file):
                    bd_directory_list_glob = os.path.join(
                        model.anchor_rootdir, 
                        bd_milestone.directory, 
                        "first_hitting_point_distribution", "lig*/")
                    bd_directory_list = glob.glob(
                        bd_directory_list_glob)
                    if len(bd_directory_list) == 0:
                        continue
                    common_analyze.combine_fhpd_results(
                        bd_milestone, bd_directory_list, bd_results_file)
                
                results_filename_list = [bd_results_file]
                transition_probabilities, transition_counts = \
                    common_analyze.browndye_parse_bd_milestone_results(
                        results_filename_list)
                data_sample.bd_transition_counts[bd_milestone.index] \
                    = transition_counts
                    
    return

def analyze_kinetics(model, analysis, max_step_list, k_on_state=None, 
                     pre_equilibrium_approx=False):
    """
    Extract relevant analysis quantities from sub-sections of the 
    data, which will later be processed for plotting.
    
    Parameters
    -----------
    model : Model()
        milestoning model object containing all milestone and 
        transition information.
        
    analysis : Analysis()
        The object which enables calculation of kinetic and 
        thermodynamic quantities.
        
    max_step_list : list
        A list of integers representing the maximum number of steps
        to analyze per anchor. Used for convergence purposes.
        
    k_on_state : int or None, default None
        If not None, then assume then this is the bound state to 
        assume for k-on calculations. A value of None will skip the
        k-on convergence..
        
    pre_equilibrium_approx : bool
        Whether to use the pre-equilibrium approximation for kinetics
        calculations.
    
    Returns
    -------
    k_on : float
        The k-on value computing using data up to the number of steps 
        in max_step_list for each milestone.
        
    k_off : float
        The k-k_off value computing using data up to the number of 
        steps in max_step_list for each milestone.
        
    N_ij : dict
        The n x n dict matrix representing how many transitions have
        occurred between milestones.
        
    R_i : dict
        An n dict representing the incubation times at each milestone.
    """
    analysis.extract_data(max_step_list=max_step_list, silence_errors=True)
    analysis.fill_out_data_samples()
    try:
        sufficient_statistics = analysis.check_extraction()
        if sufficient_statistics:
            analysis.process_data_samples(pre_equilibrium_approx)
        else:
            analyze_bd_only(model, analysis)
            
        if (k_on_state is not None) and (k_on_state in analysis.k_ons):
            k_on = analysis.k_ons[k_on_state]
        else:
            k_on = 0.0
        k_off = analysis.k_off
        main_data_sample = analysis.main_data_sample
        return k_on, k_off, main_data_sample.N_ij, main_data_sample.R_i
    except (common_analyze.MissingStatisticsError, np.linalg.LinAlgError,
            AssertionError, ValueError) as e:
        if model.get_type() == "mmvt":
            #data_sample = common_analyze.Data_sample(model)
            #data_sample = mmvt_analyze.MMVT_data_sample(model)
            pass
        elif model.get_type() == "elber":
            #data_sample = elber_analyze.Elber_data_sample(model)
            pass
            
        if model.k_on_info is not None:
            #data_sample.parse_browndye_results()
            analysis.main_data_sample.parse_browndye_results()
            
        return 0.0, 0.0, {}, {}

def get_mmvt_max_steps(model, dt):
    """
    Extract the largest simulation step number for all the sims in 
    the anchors.
    """
    max_step_list = np.zeros((model.num_anchors, DEFAULT_NUM_POINTS))
    for alpha, anchor in enumerate(model.anchors):
        max_steps = 0.0
        output_file_glob = os.path.join(
            model.anchor_rootdir, anchor.directory, 
            anchor.production_directory,
            anchor.md_output_glob)
        for output_filename in glob.glob(output_file_glob):
            if not os.path.exists(output_filename):
                continue
            with open(output_filename, "rb") as output_file:
                try:
                    if model.openmm_settings is not None:
                        output_file.seek(-2, os.SEEK_END)
                        while output_file.read(1) != b"\n":
                            output_file.seek(-2, os.SEEK_CUR)
                        last_line = output_file.readline().decode()
                        mytime = float(last_line.strip().split(",")[2])
                        step = int(mytime / dt)
                    elif model.namd_settings is not None:
                        step = 0
                        for line in output_file:
                            line = line.decode("UTF-8")
                            if not line.startswith("SEEKR") \
                                    or len(line.strip()) == 0:
                                continue
                            elif line.startswith("SEEKR: Cell"):
                                line = line.strip().split(" ")
                                step = int(line[8].strip(","))
                                
                            elif line.startswith("SEEKR: Milestone"):
                                line = line.strip().split(" ")
                                step = int(line[10].strip(","))
                                    
                except OSError:
                    step = 0
                if step > max_steps:
                    max_steps = step
                    
        max_steps = DEFAULT_NUM_POINTS * int(math.ceil(
            max_steps / DEFAULT_NUM_POINTS))
        
        if max_steps == 0:
            continue
        
        conv_stride = max_steps // DEFAULT_NUM_POINTS
        conv_intervals = np.arange(conv_stride, max_steps+conv_stride, 
                                   conv_stride)
        conv_intervals = conv_intervals + DEFAULT_SKIP
        max_step_list[alpha,:] = conv_intervals
                    
    return max_step_list
    
def get_elber_max_steps(model, dt):
    """
    Extract the largest simulation step number for all the sims in 
    the anchors.
    """
    max_step_list = np.zeros((model.num_anchors, DEFAULT_NUM_POINTS))
    for alpha, anchor in enumerate(model.anchors):
        max_steps = 0.0
        output_file_glob = os.path.join(
            model.anchor_rootdir, anchor.directory, 
            anchor.production_directory,
            anchor.md_output_glob)
        for output_filename in glob.glob(output_file_glob):
            num_lines = sum(1 for line in open(output_filename))
            max_steps += num_lines
        
        max_steps = DEFAULT_NUM_POINTS * int(math.ceil(
            max_steps / DEFAULT_NUM_POINTS))
        if max_steps == 0:
            #return max_step_list
            continue
        
        conv_stride = max_steps // DEFAULT_NUM_POINTS
        conv_intervals = np.arange(conv_stride, max_steps+conv_stride, 
                                   conv_stride)
        conv_intervals = conv_intervals + DEFAULT_SKIP
        max_step_list[alpha,:] = conv_intervals
        
    return max_step_list
    
def check_milestone_convergence(model, k_on_state=None, 
                                pre_equilibrium_approx=False, verbose=False):
    """
    Calculates the key MMVT quantities N, R, and Q as a function of 
    simulation time to estimate which milestones have been 
    sufficiently sampled. 

    Quantities are pulled from the data at step intervals determined 
    by the conv_stride value with the option to skip steps from the 
    beginning of the data with the skip variable

    Parameters
    -----------
    model : Model()
        milestoning model object containing all milestone and 
        transition information.
        
    k_on_state: int or None, default None
        If not None, then assume then this is the bound state to 
        assume for k-on calculations. A value of None will skip the
        k-on convergence.
        
    pre_equilibrium_approx : bool
        Whether to use the pre-equilibrium approximation for kinetics
        calculations.
        
    verbose : bool, Default False
        Whether to provide more verbose output information.

    Returns
    -------
    k_on_conv : list
        list of calculated on rate at each convergence interval
    
    k_off_conv : list
        list of calculated off rate at each convergence interval
    
    N_ij_conv: list
        list of transition count matrix N for each convergence interval
        
    R_i_conv : list
        list of transition time matrix R for each convergence interval
        
    max_step_list : list
        list of maximum step numbers used for each convergence sample
        
    timestep_in_ns : float
        The length of the timestep in units of nanoseconds
        
    data_sample_list : list
        A list of Data_sample objects that can be used to
        quantitatively monitor convergence.
    """
    
    data_sample_list = []
    if model.openmm_settings is not None:
        dt = model.openmm_settings.langevin_integrator.timestep
    elif model.namd_settings is not None:
        dt = model.namd_settings.langevin_integrator.timestep
    else:
        raise Exception("No OpenMM or NAMD simulation settings in model.")
    timestep_in_ns = unit.Quantity(dt, unit.picosecond).value_in_unit(
        unit.nanoseconds)
    if model.get_type() == "mmvt":
        max_step_list = get_mmvt_max_steps(model, dt)
    elif model.get_type() == "elber":
        max_step_list = get_elber_max_steps(model, dt)
    
    k_off_conv = np.zeros(DEFAULT_NUM_POINTS)
    k_on_conv = np.zeros(DEFAULT_NUM_POINTS)
    N_ij_conv = defaultdict(functools.partial(np.zeros, DEFAULT_NUM_POINTS))
    R_i_conv = defaultdict(functools.partial(np.zeros, DEFAULT_NUM_POINTS))
    analysis = analyze.Analysis(model, force_warning=False)
    for interval_index in range(DEFAULT_NUM_POINTS):
        if verbose and (interval_index % PROGRESS_UPDATE_INTERVAL == 0):
            print("Processing interval {} of {}".format(interval_index, 
                                                        DEFAULT_NUM_POINTS))
        k_on, k_off, N_ij, R_i = analyze_kinetics(
            model, analysis, max_step_list[:, interval_index], k_on_state, 
            pre_equilibrium_approx)
        data_sample_list.append(analysis.main_data_sample)
        k_on_conv[interval_index] = k_on
        k_off_conv[interval_index] = k_off
        for N_ij_key in N_ij:
            N_ij_conv[N_ij_key][interval_index] = N_ij[N_ij_key] / interval_index
        for R_i_key in R_i:
            R_i_conv[R_i_key][interval_index] = R_i[R_i_key] / interval_index
    
    return k_on_conv, k_off_conv, N_ij_conv, R_i_conv, max_step_list, \
        timestep_in_ns, data_sample_list

def plot_scalar_conv(conv_values, conv_intervals, label, title, timestep_in_ns, 
                     y_axis_logarithmic=True):
    """
    Plot convergence of off/on rate or other scalar values as a 
    function of simulation time.

    Parameters
    ----------
    conv_values : list
        list of calculated scalar values for each convergence interval
        
    conv_intervals : list
        list of convergence interval step numbers for which samples 
        are taken.
        
    label : str
        The label to give this plot.
    
    title : str
        The title of this plot.
    
    timestep_in_ns : float
        The length of the timestep in units of nanoseconds.
        
    y_axis_logarithmic : bool, default True
        Whether the y-axis will be plotted on a logarithmic scale.

    Returns
    -------
    fig : matplotlib figure
        matplotlib figure plotting N convergence for each milestone
    ax : object
        matplotlib Axes object
    """
    
    fig, ax = plt.subplots()
    ax.plot(np.multiply(conv_intervals, timestep_in_ns), conv_values, 
            linestyle="-", marker="o", markersize=1)
    plt.ylabel("$"+label+"$")
    plt.xlabel("convergence window")
    plt.title(title)
    if y_axis_logarithmic:
        plt.yscale("log", nonpositive="mask")
    return fig, ax

def plot_dict_conv(conv_dict, conv_intervals, label_base, timestep_in_ns, 
                   skip_null=True, y_axis_logarithmic=True):
    """
    Plot convergence of N_ij or R_i or other dictionary-based value 
    as a function of simulation time.

    Parameters
    ----------
    conv_dict : dict
        dict of lists of calculated off rates for each convergence 
        interval.
        
    conv_interval : list
        list of convergence interval step numbers for which samples 
        are taken.
        
    label_base : str
        The base of the label to give this plot, though the dictionary
        keys will be appended to the label.
    
    timestep_in_ns : float
        The length of the timestep in units of nanoseconds
    
    skip_null : bool, Default True
        If true, than empty convergence lists will be omitted from
        any plots.
    
    y_axis_logarithmic : bool, True
        Whether the y-axis will be plotted on a logarithmic scale.
    
    Returns
    -------
    fig_list : list
        A list of matplotlib figures plotting convergence for each 
        milestone.
        
    ax_list : list
        A list of matplotlib Axes objects.
        
    title_list : list
        A list of the plots' titles.
    
    name_list : list
        A list of the plots' file names.
    """
    
    fig_list = []
    ax_list = []
    title_list = []
    name_list = []
    for key in conv_dict:
        conv_values = conv_dict[key]
        if skip_null:
            if np.linalg.norm(conv_values) < MIN_PLOT_NORM:
                continue
        if isinstance(key, tuple):
            label = "$" + label_base + "_" + "\\rightarrow".join(
                map(str, key)) + "$"
            title = "$" + label_base + "_{" + "\\rightarrow".join(
                map(str, key)) + "}$"
            name = label_base + "_" + "".join(map(str, key))
        elif isinstance(key, int):
            label = "$" + label_base + "_" + str(key) + "$"
            title = "$" + label_base + "_" + str(key) + "$"
            name = label_base + "_" + str(key)
        else:
            raise Exception("key type not implemented: {}".format(type(key)))
        
        fig, ax = plt.subplots()
        ax.plot(np.multiply(conv_intervals, timestep_in_ns), conv_values, 
                linestyle='-', marker="o", markersize = 1)
        plt.ylabel(label)
        plt.xlabel("convergence_window")
        if y_axis_logarithmic:
            plt.yscale("log", nonpositive="mask")
        plt.title(title)
        fig_list.append(fig)
        ax_list.append(ax)
        title_list.append(title)
        name_list.append(name)
    return fig_list, ax_list, title_list, name_list

def calc_window_rmsd(conv_values):
    """
    For a window of convergence values, compute the RMSD and average
    of those values.
    """
    if len(conv_values) == 0:
        return 0.0, 0.0
    average = np.mean(conv_values)
    RMSD_sum = 0.0
    for i, conv_value in enumerate(conv_values):
        RMSD_sum += (conv_value - average)**2
    RMSD = np.sqrt(RMSD_sum / len(conv_values))
    return RMSD, average

def calc_transition_steps(model, data_sample):
    """
    For a given data_sample object, return the number of transitions
    and the minimum transitions between a pair of milestones.
    """
    transition_minima = []
    transition_prob_details = []
    transition_time_details = []
    for alpha, anchor in enumerate(model.anchors):
        transition_detail = {}
        transition_time_detail = {}
        if anchor.bulkstate:
            continue
        
        if model.get_type() == "mmvt":
            if data_sample.T_alpha is None:
                transition_minima.append(0)
                transition_prob_details.append(transition_detail)
                transition_time_details.append(transition_time_detail)
                continue
        
        elif model.get_type() == "elber":
            if data_sample.R_i_list[alpha][alpha] == 0.0:
                transition_minima.append(0)
                transition_prob_details.append(transition_detail)
                transition_time_details.append(transition_time_detail)
                continue
        
        if len(anchor.milestones) == 1:
            # if this is a dead-end milestone
            if model.get_type() == "mmvt":
                transition_dict = data_sample.N_alpha_beta
                k_rate_dict = data_sample.k_alpha_beta
                if len(transition_dict) == 0:
                    transition_quantity = 0
                    
                else:
                    lowest_value = 1e99
                    for key in transition_dict:
                        if key[0] == alpha and transition_dict[key] \
                                < lowest_value:
                            lowest_value  = transition_dict[key]
                    transition_quantity = lowest_value
                    
                    lowest_value = 1e99
                    for key in k_rate_dict:
                        if k_rate_dict[key] == 0.0:
                            lowest_value = 0.0
                            continue
                        if key[0] == alpha and k_rate_dict[key] \
                                < lowest_value:
                            lowest_value  = 1.0 / k_rate_dict[key]
                    transition_time = lowest_value
                    
                    transition_detail = {(alpha,alpha):transition_quantity}
                    transition_time_detail = {(alpha,alpha):transition_time}
                    
            elif model.get_type() == "elber":
                raise Exception("Elber simulations cannot have one milestone.")
            
        else:
            if model.get_type() == "mmvt":
                transition_dict = data_sample.N_i_j_alpha[alpha]
                time_dict = data_sample.R_i_alpha[alpha]
                
            elif model.get_type() == "elber":
                transition_dict = data_sample.N_i_j_list[alpha]
                time_dict = data_sample.R_i_list[alpha]
                
            if len(transition_dict) == 0:
                transition_quantity = 0
            else:
                lowest_value = 1e99
                highest_value = 0
                for key in transition_dict:
                    if transition_dict[key] < lowest_value:
                        lowest_value  = transition_dict[key]
                    if transition_dict[key] > highest_value:
                        highest_value  = transition_dict[key]
                    transition_detail[key] = transition_dict[key]
                transition_quantity = lowest_value
                
                for key in time_dict:
                    transition_time_detail[key] = time_dict[key] / highest_value
                
        transition_minima.append(transition_quantity)
        transition_prob_details.append(transition_detail)
        transition_time_details.append(transition_time_detail)
    return transition_minima, transition_prob_details, transition_time_details

def calc_RMSD_conv_amount(model, data_sample_list, window_size=30,
                               number_of_convergence_windows=20):
    """
    Calculate the RMSD convergence of window spanning a portion of
    a list of data samples.
    """
    convergence_results = []
    for alpha, anchor in enumerate(model.anchors):
        if anchor.bulkstate:
            continue
        conv_list = []
        RMSD_window_conv_list = []
        for window_index in range(number_of_convergence_windows):
            backwards = number_of_convergence_windows - window_index + 1
            bound1 = len(data_sample_list) - window_size - backwards
            bound2 = len(data_sample_list) - backwards
            for data_sample in data_sample_list[bound1:bound2]:
                if model.get_type() == "mmvt":
                    if data_sample.T_alpha is None:
                        RMSD_window_conv_list.append(1e99)
                        print("mark1")
                        break
                elif model.get_type() == "elber":
                    if data_sample.R_i_list[alpha][alpha] == 0.0:
                        RMSD_window_conv_list.append(1e99)
                        break
                
                if len(anchor.milestones) == 1:
                    # if this is a dead-end milestone
                    if model.get_type() == "mmvt":
                        transition_dict = data_sample.N_alpha_beta
                        T_alpha = data_sample.T_alpha[alpha]
                        lowest_value = 1e99
                        for key in transition_dict:
                            if key[0] == alpha and transition_dict[key] \
                                    < lowest_value:
                                lowest_value  = transition_dict[key]
                                conv_quantity \
                                    = transition_dict[key] / T_alpha
                        
                        if lowest_value == 1e99 or lowest_value == 0:
                            RMSD_window_conv_list.append(1e99)
                            print("mark2")
                            break
                    else:
                        raise Exception(
                            "Elber simulations cannot have one milestone.")
                else:
                    if model.get_type() == "mmvt":
                        transition_dict = data_sample.N_i_j_alpha[alpha]
                        T_alpha = data_sample.T_alpha[alpha]
                    elif model.get_type() == "elber":
                        transition_dict = data_sample.N_i_j_list[alpha]
                        T_alpha = data_sample.R_i_list[alpha][alpha]
                        
                    lowest_value = 1e99
                    
                    for key in transition_dict:
                        if transition_dict[key] < lowest_value:
                            lowest_value  = transition_dict[key]
                            conv_quantity = transition_dict[key] / T_alpha
                    
                    if lowest_value == 1e99 or lowest_value == 0:
                        RMSD_window_conv_list.append(1e99)
                        print("mark3")
                        print("alpha:", alpha)
                        print("transition_dict:", transition_dict)
                        break
                        
                conv_list.append(conv_quantity)
                
            RMSD, window_average = calc_window_rmsd(conv_list)
            if window_average == 0.0:
                RMSD_window_conv_list.append(1e99)
                print("mark4")
            else:
                fraction = RMSD / window_average
                RMSD_window_conv_list.append(fraction)
            
        max_RMSD = np.max(RMSD_window_conv_list)
        convergence_results.append(max_RMSD)
    
    return convergence_results