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
import seekr2.modules.common_analyze as common_analyze

# The default number of points to include in convergence plots
DEFAULT_NUM_POINTS = 100

# How many steps to skip before computing the convergence values
DEFAULT_SKIP = 0

# The threshold beneath which to skip plotting the convergence
MIN_PLOT_NORM = 1e-8

# The interval between which to update the user on convergence script progress
PROGRESS_UPDATE_INTERVAL = DEFAULT_NUM_POINTS // 10

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
    
    try:
        analysis.extract_data(max_step_list, silence_errors=True)
        analysis.fill_out_data_samples()
        sufficient_statistics = analysis.check_extraction()
        if sufficient_statistics:
            analysis.process_data_samples(pre_equilibrium_approx)
        if (k_on_state is not None) and (k_on_state in analysis.k_ons):
            k_on = analysis.k_ons[k_on_state]
        else:
            k_on = 0.0
        k_off = analysis.k_off
        main_data_sample = analysis.main_data_sample
        return k_on, k_off, main_data_sample.N_ij, main_data_sample.R_i
    except common_analyze.MissingStatisticsError:
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
                        for counter, line in enumerate(output_file):
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
            N_ij_conv[N_ij_key][interval_index] = N_ij[N_ij_key]
        for R_i_key in R_i:
            R_i_conv[R_i_key][interval_index] = R_i[R_i_key]
    
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
            linestyle='-', marker="o", markersize = 1)
    plt.ylabel("$"+label+"$")
    plt.xlabel("convergence window")
    plt.title("$"+title+"$")
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

def make_windows(seq, n):
    """
    Returns a sliding window (of width n) over data from the iterable
       s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   
    """
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result
        
def calc_window_rmsd(conv_values):
    """
    For a window of convergence values, compute the RMSD and average
    of those values.
    """
    average = np.mean(conv_values)
    RMSD_sum = 0.0
    for i, conv_value in enumerate(conv_values):
        RMSD_sum += (conv_value - average)**2
    RMSD = np.sqrt(RMSD_sum) / (len(conv_values)-1)
    return RMSD, average

def find_conv_min(conv_values, conv_intervals, window_size, cutoff, 
                   conv_windows):
    """
    Find the time at which the convergence values converged, if ever.
    """
    conv_times = defaultdict(float)
    for key in conv_values:
        if np.sum(conv_values[key][:]) != 0:
            conv_times[key] = np.nan
            rmsd_list = []
            windows = make_windows(conv_values[key][:], window_size)
            index = 0
            conv_counter = 0
            for w in windows:
                RMSD, window_average = calc_window_rmsd(w)
                if RMSD <= (cutoff * window_average):
                    conv_counter += 1
                    if conv_counter == conv_windows:
                        max_int = index + window_size
                        min_time = conv_intervals[max_int]
                        conv_times[key] = min_time
                        break
                else: conv_counter = 0
                index += 1
            
            """    
            if np.isnan(conv_times[key]): print(
                "Entry %s did not meet minimum convergence criteria of "\
                "%i windows" %(str(key), conv_windows))
            """
    return conv_times

def calc_transition_steps(model, data_sample):
    """
    For a given data_sample object, return the number of transitions
    and the minimum transitions between a pair of milestones.
    """
    transition_minima = []
    transition_details = []
    for alpha, anchor in enumerate(model.anchors):
        transition_detail = {}
        if anchor.bulkstate:
            continue
        
        if data_sample is None:
            transition_minima.append(0)
            continue
        
        if len(anchor.milestones) == 1:
            # if this is a dead-end milestone
            if model.get_type() == "mmvt":
                transition_dict = data_sample.N_alpha_beta
                if len(transition_dict) == 0:
                    transition_quantity = 0
                else:
                    lowest_value = 1e99
                    for key in transition_dict:
                        if key[0] == alpha and transition_dict[key] \
                                < lowest_value:
                            lowest_value  = transition_dict[key]
                    transition_quantity = lowest_value
                    transition_detail = {(alpha,alpha):transition_quantity}
                    
            elif model.get_type() == "elber":
                raise Exception("Elber simulations cannot have one milestone.")
            
        else:
            if model.get_type() == "mmvt":
                transition_dict = data_sample.N_i_j_alpha[alpha]
                
            elif model.get_type() == "elber":
                transition_dict = data_sample.N_i_j_list[alpha]
                
            if len(transition_dict) == 0:
                transition_quantity = 0
            else:
                lowest_value = 1e99
                for key in transition_dict:
                    if transition_dict[key]  < lowest_value:
                        lowest_value  = transition_dict[key]
                    transition_detail[key] = transition_dict[key]
                transition_quantity = lowest_value
        transition_minima.append(transition_quantity)
        transition_details.append(transition_detail)
    return transition_minima, transition_details

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
                if data_sample is None:
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
                        if transition_dict[key]  < lowest_value:
                            lowest_value  = transition_dict[key]
                            conv_quantity = transition_dict[key] / T_alpha
                    
                    if lowest_value == 1e99 or lowest_value == 0:
                        RMSD_window_conv_list.append(1e99)
                        break
                        
                conv_list.append(conv_quantity)
                
            RMSD, window_average = calc_window_rmsd(conv_list)
            fraction = RMSD / window_average
            if np.isnan(fraction):
                RMSD_window_conv_list.append(1e99)
            else:
                RMSD_window_conv_list.append(fraction)
        
        max_RMSD = np.max(RMSD_window_conv_list)
        convergence_results.append(max_RMSD)
    
    return convergence_results

def calc_RMSD_conv_time(model, N_conv, R_conv, conv_intervals, window_size, 
                   cutoff, conv_windows, timestep_in_ns):
    """ Estimate the convergence of sampling for each milestone using 
    a sliding window RMSD cutoff. Milestones are considered converged 
    when the values of BOTH N and R remain below the cutoff threshold 
    for a specified number of windows.

    The cutoff is a percentage of the magnutude of the corresponding 
    value (e.g. less than 5% of the magnutude of N).

    Parameters
    -----------
    model : object
        milestoning model object containing all milestone information 
        and transition statistics
    N_conv: list
        list of transition count matrix N for each convergence interval
    R_conv : list
        list of transition time matrix R for each convergence interval 
    conv_intervals : list
        list of stepnumbers for which convergence samples were 
        extracted and calculated
    window size : int
        number of convergence samples used in a window RMSD 
        calculation
    cutoff : float, Default 0.01
        RMSD must be less than the cutoff times the magnitude of the 
        quantity (N or R) to be considered converged
    conv_windows : int, Default 20
        number of consecutive windows the RMSD must remain below the 
        cutoff to be considered converged
    timestep_in_ns : float
        The length of the timestep in units of nanoseconds
    
    Return
    -------
    min_anchor_times : list
        list of times (in stepnumbers) for each voronoi cell where 
        the convergence criteria are satisfied
    """

    min_anchor_times = np.zeros((len(conv_intervals)))
    n_conv_times = find_conv_min(N_conv, conv_intervals, window_size, cutoff,
                                  conv_windows)
    r_conv_times = find_conv_min(R_conv, conv_intervals, window_size, cutoff,
                                  conv_windows)
    for anchor in model.anchors:
        n_steps = 0
        r_steps = 0
        milestones = []
        for milestone in anchor.milestones:
            milestones.append(milestone.index)
            
        for key in n_conv_times:
            if key[0] in milestones and key[1] in milestones:
                if np.isnan(n_conv_times[key]):
                    n_steps = None
                    break
                if n_conv_times[key] > n_steps:
                    n_steps = n_conv_times[key]
        
        for key in r_conv_times:
            if key in milestones:
                if np.isnan(r_conv_times[key]):
                    r_steps = None
                    break
                if r_conv_times[key] > r_steps:
                    r_steps = r_conv_times[key]
        
        if r_steps is None or n_steps is None:
            print("Anchor %i is not converged." % anchor.index)
        else:
            print("Anchor %i converged after time %f ns" \
                  % (anchor.index, timestep_in_ns * max(n_steps, r_steps)))
    return