"""
common_analyze.py

Classes, functions and other objects used for analyzing SEEKR2 
simulation outputs common to MMVT and Elber calculation types.
"""

import os
import glob
import xml.etree.ElementTree as ET
import subprocess
import random
from collections import defaultdict

import numpy as np
from scipy import linalg as la
from scipy.stats import gamma

import seekr2.modules.common_base as base

# The power to raise a matrix to when computing stationary probabilities
MATRIX_EXPONENTIAL = 9999999

GAS_CONSTANT = 0.0019872 # in kcal/mol*K
DEFAULT_IMAGE_DIR = "images_and_plots/"

class MissingStatisticsError(Exception):
    """Catch a very specific type of error in analysis stage."""
    pass

def solve_rate_matrix(Q_hat, max_iter=100, ETOL=1e-10):
    """
    Use an infinite series combined with a divide-and-conquer method
    to solve the Q-matrix since the numpy/scipy solver tends to have 
    numerical instabilities.
    """
    
    K_hat = Q_to_K(Q_hat)
    sum_of_series = np.identity(K_hat.shape[0], dtype=np.longdouble)
    for i in range(max_iter):
        n = 2**i
        K_pow_n = np.linalg.matrix_power(K_hat, n)
        section = K_pow_n @ sum_of_series
        sum_of_series += section
        error = np.linalg.norm(section)
        if error < ETOL:
            break
    T = np.zeros(Q_hat.shape, dtype=np.longdouble)
    t_vector = np.zeros((Q_hat.shape[0], 1), dtype=np.longdouble)
    for i in range(Q_hat.shape[0]):
        T[i,i] = -1.0 / Q_hat[i,i]
        t_vector[i,0] = T[i,i]
    one_vector = np.ones((Q_hat.shape[0]), dtype=np.longdouble)
    times = sum_of_series @ T @ one_vector
    return times

def Q_to_K(Q):
    """Given a rate matrix Q, compute probability transition matrix K."""
    
    K = np.zeros(Q.shape, dtype=np.longdouble)
    for i in range(Q.shape[0]):
        for j in range(Q.shape[0]):
            if i == j:
                K[i,j] = 0.0
            else:
                K[i,j] = -Q[i,j] / Q[i,i]
    return K

#def make_elber_K_matrix(oldK):
#    """Make a K matrix that is compatible with Elber milestoning."""
#    old_shape = oldK.shape
#    print("oldK:", oldK)
#    exit()
#    K_hat

def quadriture(err1, err2):
    """Add two errors in quadriture."""
    
    return float(np.sqrt(err1**2 + err2**2))

def minor1d(array1d, i):
    """Accept 1D array and return a new array with element i removed."""
    
    range1 = list(range(i)) + list(range(i+1, array1d.shape[0]))
    return array1d[np.array(range1)]
    
def minor2d(array2d, i, j):
    """
    Accept matrix array2d and return a new matrix with row i and 
    column j removed.
    """
    
    range1 = list(range(i)) + list(range(i+1, array2d.shape[0]))
    range2 = list(range(j)) + list(range(j+1,array2d.shape[1]))
    return array2d[np.array(range1)[:,np.newaxis], np.array(range2)]

def pretty_string_value_error(value, error, error_digits=2, use_unicode=True):
    """
    Returns a value/error combination of numbers in a scientifically 
    'pretty' format.
    
    Scientific quantities often come as a *value* (the actual 
    quantity) and the *error* (the uncertainty in the value).
        
    Given two floats, value and error, return the two in a 
    'pretty' formatted string: where the value and error are truncated
    at the correct precision.
    
    Parameters
    ----------
    value : float
        The quantity in question
        
    error : float
        The uncertainty of the quantity
        
    error_digits : int, default 2
        How many significant figures the error has. Scientific 
        convention holds that errors have 1 or (at most) 2 significant
        figures. The larger number of digits is chosen here by default.
        
    Returns
    -------
    new_string : str
        A new list of strings sorted numerically
    
    Examples
    --------
    
        >>> pretty_string_value_error(1.23456789e8, 4.5678e5, 
                                      error_digits=2)
        "1.2346 +/- 0.0046 * 10^+08"
    
        >>> pretty_string_value_error(5.6e-2, 2.0e-3, error_digits=1)
        "5.6 +/- 0.2 * 10^-02"
    """
    
    if value is None:
        return "None"
    if error is None or not np.isfinite(error):
        if use_unicode:
            new_string = "{:.6E} \u00B1 UNKNOWN ERROR MARGIN".format(value)
        else:
            new_string = "{:.6E} +/- UNKNOWN ERROR MARGIN".format(value)
    else:
        if not np.isfinite(value):
            return str(value)
        assert "e" in "{:e}".format(value), "Cannot convert into scientific "\
            "notation: {1}".format(value)
        value_mantissa_str, value_exponent_str = \
            "{:e}".format(value).strip().split('e')
        value_mantissa = float(value_mantissa_str)
        value_exponent = int(value_exponent_str)
        error_mantissa_str, error_exponent_str = \
            "{:e}".format(error).strip().split('e')
        error_mantissa = float(error_mantissa_str)
        error_exponent = int(error_exponent_str)
        padding = value_exponent - error_exponent + error_digits - 1
        if padding < 1: padding = 1
        exp_diff = error_exponent - value_exponent
        string_for_formatting = "{:.%df}" % padding
        new_value_mantissa = string_for_formatting.format(value_mantissa)
        new_error_mantissa = string_for_formatting.format(
            error_mantissa*10**exp_diff)
        
        if use_unicode:
            new_string = "(%s \u00B1 %s) * 10^%s" % (
                new_value_mantissa, new_error_mantissa, value_exponent_str)
        else:
            new_string = "(%s +/- %s) * 10^%s" % (
                new_value_mantissa, new_error_mantissa, value_exponent_str)
    return new_string

def browndye_run_compute_rate_constant(compute_rate_constant_program,
                                       results_filename_list, 
                                       sample_error_from_normal=False):
    """
    run the BrownDye program compute_rate_constant to find k_ons
    and the value k(b).
    
    Parameters:
    -----------
    compute_rate_constant_program : str
        The exact command to use to run the Browndye 
        compute_rate_constant program.
    
    results_filename_list : list
        A list of string paths to the results XML file, which is an 
        output from running the BrownDye program.
    
    sample_error_from_normal : bool, default False
        Add a fluctuation to the k-on value(s) sampled from a
        normal distribution with a standard deviation equal to the
        k-on's error. This is used by the error estimators to 
        incorporate the error in the k_b value obtained from BrownDye.
    
    Results:
    --------
    k_ons : defaultdict
        dictionary whose keys are various milestoning states and
        whose values are the k-ons to those states.
    
    k_on_errors : defaultdict
        dictionary whose keys are various milestoning states and
        whose values are the errors in the k-ons to those states.
    
    reaction_probabilities : defaultdict
        dictionary whose keys are various milestoning states and
        whose values are the probabilities of reaching those states
        from the b-surface.
    
    reaction_probability_errors : defaultdict
        dictionary whose keys are various milestoning states and
        whose values are the errors in the probabilities of 
        reaching those states from the b-surface.
        
    transition_counts : defaultdict
        dictionary whose keys are the various milestone states and 
        also the "escaped" and "stuck" states. The values are the 
        counts of crossings of these various states.
    """
    
    k_ons = defaultdict(float)
    k_on_errors = defaultdict(float)
    reaction_probabilities = defaultdict(float)
    reaction_probability_errors = defaultdict(float)
    transition_counts = defaultdict(int)
    k_b = None
    total_n_trajectories = 0
    
    for results_filename in results_filename_list:
        cmd = "%s < %s" % (compute_rate_constant_program, 
                           results_filename)
        output_string = subprocess.check_output(cmd, shell=True)
        root = ET.fromstring(output_string)
        
        # first, obtain counts directly from the results file
        results_root = ET.parse(results_filename)
        reactions = results_root.find("reactions")
        n_trajectories = int(reactions.find("n_trajectories").text.strip())
        total_n_trajectories += n_trajectories
        stuck = int(reactions.find("stuck").text.strip())
        escaped = int(reactions.find("escaped").text.strip())
        for completed in reactions.iter("completed"):
            name = int(completed.find("name").text.strip())
            n = int(completed.find("n").text.strip())
            transition_counts[name] += n
        transition_counts["escaped"] += escaped
        transition_counts["stuck"] += stuck
        
        # now obtain probabilities and rates
        for rate in root.iter("rate"):
            name = int(rate.find("name").text.strip())
            rate_constant = rate.find("rate_constant")
            k_on_tag = rate_constant.find("mean")
            k_on_tag_high = rate_constant.find("high")
            assert k_on_tag is not None
            assert k_on_tag_high is not None
            k_on = float(k_on_tag.text)
            k_on_error = 0.5*(float(k_on_tag_high.text) - k_on)
            reaction_probability = rate.find("reaction_probability")
            beta_tag = reaction_probability.find("mean")
            beta_tag_high = reaction_probability.find("high")
            assert beta_tag is not None
            assert beta_tag_high is not None
            beta = float(beta_tag.text)
            beta_error = 0.5*(float(beta_tag_high.text) - beta)
            if beta != 0.0:
                k_b = k_on / beta
            
            if sample_error_from_normal:
                k_on = np.random.normal(loc=k_on, scale=k_on_error)
            
    assert k_b is not None, "No BD reactions from the b-surface " \
        "successfully reached any of the milestone surfaces."
    
    for key in transition_counts:
        reaction_probabilities[key] \
            = transition_counts[key] / total_n_trajectories
        k_ons[key] = k_b * reaction_probabilities[key]
        
        if transition_counts[key] > 1:
            reaction_probability_errors[key] = reaction_probabilities[key] \
                / np.sqrt(transition_counts[key]-1)
            k_on_errors[key] = k_b * reaction_probability_errors[key]
        else:
            reaction_probability_errors[key] = 1e99
            k_on_errors[key] = 1e99
        
        if sample_error_from_normal:
            k_ons[key] = np.random.normal(loc=k_ons[key], 
                                          scale=k_on_errors[key])
    
    assert k_ons.keys() == transition_counts.keys()
    return k_ons, k_on_errors, reaction_probabilities, \
        reaction_probability_errors, transition_counts
        
def browndye_parse_bd_milestone_results(results_filename_list, 
                                        sample_error_from_normal=False):
    """
    Read and extract transition probabilities for a BD milestone.
    
    Parameters:
    -----------
    results_filename_list : list
        A list of paths to the results XML files for a BD milestone,
        which is an output from running the BrownDye program.
        
    sample_error_from_normal : bool, default False
        Add a fluctuation to the probabilities sampled from a
        normal distribution with a standard deviation equal to
        1/sqrt(n-1). This is used by the error estimators to 
        incorporate the error in the probabilities.
    
    Results:
    --------
    transition_probabilities : defaultdict
        The probabilities of transitions.
    
    transition_counts : defaultdict
        The counts of transitions to all visited states.    
    """
    
    transition_counts = defaultdict(int)
    transition_probabilities = defaultdict(float)
    completed_prob = 0.0
    completed_count = 0
    
    for results_filename in results_filename_list:
        assert os.path.exists(results_filename), "You must perform successful "\
            "extractions and runs of all BD milestones if k-on settings are "\
            "provided in the model XML. Missing file: " \
            + results_filename
        root = ET.parse(results_filename)
        reactions = root.find("reactions")
        n_trajectories = int(reactions.find("n_trajectories").text.strip())
        stuck = int(reactions.find("stuck").text.strip())
        escaped = int(reactions.find("escaped").text.strip())
        
        for completed in reactions.iter("completed"):
            name = int(completed.find("name").text.strip())
            n = int(completed.find("n").text.strip())
            transition_counts[name] += n
            completed_count += n
            
        transition_counts["escaped"] += escaped
        transition_counts["stuck"] += stuck
    
    for key in transition_counts:
        n = transition_counts[key]
        avg = n / (completed_count + escaped)
        if sample_error_from_normal:
            assert n > 1, "Too few transitions to compute error."
            std_dev = avg / np.sqrt(n-1)
            transition_probabilities[key] = np.random.normal(
                loc=avg, scale=std_dev)
        else:
            transition_probabilities[key] = avg
            
        if key not in ["escaped", "stuck"]:
            completed_prob += transition_probabilities[key]
    
    transition_probabilities["escaped"] = 1.0 - completed_prob
    return transition_probabilities, transition_counts

def combine_fhpd_results(bd_milestone, fhpd_directories, 
                         combined_results_filename):
    """
    Read all results files in the first-hitting-point-distribution 
    (FHPD) directories and combine them into a single results file
    to be placed in the bd_milestone directory.
    """
    
    reaction_dict = defaultdict(float)
    number_escaped = 0
    number_stuck = 0
    number_total = 0
    number_total_check = 0
    results_filename_list = []
    if len(fhpd_directories) == 0:
        return
    
    for fhpd_directory in fhpd_directories:
        results_glob = os.path.join(fhpd_directory, 
                                    "results*.xml")
        results_filename_list += glob.glob(results_glob)
    
    if len(results_filename_list) == 0:
        print("No BD output files found.")
        return
    
    for results_filename in results_filename_list:
        tree = ET.parse(results_filename)
        root = tree.getroot()
        reactions_XML = root.find("reactions")
        number_total += int(reactions_XML.find("n_trajectories").text.strip())
        number_stuck += int(reactions_XML.find("stuck").text.strip())
        number_escaped += int(reactions_XML.find("escaped").text.strip())
        for completed_XML in reactions_XML.iter("completed"):
            name = completed_XML.find("name").text.strip()
            n = int(completed_XML.find("n").text.strip())
            reaction_dict[name] += n
            number_total_check += n
        
    assert number_total == number_total_check + number_stuck + number_escaped
    for completed_XML in reactions_XML.iter("completed"):
        reactions_XML.remove(completed_XML)
    
    reactions_XML.find("n_trajectories").text = str(number_total)
    reactions_XML.find("stuck").text = str(number_stuck)
    reactions_XML.find("escaped").text = str(number_escaped)
    for key in reaction_dict:
        completed_XML = ET.SubElement(reactions_XML, "completed")
        completed_XML.text = "\n      "
        completed_XML.tail = "\n  "
        name_XML = ET.SubElement(completed_XML, "name")
        name_XML.text = key
        name_XML.tail = "\n      "
        n_XML = ET.SubElement(completed_XML, "n")
        n_XML.text = str(int(reaction_dict[key]))
        n_XML.tail = "\n    "
        
    xmlstr = ET.tostring(root).decode("utf-8")
    with open(combined_results_filename, 'w') as f:
        f.write(xmlstr)
        
    return

class Data_sample():
    """
    Represent a set of data needed to compute kinetic or thermodynamic
    quantities of interest using the SEEKR method. The Data_sample may
    be constructed directly from average times and rates, or may be 
    generated from a random distribution.
    
    Attributes:
    -----------
    model : Model()
        The model contains all the settings, parameters, directories,
        and file names used to perform a SEEKR2 simulation.
        
    N_ij : numpy.array()
        2-dimensional array of transition counts from milestone i to 
        milestone j.
    
    R_i : numpy.array()
        1-dimensional array of incubation times for each milestone i.
    
    Q : numpy.array()
        The rate matrix constructed from transition counts and times.
        No sink state is defined - probability is conserved for 
        thermodynamics calculations.
    
    Q_hat : numpy.array()
        Same as attribute Q, but one or more sink states are defined - 
        probability is not conserved, allowing kinetics calculations.
    
    K : numpy.array()
        The probability transition matrix made only from transition
        counts - it is equivalent to a Markov transition matrix.
        Probability is conserved for thermodynamics calculations.
    
    K_hat : numpy.array()
        Same as attribute Q, except one or more sink states are defined
        and probability is not conserved, allowing kinetics 
        calculations.
    
    p_i : numpy.array()
        Stationary probability distribution along the milestones, 
        which can be used to compute other thermodynamic quantities, 
        such as the attribute free_energy_profile.
    
    free_energy_profile : numpy.array()
        The stationary free energy profile along the milestones.
        
    MFPTs : dict
        A dictionary whose keys are 2-tuples of milestone indices
        (source milestone, destination milestone) and whose values
        are mean first passage times (MFPTs) between these states.
        States can be added to this dictionary by being either end
        states or the bulk state.
        
    k_off : float
        The off-rate constant from a stationary distribution across
        all the milestones to the bulk state.
        
    k_on : dict
        A dictionary whose keys are milestone indices (which are end
        states) and whose values represent the 2nd-order on-rate 
        constant.
        
    bd_transition_counts : dict
        Keep track of transition counts for each BD state. This is a
        dictionary of dictionaries, with keys of "b_surface" or 
        integers representing BD milestones. The values are 
        dictionaries whose keys are states such as milestone indices
        or the string "escaped", "stuck", and whose values are counts
        of transitions. This attribute is primarily used for 
        convergence estimates.
    """
    
    def __init__(self, model):
        self.model = model
        self.N_ij = None
        self.R_i = None
        self.Q = None
        self.Q_hat = None
        self.K = None
        self.K_hat = None
        self.p_i = None
        #self.p_i_hat = None
        self.free_energy_profile = None
        self.MFPTs = None
        self.k_off = None
        self.k_ons = {}
        self.b_surface_k_ons_src = None
        self.b_surface_k_on_errors_src = None
        #self.b_surface_reaction_probabilities = None # REMOVE?
        #self.b_surface_reaction_probability_errors = None # REMOVE?
        #self.b_surface_transition_counts = None # REMOVE?
        self.bd_transition_counts = {}
        self.bd_transition_probabilities = {}
    
    def parse_browndye_results(self, bd_sample_from_normal=False):
        """
        Parse Browndye2 output files to fill out the milestoning model.
        
        Parameters:
        -----------
        
        bd_sample_from_normal : bool, default False
            If set to True, then k-on quantities will have a random
            fluctuation introduced in a magnitude proportional to k-on
            errors. This is used only for error estimations.
        """
        
        b_surface_output_file_glob = os.path.join(
            self.model.anchor_rootdir,
            self.model.k_on_info.b_surface_directory, 
            self.model.k_on_info.bd_output_glob)
        output_file_list = glob.glob(b_surface_output_file_glob)
        output_file_list = base.order_files_numerically(output_file_list)
        if len(output_file_list) > 0:
            if self.model.browndye_settings is not None:
                k_ons_src, k_on_errors_src, reaction_probabilities, \
                    reaction_probability_errors, transition_counts = \
                    browndye_run_compute_rate_constant(os.path.join(
                        self.model.browndye_settings.browndye_bin_dir,
                        "compute_rate_constant"), output_file_list, 
                        sample_error_from_normal=bd_sample_from_normal)
                self.bd_transition_counts["b_surface"] = transition_counts
                self.bd_transition_probabilities["b_surface"] = reaction_probabilities
                self.b_surface_k_ons_src = k_ons_src
                self.b_surface_b_surface_k_on_errors_src = k_on_errors_src
                #self.b_surface_reaction_probabilities = reaction_probabilities # REMOVE?
                #self.b_surface_reaction_probability_errors = reaction_probability_errors # REMOVE?
            else:
                raise Exception("No valid BD program settings provided.")
            
        if len(self.model.k_on_info.bd_milestones) > 0:
            for bd_milestone in self.model.k_on_info.bd_milestones:
                bd_milestone_results_file = os.path.join(
                    self.model.anchor_rootdir, bd_milestone.directory, 
                    "results.xml")
                
                if not os.path.exists(bd_milestone_results_file):
                    bd_directory_list_glob = os.path.join(
                        self.model.anchor_rootdir, 
                        bd_milestone.directory, 
                        "first_hitting_point_distribution", "lig*/")
                    bd_directory_list = glob.glob(
                        bd_directory_list_glob)
                    if len(bd_directory_list) == 0:
                        continue
                    combine_fhpd_results(bd_milestone, 
                                         bd_directory_list, 
                                         bd_milestone_results_file)
                    
                transition_probabilities, transition_counts = \
                    browndye_parse_bd_milestone_results(
                        [bd_milestone_results_file])
                self.bd_transition_counts[bd_milestone.index] \
                    = transition_counts
                self.bd_transition_probabilities[bd_milestone.index] \
                        = transition_probabilities
                                                
        return
    
    def compute_rate_matrix(self):
        """
        Compute Q and K from N_ij and R_i.
        """
        self.Q = np.zeros((self.model.num_milestones, 
                           self.model.num_milestones), dtype=np.longdouble)
        for i in range(self.model.num_milestones):
            for j in range(self.model.num_milestones):
                if self.R_i[i] == 0.0:
                    self.Q[i,j] = 0.0
                else:
                    self.Q[i,j] = self.N_ij[i,j] / self.R_i[i]
                assert self.Q[i,j] >= 0.0, "self.Q[i,j]: {}".format(self.Q[i,j])
                    
        for i in range(self.model.num_milestones):
            row_sum = np.sum(self.Q[i])
            if row_sum == 0:
                new_row_sum = 0.0
                for j in range(self.model.num_milestones):
                    self.Q[i,j] = self.Q[j,i]
                    new_row_sum += self.Q[j,i]
                    
                self.Q[i,i] = -new_row_sum
            else:
                self.Q[i,i] = -row_sum
        
        self.K = Q_to_K(self.Q)
        return
    
    def calculate_thermodynamics(self):
        """
        Use this data sample's statistics to construct the 
        thermodynamic quantities.
        """
        eigvals, eigvecs = la.eig(self.K.T)
        closest_eigval = -1e9
        closest_eigvec = None
        for i, eigval, in enumerate(eigvals):
            if (eigval - 1.0)**2 < (closest_eigval - 1.0)**2:
                closest_eigval = eigval
                closest_eigvec = eigvecs[:,i]
        
        K = abs(self.K)
        q_i = np.real(closest_eigvec / np.sum(closest_eigvec))
        q_i = q_i @ np.linalg.matrix_power(K, MATRIX_EXPONENTIAL)
        q_i = q_i / np.sum(q_i)
        incubation_times = np.zeros(self.model.num_milestones)
        p_i = np.zeros(self.model.num_milestones)
        for i in range(self.model.num_milestones):
            N_i = 0
            for j in range(self.model.num_milestones):
                N_i += self.N_ij[i,j]
            
            if N_i == 0:
                incubation_times[i] = -1.0/self.Q[i,i]
            else:
                incubation_times[i] = self.R_i[i] / N_i

        for i, q_i_val in enumerate(q_i):
            p_i[i] = abs(q_i_val * incubation_times[i])
        
        sum_p_i = np.sum(p_i)
        for i in range(self.model.num_milestones):
            p_i[i] /= sum_p_i
        
        self.q_i = q_i
        self.p_i = p_i
        self.free_energy_profile = np.zeros(self.p_i.shape)
        highest_p_i = max(self.p_i)
        for i, p_i_val in enumerate(self.p_i):
            free_energy = -GAS_CONSTANT*self.model.temperature*np.log(
                p_i_val / highest_p_i)
            self.free_energy_profile[i] = free_energy
        
        return
        
    def calculate_kinetics(self, pre_equilibrium_approx=False):
        """
        Once the rate matrix Q is computed, determine the timescales 
        and probabilities of transfers between different states. Fill
        out all kinetics quantities.
        
        Parameters:
        -----------
        
        pre_equilibrium_approx : bool, default False
            Whether to use the pre-equilibrium approximation for
            computing kinetics.
            
        
        """
        
        end_milestones = []
        bulk_milestones = []
        MFPTs = {}
        k_off = 0.0
        k_ons = {}
        for alpha, anchor in enumerate(self.model.anchors):
            if anchor.endstate:
                for milestone_id in anchor.get_ids():
                    if self.model.get_type() == "elber":
                        if anchor.alias_from_id(milestone_id) == 3: 
                            # TODO: hacky
                            continue
                    end_milestones.append(milestone_id)
            if anchor.bulkstate:
                for milestone_id in anchor.get_ids():
                    bulk_milestones.append(milestone_id)

        # first, make the bulk state the sink state to compute k_offs
        Q_hat = self.Q[:,:]
        p_i_hat = self.p_i[:]
        #if self.model.k_on_info:
        #    K_hat = self.K[:,:]
        
        n = len(self.Q)
        for bulk_milestone in sorted(bulk_milestones, reverse=True):
            Q_hat = minor2d(Q_hat, bulk_milestone, bulk_milestone)
            p_i_hat = minor1d(p_i_hat, bulk_milestone)
        
        Q_hat = Q_hat.astype(dtype=np.longdouble)
        if pre_equilibrium_approx:
            lowest_p_i = np.min(self.p_i)
            lowest_i = np.argmin(self.p_i)
            assert lowest_p_i >= 0.0, \
                "Negative stationary probability detected."
            if lowest_i == n-1:
                k_off = lowest_p_i * Q_hat[lowest_i-1,lowest_i]
            else: 
                k_off = lowest_p_i * Q_hat[lowest_i,lowest_i+1]
            bulk_times = np.ones(p_i_hat.shape) / k_off
            
        else:
            #negative_unity = np.zeros((len(Q_hat)), dtype=np.longdouble)    
            #negative_unity[:] = -1.0
            #bulk_times = la.solve(Q_hat, negative_unity)
            bulk_times = solve_rate_matrix(Q_hat)
        
        for end_milestone in end_milestones:
            if end_milestone in bulk_milestones:
                continue
            # must account for the removal of bulk state to matrix indices
            no_bulk_index = end_milestone
            for bulk_milestone in bulk_milestones:
                if end_milestone > bulk_milestone: 
                    no_bulk_index -= 1
                
            mfpt = bulk_times[no_bulk_index]    
            MFPTs[(end_milestone, "bulk")] = mfpt
        
        MFPT_to_bulk = 0
        assert bulk_times.shape == p_i_hat.shape
        for i, bulk_time in enumerate(bulk_times):
            MFPT_to_bulk += bulk_time * p_i_hat[i]
        
        # convert to 1/s
        k_off = 1.0e12 / MFPT_to_bulk
        
        # Next, compute the MFPTs between different states
        for end_milestone_dest in end_milestones:
            if end_milestone_dest in bulk_milestones:
                continue
            Q_hat = minor2d(self.Q[:], end_milestone_dest, end_milestone_dest)
            #I = np.zeros((len(Q_hat)), dtype = float)    
            #I[:] = 1.0
            #end_state_times = la.solve(Q_hat, -I)
            end_state_times = solve_rate_matrix(Q_hat)
            for end_milestone_src in end_milestones:
                if end_milestone_dest == end_milestone_src:
                    # don't get the MFPT from a milestone to itself
                    continue
                if end_milestone_src in bulk_milestones:
                    # a bulk milestone will never be a source
                    continue
                mfpt = end_state_times[end_milestone_src]
                MFPTs[(end_milestone_src, end_milestone_dest)] = mfpt
        
        if self.model.k_on_info:
            K_hat = self.K[:,:]
            for end_milestone in end_milestones:
                K_hat[end_milestone, :] = 0.0
                K_hat[end_milestone, end_milestone] = 1.0
            
            p_i_hat = self.p_i[:]
            n = K_hat.shape[0]
            source_vec = np.zeros((n,1))
            
            if len(self.bd_transition_counts) > 0:
                if len(bulk_milestones) > 0:
                    bulk_milestone = bulk_milestones[0]
                    for bd_milestone in self.model.k_on_info.bd_milestones:
                        source_index = bd_milestone.outer_milestone.index
                        if self.b_surface_k_ons_src is None:
                            raise Exception("Missing b-surface simulations.")
                            
                        source_vec[source_index] \
                            = self.b_surface_k_ons_src[source_index]
                            
                        transition_probabilities \
                            = self.bd_transition_probabilities[
                                bd_milestone.index]
                        K_hat[source_index, :] = 0.0
                        for key in transition_probabilities:
                            value = transition_probabilities[key]
                            if key in ["escaped", "stuck"]:
                                pass
                            else:
                                K_hat[source_index, key] = value

                K_hat_inf = np.linalg.matrix_power(K_hat, MATRIX_EXPONENTIAL)
                end_k_ons = np.dot(K_hat_inf.T, source_vec)
                for end_milestone in end_milestones:
                    k_ons[end_milestone] = end_k_ons[end_milestone]
                self.K_hat = K_hat
                self.k_ons = k_ons
        
        self.Q_hat = Q_hat
        #self.p_i_hat = p_i_hat # TODO: remove after successful CI test
        self.MFPTs = MFPTs
        self.k_off = k_off
        return
    
    def monte_carlo_milestoning_error(self, num=1000, skip=100, stride=1, 
                                      verbose=False, 
                                      pre_equilibrium_approx=False):
        """
        Calculates an error estimate by sampling a distribution of rate 
        matrices assumming a Poisson (gamma) distribution with 
        parameters Nij and Ri using Markov chain Monte Carlo.
            
        Enforces detailed Balance-- using a modified version of 
        Algorithm 4 from Noe 2008 for rate matrices.-- 
        Distribution is:  p(Q|N) = p(Q)p(N|Q)/p(N) = 
            p(Q) PI(q_ij**N_ij * exp(-q_ij * Ri))
            
        Parameters
        ----------
        
        num : int, default 1000
            number of rate matrix (Q) samples to be generated
            
        skip : int, default 100
            number of inital rate matrix samples to skip for "burn in"
            
        stride : int, default 1
            frequency at which rate matrix samples are recorded- larger
            frequency reduces correlation between samples
            
        verbose : bool, default False
            allow additional verbosity/printing
        
        pre_equilibrium_approx : bool, default False
            Whether to use the pre-equilibrium approximation for
            computing kinetics.
        """
        N = self.N_ij
        R = self.R_i
        Q = self.Q[:]
        m = Q.shape[0] #get size of count matrix
        Q_mats = []
        MFPTs_list = defaultdict(list)
        MFPTs_error = defaultdict(float)
        k_off_list = []
        k_ons_list = defaultdict(list)
        k_ons_error = defaultdict(float)
        p_i_list = []
        free_energy_profile_list = []
        Qnew = Q[:,:]
        if verbose: print("collecting ", num, " MCMC samples from ", 
                          num*(stride) + skip, " total moves")  
        for counter in range(num * (stride) + skip):
            #if verbose: print("MCMC stepnum: ", counter)
            Qnew = Q[:,:]
            for i in range(m): # rows
                for j in range(m): # columns
                    if i == j: continue
                    if Qnew[i,j] == 0.0: continue
                    if Qnew[j,j] == 0.0: continue
                    if N[i,j] == 0: continue
                    if R[i] == 0: continue
                    Q_gamma = 0
                    delta = Qnew[i,j]
                    while ((delta) >= (Qnew[i,j])):
                        Q_gamma = gamma.rvs(a=N[i,j], scale = 1/R[i],)
                        delta =  Qnew[i,j] - Q_gamma
        
                    log_p_Q_old = N[i,j] * np.log(Qnew[i,j])  - Qnew[i,j] * R[i] 
                    log_p_Q_new = N[i,j] * np.log(Qnew[i,j] - delta) - \
                        (Qnew[i,j] - delta) * R[i]     
                    if verbose: print("log P(Q_new)", log_p_Q_new)
                    if verbose: print("log P(Q_old)", log_p_Q_old)
    
                    r2 = random.random()  
                    p_acc =  log_p_Q_new - log_p_Q_old
                    if verbose: print("p_acc", p_acc, "r", np.log(r2))
                        
                    if np.log(r2) <= p_acc: 
                        #log(r) can be directly compared to 
                        # log-likelihood acceptance, p_acc
                        if verbose: print("performing non-reversible element "\
                                          "shift...")
                        Qnew[i,i] = (Qnew[i,i]) + delta
                        Qnew[i,j] = Qnew[i,j] - delta
                        if verbose: print(Qnew)
    
            if counter > skip and counter % stride == 0:
                self.Q = Qnew
                self.K = Q_to_K(self.Q)
                self.calculate_thermodynamics()
                self.calculate_kinetics(pre_equilibrium_approx)
                p_i_list.append(self.p_i)
                free_energy_profile_list.append(self.free_energy_profile)
                for key in self.MFPTs:
                    MFPTs_list[key].append(self.MFPTs[key])
                k_off_list.append(self.k_off)
                if self.model.k_on_info:
                    for key in self.k_ons:
                        k_ons_list[key].append(self.k_ons[key])  
     
            Q = Qnew
        if verbose: print("final MCMC matrix", Q)
        p_i_error = np.zeros(self.p_i.shape)
        free_energy_profile_err = np.zeros(self.free_energy_profile.shape)
        for i in range(p_i_error.shape[0]):
            p_i_val_list = []
            for j in range(len(p_i_list)):
                p_i_val_list.append(p_i_list[j][i])
            p_i_error[i] = np.std(p_i_val_list)
        
        for i in range(free_energy_profile_err.shape[0]):
            free_energy_profile_val_list = []
            for j in range(len(free_energy_profile_list)):
                free_energy_profile_val_list.append(free_energy_profile_list[j][i])
            free_energy_profile_err[i] = np.std(free_energy_profile_val_list)
        for key in self.MFPTs:
            MFPTs_error[key] = np.std(MFPTs_list[key])
        k_off_error = np.std(k_off_list)
        if self.model.k_on_info:
            for key in self.k_ons:
                k_ons_error[key] = np.std(k_ons_list[key])
        
        return p_i_error, free_energy_profile_err, MFPTs_error, k_off_error, \
            k_ons_error
