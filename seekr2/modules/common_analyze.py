"""
Functions and objects for analyzing SEEKR2 simulation outputs
"""

import os
import glob
import xml.etree.ElementTree as ET
import subprocess
import random
from collections import defaultdict

import numpy as np
#from numpy import log, exp
from scipy import linalg as la
from scipy.stats import gamma

import seekr2.modules.common_base as base

MATRIX_EXPONENTIAL = 9999999
GAS_CONSTANT = 0.0019872 # in kcal/mol*K
BOLTZMANN = 1.3806488e-23
MAX_MATRIX_CONDITION = 1e16
DEFAULT_IMAGE_DIR = "images_and_plots/"

class MissingStatisticsError(Exception):
    # TODO: replace with generic error
    pass

class IllCondMatrixError(Exception):
    # TODO: replace with generic error
    pass

def solve_rate_matrix(Q_hat, max_iter=100, ETOL=1e-10):
    """
    Use an infinite series combined with a divide-and-conquer method
    to solve the Q-matrix in case the numpy/scipy solver has numerical
    instabilities.
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
    """
    Given a rate matrix Q, compute the probability transition matrix
    K.
    """
    K = np.zeros(Q.shape, dtype=np.longdouble)
    for i in range(Q.shape[0]):
        for j in range(Q.shape[0]):
            if i == j:
                K[i,j] = 0.0
            else:
                K[i,j] = -Q[i,j] / Q[i,i]
    return K

def quadriture(err1, err2):
    """
    Add two errors in quadriture.
    """
    return float(np.sqrt(err1**2 + err2**2))

def minor1d(array1d, i):
    """
    Accept 1D array array1d and return a new array with element i 
    removed.
    """
    range1 = list(range(i)) + list(range(i+1, array1d.shape[0]))
    return array1d[np.array(range1)]
    
def minor2d(array2d, i, j):
    """
    Accept matrix array2d and returns a new matrix with row i and 
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
            new_string = "%s \u00B1 %s * 10^%s" % (
                new_value_mantissa, new_error_mantissa, value_exponent_str)
        else:
            new_string = "%s +/- %s * 10^%s" % (
                new_value_mantissa, new_error_mantissa, value_exponent_str)
    return new_string

def browndye_run_compute_rate_constant(compute_rate_constant_program,
                                       results_filename, 
                                       sample_error_from_normal=False):
    # common
    """
    run the BrownDye program compute_rate_constant to find k_ons
    and the value k(b).
    
    Parameters:
    -----------
    compute_rate_constant_program : str
        The exact command to use to run the Browndye 
        compute_rate_constant program.
    
    results_filename : str
        A path to the results XML file, which is an output from
        running the BrownDye program.
    
    sample_error_from_normal : bool, default False
        Add a fluctuation to the k-on value(s) sampled from a
        normal distribution with a standard deviation equal to the
        k-on's error. This is used by the Monte Carlo error
        estimator to incorporate the error in the k_b value
        obtained from BrownDye.
    
    Results:
    --------
    k_ons : dict
        dictionary whose keys are various milestoning states and
        whose values are the k-ons to those states.
    
    k_on_errors : dict
        dictionary whose keys are various milestoning states and
        whose values are the errors in the k-ons to those states.
    
    reaction_probabilities : dict
        dictionary whose keys are various milestoning states and
        whose values are the probabilities of reaching those states
        from the b-surface.
    
    reaction_probability_errors : dict
        dictionary whose keys are various milestoning states and
        whose values are the errors in the probabilities of 
        reaching those states from the b-surface.
    
    """
    cmd = "%s < %s" % (compute_rate_constant_program, 
                       results_filename)
    #print("running command:", cmd)
    output_string = subprocess.check_output(cmd, shell=True)
    root = ET.fromstring(output_string)
    k_ons = {}
    k_on_errors = {}
    reaction_probabilities = {}
    reaction_probability_errors = {}
    k_b = None
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
        beta_error = float(beta_tag_high.text) - beta
        reaction_probabilities[name] = beta
        reaction_probability_errors[name] = beta_error
        if beta != 0.0:
            k_b = k_on / beta
        
        if sample_error_from_normal:
            k_on = np.random.normal(loc=k_on, scale=k_on_error)
            
        k_ons[name] = k_on
        k_on_errors[name] = k_on_error
    
    assert k_b is not None, "No BD reactions from the b-surface " \
        "successfully reached any of the milestone surfaces."
    #k_ons["b_surface"] = k_b
    #k_on_errors["b_surface"] = 0.0
    #reaction_probabilities["b_surface"] = 1.0
    #reaction_probability_errors["b_surface"] = 0.0
    return k_ons, k_on_errors, reaction_probabilities, \
        reaction_probability_errors
        
def browndye_parse_bd_milestone_results(results_filename, 
                                        sample_error_from_normal=False):
    # common
    """
    Read and extract transition probabilities for a BD milestone.
    
    Parameters:
    -----------
    results_filename : str
        A path to the results XML file for a BD milestone, which 
        is an output from running the BrownDye program.
        
    sample_error_from_normal : bool, default False
        Add a fluctuation to the probabilities sampled from a
        normal distribution with a standard deviation equal to
        1/sqrt(n-1). This is used by the Monte Carlo error
        estimator to incorporate the error in the probabilities.
    
    Results:
    --------
    transition_probabilities : dict
        The probabilities of transitions.
    """
    transition_counts = {}
    transition_probabilities = {}
    assert os.path.exists(results_filename), "You must perform successful "\
        "extractions and runs of all BD milestones if k-on settings are "\
        "provided in the model XML. Missing file: " \
        + results_filename
    root = ET.parse(results_filename)
    reactions = root.find("reactions")
    n_trajectories = int(reactions.find("n_trajectories").text.strip())
    stuck = int(reactions.find("stuck").text.strip())
    escaped = int(reactions.find("escaped").text.strip())
    completed_prob = 0.0
    completed_count = 0
    for completed in root.iter("completed"):
        name = int(completed.find("name").text.strip())
        n = int(completed.find("n").text.strip())
        transition_counts[name] = n
        completed_count += n
        
    transition_counts["escaped"] = escaped
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
            
        completed_prob += transition_probabilities[key]
        
    transition_probabilities["escaped"] = 1.0 - completed_prob
    return transition_probabilities

class Data_sample():
    # common
    """
    Represent a set of data needed to compute kinetic or thermodynamic
    quantities of interest using the SEEKR method. The Data_sample may
    be constructed directly from average times and rates, or may be 
    generated from a random distribution.
    
    Attributes:
    -----------
    model : Model()
    
    N_ij :
    
    R_i :
    
    Q : 
    
    Q_hat : 
    
    K : 
    
    K_hat :
    
    p_i : 
    
    p_i_hat : 
    
    free_energy_profile : 
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
        self.p_i_hat = None
        self.free_energy_profile = None
        self.MFPTs = None
        self.k_off = None
        self.k_ons = {}
    
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
                    
        for i in range(self.model.num_milestones):
            self.Q[i][i] = -np.sum(self.Q[i])
        
        self.K = Q_to_K(self.Q)
        return
    
    def calculate_thermodynamics(self):
        """
        With the rate matrix computed, make a probability matrix
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
        
    def calculate_kinetics(self, pre_equilibrium_approx=False, 
                           bd_sample_from_normal=False):
        """
        Once the rate matrix Q is computed, determine the timescales 
        and probabilities of transfers between different states.
        
        Parameters:
        -----------
        
        assign_to_self : bool, default True
            Assign all results from this calculation, including MFPT
            values, k-off, and k-on values to this Analyze object. This
            parameter would only be False if a different Q, p_i, or K
            are needed to, say, sample a distribution of those 
            quantities in an error estimation.
            
        Q : numpy.array or None, default None
            If provided, then a different rate matrix Q will be used
            to compute kinetics and thermodynamics than self.Q.
            
        p_i : numpy.array or None, default None
            If provided, then a different probability distribution p_i
            than self.p_i will be used to compute thermodynamics.
            
        K : numpy.array or None, default None
            If provided, then a different transition probability matrix
            K will be used instead of self.K to compute thermodynamics
            and kinetics.
            
        bd_sample_from_normal : bool, default False
            If set to True, then k-on quantities will have a random
            fluctuation introduced in a magnitude proportional to k-on
            errors. This is used only for error estimations.
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
                            continue
                    end_milestones.append(milestone_id)
            if anchor.bulkstate:
                for milestone_id in anchor.get_ids():
                    bulk_milestones.append(milestone_id)

        # first, make the bulk state the sink state to compute k_offs
        Q_hat = self.Q[:,:]
        p_i_hat = self.p_i[:]
        if self.model.k_on_info:
            K_hat = self.K[:,:]
        
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
            """
            if np.linalg.cond(Q_hat) > MAX_MATRIX_CONDITION:
                warnstr="The MMVT rate matrix Q_hat is ill-conditioned. This "\
                      "is most likely caused by very long-timescale kinetics "\
                      "within this system. It is now recommended to use the "\
                      "--pre_equilibrium_approx (-p) option, which does not "\
                      "use the Q_mat matrix. Alternatively, you can use the "\
                      "--force_warning (-f) option to forge ahead, at your "\
                      "own risk, with normal MMVT. When choosing the latter "\
                      "option, keep in mind that ill-conditioned matrices "\
                      "have been observed to produce inaccurate and/or "\
                      "nonphysical kinetics."
                if self.force_warning:
                    warnings.warn(warnstr)
                else:
                    raise IllCondMatrixError(warnstr)
            """
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
                    # doesn't make sense to get the MFPT from a milestone to
                    #  itself
                    continue
                if end_milestone_src in bulk_milestones:
                    # a bulk milestone will never be a source
                    continue
                mfpt = end_state_times[end_milestone_src]
                MFPTs[(end_milestone_src, end_milestone_dest)] = mfpt
        
        if self.model.k_on_info:
            #K_hat = self.K[:,:]
            p_i_hat = self.p_i[:]
            for end_milestone in end_milestones:
                K_hat[end_milestone, :] = 0.0
                K_hat[end_milestone, end_milestone] = 1.0
            n = K_hat.shape[0]
            source_vec = np.zeros((n,1))
            output_file_glob = os.path.join(
                self.model.anchor_rootdir,
                self.model.k_on_info.b_surface_directory, 
                self.model.k_on_info.bd_mmvt_output_glob)
            output_file_list = glob.glob(output_file_glob)
            output_file_list = base.order_files_numerically(output_file_list)
            #assert len(output_file_list) > 0, "No BD output file found"
            if len(output_file_list) > 0:
                assert len(output_file_list) == 1, "More than one BD output file "\
                    "in b_surface folder not yet supported"
                results_filename = output_file_list[0]
                if self.model.browndye_settings is not None:
                    k_ons_src, k_on_errors_src, reaction_probabilities, \
                        reaction_probability_errors = \
                        browndye_run_compute_rate_constant(os.path.join(
                            self.model.browndye_settings.browndye_bin_dir,
                            "compute_rate_constant"), results_filename, 
                            sample_error_from_normal=bd_sample_from_normal)
                else:
                    raise Exception("No valid BD program settings provided.")
                
                if len(bulk_milestones) > 0:
                    bulk_milestone = bulk_milestones[0]
                    for bd_milestone in self.model.k_on_info.bd_milestones:
                        source_vec[bd_milestone.outer_milestone.index] = \
                            k_ons_src[bd_milestone.outer_milestone.index]
                        results_filename = os.path.join(
                            self.model.anchor_rootdir, bd_milestone.directory,
                            bd_milestone.bd_mmvt_output_glob)
                        transition_probabilities = \
                            browndye_parse_bd_milestone_results(results_filename)
                        src_index = bd_milestone.outer_milestone.index
                        K_hat[src_index, :] = 0.0
                        for key in transition_probabilities:
                            value = transition_probabilities[key]
                            if key == "escaped":
                                pass
                            else:
                                K_hat[src_index, key] = value
                
                K_hat_inf = np.linalg.matrix_power(K_hat, MATRIX_EXPONENTIAL)
                
                end_k_ons = np.dot(K_hat_inf.T, source_vec)
                for end_milestone in end_milestones:
                    k_ons[end_milestone] = end_k_ons[end_milestone]
                self.K_hat = K_hat
                self.k_ons = k_ons
        
        self.Q_hat = Q_hat
        self.p_i_hat = p_i_hat
        self.MFPTs = MFPTs
        self.k_off = k_off
        
        """
        self.p_i_error = np.zeros(p_i.shape)
        for i in range(self.p_i_error.shape[0]):
            p_i_val_list = []
            for j in range(len(p_i_list)):
                p_i_val_list.append(p_i_list[j][i])
            self.p_i_error[i] = np.std(p_i_val_list)
        
        self.free_energy_profile_err = np.zeros(self.p_i.shape)
        highest_p_i = max(self.p_i)
        for i, p_i_val in enumerate(self.p_i):
            p_i_val_err = self.p_i_error[i]
            free_energy_err = GAS_CONSTANT*self.model.temperature*highest_p_i\
                *p_i_val_err/p_i_val
            self.free_energy_profile_err[i] = free_energy_err
        """
        
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
                self.calculate_kinetics(pre_equilibrium_approx,
                                        bd_sample_from_normal=True)
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
