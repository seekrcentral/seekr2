"""
common_analyze.py

Classes, functions and other objects used for analyzing SEEKR2 
simulation outputs common to MMVT and Elber calculation types.
"""

import os
import glob
import xml.etree.ElementTree as ET
import subprocess
from collections import defaultdict
from copy import deepcopy

import numpy as np
from scipy import linalg as la
import PyGT

import seekr2.modules.common_base as base

# The power to raise a matrix to when computing stationary probabilities
MATRIX_EXPONENTIAL = 9999999

GAS_CONSTANT = 0.0019872 # in kcal/mol*K
DEFAULT_IMAGE_DIR = "images_and_plots/"

N_ALPHA_BETA_DIR = "N_alpha_beta/"
T_ALPHA_DIR = "T_alpha/"
K_ALPHA_BETA_DIR = "k_alpha_beta/"
PI_ALPHA_DIR = "pi_alpha/"
N_IJ_ALPHA_DIR = "N_ij_alpha/"
R_I_ALPHA_DIR = "R_i_alpha/"
N_IJ_DIR = "N_ij/"
R_I_DIR = "R_i/"

class MissingStatisticsError(Exception):
    """Catch a very specific type of error in analysis stage."""
    pass

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

def combine_dest_states(Q, src_indices, dest_indices):
    """
    When multiple milestones represent a given destination state, then
    adjust the rate matrix to consider all states, with fluxes into them,
    as a single state (single row and column in the matrix).
    """
    newQ = deepcopy(Q)
    new_src_indices = deepcopy(src_indices)
    new_dest_indices = deepcopy(dest_indices)
    n = newQ.shape[0]
    # insert a new row and column at the very end
    newQ = np.insert(newQ, n, np.zeros(n), 0)
    newQ = np.insert(newQ, n, np.zeros(n+1), 1)
    for i in range(len(dest_indices)):
        dest_index = new_dest_indices[i]
        dest_row = np.delete(newQ[dest_index,:], dest_index, 0)
        dest_column = np.delete(newQ[:,dest_index], dest_index, 0)
        newQ = np.delete(newQ, dest_index, 0)
        newQ = np.delete(newQ, dest_index, 1)
        newQ[-1,:] += dest_row
        newQ[:,-1] += dest_column
        
        for j in range(len(new_src_indices)):
            old_src_index = new_src_indices[j]
            if old_src_index > dest_index:
                new_src_indices[j] = old_src_index - 1
            else:
                new_src_indices[j] = old_src_index
        
        for j in range(len(new_dest_indices)):
            old_dest_index = new_dest_indices[j]
            if old_dest_index > dest_index:
                new_dest_indices[j] = old_dest_index - 1
            else:
                new_dest_indices[j] = old_dest_index
        
    n = newQ.shape[0]
    for i in range(n):
        newQ[i,i] = 0.0
        newQ[i,i] = -np.sum(newQ[i,:])
    
    return newQ, new_src_indices, n-1

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

def make_image_directory(model, image_directory):
    """
    Create the directory and sub-directories for images to be written.
    """
    if image_directory is None:
        image_directory = os.path.join(model.anchor_rootdir, DEFAULT_IMAGE_DIR)
    
    if image_directory != "" and not os.path.exists(image_directory):
        os.mkdir(image_directory)
    
    if not os.path.exists(os.path.join(image_directory, N_ALPHA_BETA_DIR)):
        os.mkdir(os.path.join(image_directory, N_ALPHA_BETA_DIR))
        
    if not os.path.exists(os.path.join(image_directory, T_ALPHA_DIR)):
        os.mkdir(os.path.join(image_directory, T_ALPHA_DIR))
        
    if not os.path.exists(os.path.join(image_directory, K_ALPHA_BETA_DIR)):
        os.mkdir(os.path.join(image_directory, K_ALPHA_BETA_DIR))
        
    if not os.path.exists(os.path.join(image_directory, PI_ALPHA_DIR)):
        os.mkdir(os.path.join(image_directory, PI_ALPHA_DIR))
        
    if not os.path.exists(os.path.join(image_directory, N_IJ_ALPHA_DIR)):
        os.mkdir(os.path.join(image_directory, N_IJ_ALPHA_DIR))
        
    if not os.path.exists(os.path.join(image_directory, R_I_ALPHA_DIR)):
        os.mkdir(os.path.join(image_directory, R_I_ALPHA_DIR))
    
    if not os.path.exists(os.path.join(image_directory, N_IJ_DIR)):
        os.mkdir(os.path.join(image_directory, N_IJ_DIR))
        
    if not os.path.exists(os.path.join(image_directory, R_I_DIR)):
        os.mkdir(os.path.join(image_directory, R_I_DIR))
    
    return image_directory

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
    k_b = 0.0
    #total_n_trajectories = 0
    
    for results_filename in results_filename_list:
        cmd = "%s < %s" % (compute_rate_constant_program, 
                           results_filename)
        output_string = subprocess.check_output(cmd, shell=True)
        root = ET.fromstring(output_string)
        
        # first, obtain counts directly from the results file
        results_root = ET.parse(results_filename)
        reactions = results_root.find("reactions")
        n_trajectories = int(reactions.find("n_trajectories").text.strip())
        #total_n_trajectories += n_trajectories
        stuck = int(reactions.find("stuck").text.strip())
        escaped = int(reactions.find("escaped").text.strip())
        for completed in reactions.iter("completed"):
            name = completed.find("name").text.strip().split("_")
            dest = int(name[1])
            n = int(completed.find("n").text.strip())
            transition_counts[dest] += n
        transition_counts["escaped"] += escaped
        transition_counts["stuck"] += stuck
        transition_counts["total"] += n_trajectories
        
        # now obtain probabilities and rates
        for rate in root.iter("rate"):
            rate_constant = rate.find("rate_constant")
            k_on_tag = rate_constant.find("mean")
            k_on_tag_high = rate_constant.find("high")
            assert k_on_tag is not None
            assert k_on_tag_high is not None
            k_on = float(k_on_tag.text)
            k_on_error = 0.5*(float(k_on_tag_high.text) - k_on)
            reaction_probability = rate.find("reaction_probability")
            beta_tag = reaction_probability.find("mean")
            assert beta_tag is not None
            beta = float(beta_tag.text)
            if beta != 0.0:
                k_b = k_on / beta
            
            if sample_error_from_normal:
                k_on = np.random.normal(loc=k_on, scale=k_on_error)
    
    for key in transition_counts:
        reaction_probabilities[key] \
            = transition_counts[key] / transition_counts["total"]
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
        self.K = None
        self.K_hat = None
        self.p_i = None
        self.free_energy_profile = None
        self.MFPTs = None
        self.k_off = None
        self.k_ons = {}
        self.b_surface_k_ons_src = None
        self.b_surface_k_on_errors_src = None
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
                self.bd_transition_probabilities["b_surface"] \
                    = reaction_probabilities
                self.b_surface_k_ons_src = k_ons_src
                self.b_surface_b_surface_k_on_errors_src = k_on_errors_src
            else:
                raise Exception("No valid BD program settings provided.")
        
        if len(self.model.k_on_info.bd_milestones) > 0 \
                and len(output_file_list) > 0:
            for bd_milestone in self.model.k_on_info.bd_milestones:
                transition_counts_bd_milestone = defaultdict(int)
                transition_probabilities_bd_milestone = defaultdict(float)
                inner_milestone_index = bd_milestone.inner_milestone.index
                outer_milestone_index = bd_milestone.outer_milestone.index
                assert inner_milestone_index in transition_counts
                assert outer_milestone_index in transition_counts
                transition_counts_bd_milestone[inner_milestone_index] \
                    = transition_counts[inner_milestone_index]
                transition_counts_bd_milestone["escaped"] \
                    = transition_counts[outer_milestone_index] \
                    - transition_counts[inner_milestone_index]
                transition_counts_bd_milestone["total"] \
                    = transition_counts[outer_milestone_index]
                if transition_counts_bd_milestone["escaped"] == 0:
                    transition_probabilities_bd_milestone[
                        inner_milestone_index] = 1.0
                    transition_probabilities_bd_milestone["escaped"] \
                        = 0
                else:
                    transition_probabilities_bd_milestone[
                        inner_milestone_index] \
                        = transition_counts_bd_milestone[
                            inner_milestone_index] \
                        / transition_counts_bd_milestone["total"]
                    transition_probabilities_bd_milestone["escaped"] \
                        = transition_counts_bd_milestone["escaped"] \
                        / transition_counts_bd_milestone["total"]
                
                self.bd_transition_counts[bd_milestone.index] \
                    = transition_counts_bd_milestone
                self.bd_transition_probabilities[bd_milestone.index] \
                        = transition_probabilities_bd_milestone
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
            row_sum = np.sum(self.Q[i], dtype=np.longdouble)
            if row_sum == 0:
                new_row_sum = np.longdouble(0.0)
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
        
    def calculate_kinetics(self):
        """
        Once the rate matrix Q is computed, determine the timescales 
        and probabilities of transfers between different states. Fill
        out all kinetics quantities.
        """
        end_milestones = {}
        end_states = defaultdict(list)
        end_milestone_list = []
        bulk_milestones = []
        MFPTs = defaultdict(float)
        MFPTs_anchors_to_bulk_milestones = defaultdict(float)
        k_off = 0.0
        k_ons = {}
        
        for anchor in self.model.anchors:
            if anchor.bulkstate:
                for milestone_id in anchor.get_ids():
                    if self.model.get_type() == "elber":
                        if anchor.alias_from_id(milestone_id) == 1: 
                            # TODO: hacky
                            continue
                    bulk_milestones.append(milestone_id)
                    #bulk_milestone = milestone_id
                    
        for anchor in self.model.anchors:
            if anchor.endstate:
                for milestone_id in anchor.get_ids():
                    if self.model.get_type() == "elber":
                        if anchor.alias_from_id(milestone_id) == 3: 
                            # TODO: hacky
                            continue
                        
                    if milestone_id not in bulk_milestones:
                        end_milestones[milestone_id] = anchor.name
                        end_states[anchor.name].append(milestone_id)
                        end_milestone_list.append(milestone_id)
        
        assert len(end_milestones) > 0, "No end (bound or otherwise) states "\
            "defined in this model. Kinetics calculations will not work."
        
        if np.any(self.Q.sum(axis=1) > 1.E-10):
            problem_milestone = np.argmin(self.Q.T.sum(axis=1))
            error_msg = """The rate matrix Q has a numerically overflowed row 
            sum placed into the main diagonal of milestone {}. This can be
            caused by large discrepancies of rates between adjacent milestones.
            You can try additional sampling of whichever anchor(s) have this 
            milestone. If possible, you may consider adding more milestones in
            this area, or omitting calculations in this portion of the CV(s). 
            If this problem persists, please contact the developers.""".format(
                problem_milestone)
            raise Exception(error_msg)
        
        prob_marginalized_per_anchor = defaultdict(float)
        prob_marginalized_all_anchors = 0.0
        for end_milestone_src in end_milestones:
            src_anchor_name = end_milestones[end_milestone_src]
            prob_marginalized_per_anchor[src_anchor_name] \
                    += self.p_i[end_milestone_src]
            prob_marginalized_all_anchors += self.p_i[end_milestone_src]
            
        if len(bulk_milestones) > 0:
            bulk_combined_Q, new_end_milestone_list, bulk_index \
                = combine_dest_states(self.Q, end_milestone_list, bulk_milestones)
            B, tau = PyGT.tools.load_CTMC(bulk_combined_Q.T)
            MFPTs_end_milestones_to_bulk = {}
            for i, end_milestone_src in enumerate(end_milestones):
                new_end_milestone_index = new_end_milestone_list[i]
                mfpt_ij, mfpt_ji = PyGT.mfpt.compute_MFPT(
                    new_end_milestone_index, bulk_index, B, tau)
                mfpt = mfpt_ji
                MFPTs_end_milestones_to_bulk[end_milestone_src] = mfpt
            
            MFPT_to_bulk = 0.0
            for end_milestone_src in end_milestones:
                if end_milestone_src in bulk_milestones:
                    continue
                src_anchor_name = end_milestones[end_milestone_src]
                mftp_milestone_to_bulk = MFPTs_end_milestones_to_bulk[end_milestone_src]
                prob_milestone_in_anchor = self.p_i[end_milestone_src] \
                    / prob_marginalized_per_anchor[src_anchor_name]
                prob_milestone_all_anchors = self.p_i[end_milestone_src] \
                    / prob_marginalized_all_anchors
                MFPTs[(src_anchor_name, "bulk")] += prob_milestone_in_anchor \
                    * mftp_milestone_to_bulk
                MFPT_to_bulk += prob_milestone_all_anchors * mftp_milestone_to_bulk
            
            k_off = 1.0e12 / MFPT_to_bulk
            self.k_off = k_off
        
        # Next, compute the MFPTs between different states
        
        milestone_src_to_anchor_dest_MFPTs = {}
        for dest_state_name, dest_milestones in end_states.items():
            for dest_milestone in dest_milestones:
                assert dest_milestone not in bulk_milestones
            
            for src_state_name, src_milestone_list in end_states.items():
                if src_state_name == dest_state_name:
                    continue
                dest_combined_Q, new_src_milestone_list, dest_index \
                    = combine_dest_states(self.Q, src_milestone_list, dest_milestones)
                B, tau = PyGT.tools.load_CTMC(dest_combined_Q.T)
                for i, src_milestone_index in enumerate(src_milestone_list):
                    new_src_milestone_index = new_src_milestone_list[i]
                    mfpt_ij, mfpt_ji = PyGT.mfpt.compute_MFPT(
                        new_src_milestone_index, dest_index, B, tau)
                    mfpt = mfpt_ji
                    milestone_src_to_anchor_dest_MFPTs[
                        (src_milestone_index, dest_state_name)] = mfpt
            
            for src_state_name, src_milestone_list in end_states.items():
                if src_state_name == dest_state_name:
                    continue
                for src_milestone_index in src_milestone_list:
                    mftp_milestone_to_dest \
                        = milestone_src_to_anchor_dest_MFPTs[
                            (src_milestone_index, dest_state_name)]
                    prob_milestone_in_anchor = self.p_i[src_milestone_index] \
                        / prob_marginalized_per_anchor[src_state_name]
                    MFPTs[(src_state_name, dest_state_name)] \
                        += prob_milestone_in_anchor * mftp_milestone_to_dest
            
        if self.model.using_bd():
            K_hat = self.K[:,:]
            for end_milestone in end_milestones:
                K_hat[end_milestone, :] = 0.0
                K_hat[end_milestone, end_milestone] = 1.0
            
            n = K_hat.shape[0]
            source_vec = np.zeros((n,1))
            
            if len(self.bd_transition_counts) > 0:
                if len(bulk_milestones) > 0:
                    for bd_milestone in self.model.k_on_info.bd_milestones:
                        source_index = bd_milestone.outer_milestone.index
                        if self.b_surface_k_ons_src is None:
                            raise Exception("Missing b-surface simulations.")
                            
                        source_vec[source_index] \
                            = self.b_surface_k_ons_src[source_index]
                        
                        if bd_milestone.index \
                                not in self.bd_transition_probabilities:
                            break
                        
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
        
        self.MFPTs = MFPTs
        return