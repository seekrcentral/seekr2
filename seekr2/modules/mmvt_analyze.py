"""
Routines for the analysis stage of an MMVT calculation
"""
import os
import re
from collections import defaultdict
from copy import deepcopy

import numpy as np
import scipy.linalg as la

import seekr2.modules.common_analyze as common_analyze
import seekr2.modules.markov_chain_monte_carlo as markov_chain_monte_carlo

FLUX_MATRIX_K_EXPONENT = 1000
LOW_N_IJ = 1e-10
HIGH_N_IJ = 1e3
DEFAULT_R_I = 1.0
DEFAULT_Q_IJ = 1.0
LOW_Q_IJ = 1e-15
HIGH_FLUX = 1e10
LOW_PROBABILITY = 1e-25
from seekr2.modules.common_analyze import GAS_CONSTANT

def flux_matrix_to_K(M):
    """Given a flux matrix M, compute probability transition matrix K."""
    K = np.zeros((M.shape[0]-1, M.shape[0]-1))
    for i in range(M.shape[0]-1):
        for j in range(M.shape[0]-1):
            if i == j:
                K[i,j] = 0.0
            else:
                K[i,j] = -M[i,j] / M[i,i]
                
            assert K[i,j] >= 0.0, "Negative values for K matrix not allowed."
    
    return K

def make_new_Nij_alpha(mmvt_Qij_alpha, mmvt_Ri_alpha):
    """
    After a new matrix Q has been sampled, the "data" N_ij must be
    adjusted to be consistent with with newly sampled matrix Q and 
    R_i.
    """
    n_milestones = mmvt_Qij_alpha.shape[0]
    new_mmvt_Nij_alpha = np.zeros((n_milestones, n_milestones))
    for i in range(n_milestones):
        for j in range(n_milestones):
            if i == j: continue
            new_mmvt_Nij_alpha[i,j] = mmvt_Qij_alpha[i,j] * mmvt_Ri_alpha[i,0]
            
    return new_mmvt_Nij_alpha

def openmm_read_output_file_list(output_file_list, min_time=None, max_time=None, 
                                 existing_lines=None, skip_restart_check=False):
    """
    Read the output files produced by the plugin (backend) of 
    SEEKR2 and extract transition statistics and times
    """
    
    MAX_ITER = 1000000000
    NEW_SWARM = "NEW_SWARM"
    swarm_index = 0
    no_checkpoints = True
    if existing_lines is None:
        existing_lines = []
        files_lines = []
        start_times = []
        checkpoint_times = []
        for i, output_file_name in enumerate(output_file_list):
            start_time = None
            checkpoint_time = None
            file_lines = []
            outfile_basename = os.path.basename(output_file_name)
            if re.search("swarm", outfile_basename):
                old_swarm_index = swarm_index
                swarm_index = int(outfile_basename.split(".")[1].split("_")[1])
                if swarm_index != old_swarm_index:
                    file_lines.append(NEW_SWARM)
                
            with open(output_file_name, "r") as output_file:
                for counter, line in enumerate(output_file):
                    if line.startswith("#") or len(line.strip()) == 0:
                        continue
                    line_list = line.strip().split(",")
                    if line.startswith("CHECKPOINT"):
                        no_checkpoints = False
                        if len(line_list) != 2:
                            continue
                        if re.match(r"^-?\d+\.\d{3,20}$", line_list[1]):
                            if checkpoint_time is None:
                                checkpoint_time = float(line_list[1])
                        else:
                            continue
                        continue
                    if len(line_list) != 3:
                        continue
                    dest_boundary = int(line_list[0])
                    bounce_index = int(line_list[1])
                    if re.match(r"^-?\d+\.\d{3,20}$", line_list[2]):
                        dest_time = float(line_list[2]) 
                    else:
                        continue
                    
                    file_lines.append(line_list)
                    if start_time is None:
                        start_time = dest_time
                    
                
                if len(file_lines) == 0:
                    continue
                
                if start_time is None:
                    start_times.append(0.0)
                else:
                    start_times.append(start_time)
                
                checkpoint_times.append(checkpoint_time)
            
            files_lines.append(file_lines)
                    
        lines = []
        
        checkpoint_times.append(1e99)
        
        for i, file_lines in enumerate(files_lines):
            if not skip_restart_check and not len(file_lines) == 0:
                # make sure that we delete lines that occurred after a
                # checkpoint backup
                if i < len(files_lines)-1:
                    # then it's not the last file
                    if files_lines[i+1][0] == NEW_SWARM:
                        # If the current file is the beginning of a swarm
                        next_start_time = 1e99
                    else:
                        
                        if no_checkpoints:
                            next_start_time = start_times[i+1]
                        else:
                            if checkpoint_times[i] is not None:
                                next_start_time = checkpoint_times[i+1]
                            else:
                                next_start_time = start_times[i+1]
                else:
                    next_start_time = 1e99
                counter = 0
                
                if file_lines[-1]==NEW_SWARM:
                    lines.append(NEW_SWARM)
                    file_lines.pop()                    
                
                if len(file_lines) == 0:
                    continue
                
                while float(file_lines[-1][2]) >= next_start_time:
                    file_lines.pop()
                    if counter > MAX_ITER or len(file_lines)==0:
                        break
                    if file_lines[-1]==NEW_SWARM:
                        file_lines.pop()
                        lines.append(NEW_SWARM)
                        break
                        
                    counter += 1
                                    
            lines += file_lines
        
    else:
        lines = existing_lines
        
    N_i_j_alpha = defaultdict(int)
    R_i_alpha_list = defaultdict(list)
    N_alpha_beta = defaultdict(int)
    T_alpha_list = []
    
    last_bounce_time = -1.0
    src_boundary = None
    src_time = None
    for counter, line in enumerate(lines):
        if line == NEW_SWARM:
            last_bounce_time = -1.0
            src_boundary = None
            src_time = None
            continue 
        
        dest_boundary = int(line[0])
        bounce_index = int(line[1])
        dest_time = float(line[2]) 
        # This is used to cut out early transitions for analysis
        if min_time is not None:
            if dest_time < min_time:
                continue
            
        # This is used in convergence
        if max_time is not None:
            if dest_time > max_time:
                break
            
        if src_boundary is None:
            src_boundary = dest_boundary
            src_time = dest_time
            last_bounce_time = dest_time
            continue
        
        directory = None
        if len(output_file_list) > 0:
            directory = os.path.dirname(output_file_list[0])
                
        if src_boundary != dest_boundary:
            
            time_diff = dest_time - src_time
            if not skip_restart_check:
                assert time_diff >= 0.0, "incubation times cannot be "\
                    "negative. Has an output file been concatenated "\
                    "incorrectly? directory: %s, line number: %d" % (
                    directory, counter)
            N_i_j_alpha[(src_boundary, dest_boundary)] += 1
            R_i_alpha_list[src_boundary].append(time_diff)
            src_boundary = dest_boundary
            src_time = dest_time
            
        # dest_boundary is beta
        N_alpha_beta[dest_boundary] += 1
        assert dest_time - last_bounce_time >= 0.0, "times between bounces "\
            "cannot be negative. bounce index: "\
            "{}, Dest_time: {}, last_bounce_time: {}, directory: {}".format(
                bounce_index, dest_time, last_bounce_time, directory)
        T_alpha_list.append(dest_time - last_bounce_time)
        last_bounce_time = dest_time
    
    R_i_alpha_total = defaultdict(float)
    
    for key in R_i_alpha_list:
        assert key is not None
        R_i_alpha_total[key] = np.sum(R_i_alpha_list[key])
        
    T_alpha_total = np.sum(T_alpha_list)
    
    if len(R_i_alpha_list) == 0 and src_boundary is not None:
        R_i_alpha_total[src_boundary] = T_alpha_total
        R_i_alpha_list[src_boundary] = T_alpha_list
        
    return N_i_j_alpha, R_i_alpha_total, N_alpha_beta, T_alpha_total, lines

def openmm_read_statistics_file(statistics_file_name):
    """
    Read the statistics file produced by the plugin (backend) of
    SEEKR2 and extract transition statistics and times.

    Parameters
    ----------
    statistics_file_name : string
        The name of the openmmvt output file to read for extracting
        transition statistics

    Returns
    -------
    N_i_j_alpha_dict : dict
        dictionary of counts of transitions within cell alpha 
        between boundaries i and j.
        
    R_i_alpha_dict : dict
        dictionary of transition incubation times spent before 
        touching surface i
    
    N_alpha_beta_dict : dict
        dictionary array of counts that a transition between 
        cell alpha and beta
        
    T_alpha : float
        Total time spent within this cell.
    """
    raise Exception("Statistics needed for error analysis are not yet "\
                    "extracted in this function.")
    N_i_j_alpha_dict = defaultdict(int)
    R_i_alpha_dict = defaultdict(float)
    N_alpha_beta_dict = defaultdict(int)
    T_alpha = 0.0
    with open(statistics_file_name, 'r') as output_file:
        for line in output_file:
            linesplit = line.strip().split(' ')
            if line.startswith("N_alpha_"):
                keysplit = linesplit[0].split('_')
                dest_boundary = int(keysplit[2].strip(':'))
                value = int(linesplit[1])
                N_alpha_beta_dict[dest_boundary] += value
                
            elif line.startswith("N_"):
                keysplit = linesplit[0].split('_')
                src_boundary = int(keysplit[1])
                dest_boundary = int(keysplit[2])
                if src_boundary == dest_boundary:
                    continue
                value = int(linesplit[1])
                N_i_j_alpha_dict[(src_boundary, dest_boundary)] += value
                
            elif line.startswith("R_"):
                keysplit = linesplit[0].split('_')
                src_boundary = int(keysplit[1])
                value = float(linesplit[1])
                R_i_alpha_dict[src_boundary] += value
                
            elif line.startswith("T_alpha"):
                value = float(linesplit[1])
                T_alpha = value
                
            else:
                raise Exception("Unexpected line: %s in file %s" % \
                                (line, statistics_file_name))
            
    return N_i_j_alpha_dict, R_i_alpha_dict, N_alpha_beta_dict, T_alpha

def namd_read_output_file_list(output_file_list, anchor, timestep, 
                               min_time=None, max_time=None, existing_lines=None, 
                               skip_restart_check=False):
    """
    Read the output files produced by the plugin (backend) of 
    NAMD and extract transition statistics and times

    Parameters
    ----------
    output_file_list : list
        The names of the NAMD output files to read for 
        extracting transition statistics

    Returns
    -------
    N_i_j_alpha_dict : dict
        dictionary of counts of transitions within cell alpha 
        between boundaries i and j.
        
    R_i_alpha_dict : dict
        dictionary of transition incubation times spent before 
        touching surface i
    
    N_alpha_beta_dict : dict
        dictionary array of counts that a transition between cell 
        alpha and beta
        
    T_alpha : float
        Total time spent within this cell.
    """
    MAX_ITER = 1000000000
    no_checkpoints = True
    if existing_lines is None:
        existing_lines = []
        files_lines = []
        start_times = []
        checkpoint_times = []
        for i, output_file_name in enumerate(output_file_list):
            start_time = None
            checkpoint_time = None
            file_lines = []            
            with open(output_file_name, "r") as output_file:
                for counter, line in enumerate(output_file):
                    if line.startswith("CHECKPOINT"):
                        no_checkpoints = False
                        checkpoint_time = float(line.strip().split(" ")[3])
                        continue
                    elif not line.startswith("SEEKR: Cell") \
                            or len(line.strip()) == 0:
                        continue
                    #elif line.startswith("SEEKR: Cell"):
                    line = line.strip().split(" ")
                    file_lines.append(line)
                    current_stepnum = int(line[8].strip(","))
                    if start_time is None:
                        start_time = current_stepnum
                    
                    # Not used anymore
                    #elif line.startswith("SEEKR: Milestone"):
                    #    line = line.strip().split(" ")
                    #    file_lines.append(line)
                    #    current_stepnum = int(line[10].strip(","))
                    #    if start_time is None:
                    #        start_time = current_stepnum
                    
                if len(file_lines) == 0:
                    continue
                
                if start_time is None:
                    start_times.append(0.0)
                else:
                    start_times.append(start_time)
                
                checkpoint_times.append(checkpoint_time)
            
            files_lines.append(file_lines)
            
        lines = []
        for i, file_lines in enumerate(files_lines):
            if not skip_restart_check and not len(file_lines) == 0:
                # make sure that we delete lines that occurred after a
                # checkpoint backup
                if i < len(files_lines)-1:
                    # then it's not the last file
                    if no_checkpoints:
                        next_start_time = start_times[i+1]
                    else:
                        if checkpoint_times[i] is not None:
                            next_start_time = checkpoint_times[i]
                else:
                    next_start_time = 1e99
                counter = 0
                
                if len(file_lines) == 0:
                    continue
                
                while float(file_lines[-1][8]) >= next_start_time:
                    file_lines.pop()
                    if counter > MAX_ITER or len(file_lines)==0:
                        break
                    
                    counter += 1
                
            lines += file_lines
        
    else:
        lines = existing_lines
    
    N_i_j_alpha = defaultdict(int)
    R_i_alpha_list = defaultdict(list)
    N_alpha_beta = defaultdict(int)
    T_alpha_list = []
    
    last_bounce_time = -1.0
    src_boundary = None
    src_time = None
    for counter, line in enumerate(lines):
        next_anchor_raw = int(line[6].strip(","))
        current_stepnum = int(line[8].strip(","))
                
        dest_time = current_stepnum * timestep
        dest_boundary = None
        for milestone in anchor.milestones:
            if milestone.neighbor_anchor_index == next_anchor_raw:
                dest_boundary = milestone.alias_index
                
        # This is used to cut out early transitions for analysis
        if min_time is not None:
            if dest_time < min_time:
                continue
            
        # This is used in convergence
        if max_time is not None:
            if dest_time > max_time:
                break
            
        if src_boundary is None:
            src_boundary = dest_boundary
            src_time = dest_time
            last_bounce_time = dest_time
            continue
        
        if src_boundary != dest_boundary:
            time_diff = dest_time - src_time
            if not skip_restart_check:
                assert time_diff >= 0.0, "incubation times cannot be "\
                    "negative. Has an output file been concatenated "\
                    "incorrectly? file name(s): %s, line number: %d" % (
                    ",".join(output_file_list), counter)
            N_i_j_alpha[(src_boundary, dest_boundary)] += 1
            R_i_alpha_list[src_boundary].append(time_diff)
            src_boundary = dest_boundary
            src_time = dest_time
            
        # dest_boundary is beta
        N_alpha_beta[dest_boundary] += 1
        assert dest_time - last_bounce_time >= 0.0, "times between bounces "\
            "cannot be negative. bounce index: "\
            "{}, Dest_time: {}, last_bounce_time: {}".format(
                counter, dest_time, last_bounce_time)
        T_alpha_list.append(dest_time - last_bounce_time)
        last_bounce_time = dest_time
    
    R_i_alpha_total = defaultdict(float)
    
    for key in R_i_alpha_list:
        assert key is not None
        R_i_alpha_total[key] = np.sum(R_i_alpha_list[key])
        
    T_alpha_total = np.sum(T_alpha_list)
    
    if len(R_i_alpha_list) == 0 and src_boundary is not None:
        R_i_alpha_total[src_boundary] = T_alpha_total
        R_i_alpha_list[src_boundary] = T_alpha_list
        
    return N_i_j_alpha, R_i_alpha_total, N_alpha_beta, T_alpha_total, lines

class MMVT_anchor_statistics():
    """
    Contain all the transition statistics for an anchor within
     an MMVT calculation.
    
    Attributes:
    -----------
    N_i_j_alpha_dict : dict
        Represents the N_i_j_alpha quantities in the MMVT calculation.
        The attribute represents the number of times a transition was
        observed within cell alpha going from surface i to surface j.
        The keys of the dict are a tuple (i,j) composed of two ints,
        and the value is an int count of transitions.
        
    R_i_alpha_dict : dict
        Represents the R_i_alpha quantities in the MMVT calculation.
        The attribute represents the total time a transition spent
        after bouncing off of surface i before touching another 
        surface, and was observed within cell alpha. The keys of the 
        dict are an int i, and the value is a float representing time.
    
    R_i_average : dict
        The averages of the values in R_i_list.
        
    R_i_std_dev : dict
        The standard deviations of the values in R_i_list.
        
    R_i_total : dict
        The sum total of the values in R_i_list.
    
    N_alpha_beta_dict : dict
        Represents the N_alpha_beta quantities in the MMVT calculation.
        The attribute represents the number of times the system bounced
        off of the surface between the current cell alpha and another
        cell beta.
    
    T_alpha_list : list
        A list of times between bounces within anchor alpha, regardless
        of the identity of the source or destination milestones.
        
    T_alpha_average : float
        The average of the values in T_alpha_list.
        
    T_alpha_std_dev : float
        The standard deviation of the values in in T_alpha_list.
    
    T_alpha_total : float
        Represents the quantity T_alpha in the MMVT calculation. The
        attribute represents the total time spent within the simulation
        not counting the time at the beginning spent before any bounces
        occurred, and also not counting the time spent at the end after
        the last bounce occurred. This is equivalent to the sum total
        of the values in T_alpha_list.
    
    N_alpha_beta : dict
        A dictionary whose keys are the indices of adjacent anchors and
        whose values are the counts of bounces against those anchors.
        
    k_alpha_beta : dict
        A dictionary whose keys are the indices of adjacent anchors and
        whose values are the rate of transitions between those anchors.
    
    alpha : int
        The index of the anchor whose transition statistics are 
        represented by this object.
    
    existing_lines : list or None
        If the output files for this anchor has already been read, then
        those lines will be stored here to save on I/O.
    """
    
    def __init__(self, alpha):
        self.N_i_j_alpha = None
        self.R_i_alpha_total = None
        self.T_alpha_total = None
        self.N_alpha_beta = None
        self.k_alpha_beta = {}
        self.alpha = alpha
        self.existing_lines = None
        return
    
    def read_output_file_list(self, engine, output_file_list, min_time, 
                              max_time, anchor, timestep):
        """
        Depending on the engine and other settings, read the SEEKR2 
        output files to fill out transition statistics.
        """
        if engine == "openmm":
            self.N_i_j_alpha, self.R_i_alpha_total, self.N_alpha_beta, \
            self.T_alpha_total, self.existing_lines \
                = openmm_read_output_file_list(
                    output_file_list, min_time, max_time, self.existing_lines)
        elif engine == "namd":
            self.N_i_j_alpha, self.R_i_alpha_total, self.N_alpha_beta, \
            self.T_alpha_total, self.existing_lines \
                = namd_read_output_file_list(
                    output_file_list, anchor, timestep, min_time, max_time, 
                    self.existing_lines)
        elif engine == "browndye2":
            pass
        else:
            raise Exception("Engine not allowed: {}. Must be 'openmm', "\
                            "'namd', or 'browndye'.")
        for key in self.N_alpha_beta:
            self.k_alpha_beta[key] = self.N_alpha_beta[key] / self.T_alpha_total
            
        assert self.T_alpha_total >= 0.0
        return
    
    def print_stats(self):
        """ Print statistics for easy viewing."""
        print("self.alpha:", self.alpha)
        print("self.N_i_j_alpha:", self.N_i_j_alpha)
        print("self.R_i_alpha_total:", self.R_i_alpha_total)
        print("self.T_alpha_total:", self.T_alpha_total)
        print("self.N_alpha_beta:", self.N_alpha_beta)
        print("self.k_alpha_beta:", self.k_alpha_beta)

class MMVT_data_sample(common_analyze.Data_sample):
    """
    Represent a data sample specific to an MMVT calculation.
    
    Attributes:
    -----------
    model : Model()
        The SEEKR2 model the calculation is being performed on.
        
    k_alpha_beta : defaultdict
        A dictionary whose keys are 2-tuples of source/destination
        anchor indices and whose values are the rates between those
        anchors.
    
    N_i_j_alpha : list
        A list where each entry corresponds to anchor index alpha, and
        are defaultdicts whose keys are 2-tuples of source/destination
        milestone indices, and whose values are counts of transitions
        between those milestones.
    
    R_i_alpha : list
        A list where each entry corresponds to anchor index alpha, and
        are defaultdicts whose keys are integers of source
        milestone indices, and whose values are incubation times of
        transitions between milestones.
        
    T_alpha : list
        A list of total times spent in each anchor.
    
    pi_alpha : numpy.array()
        For each anchor alpha, contains the stationary distribution
        within that anchor.
    
    N_ij : dict
        A dictionary with keys of 2-tupes of source and destination
        milestone, and whose values are the counts of transitions.
        
    R_i : dict
        A dictionary with keys of source milestone indices, and whose
        values are incubation times spent at those milestones.
        
    T : float
        Total time: nothing more than a normalization factor for the 
        calculations.
    
    Q : numpy.array()
        A 2x2 rate matrix for a milestoning calculation. No sink states.
        
    K : numpy.array()
        A Markov transition probability matrix constructed from Q. No
        sink states.
        
    K_hat : numpy.array()
        A Markov transition probability matrix constructed from Q. 
        Sink state(s) included.
        
    p_i : numpy.array()
        Stationary probabilities along the milestones.
        
    free_energy_profile : numpy.array()
        The free energy profile as computed by the stationary 
        probabilities.
        
    MFPTs : dict
        A dictionary whose keys are 2-tuples of milestone end states, 
        and whose values are mean first passage times (MFPTs) between
        those states.
        
    k_off : float
        The overall k-off (in s^-1) from a stationary distribution to 
        the bulk state.
        
    k_ons : dict
        A set of k-ons (in M^-1 * s^-1) from the bulk state to each of
        the end states.
        
    bd_transition_counts : dict
        A dictionary whose keys are bd_milestones, and whose values are
        more dictionaries, whose keys are states the BD simulations 
        can end at, and whose values are the counts of encountering
        those states.
    """
    def __init__(self, model, N_alpha_beta=None, k_alpha_beta=None, 
                 N_i_j_alpha=None, R_i_alpha=None, T_alpha=None):
        self.model = model
        self.N_alpha_beta = N_alpha_beta
        self.k_alpha_beta = k_alpha_beta
        self.k_alpha = []
        self.N_i_j_alpha = N_i_j_alpha
        self.N_alpha = []
        self.R_i_alpha = R_i_alpha
        self.T_alpha = T_alpha
        self.pi_alpha = None
        self.N_ij = {}
        self.R_i = {}
        self.N_ij_unmodified = {}
        self.R_i_unmodified = {}
        self.T = None
        self.Q = None
        self.K = None
        self.K_hat = None
        self.p_i = None
        self.p_i_hat = None
        self.free_energy_profile = None
        self.free_energy_anchors = None
        self.MFPTs = {}
        self.k_off = None
        self.k_ons = {}
        self.b_surface_k_ons_src = None
        self.b_surface_k_on_errors_src = None
        self.bd_transition_counts = {}
        self.bd_transition_probabilities = {}
        
        if self.N_alpha_beta is None or self.k_alpha_beta is None \
                or self.N_i_j_alpha is None or self.R_i_alpha is None \
                or self.T_alpha is None:
            return
        
        # Fill out N_alpha
        for alpha, anchor in enumerate(model.anchors):
            if anchor.bulkstate:
                continue
            self.N_alpha.append(0)
            for key in self.N_i_j_alpha[alpha]:
                self.N_alpha[alpha] += self.N_i_j_alpha[alpha][key]
            
            k_alpha = 0.0
            for beta, anchor2 in enumerate(model.anchors):
                key = (alpha, beta)
                if key in self.k_alpha_beta:
                    k_alpha += self.k_alpha_beta[key]
            self.k_alpha.append(k_alpha)
            
        return
    
    def calculate_pi_alpha(self):
        """
        Computes the equilibrium probability quantities "pi_alpha" used
        in MMVT theory. The value self.pi_alpha gets set by this 
        function.
        """
        if self.k_alpha_beta is None:
            raise Exception("Unable to call calculate_pi_alpha(): "\
                            "No statistics present in Data Sample.")
        # First, determine if there is a bulk anchor, and if so, use it as
        # the "pivot", if not, make our own pivot.
        # The "pivot" is the entry in the flux_matrix that 
        
        # TODO: replace with model.get_bulk_index() function
        bulk_indices = []
        found_bulk = False
        for alpha, anchor in enumerate(self.model.anchors):
            if anchor.bulkstate:
                bulk_indices.append(alpha)
                found_bulk = True
            else:
                assert not found_bulk, \
                "The bulk anchors must be the last ones in the model."
        
        if len(bulk_indices) == 0:
            # Then we have a model without a bulk anchor: we need to make our
            # own pivot
            pivot_index = self.model.num_anchors
            flux_matrix_dimension = self.model.num_anchors + 1
        else:
            # If bulk anchors exist, we can use them as the pivot
            pivot_index = bulk_indices[0]
            flux_matrix_dimension = self.model.num_anchors-len(bulk_indices)+1
        
        self.pi_alpha = np.zeros(flux_matrix_dimension)
        flux_matrix = np.zeros((flux_matrix_dimension, flux_matrix_dimension))
        column_sum = np.zeros(flux_matrix_dimension)
        flux_matrix[pivot_index, pivot_index] = 1.0
        for alpha, anchor1 in enumerate(self.model.anchors):
            flux_matrix[alpha, pivot_index] = 1.0
            if anchor1.bulkstate:
                break
            id_alias = anchor1.alias_from_neighbor_id(pivot_index)
            flux_matrix[pivot_index, alpha] = 0.0
            dead_end_anchor = False
            if len(anchor1.milestones) == 1:
                dead_end_anchor = True
            for beta, anchor2 in enumerate(self.model.anchors):
                if beta == pivot_index:
                    break
                if alpha == beta:
                    pass 
                else:
                    id_alias = anchor1.alias_from_neighbor_id(
                        anchor2.index)
                    if id_alias is None:
                        flux_matrix[alpha, beta] = 0.0
                    else:
                        if dead_end_anchor:
                            # This line was supposed to work for a 1D
                            # Smoluchowski system, but with a 3D
                            # spherical system, the 2.0 needs to be 1.0.
                            #flux_matrix[alpha, beta] = 2.0 *\
                            #     self.k_alpha_beta[(alpha, beta)]
                            flux_matrix[alpha, beta] = 1.0 *\
                                 self.k_alpha_beta[(alpha, beta)]
                        else:
                            flux_matrix[alpha, beta] = \
                                self.k_alpha_beta[(alpha, beta)]
                        column_sum[alpha] += flux_matrix[alpha, beta]
                
            flux_matrix[alpha, alpha] = -column_sum[alpha]
            
        flux_matrix[pivot_index, pivot_index-1] = HIGH_FLUX
        prob_equil = np.zeros((flux_matrix_dimension,1))
        prob_equil[pivot_index] = 1.0
        pi_alpha = la.solve(flux_matrix.T, prob_equil)
        self.pi_alpha = abs(pi_alpha)
        
        # refine pi_alpha
        pi_alpha_slice = np.zeros((flux_matrix_dimension-1, 1))
        for i in range(flux_matrix_dimension-1):
            pi_alpha_slice[i,0] = -self.pi_alpha[i,0] * flux_matrix[i,i]
                    
        K = flux_matrix_to_K(flux_matrix)
        K_inf = np.linalg.matrix_power(K, FLUX_MATRIX_K_EXPONENT)
        stationary_dist = K_inf.T @ pi_alpha_slice
        for i in range(flux_matrix_dimension-1):
            self.pi_alpha[i,0] = -stationary_dist[i,0] / flux_matrix[i,i]
        for alpha, anchor1 in enumerate(self.model.anchors):
            if anchor1.bulkstate:
                break
            # Set low probabilities to anchors with no transitions observed.
            if self.k_alpha[alpha] == 0:
                self.pi_alpha[alpha,0] = LOW_PROBABILITY
        self.pi_alpha[-1,0] = 0.0
        self.pi_alpha = self.pi_alpha / np.sum(self.pi_alpha)
        return
    
    def fill_out_data_quantities(self):
        """
        Compute quantities such as N_ij, R_i, and T for eventual 
        construction of rate matrix Q.
        """
        if self.N_alpha_beta is None or self.k_alpha_beta is None \
            or self.N_i_j_alpha is None or self.R_i_alpha is None \
            or self.T_alpha is None:
            raise Exception("Unable to call fill_out_data_quantities(): "\
                        "No statistics present in Data Sample.")
        
        self.N_ij = defaultdict(float)
        self.R_i = defaultdict(float)
        self.T = 0.0
        sum_pi_alpha_over_T_alpha = 0.0
        for alpha, anchor in enumerate(self.model.anchors):
            if anchor.bulkstate:
                break
            if self.T_alpha[alpha] == 0.0:
                continue
            this_anchor_pi_alpha = self.pi_alpha[alpha]
            sum_pi_alpha_over_T_alpha += this_anchor_pi_alpha \
                / self.T_alpha[alpha]
                
        self.T = float(sum_pi_alpha_over_T_alpha ** -1)
        assert self.T >= 0.0, "self.T should be positive: {}".format(self.T)
        for alpha, anchor in enumerate(self.model.anchors):
            if anchor.bulkstate:
                break
            this_anchor_pi_alpha = float(self.pi_alpha[alpha])
            T_alpha = self.T_alpha[alpha]
            assert T_alpha >= 0.0, \
                "T_alpha should be positive: {}".format(T_alpha)
            if T_alpha > 0.0:
                time_fraction = self.T / T_alpha
            else:
                time_fraction = 0.0
            N_i_j_alpha = self.N_i_j_alpha[alpha]
            for key in N_i_j_alpha:
                assert time_fraction >= 0.0, \
                    "time_fraction should be positive: {}".format(time_fraction)
                assert this_anchor_pi_alpha >= 0.0, \
                    "this_anchor_pi_alpha should be positive: {}".format(
                        this_anchor_pi_alpha)
                if N_i_j_alpha[key] > 0.0:
                    self.N_ij[key] += time_fraction * \
                        this_anchor_pi_alpha * N_i_j_alpha[key]
                else:
                    # Bring the entry into existence
                    self.N_ij[key] += 0.0
            
            R_i_alpha = self.R_i_alpha[alpha]
            if len(R_i_alpha) > 0:
                for key in R_i_alpha:
                    assert time_fraction >= 0.0, \
                        "time_fraction should be positive: {}".format(
                            time_fraction)
                    assert this_anchor_pi_alpha >= 0.0, \
                        "this_anchor_pi_alpha should be positive: {}".format(
                            this_anchor_pi_alpha)
                    if R_i_alpha[key] > 0.0:
                        self.R_i[key] += time_fraction * this_anchor_pi_alpha \
                            * R_i_alpha[key]
                    else:
                        # Bring the entry into existence
                        self.R_i[key] += 0.0
            else:
                raise Exception("R_i_alpha should always be set.")
        
        # Now we need to try and correct the zero entries
        self.N_ij_unmodified = deepcopy(self.N_ij)
        self.R_i_unmodified = deepcopy(self.R_i)
        for key in self.N_ij:
            if self.N_ij[key] == 0.0:
                i = key[0]
                if self.R_i[i] == 0.0:
                    # Then this milestone has no transitions out: give
                    #  it a high N_ij
                    self.N_ij[key] = HIGH_N_IJ
                else:
                    # Then this milestone has other transitions out: give it
                    #  a low N_IJ
                    self.N_ij[key] = LOW_N_IJ
                    
        for i in self.R_i:
            if self.R_i[i] == 0.0:
                self.R_i[i] = DEFAULT_R_I
            
        return
        
    def make_mcmc_quantities(self):
        """
        The MCMC matrix sampler requires inputs to be in Numpy matrix
        form. The Data_sample objects do not have their inputs in matrix
        form, so convert and return matrix objects that may be used in
        the MCMC matrix sampler.
        
        Finds k_alpha_beta_matrix, N_alpha_beta_matrix, T_alpha_matrix,
        mmvt_Nij_alpha, mmvt_Ri_alpha, mmvt_Qij_alpha, T
        """
        n_anchors = self.model.num_anchors
        n_milestones = self.model.num_milestones
        k_alpha_beta_matrix = np.zeros((n_anchors, n_anchors))
        N_alpha_beta_matrix = np.zeros((n_anchors, n_anchors))
        T_alpha_matrix = np.zeros((n_anchors, 1))
        for alpha in range(n_anchors):
            if self.model.anchors[alpha].bulkstate:
                continue
            for beta in range(n_anchors):
                if alpha == beta: continue
                if (alpha,beta) in self.k_alpha_beta:
                    k_alpha_beta_matrix[alpha,beta] \
                        = self.k_alpha_beta[(alpha,beta)]
                    k_alpha_beta_matrix[alpha,alpha] \
                        -= self.k_alpha_beta[(alpha,beta)]
                    N_alpha_beta_matrix[alpha,beta] \
                        = self.N_alpha_beta[alpha,beta]
            
            T_alpha_matrix[alpha,0] = self.T_alpha[alpha]

        invT = 0.0
        mmvt_Nij_alpha = []
        mmvt_Ri_alpha = []
        mmvt_Qij_alpha = []
        for alpha in range(n_anchors):
            if self.model.anchors[alpha].bulkstate:
                continue
            if self.T_alpha[alpha] > 0.0:
                invT += self.pi_alpha[alpha] / self.T_alpha[alpha]
            this_mmvt_Nij_alpha = np.zeros((n_milestones, n_milestones))
            this_mmvt_Ri_alpha = np.zeros((n_milestones, 1))
            this_mmvt_Qij_alpha = np.zeros((n_milestones, n_milestones))
            for i in range(n_milestones):
                for j in range(n_milestones):
                    key = (i,j)
                    if key in self.N_i_j_alpha[alpha]:
                        this_mmvt_Nij_alpha[i,j] = self.N_i_j_alpha[alpha][key]
                
                if i in self.R_i_alpha[alpha]:
                    this_mmvt_Ri_alpha[i,0] = self.R_i_alpha[alpha][i]
            
            for i in range(n_milestones):
                for j in range(n_milestones):
                    if i == j: continue
                    if this_mmvt_Ri_alpha[i,0] > 0.0:
                        this_mmvt_Qij_alpha[i,j] = this_mmvt_Nij_alpha[i,j] \
                            / this_mmvt_Ri_alpha[i,0]
                        this_mmvt_Qij_alpha[i,i] -= this_mmvt_Nij_alpha[i,j] \
                            / this_mmvt_Ri_alpha[i,0]
                        
            mmvt_Nij_alpha.append(this_mmvt_Nij_alpha)
            mmvt_Ri_alpha.append(this_mmvt_Ri_alpha)
            mmvt_Qij_alpha.append(this_mmvt_Qij_alpha)
            
        T = 1.0 / invT 
        
        return k_alpha_beta_matrix, N_alpha_beta_matrix, T_alpha_matrix,\
            mmvt_Nij_alpha, mmvt_Ri_alpha, mmvt_Qij_alpha, T
            
    def fill_k_from_matrices(self, k_alpha_beta_matrix):
        """
        When the k-matrix is sampled by the MCMC procedure, the entries
        need to be converted back into a dictionary form in order to
        use the existing data_sample methods.
        """
        n_anchors = self.model.num_anchors
        self.k_alpha_beta = defaultdict(float)
        self.k_alpha = []
        for alpha in range(n_anchors):
            if self.model.anchors[alpha].bulkstate:
                continue
            
            k_alpha = 0.0
            for beta in range(n_anchors):
                if alpha == beta: continue
                if k_alpha_beta_matrix[alpha,beta] > 0.0:
                    self.k_alpha_beta[(alpha,beta)] \
                        = k_alpha_beta_matrix[alpha,beta]
                    k_alpha += self.k_alpha_beta[(alpha,beta)]
                self.k_alpha.append(k_alpha)
            
        return
    
    def fill_N_R_alpha_from_matrices(self, mmvt_Nij_alpha_list, 
                                     mmvt_Ri_alpha_list):
        """
        When the N matrix and R vector is sampled by the MCMC 
        procedure, the entries need to be converted back into a 
        dictionary form in order to use the existing data_sample 
        methods.
        """
        for alpha, anchor in enumerate(self.model.anchors):
            if self.N_i_j_alpha is None or self.R_i_alpha is None:
                continue
            if anchor.bulkstate:
                break
            new_T_alpha = 0.0
            N_i_j_alpha_keys = list(self.N_i_j_alpha[alpha].keys())
            R_i_alpha_keys = list(self.R_i_alpha[alpha].keys())
            for key in N_i_j_alpha_keys:
                i, j = key
                self.N_i_j_alpha[alpha][key] = float(
                    mmvt_Nij_alpha_list[alpha][i,j])
            
            for key in R_i_alpha_keys:
                self.R_i_alpha[alpha][key] = float(
                    mmvt_Ri_alpha_list[alpha][key,0])
                new_T_alpha += self.R_i_alpha[alpha][key]
            self.T_alpha[alpha] = new_T_alpha
        
        return
    
    def calculate_extra_thermodynamics(self):
        """
        Use this data sample's statistics to construct the 
        anchor free energy profile.
        """
        self.free_energy_anchors = np.zeros(self.model.num_anchors)
        highest_pi_alpha = max(self.pi_alpha)
        for alpha in range(self.model.num_anchors):
            if self.model.anchors[alpha].bulkstate:
                pi_alpha_val = self.pi_alpha[-1, 0]
            else:
                pi_alpha_val = self.pi_alpha[alpha, 0]
            free_energy = -GAS_CONSTANT*self.model.temperature*np.log(
                pi_alpha_val / highest_pi_alpha)
            self.free_energy_anchors[alpha] = free_energy
            

def find_nonzero_matrix_entries(M):
    """
    For a given matrix M find all the nonzero entries. Return the index
    tuples in a list.
    """
    assert len(M.shape) == 2, "Only 2D matrices allowed."
    nonzero_matrix_indices = []
    n = M.shape[0]
    m = M.shape[1]
    for i in range(n):
        for j in range(m):
            if i == j:
                continue
            if M[i,j] != 0:
                nonzero_matrix_indices.append((i,j))
    return nonzero_matrix_indices

def monte_carlo_milestoning_error(
        main_data_sample, num=100, stride=None, skip=None, verbose=False):
    """
    Calculates an error estimate by sampling a distribution of rate 
    matrices assumming a Poisson (gamma) distribution with 
    parameters Nij and Ri using Markov chain Monte Carlo.
        
    Enforces detailed Balance-- using a modified version of 
    Algorithm 4 from Noe 2008 for rate matrices.-- 
    Citation: Noe, F. "Probability Distributions of molecular observables
    computed from Markov models." J. Chem. Phys. 2008, 128, No. 244103.
    Distribution is:  p(Q|N) = p(Q)p(N|Q)/p(N) = 
        p(Q) PI(q_ij**N_ij * exp(-q_ij * Ri))
        
    Parameters
    ----------
    
    main_data_sample : MMVT_data_sample()
        The data sample contains all data extracted from simulations as
        well as model information.
    
    num : int, default 1000
        number of rate matrix (Q) samples to be generated
        
    skip : int, default None
        number of inital rate matrix samples to skip for "burn in". If
        None, it will be assigned 10 * N**2, where N is the size of the
        rate matrix Q
        
    stride : int, default None
        frequency at which rate matrix samples are recorded- larger
        frequency reduces correlation between samples. If None, it will
        be assigned N**2, where N is the size of the rate matrix Q.
        
    verbose : bool, default False
        allow additional verbosity/printing
    """
    MAX_SKIP = 1000
    MAX_STRIDE = 100
    model = main_data_sample.model
    Q_ij = main_data_sample.Q
    Q = deepcopy(Q_ij)
    n_milestones = Q.shape[0] #get size of count matrix
    if skip is None:
        skip = 40*n_milestones
        if skip > MAX_SKIP:
            skip = MAX_SKIP
    if stride is None:
        stride = 4*n_milestones
        if stride > MAX_STRIDE:
            stride = MAX_STRIDE
    
    k_alpha_beta_matrix, N_alpha_beta_matrix, T_alpha_matrix, mmvt_Nij_alpha, \
        mmvt_Ri_alpha, mmvt_Qij_alpha, T \
        = main_data_sample.make_mcmc_quantities()
    MFPTs_list = defaultdict(list)
    MFPTs_error = defaultdict(float)
    k_off_list = []
    k_ons_list = defaultdict(list)
    k_ons_error = defaultdict(float)
    k_off_error = None
    p_i_list = []
    free_energy_profile_list = []
    free_energy_anchors_list = []
    data_sample_list = []
    if verbose: print("collecting ", num, " MCMC samples from ", 
                      num*(stride) + skip, " total moves")  
    new_data_sample = None
    nonzero_indices_k_alpha_beta_matrix = find_nonzero_matrix_entries(
        k_alpha_beta_matrix)
    nonzero_indices_Q_alpha_list = []
    for alpha, anchor in enumerate(model.anchors):
        if anchor.bulkstate:
            continue
        nonzero_indices_Q_alpha = find_nonzero_matrix_entries(
            mmvt_Qij_alpha[alpha])
        nonzero_indices_Q_alpha_list.append(nonzero_indices_Q_alpha)
    
    for counter in range(num * (stride) + skip + 1):
        new_data_sample = MMVT_data_sample(model)
        new_data_sample.N_alpha_beta = main_data_sample.N_alpha_beta
        new_data_sample.N_i_j_alpha = main_data_sample.N_i_j_alpha
        new_data_sample.R_i_alpha = main_data_sample.R_i_alpha
        new_data_sample.T_alpha = main_data_sample.T_alpha
        k_alpha_beta_matrix_new = markov_chain_monte_carlo\
            .irreversible_stochastic_matrix_algorithm_sample(
                k_alpha_beta_matrix, N_alpha_beta_matrix, T_alpha_matrix,
                nonzero_indices_k_alpha_beta_matrix)
        mmvt_Qnew_list = []
        for alpha, anchor in enumerate(model.anchors):
            if anchor.bulkstate:
                continue
            mmvt_Qnew_alpha = markov_chain_monte_carlo\
                .irreversible_stochastic_matrix_algorithm_sample(
                    mmvt_Qij_alpha[alpha], mmvt_Nij_alpha[alpha], 
                    mmvt_Ri_alpha[alpha], nonzero_indices_Q_alpha_list[alpha])
            mmvt_Qnew_list.append(mmvt_Qnew_alpha)
            
        if counter > skip and counter % stride == 0:
            mmvt_Nij_list = []
            for alpha, anchor in enumerate(model.anchors):
                if anchor.bulkstate:
                    continue
                new_mmvt_Nij_alpha = make_new_Nij_alpha(
                    mmvt_Qnew_list[alpha], mmvt_Ri_alpha[alpha])
                mmvt_Nij_list.append(new_mmvt_Nij_alpha)
                
            new_data_sample.fill_k_from_matrices(k_alpha_beta_matrix_new)
            new_data_sample.calculate_pi_alpha()
            new_data_sample.fill_N_R_alpha_from_matrices(
                mmvt_Nij_list, mmvt_Ri_alpha)
            new_data_sample.fill_out_data_quantities()
            new_data_sample.compute_rate_matrix()
            new_data_sample.N_ij = main_data_sample.N_ij
            new_data_sample.R_i = main_data_sample.R_i
            new_data_sample.K = common_analyze.Q_to_K(new_data_sample.Q)
            if model.using_bd():
                new_data_sample.parse_browndye_results(
                    bd_sample_from_normal=True)
            new_data_sample.calculate_thermodynamics()
            new_data_sample.calculate_extra_thermodynamics()
            new_data_sample.calculate_kinetics()
            p_i_list.append(new_data_sample.p_i)
            free_energy_profile_list.append(new_data_sample.free_energy_profile)
            free_energy_anchors_list.append(new_data_sample.free_energy_anchors)
            for key in new_data_sample.MFPTs:
                MFPTs_list[key].append(new_data_sample.MFPTs[key])
            if new_data_sample.k_off is not None:
                k_off_list.append(new_data_sample.k_off)
            if model.using_bd():
                for key in new_data_sample.k_ons:
                    k_ons_list[key].append(new_data_sample.k_ons[key])  
            data_sample_list.append(new_data_sample)
    
    assert new_data_sample is not None, "Nothing sampled in Monte Carlo "\
        "milestoning procedure. Please choose different arrangement of num, "\
        "skip, and stride."
    p_i_error = np.zeros(main_data_sample.p_i.shape)
    free_energy_profile_err = np.zeros(
        main_data_sample.free_energy_profile.shape)
    free_energy_anchors_err = np.zeros(
        main_data_sample.free_energy_anchors.shape)
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
    for i in range(free_energy_anchors_err.shape[0]):
        free_energy_anchors_val_list = []
        for j in range(len(free_energy_anchors_list)):
            free_energy_anchors_val_list.append(free_energy_anchors_list[j][i])
        free_energy_anchors_err[i] = np.std(free_energy_anchors_val_list)
    for key in main_data_sample.MFPTs:
        MFPTs_error[key] = np.std(MFPTs_list[key])    
    if len(k_off_list) > 0:
        k_off_error = np.std(k_off_list)
    if main_data_sample.model.using_bd():
        for key in main_data_sample.k_ons:
            k_ons_error[key] = np.std(k_ons_list[key])
    
    return data_sample_list, p_i_error, free_energy_profile_err, \
        free_energy_anchors_err, MFPTs_error, k_off_error, k_ons_error
