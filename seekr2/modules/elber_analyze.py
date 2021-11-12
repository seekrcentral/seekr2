"""
elber_analyze.py

Routines for the analysis stage of an Elber milestoning calculation.
"""

from collections import defaultdict
from copy import deepcopy

import numpy as np

import seekr2.modules.common_analyze as common_analyze
import seekr2.modules.markov_chain_monte_carlo as markov_chain_monte_carlo

def openmm_read_output_file_list(output_file_list, min_time=None, max_time=None, 
                                 existing_lines=[], skip_restart_check=False):
    """
    Read the output files produced by the plugin (backend) of 
    openmmvt and extract transition statistics and times
    """
    MAX_ITER = 1000000000
    if len(existing_lines) == 0:
        files_lines = []
        for i, output_file_name in enumerate(output_file_list):
            file_lines = []
            with open(output_file_name, "r") as output_file:
                for counter, line in enumerate(output_file):
                    if line.startswith("#") or len(line.strip()) == 0:
                        continue
                    line_list = line.strip().split(",")
                    dest_boundary = line_list[0]
                    transition_counter = line_list[1]
                    incubation_time = float(line_list[2])
                    if dest_boundary.endswith("*"):
                        continue
                    line_list = (int(dest_boundary), incubation_time)
                    file_lines.append(line_list)
            
            files_lines.append(file_lines)
                
        lines = []
        for i, file_lines in enumerate(files_lines):
            lines += file_lines
        
    else:
        lines = existing_lines
    
    N_i_j = defaultdict(int)
    R_i_list = []
    
    for counter, line in enumerate(lines):
        dest_boundary = int(line[0])
        incubation_time = float(line[1])
        
        # TODO: fix this once the new Elber plugin is created.
        # Absolute time will need to be written by the plugin
        #if min_time is not None:
        #    if counter > min_time:
        #        break
        #if max_time is not None:
        #    if counter > max_time:
        #        break
            
        if not skip_restart_check:
            assert incubation_time >= 0.0, "incubation times cannot be "\
                "negative. Has an output file been concatenated "\
                "incorrectly? file name: %s, line number: %d" % (
                output_file_name, counter)
        N_i_j[dest_boundary] += 1
        R_i_list.append(incubation_time)
            
    
    R_i_average = np.mean(R_i_list)
    R_i_std_dev = np.std(R_i_list)
    R_i_total = np.sum(R_i_list)
        
    return N_i_j, R_i_list, R_i_average, R_i_std_dev, R_i_total, lines
        
# def openmm_read_statistics_file(statistics_file_name): ?

class Elber_anchor_statistics():
    """
    Contain all the transition statistics for an anchor within
     an Elber calculation.
    
    Attributes:
    -----------
    N_i_j : dict
        Represents the N_i_j_alpha quantities in the Elber calculation.
        The attribute represents the number of times a transition was
        observed within cell alpha going from surface i to surface j.
        The keys of the dict are a tuple (i,j) composed of two ints,
        and the value is an int count of transitions. The int i and j
        are alias_id's.
        
    R_i_list : list
        Represents the R_i_alpha quantities in the Elber calculation.
        The attribute represents the total time a transition spent
        after bouncing off of surface i before touching another 
        surface, and was observed within cell alpha. The keys of the 
        dict are an int i, and the value is a float representing time.
    
    R_i_average : float
        The average of the values in R_i_list.
        
    R_i_std_dev : float
        The standard deviation of the values in R_i_list.
        
    R_i_total : float
        The sum total of the values in R_i_list.
    
    i : int
        The index of the milestone whose transition statistics are 
        represented by this object.
    
    existing_lines : list
        If the output files for this anchor has already been read, then
        those lines will be stored here to save on I/O.
    """
    
    def __init__(self, i):
        self.N_i_j = None
        self.R_i_list = None
        self.R_i_average = None
        self.R_i_std_dev = None
        self.R_i_total = None
        self.i = i
        self.existing_lines = []
        return
    
    def read_output_file_list(self, engine, output_file_list, min_time, 
                              max_time, anchor, timestep):
        """Parse the statistics from the plugin's output file."""
        if engine == "openmm":
            self.N_i_j, self.R_i_list, self.R_i_average, self.R_i_std_dev, \
            self.R_i_total, self.existing_lines \
                = openmm_read_output_file_list(
                    output_file_list, min_time, max_time, self.existing_lines)
        elif engine == "browndye2":
            pass
        else:
            raise Exception("Engine not allowed: {}. Must be 'openmm' "\
                            "or 'browndye'.")
                
        #assert len(self.N_i_j) > 0, \
        #    "Missing statistics for anchor %d" % anchor.index
        return
    
    def print_stats(self):
        """ Print statistics for easy viewing."""
        print("self.i:", self.i)
        print("self.N_i_j:", self.N_i_j)
        print("self.R_i_average:", self.R_i_average)
        print("self.R_i_std_dev:", self.R_i_std_dev)
        print("self.R_i_total:", self.R_i_total)
        print("len(self.R_i_list):", len(self.R_i_list))
        return

class Elber_data_sample(common_analyze.Data_sample):
    """
    Represent a data sample specific to an Elber calculation.
    
    Attributes:
    -----------
    model : Model()
        The SEEKR2 model the calculation is being performed on.
    N_i_j_list : list
        A list whose indices are milestones, and each value is a dictionary
        of transitions to adjacent milestones.
    R_i_list : list
        A list whose indices are milestones, and each value is the 
        average incubation time spent in that milestone.
    N_ij : dict
        A dictionary with keys of 2-tupes of source and destination
        milestone, and whose values are the counts of transitions.
    R_i : dict
        A dictionary with keys of source milestone indices, and whose
        values are incubation times spent at those milestones.
    Q : numpy.array()
        A 2x2 rate matrix for a milestoning calculation. No sink states.
    Q_hat : numpy.array()
        A 2x2 rate matrix for a milestoning calculation. Sink state(s) 
        included.
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
    """
    
    def __init__(self, model, N_i_j_list=None, R_i_list=None):
        self.model = model
        self.N_i_j_list = N_i_j_list
        self.R_i_list = R_i_list
        self.N_ij = None
        self.R_i = None
        self.Q = None
        self.Q_hat = None
        self.K = None
        self.K_hat = None
        self.p_i = None
        self.free_energy_profile = None
        self.MFPTs = {}
        self.k_off = None
        self.k_ons = {}
        self.b_surface_k_ons_src = None
        self.b_surface_k_on_errors_src = None
        self.bd_transition_counts = {}
        self.bd_transition_probabilities = {}
        return
    
    def fill_out_data_quantities(self):
        """
        Compute quantities such as N_ij and R_i for eventual 
        construction of rate matrix Q.
        """
        if self.N_i_j_list is None or self.R_i_list is None:
            raise Exception("Unable to call fill_out_data_quantities(): "\
                            "No statistics present in Data Sample.")
            
        self.N_ij = defaultdict(int)
        self.R_i = defaultdict(float)
        for i, anchor in enumerate(self.model.anchors):
            if anchor.bulkstate:
                continue
            
            N_i_j_item = self.N_i_j_list[i]
            for key in N_i_j_item:
                self.N_ij[key] += N_i_j_item[key]
            
            R_i_item = self.R_i_list[i]
            for key in R_i_item:
                self.R_i[key] += R_i_item[key]
            
        return
        
def monte_carlo_milestoning_error(
        main_data_sample, num=1000, skip=None, stride=None, verbose=False, 
        pre_equilibrium_approx=False):
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
    
    pre_equilibrium_approx : bool, default False
        Whether to use the pre-equilibrium approximation for
        computing kinetics.
    """
    model = main_data_sample.model
    Q_ij = main_data_sample.Q
    N_ij = main_data_sample.N_ij
    R_i = main_data_sample.R_i
    Q = deepcopy(Q_ij)
    m = Q.shape[0] #get size of count matrix
    N = np.zeros((m, m))
    R = np.zeros((m, 1))
    if skip is None:
        skip = 10 * Q.shape[0]**2
    if stride is None:
        stride = Q.shape[0]**2
        
    for i in range(m):
        for j in range(m):
            N[i,j] = N_ij[i,j]
        
        R[i,0] = R_i[i]
        
    MFPTs_list = defaultdict(list)
    MFPTs_error = defaultdict(float)
    k_off_list = []
    k_ons_list = defaultdict(list)
    k_ons_error = defaultdict(float)
    p_i_list = []
    free_energy_profile_list = []
    data_sample_list = []
    if verbose: print("collecting ", num, " MCMC samples from ", 
                      num*(stride) + skip, " total moves")  
    new_data_sample = None
    for counter in range(num * (stride) + skip):
        #if verbose: print("MCMC stepnum: ", counter)
        Qnew = markov_chain_monte_carlo\
            .irreversible_stochastic_matrix_algorithm_sample(Q, N, R)
        
        if counter > skip and counter % stride == 0:
            new_data_sample = Elber_data_sample(model)
            new_data_sample.N_ij = N_ij
            new_data_sample.R_i = R_i
            new_data_sample.Q = Qnew
            new_data_sample.K = common_analyze.Q_to_K(Qnew)
            if model.k_on_info:
                new_data_sample.parse_browndye_results(
                    bd_sample_from_normal=True)
            new_data_sample.calculate_thermodynamics()
            new_data_sample.calculate_kinetics(pre_equilibrium_approx)
            p_i_list.append(new_data_sample.p_i)
            free_energy_profile_list.append(new_data_sample.free_energy_profile)
            for key in new_data_sample.MFPTs:
                MFPTs_list[key].append(new_data_sample.MFPTs[key])
            k_off_list.append(new_data_sample.k_off)
            if model.k_on_info:
                for key in new_data_sample.k_ons:
                    k_ons_list[key].append(new_data_sample.k_ons[key])  
            data_sample_list.append(new_data_sample)
            
        Q = Qnew
    
    assert new_data_sample is not None, "Nothing sampled in Monte Carlo "\
        "milestoning procedure. Please choose different arrangement of num, "\
        "skip, and stride."
    if verbose: print("final MCMC matrix", Q)
    p_i_error = np.zeros(new_data_sample.p_i.shape)
    free_energy_profile_err = np.zeros(
        new_data_sample.free_energy_profile.shape)
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
    for key in new_data_sample.MFPTs:
        MFPTs_error[key] = np.std(MFPTs_list[key])
    k_off_error = np.std(k_off_list)
    if new_data_sample.model.k_on_info:
        for key in new_data_sample.k_ons:
            k_ons_error[key] = np.std(k_ons_list[key])
    
    return data_sample_list, p_i_error, free_energy_profile_err, MFPTs_error, \
        k_off_error, k_ons_error
