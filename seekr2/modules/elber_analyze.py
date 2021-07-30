"""
elber_analyze.py

Routines for the analysis stage of an Elber milestoning calculation.
"""

from collections import defaultdict

import numpy as np

import seekr2.modules.common_analyze as common_analyze

def openmm_read_output_file_list(output_file_list, max_time=None, 
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
    
    def read_output_file_list(self, engine, output_file_list, max_time, anchor, 
                              timestep):
        """Parse the statistics from the plugin's output file."""
        if engine == "openmm":
            self.N_i_j, self.R_i_list, self.R_i_average, self.R_i_std_dev, \
            self.R_i_total, self.existing_lines \
                = openmm_read_output_file_list(
                    output_file_list, max_time, self.existing_lines)
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
    
    def __init__(self, model, N_i_j_list, R_i_list):
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
        #self.b_surface_reaction_probabilities = None # REMOVE?
        #self.b_surface_reaction_probability_errors = None # REMOVE?
        #self.b_surface_transition_counts = None # REMOVE?
        self.bd_transition_counts = {}
        self.bd_transition_probabilities = {}
        return
    
    def fill_out_data_quantities(self):
        """
        Compute quantities such as N_ij and R_i for eventual 
        construction of rate matrix Q.
        """
        
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
        