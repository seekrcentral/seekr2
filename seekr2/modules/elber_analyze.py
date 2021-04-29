"""
Routines for the analysis stage of an Elber milestoning calculation
"""

from collections import defaultdict

import numpy as np

import seekr2.modules.common_analyze as common_analyze

def openmm_read_output_file_list(output_file_list, max_time=None, 
                                 existing_lines=[], skip_restart_check=False):
    """
    Read the output files produced by the plugin (backend) of 
    openmmvt and extract transition statistics and times

    Parameters
    ----------
    output_file_list : list
        The names of the openmmvt output file to read for 
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
    MAX_ITER = 200000
    if len(existing_lines) == 0:
        files_lines = []
        for i, output_file_name in enumerate(output_file_list):
            file_lines = []
            with open(output_file_name, "r") as output_file:
                for counter, line in enumerate(output_file):
                    if line.startswith("#") or len(line.strip()) == 0:
                        continue
                    line_list = line.strip().split(" ")
                    dest_boundary = line_list[0]
                    incubation_time = float(line_list[1])
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
        Represents the N_i_j_alpha quantities in the MMVT calculation.
        The attribute represents the number of times a transition was
        observed within cell alpha going from surface i to surface j.
        The keys of the dict are a tuple (i,j) composed of two ints,
        and the value is an int count of transitions.
        
    R_i_list : dict
        Represents the R_i_alpha quantities in the MMVT calculation.
        The attribute represents the total time a transition spent
        after bouncing off of surface i before touching another 
        surface, and was observed within cell alpha. The keys of the 
        dict are an int i, and the value is a float representing time.
        
    i : int
        The index of the milestone whose transition statistics are 
        represented by this object.
    
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
        """
        
        """
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
    Represent a data sample specific to an MMVT calculation.
    
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
    
    free_energy_profile : 
    
    MFPTs : 
    
    k_off : 
    
    k_ons : 
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
        self.p_i_hat = None
        self.free_energy_profile = None
        self.MFPTs = {}
        self.k_off = None
        self.k_ons = {}
        self.bd_transition_counts = {}
        return
    
    
    def fill_out_data_quantities(self):
        """
        Compute quantities such as N_ij, R_i, and T for eventual 
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
        