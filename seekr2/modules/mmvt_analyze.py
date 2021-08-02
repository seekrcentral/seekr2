"""
Routines for the analysis stage of an MMVT calculation
"""

from collections import defaultdict

import numpy as np
import scipy.linalg as la

import seekr2.modules.common_analyze as common_analyze

def openmm_read_output_file_list(output_file_list, max_time=None, 
                                 existing_lines=[], skip_restart_check=False):
    """
    Read the output files produced by the plugin (backend) of 
    SEEKR2 and extract transition statistics and times
    """
    
    MAX_ITER = 1000000000
    if len(existing_lines) == 0:
        files_lines = []
        start_times = []
        for i, output_file_name in enumerate(output_file_list):
            start_time = None
            file_lines = []
            with open(output_file_name, "r") as output_file:
                for counter, line in enumerate(output_file):
                    if line.startswith("#") or len(line.strip()) == 0:
                        continue
                    line_list = line.strip().split(',')
                    file_lines.append(line_list)
                    dest_boundary = int(line_list[0])
                    bounce_index = int(line_list[1])
                    # TODO: change the following value to interval time
                    dest_time = float(line_list[2]) 
                    if start_time is None:
                        start_time = dest_time
                
                if len(file_lines) == 0:
                    continue
                
                if start_time is None:
                    start_times.append(0.0)
                else:
                    start_times.append(start_time)
            
            files_lines.append(file_lines)
                
        lines = []
        for i, file_lines in enumerate(files_lines):
            if not skip_restart_check and not len(file_lines) == 0:
                # make sure that we delete lines that occurred after a restart 
                # backup
                if i < len(files_lines)-1:
                    # then it's not the last file
                    next_start_time = start_times[i+1]
                else:
                    next_start_time = 1e99
                counter = 0
                while float(file_lines[-1][2]) > next_start_time:
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
        
        dest_boundary = int(line[0])
        bounce_index = int(line[1])
        # TODO: change the following value to interval time
        dest_time = float(line[2]) 
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
            """
            Problem situation: If a simulation fails, but bounces
            have occurred after the checkpoint, then when the 
            restart is restored, then those bounces will be double-
            counted in the restart as well.
            Solution: look at all starting times in the files, and 
            disregard all bounces that occur after the next files'
            starting time.
            """
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
        T_alpha_list.append(dest_time - last_bounce_time)
        #R_i_alpha_dict[src_boundary] += dest_time - last_bounce_time
        #assert dest_time - last_bounce_time != 0, \
        #    "Two bounces occurred at the same time and at the same "\
        #    "milestone in file name: %s, line number: %d" % (
        #        output_file_name, counter)
        last_bounce_time = dest_time
    
    R_i_alpha_average = defaultdict(float)
    R_i_alpha_std_dev = defaultdict(float)
    R_i_alpha_total = defaultdict(float)
    
    for key in R_i_alpha_list:
        R_i_alpha_average[key] = np.mean(R_i_alpha_list[key])
        R_i_alpha_std_dev[key] = np.std(R_i_alpha_list[key])
        R_i_alpha_total[key] = np.sum(R_i_alpha_list[key])
        
    T_alpha_average = np.mean(T_alpha_list)
    T_alpha_std_dev = np.std(T_alpha_list)
    T_alpha_total = np.sum(T_alpha_list)
    
    if len(R_i_alpha_list) == 0:
        R_i_alpha_average[src_boundary] = T_alpha_average
        R_i_alpha_std_dev[src_boundary] = T_alpha_std_dev
        R_i_alpha_total[src_boundary] = T_alpha_total
        R_i_alpha_list[src_boundary] = T_alpha_list
        
    return N_i_j_alpha, R_i_alpha_list, R_i_alpha_average, \
        R_i_alpha_std_dev, R_i_alpha_total, N_alpha_beta, T_alpha_list, \
        T_alpha_average, T_alpha_std_dev, T_alpha_total, lines

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
                               max_time=None, existing_lines=[], 
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
    if len(existing_lines) == 0:
        files_lines = []
        start_times = []
        for i, output_file_name in enumerate(output_file_list):
            start_time = None
            file_lines = []
            with open(output_file_name, "r") as output_file:
                for counter, line in enumerate(output_file):
                    if not line.startswith("SEEKR") or len(line.strip()) == 0:
                        continue
                    elif line.startswith("SEEKR: Cell"):
                        line = line.strip().split(" ")
                        file_lines.append(line)
                        current_stepnum = int(line[8].strip(","))
                        if start_time is None:
                            start_time = current_stepnum
                    
                    elif line.startswith("SEEKR: Milestone"):
                        line = line.strip().split(" ")
                        file_lines.append(line)
                        current_stepnum = int(line[10].strip(","))
                        if start_time is None:
                            start_time = current_stepnum
                
                start_times.append(start_time)
            
            files_lines.append(file_lines)
                
        lines = []
        for i, file_lines in enumerate(files_lines):
            if not skip_restart_check:
                # make sure that we delete lines that occurred after a restart 
                # backup
                if i < len(files_lines)-1:
                    # then it's not the last file
                    # TODO: this feature is broken
                    #next_start_time = start_times[i+1]
                    next_start_time = 1e99
                else:
                    next_start_time = 1e99
                counter = 0
                last_line = file_lines[-1]
                if last_line[1].startswith("Cell"):
                    last_stepnum = int(last_line[8].strip(","))
                elif last_line[1].startswith("Milestone"):
                    last_stepnum = int(last_line[10].strip(","))
                else:
                    raise Exception("Problem: last_line: {}".format(last_line))
                    
                while last_stepnum > next_start_time:
                    file_lines.pop()
                    last_line = file_lines[-1]
                    if last_line[1].startswith("SEEKR: Cell"):
                        last_stepnum = int(last_line[8].strip(","))
                    elif last_line[1].startswith("SEEKR: Milestone"):
                        last_stepnum = int(last_line[10].strip(","))
                    if counter > MAX_ITER:
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
    alias_src = None
        
    for counter, line in enumerate(lines):
        if line[1].startswith("Cell"):
            current_anchor_raw = int(line[4].strip(","))
            next_anchor_raw = int(line[6].strip(","))
            current_stepnum = int(line[8].strip(","))
            diff = next_anchor_raw - current_anchor_raw
            next_anchor = anchor.index + diff
            dest_boundary = None
            dest_time = current_stepnum * timestep
            if max_time is not None:
                if dest_time > max_time:
                    break
                
            for milestone in anchor.milestones:
                if milestone.neighbor_anchor_index == next_anchor:
                    dest_boundary = milestone.index
                    alias_dest = milestone.alias_index
                    alias_src = alias_dest
                    
            assert dest_boundary is not None
            
        elif line[1].startswith("Milestone"):
            src_boundary = line[6].strip(",")
            dest_boundary = int(line[8].strip(",")) # NOTE: make it alias_id, not id
            current_stepnum = int(line[10].strip(","))
            dest_time = current_stepnum * timestep
            if max_time is not None:
                if dest_time > max_time:
                    break
                
            if src_boundary == "none":
                src_boundary = dest_boundary
                src_time = dest_time
                # TODO: the situation will need to be eventually resolved
                # when a restart backup is restored and the later files
                # have bounces which occur at a later time than the first
                # bounces of the open file. Those overwritten statistics
                # will have to be discarded somehow. The following
                # assertion should, in the meantime, prevent those
                # statistics from being inadvertently counted
                """
                if dest_time < last_bounce_time:
                    print("dest_time:", dest_time)
                    print("last_bounce_time:", last_bounce_time)
                assert dest_time > last_bounce_time
                """
                # Handle this if 'restart' is in the file name.
                
                continue
            
            else:
                src_boundary = int(src_boundary)
                incubation_steps = int(line[13].strip(","))
                incubation_time = incubation_steps * timestep
                time_diff = dest_time - src_time
                assert time_diff >= 0.0, "incubation times cannot be "\
                    "negative. Has an output file been concatenated "\
                    "incorrectly? file name: %s, line number: %d" % (
                        output_file_name, counter)
                alias_src = None
                alias_dest = None
                for milestone in anchor.milestones:
                    if milestone.index == src_boundary:
                        alias_src = milestone.alias_index
                    elif milestone.index == dest_boundary:
                        alias_dest = milestone.alias_index
                
                N_i_j_alpha[(alias_src, alias_dest)] += 1
                R_i_alpha_list[alias_src].append(incubation_time)
            src_boundary = dest_boundary
            src_time = dest_time
            
        else:
            raise Exception("Unable to handle line: %s" % line)
        
        N_alpha_beta[alias_dest] += 1
        T_alpha_list.append(dest_time - last_bounce_time)
        """
        assert dest_time - last_bounce_time != 0, \
            "Two bounces occurred at the same time and at the same "\
            "milestone in file name: %s, line number: %d" % (
                output_file_name, counter)
        """
        last_bounce_time = dest_time
    
    R_i_alpha_average = defaultdict(float)
    R_i_alpha_std_dev = defaultdict(float)
    R_i_alpha_total = defaultdict(float)
    for key in R_i_alpha_list:
        R_i_alpha_average[key] = np.mean(R_i_alpha_list[key])
        R_i_alpha_std_dev[key] = np.std(R_i_alpha_list[key])
        R_i_alpha_total[key] = np.sum(R_i_alpha_list[key])
        
    T_alpha_average = np.mean(T_alpha_list)
    T_alpha_std_dev = np.std(T_alpha_list)
    T_alpha_total = np.sum(T_alpha_list)
    
    if len(R_i_alpha_list) == 0 and alias_src is not None:
        R_i_alpha_average[alias_src] = T_alpha_average
        R_i_alpha_std_dev[alias_src] = T_alpha_std_dev
        R_i_alpha_total[alias_src] = T_alpha_total
        R_i_alpha_list[alias_src] = T_alpha_list
        
        
    return N_i_j_alpha, R_i_alpha_list, R_i_alpha_average, \
        R_i_alpha_std_dev, R_i_alpha_total, N_alpha_beta, T_alpha_list, \
        T_alpha_average, T_alpha_std_dev, T_alpha_total, lines

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
    
    existing_lines : list
        If the output files for this anchor has already been read, then
        those lines will be stored here to save on I/O.
    """
    
    def __init__(self, alpha):
        
        self.N_i_j_alpha = None
        self.R_i_alpha_list = None
        self.R_i_alpha_average = None
        self.R_i_alpha_std_dev = None
        self.R_i_alpha_total = None
        self.T_alpha_list = []
        self.T_alpha_average = None
        self.T_alpha_std_dev = None
        self.T_alpha_total = None
        self.N_alpha_beta = None
        self.k_alpha_beta = {}
        self.alpha = alpha
        self.existing_lines = []
        return
    
    def read_output_file_list(self, engine, output_file_list, max_time, anchor, 
                              timestep):
        """
        Depending on the engine and other settings, read the SEEKR2 
        output files to fill out transition statistics.
        """
        if engine == "openmm":
            self.N_i_j_alpha, self.R_i_alpha_list, self.R_i_alpha_average, \
            self.R_i_alpha_std_dev, self.R_i_alpha_total, self.N_alpha_beta, \
            self.T_alpha_list, self.T_alpha_average, self.T_alpha_std_dev, \
            self.T_alpha_total, self.existing_lines \
                = openmm_read_output_file_list(
                    output_file_list, max_time, self.existing_lines)
        elif engine == "namd":
            self.N_i_j_alpha, self.R_i_alpha_list, self.R_i_alpha_average, \
            self.R_i_alpha_std_dev, self.R_i_alpha_total, self.N_alpha_beta, \
            self.T_alpha_list, self.T_alpha_average, self.T_alpha_std_dev, \
            self.T_alpha_total, self.existing_lines \
                = namd_read_output_file_list(
                    output_file_list, anchor, timestep, max_time, 
                    self.existing_lines)
        elif engine == "browndye2":
            pass
        else:
            raise Exception("Engine not allowed: {}. Must be 'openmm', "\
                            "'namd', or 'browndye'.")
        for key in self.N_alpha_beta:
            self.k_alpha_beta[key] = self.N_alpha_beta[key] / self.T_alpha_total
        assert self.T_alpha_total >= 0.0
        #assert len(self.k_alpha_beta) > 0, \
        #    "Missing statistics for anchor %d" % anchor.index
        #self.print_stats()
        return
    
    def print_stats(self):
        """ Print statistics for easy viewing."""
        print("self.alpha:", self.alpha)
        print("self.N_i_j_alpha:", self.N_i_j_alpha)
        print("self.R_i_alpha_average:", self.R_i_alpha_average)
        print("self.R_i_alpha_std_dev:", self.R_i_alpha_std_dev)
        print("self.R_i_alpha_total:", self.R_i_alpha_total)
        if self.R_i_alpha_list is not None:
            print("len(self.R_i_alpha_list):", len(self.R_i_alpha_list))
        print("self.T_alpha_average:", self.T_alpha_average)
        print("self.T_alpha_std_dev:", self.T_alpha_std_dev)
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
        self.N_ij = None
        self.R_i = None
        self.T = None
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
        self.b_surface_k_ons_src = None
        self.b_surface_k_on_errors_src = None
        #self.b_surface_reaction_probabilities = None # REMOVE?
        #self.b_surface_reaction_probability_errors = None # REMOVE?
        #self.b_surface_transition_counts = None # REMOVE?
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
        if self.N_alpha_beta is None or self.k_alpha_beta is None \
                or self.N_i_j_alpha is None or self.R_i_alpha is None \
                or self.T_alpha is None:
            raise Exception("Unable to call calculate_pi_alpha(): "\
                            "No statistics present in Data Sample.")
        
        flux_matrix_dimension = self.model.num_anchors
        self.pi_alpha = np.zeros(flux_matrix_dimension)
        flux_matrix = np.zeros((flux_matrix_dimension, flux_matrix_dimension))
        column_sum = np.zeros(flux_matrix_dimension)
        bulk_index = None
        for alpha in range(flux_matrix_dimension):
            anchor1 = self.model.anchors[alpha] # source anchor
            if anchor1.bulkstate:
                assert bulk_index is None, "Only one bulk state is allowed "\
                    "in model"
                flux_matrix[alpha,alpha] = 1.0
                bulk_index = alpha
                continue
                
            dead_end_anchor = False
            if len(anchor1.milestones) == 1:
                dead_end_anchor = True
            for beta in range(flux_matrix_dimension):
                anchor2 = self.model.anchors[beta] # destination anchor
                if anchor2.bulkstate :
                    flux_matrix[alpha, beta] = 1.0
                    id_alias = anchor1.alias_from_neighbor_id(
                        anchor2.index)
                    if id_alias is None:
                        flux_matrix[beta, alpha] = 0.0
                    else:
                        flux_matrix[beta, alpha] = 1.0
                else:
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
            
        prob_equil = np.zeros((flux_matrix_dimension,1))
        prob_equil[bulk_index] = 1.0
        self.pi_alpha = abs(la.solve(flux_matrix.T, prob_equil))
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
                continue
            this_anchor_pi_alpha = float(self.pi_alpha[alpha])
            sum_pi_alpha_over_T_alpha += this_anchor_pi_alpha \
                / self.T_alpha[alpha]
        
        self.T = 1.0 / sum_pi_alpha_over_T_alpha
        assert self.T >= 0.0, "self.T should be positive: {}".format(self.T)
        for alpha, anchor in enumerate(self.model.anchors):
            if anchor.bulkstate:
                continue
            this_anchor_pi_alpha = float(self.pi_alpha[alpha])
            T_alpha = self.T_alpha[alpha]
            assert T_alpha >= 0.0, \
                "T_alpha should be positive: {}".format(T_alpha)
            time_fraction = self.T / T_alpha
            N_i_j_alpha = self.N_i_j_alpha[alpha]
            
            for key in N_i_j_alpha:
                assert time_fraction >= 0.0, \
                    "time_fraction should be positive: {}".format(time_fraction)
                assert this_anchor_pi_alpha >= 0.0, \
                    "this_anchor_pi_alpha should be positive: {}".format(
                        this_anchor_pi_alpha)
                assert N_i_j_alpha[key] >= 0.0, \
                    "N_i_j_alpha[key] should be positive: {}".format(
                        N_i_j_alpha[key])
                self.N_ij[key] += time_fraction * \
                    this_anchor_pi_alpha * N_i_j_alpha[key]
            
            R_i_alpha = self.R_i_alpha[alpha]
            if len(R_i_alpha) > 0:
                for key in R_i_alpha:
                    assert time_fraction >= 0.0, \
                        "time_fraction should be positive: {}".format(time_fraction)
                    assert this_anchor_pi_alpha >= 0.0, \
                        "this_anchor_pi_alpha should be positive: {}".format(
                            this_anchor_pi_alpha)
                    assert R_i_alpha[key] >= 0.0, \
                        "R_i_alpha[key] should be positive: {}".format(
                            R_i_alpha[key])
                    self.R_i[key] += time_fraction * this_anchor_pi_alpha * \
                        R_i_alpha[key]
                        
            else:
                raise Exception("R_i_alpha should always be set.")
                
        return
