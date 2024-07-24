"""
analyze.py

Functions and objects for analyzing MMVT simulation outputs
"""

import os
import argparse
import warnings
import glob
from collections import defaultdict
import math

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

import seekr2.modules.common_base as base
import seekr2.modules.common_analyze as common_analyze
import seekr2.modules.mmvt_analyze as mmvt_analyze
import seekr2.modules.elber_analyze as elber_analyze
import seekr2.modules.check as check

LOW_FLUX = 1e-10
FAST_FLUX = 1e6

def check_graph_connected(graph, state_index1, state_index2):
    """
    A graph is made of anchor indices that have recorded bounces against
    neighboring anchors. If pathways can be found, then sufficient statistics
    likely exist to compute kinetics from the resulting bounces.
    
    Parameters:
    -----------
    graph : dict
        A dictionary whose keys are node incides and whose values are
        lists of connected node indices
    state_index1 : int
        The index of the "source" node.
    state_index2 : int
        The index o the "sink" node.
        
    Returns:
    --------
    connected : bool
        If a valid path can be found between state_index1 and 
        state_index2, then return True. Otherwise return False.
    """
    max_iter = len(graph) + 5
    if state_index1 == state_index2:
        return True
    explored_nodes = [state_index1]
    edge_nodes = graph[state_index1]
    counter = 0
    while len(edge_nodes) > 0:
        edge_node = edge_nodes.pop()
        if state_index2 == edge_node:
            return True
        explored_nodes.append(edge_node)
        for possible_edge_node in graph[edge_node]:
            if possible_edge_node not in edge_nodes \
                    and possible_edge_node not in explored_nodes:
                edge_nodes.append(possible_edge_node)
        
        counter += 1
        assert counter < max_iter, "Maximum iterations exceeded."
        
    return False

class Analysis:
    """
    Compute MMVT thermodynamics and kinetics from model data and 
    parameters. Utilizes multiple Data_sample objects to get average
    and uncertainties in kinetics/thermo quantities.
    
    Attributes
    ----------
    model : Model()
        The Model() object that will be used to read transition 
        statistics and also provides relevant parameters for 
        calculations.
        
    pi_alpha : numpy.array
        The pi_alpha quantity in MMVT theory - it represents the
        relative probability of finding the system within anchor alpha.
        It is an array of floats.
        
    pi_alpha_error : numpy.array
        The uncertainty for each value in pi_alpha.
        
        
    K : numpy.array
        The transition matrix K, which appears in MMVT theory as well
        as classical milestoning theory. K represents probabilities of
        transitions (whereas Q represents rates of transitions). Each
        element in row i and column j represent the probability of 
        transition between i and j.
        
    p_i : numpy.array
        The probability vector p_i, which appears in MMVT theory as 
        well as classical milestoning theory. p_i represents 
        the probabilities of finding the system within the vicinity of
        milestone i. Therefore, p_i can be used to obtain thermodynamic
        quantities such as the free energy of milestone i.
        
    p_i_error : numpy.array
        The uncertainty for each value in p_i.
        
    free_energy_profile : numpy.array
        The free energy profile as computed from the Boltzmann 
        distribution of p_i.
        
    free_energy_profile_err : numpy.array
        The uncertainty of each value in p_i.
        
    MFPTs : dict
        The mean first passage times (MFPTs) between various end states
        as indicated by the model. The MFPT dict has keys (i,j) where
        i and j are the indices of system end states and the dict has
        values which represent the expectation values of time spent 
        transitioning between states i and j.
        
    MFPT_errors : dict
        The errors in the MFPTs values. THe dictionary keys are 
        structured in an identical manner as MFPTs.
        
    k_off : float
        The k-off (off-rate) value. That is, the expected rate of
        unbinding from all states indicated as an end state (weighted
        by their relative probabilities defined by p_i.
        
    k_off_error : float
        The error (uncertainty) in the k-off value.
        
    k_ons : dict
        The k-ons predicted to the various end states within the
        system. The keys are the indices of the end state milestones
        and the values are the k-ons.
        
    k_on_errors : dict
        The errors in the k-on values. The keys are structured in an
        identical manner as the k_ons dictionary.
    
    force_warning : bool
        Whether to bypass certain errors with a mere warning.
        
    num_error_samples : int
        The number of data samples to create for generating uncertainty
        values.
    """
    
    def __init__(self, model, force_warning=False, num_error_samples=0,
                 stride_error_samples=None, skip_error_samples=None):
        """
        Creates the Analyze() object, which applies transition 
        statistics and times, as well as MMVT theory, to compute 
        kinetics and thermodynamics quantities.
        """
        self.model = model
        self.anchor_stats_list = []
        self.main_data_sample = None
        self.data_sample_list = []
        self.pi_alpha = None
        self.pi_alpha_error = None
        self.p_i = None
        self.p_i_error = None
        self.free_energy_profile = None
        self.free_energy_profile_err = None
        self.free_energy_anchors = None
        self.free_energy_anchors_err = None
        self.MFPTs = {}
        self.MFPTs_error = {}
        self.k_off = None
        self.k_off_error = None
        self.k_ons = {}
        self.k_ons_error = {}
        self.force_warning = force_warning
        self.num_error_samples = num_error_samples
        self.stride_error_samples = stride_error_samples
        self.skip_error_samples = skip_error_samples
        return
    
    def elber_check_anchor_stats(self, silent=False):
        """
        Check the anchor statistics to make sure that enough bounces
        have been observed to perform the analysis
        """
        anchors_missing_statistics = []
        for i, anchor in enumerate(self.model.anchors):
            if anchor.bulkstate:
                continue
            anchor_stats = self.anchor_stats_list[i]
            existing_alias_ids = []
            existing_alias_transitions = []
            for key in anchor_stats.N_i_j:
                existing_alias_transitions.append(key)
            
            # Hacky!
            existing_alias_transitions.append(2)
            for milestone in anchor.milestones:
                found_problem = False
                if milestone.alias_index not in existing_alias_transitions:
                    anchors_missing_statistics.append(anchor.index)
                    break
                
                if found_problem:
                    break
                
        if len(anchors_missing_statistics) > 0:
            if silent:
                return False
            else:
                error_warning_string = "Anchor(s) {0} are missing sufficient "\
                    "statistics. Consider running simulations of anchor(s) {0} "\
                    "for longer time scales or readjust anchor locations to "\
                    "make transitions more frequent. You may skip this check "\
                    "with the --skip_checks (-s) option.".format(
                        anchors_missing_statistics)
                if self.force_warning:
                    warnings.warn(error_warning_string)
                else:
                    raise common_analyze.MissingStatisticsError(
                        error_warning_string)
        return True
    
    def mmvt_check_anchor_stats(self, silent=False):
        """
        Check the anchor statistics to make sure that enough transitions
        have been observed to perform the analysis
        """
        error_warning_base = "Anchor(s) {0} are missing some "\
            "statistics. Consider running simulations of anchor(s) {0} "\
            "for longer time scales or readjust anchor locations to "\
            "make transitions more frequent. You may skip this check "\
            "with the --skip_checks (-s) option."
        path_not_connected_string = "Insufficient statistics were found "\
            "to create a path between end states \"{0}\" and \"{1}\"."
        error_warning_string = ""
        anchors_missing_alpha_beta_statistics = []
        anchors_missing_i_j_statistics = []
        connection_graph = defaultdict(list)
        endstate_list = []
        bulk_anchor_index = None
        for i, anchor in enumerate(self.model.anchors):
            if anchor.bulkstate:
                # we want the bulk state to indicate all connections because
                # they will be handled by BD.
                for milestone in anchor.milestones:
                    connection_graph[anchor.index].append(
                        milestone.neighbor_anchor_index)
                bulk_anchor_index = anchor.index
                continue
            if anchor.endstate:
                endstate_list.append(anchor.index)
            anchor_stats = self.anchor_stats_list[i]
            existing_alias_ids = []
            existing_alias_transitions = []
            for key in anchor_stats.k_alpha_beta:
                existing_alias_ids.append(key)
                assert anchor_stats.k_alpha_beta[key] >= 0.0, \
                    "negative k_alpha_beta values not allowed."
            for key in anchor_stats.N_i_j_alpha:
                existing_alias_transitions.append(key)
            
            for milestone in anchor.milestones:
                if milestone.alias_index in existing_alias_ids:
                    connection_graph[anchor.index].append(
                        milestone.neighbor_anchor_index)
                else:
                    if anchor.index not in anchors_missing_alpha_beta_statistics:
                        anchors_missing_alpha_beta_statistics.append(anchor.index)
                
                for milestone2 in anchor.milestones:
                    if milestone.alias_index == milestone2.alias_index:
                        continue
                    if (milestone.alias_index, milestone2.alias_index) \
                            not in existing_alias_transitions:
                        if anchor.index not in anchors_missing_i_j_statistics:
                            anchors_missing_i_j_statistics.append(anchor.index)
                        break
        
        # Check a graph connectivity
        missing_path = False
        for endstate_index in endstate_list:
            if bulk_anchor_index is not None:
                states_connected = check_graph_connected(connection_graph, endstate_index, bulk_anchor_index)
                if not states_connected:
                    if not silent:
                        print(path_not_connected_string.format(self.model.anchors[endstate_index].name, "Bulk"))
                    missing_path = True
                
                states_connected = check_graph_connected(connection_graph, bulk_anchor_index, endstate_index)
                if not states_connected:
                    if not silent:
                        print(path_not_connected_string.format("Bulk", self.model.anchors[endstate_index].name))
                    missing_path = True
            
            for endstate_index2 in endstate_list:
                if endstate_index == endstate_index2:
                    continue
                states_connected = check_graph_connected(connection_graph, endstate_index, endstate_index2)
                if not states_connected:
                    if not silent:
                        print(path_not_connected_string.format(self.model.anchors[endstate_index].name, self.model.anchors[endstate_index2].name))
                    missing_path = True
            
        if missing_path:
            error_warning_string = error_warning_base.format(
                anchors_missing_alpha_beta_statistics)
        
        if error_warning_string != "":
            if silent:
                return False
            else:
                if self.force_warning:
                    warnings.warn(error_warning_string)
                else:
                    raise common_analyze.MissingStatisticsError(
                        error_warning_string)
        return True
    
    def extract_data(self, min_time=None, max_time=None, max_step=None):
        """
        Extract the data from simulations used in this analysis.
        """
        # If possible, avoid expensive I/O
        files_already_read = False
        if len(self.anchor_stats_list) > 0:
            files_already_read = True
            
        timestep = self.model.get_timestep()
        for alpha, anchor in enumerate(self.model.anchors):
            if anchor.bulkstate:
                continue
            
            if self.model.get_type() == "mmvt":
                if max_step is not None:
                    max_time = max_step * timestep
            else:
                if max_step is not None:
                    max_time = max_step
            
            # These contain only alias_id keys, not the true id values
            if not files_already_read:
                if self.model.get_type() == "mmvt":
                    anchor_stats = mmvt_analyze.MMVT_anchor_statistics(alpha)
                    
                elif self.model.get_type() == "elber":
                    anchor_stats = elber_analyze.Elber_anchor_statistics(alpha)
            else:
                anchor_stats = self.anchor_stats_list[alpha]
                
            if anchor.md:
                output_file_glob = os.path.join(
                    self.model.anchor_rootdir, anchor.directory, 
                    anchor.production_directory, anchor.md_output_glob)
                output_file_list = glob.glob(output_file_glob)
                output_file_list = base.order_files_numerically(
                    output_file_list)
                if self.model.openmm_settings is not None:
                    anchor_stats.read_output_file_list(
                        "openmm", output_file_list, min_time, max_time, anchor,
                        timestep)
                elif self.model.namd_settings is not None:
                    anchor_stats.read_output_file_list(
                        "namd", output_file_list, min_time, max_time, anchor, 
                        timestep)
                else:
                    raise Exception("Both OpenMM and NAMD settings missing. "\
                                    "One of these must be present in the "\
                                    "model XML.")
            else:
                pass    
            
            if not files_already_read:
                self.anchor_stats_list.append(anchor_stats)
        return
        
    def check_extraction(self, silent=False):
        """
        Check whether sufficient and correct anchor statistics can 
        be used for analysis.
        """
        if self.model.get_type() == "mmvt":
            result = self.mmvt_check_anchor_stats(silent)
        if self.model.get_type() == "elber":
            result = self.elber_check_anchor_stats(silent)
        return result
        
    def fill_out_data_samples_mmvt(self):
        """
        Now that the statistics for each anchor have been extracted
        from the output files, construct the global transition
        statistics objects. Applies to systems using MMVT milestoning.
        """
        N_alpha_beta = defaultdict(int)
        k_alpha_beta = defaultdict(float)
        N_i_j_alpha = []
        R_i_alpha_total = []
        T_alpha_total = []
        
        for alpha, anchor1 in enumerate(self.model.anchors):
            if anchor1.bulkstate:
                continue
            anchor_N_alpha_beta = self.anchor_stats_list[alpha].N_alpha_beta
            anchor_k_alpha_beta = self.anchor_stats_list[alpha].k_alpha_beta
            for beta, anchor2 in enumerate(self.model.anchors):
                if anchor2.bulkstate:
                    continue
                if alpha == beta:
                    continue
                id_alias = anchor1.alias_from_neighbor_id(anchor2.index)
                if id_alias is None:
                    continue
                if len(anchor_N_alpha_beta) == 0:
                    # Then no transitions were observed in this anchor
                    default_flux = FAST_FLUX
                else:
                    default_flux = LOW_FLUX
                
                if id_alias in anchor_N_alpha_beta:
                    N_alpha_beta[(alpha, beta)] = anchor_N_alpha_beta[id_alias]
                    k_alpha_beta[(alpha, beta)] = anchor_k_alpha_beta[id_alias]
                else:
                    N_alpha_beta[(alpha, beta)] = 0
                    k_alpha_beta[(alpha, beta)] = default_flux #0.0
                
            anchor_N_i_j_alpha = self.anchor_stats_list[alpha].N_i_j_alpha
            N_i_j_alpha_element = defaultdict(int)
            R_i_alpha_element = defaultdict(float)
            
            # Fill out empty milestone statistics
            for milestone1 in anchor1.milestones:
                R_i_alpha_element[milestone1.index] = 0.0
                for milestone2 in anchor1.milestones:
                    if milestone1.index == milestone2.index:
                        continue
                    new_key = (milestone1.index, milestone2.index)
                    N_i_j_alpha_element[new_key] = 0
            
            for key in anchor_N_i_j_alpha:
                (alias_id_i, alias_id_j) = key
                id_i = anchor1.id_from_alias(alias_id_i)
                id_j = anchor1.id_from_alias(alias_id_j)
                new_key = (id_i, id_j)
                N_i_j_alpha_element[new_key] = anchor_N_i_j_alpha[key]
            N_i_j_alpha.append(N_i_j_alpha_element)
            
            anchor_R_i_alpha = self.anchor_stats_list[alpha].R_i_alpha_total
            
            for key in anchor_R_i_alpha:
                alias_id_i = key
                id_i = anchor1.id_from_alias(alias_id_i)
                assert id_i is not None, \
                    "Anchor {} missing alias {}.".format(anchor1.index, 
                                                         alias_id_i)
                R_i_alpha_element[id_i] = anchor_R_i_alpha[key]
            
            R_i_alpha_total.append(R_i_alpha_element)
            anchor_T_alpha = self.anchor_stats_list[alpha].T_alpha_total
            T_alpha_total.append(anchor_T_alpha)
            
        self.main_data_sample = mmvt_analyze.MMVT_data_sample(
            self.model, N_alpha_beta, k_alpha_beta, N_i_j_alpha, 
            R_i_alpha_total, T_alpha_total)
        
        return

    def process_data_samples_mmvt(self):
        """
        Since the global, system-side statistics have been gathered, 
        compute the thermodynamic and kinetic quantities and their
        uncertainties. Applies to systems using MMVT milestoning.
        """
        self.main_data_sample.calculate_pi_alpha()
        self.main_data_sample.fill_out_data_quantities()
        if self.model.using_bd():
            self.main_data_sample.parse_browndye_results()
        self.main_data_sample.compute_rate_matrix()
        self.main_data_sample.calculate_thermodynamics()
        self.main_data_sample.calculate_extra_thermodynamics()
        self.main_data_sample.calculate_kinetics()
        
        self.pi_alpha = self.main_data_sample.pi_alpha
        self.p_i = self.main_data_sample.p_i
        self.free_energy_profile = self.main_data_sample.free_energy_profile
        self.free_energy_anchors = self.main_data_sample.free_energy_anchors
        self.MFPTs = self.main_data_sample.MFPTs
        self.k_off = self.main_data_sample.k_off
        self.k_ons = self.main_data_sample.k_ons
        if self.num_error_samples > 0:
            data_sample_list, p_i_error, free_energy_profile_err, \
                free_energy_anchors_err, MFPTs_error, k_off_error, k_ons_error\
                = mmvt_analyze.monte_carlo_milestoning_error(
                    self.main_data_sample, num=self.num_error_samples, 
                    stride=self.stride_error_samples,
                    skip=self.skip_error_samples)
            self.data_sample_list = data_sample_list
            self.pi_alpha_error = None #pi_alpha_error
            self.p_i_error = p_i_error
            self.free_energy_profile_err = free_energy_profile_err
            self.free_energy_anchors_err = free_energy_anchors_err
            self.MFPTs_error = MFPTs_error
            self.k_off_error = k_off_error
            self.k_ons_error = k_ons_error
            
        return
    
    def fill_out_data_samples_elber(self):
        """
        Now that the statistics for each anchor have been extracted
        from the output files, construct the global transition
        statistics objects. Applies to systems using Elber milestoning.
        """
        N_i_j_list = []
        R_i_total = []
        for i, anchor1 in enumerate(self.model.anchors):
            if anchor1.bulkstate:
                continue
            
            anchor_N_i_j = self.anchor_stats_list[i].N_i_j
            N_i_j_element = defaultdict(int)
            for key in anchor_N_i_j:
                alias_id_j = key
                id_j = anchor1.id_from_alias(alias_id_j)
                new_key = (i, id_j)
                N_i_j_element[new_key] = anchor_N_i_j[key]
            N_i_j_list.append(N_i_j_element)
            
            anchor_R_i = self.anchor_stats_list[i].R_i_total
            R_i_total.append(anchor_R_i)
        
        self.main_data_sample = elber_analyze.Elber_data_sample(
            self.model, N_i_j_list, R_i_total)
        self.main_data_sample.fill_out_data_quantities()
        error_sample = elber_analyze.Elber_data_sample(
            self.model, N_i_j_list, R_i_total)
        error_sample.fill_out_data_quantities()
        self.data_sample_list.append(error_sample)
        return
        
        
    def process_data_samples_elber(self):
        """
        Since the global, system-side statistics have been gathered, 
        compute the thermodynamic and kinetic quantities and their
        uncertainties. Applies to systems using Elber milestoning.
        """
        if self.model.using_bd():
            self.main_data_sample.parse_browndye_results()
        self.main_data_sample.compute_rate_matrix()
        self.main_data_sample.calculate_thermodynamics()
        self.main_data_sample.calculate_kinetics()
        self.p_i = self.main_data_sample.p_i
        self.free_energy_profile = self.main_data_sample.free_energy_profile
        self.MFPTs = self.main_data_sample.MFPTs
        self.k_off = self.main_data_sample.k_off
        self.k_ons = self.main_data_sample.k_ons
        if self.num_error_samples > 0:
            data_sample_list, p_i_error, free_energy_profile_err, MFPTs_error, \
                k_off_error, k_ons_error \
                = elber_analyze.monte_carlo_milestoning_error(
                    self.main_data_sample, num=self.num_error_samples,
                    stride=self.stride_error_samples,
                    skip=self.skip_error_samples)
            self.data_sample_list = data_sample_list
            self.p_i_error = p_i_error
            self.free_energy_profile_err = free_energy_profile_err
            self.MFPTs_error = MFPTs_error
            self.k_off_error = k_off_error
            self.k_ons_error = k_ons_error
        
        return
        
    def fill_out_data_samples(self):
        """
        Based on the type of milestoning, construct the data samples
        and fill out their statistics.
        """
        if self.model.get_type() == "mmvt":
            self.fill_out_data_samples_mmvt()
        elif self.model.get_type() == "elber":
            self.fill_out_data_samples_elber()
        return
    
    def process_data_samples(self):
        """
        Based on the type of milestoning, use the data samples to 
        compute thermo and kinetics quantities and their uncertainties.
        """
        if self.model.get_type() == "mmvt":
            self.process_data_samples_mmvt()
        elif self.model.get_type() == "elber":
            self.process_data_samples_elber()
        return
    
    def print_results(self):
        """Print all results of the analysis calculation."""
        print("Printing results from MMVT SEEKR calculation")
        if self.k_off is not None:
            print("k_off (1/s):", common_analyze.pretty_string_value_error(
                self.k_off, self.k_off_error))
        if len(self.k_ons) > 0:
            print("k_ons :")
        for key in self.k_ons:
            k_on = float(self.k_ons[key])
            diss_constant = self.k_off / k_on
            delta_G = common_analyze.GAS_CONSTANT*self.model.temperature\
                *math.log(diss_constant)
            if key in self.k_ons_error:
                k_on_err = float(self.k_ons_error[key])
                print("  k_on (1/s * 1/M) to state", key, ":", 
                      common_analyze.pretty_string_value_error(k_on, k_on_err))
                if k_on > 0.0 and self.k_off > 0.0:
                    diss_constant_err = diss_constant \
                        * common_analyze.quadriture(
                        k_on_err/k_on, self.k_off_error/self.k_off)
                    delta_G_err = diss_constant_err*common_analyze.GAS_CONSTANT\
                    *self.model.temperature/diss_constant
                else:
                    diss_constant_err = None
                    delta_G_err = None
            else:
                print("  k_on (1/s * 1/M) to state", key, ":", 
                      common_analyze.pretty_string_value_error(k_on, None))
                diss_constant_err = None
                delta_G_err = None
            
            print("  Dissociation constant (M) to state", key, ":", 
                  common_analyze.pretty_string_value_error(
                      diss_constant, diss_constant_err))
            if key in self.k_ons_error:
                print("  \u0394G (kcal/mol) to state", key, ":", 
                      common_analyze.pretty_string_value_error(
                          delta_G, delta_G_err))
        
        print("Mean first passage times (s):")
        for key in self.MFPTs:
            state1 = key[0]
            state2 = key[1]
            if key in self.MFPTs_error:
                print("  MFPT from", state1, "to", state2, ":",
                      common_analyze.pretty_string_value_error(
                          float(self.MFPTs[key]*1.0e-12), 
                          float(self.MFPTs_error[key]*1.0e-12)))
            else:
                print("  MFPT from", state1, "to", state2, ":",
                      common_analyze.pretty_string_value_error(
                          float(self.MFPTs[key]*1.0e-12), None))
        return
    
    def make_milestone_value_dict(self):
        """
        Make a dictionary with milestone indices as keys, and
        milestone values as dictionary values.
        """
        milestone_value_dict = {}
        for i, anchor in enumerate(self.model.anchors):
            for milestone in anchor.milestones:
                if "radius" in milestone.variables:
                    value = milestone.variables["radius"]
                else:
                    value = milestone.variables["value"]
                if milestone.index not in milestone_value_dict:
                    milestone_value_dict[milestone.index] = value
                    
        return milestone_value_dict
    
    def save_plots(self, image_directory):
        """
        Save a potentially useful series of plots of some quantities
        obtained during the analysis.
        
        TODO: interact with model, because the way these plots are saved
        depends on the structure of the CVs.
        """
        anchor_indices = np.zeros(len(self.model.anchors), dtype=np.int8)
        anchor_values = np.zeros(len(self.model.anchors), dtype=np.float64)
        if len(self.model.anchors[0].variables) > 1:
            # cannot make plots for multidimensional CVs
            return
        for i, anchor in enumerate(self.model.anchors):
            anchor_indices[i] = anchor.index
            anchor_values[i] = list(anchor.variables.values())[0]
        milestone_indices = np.zeros(self.p_i.shape[0], dtype=np.int8)
        milestone_values = np.zeros(self.p_i.shape[0], dtype=np.float64)
        milestone_value_dict = self.make_milestone_value_dict()
        for i in range(self.p_i.shape[0]):
            milestone_indices[i] = i
            milestone_values[i] = milestone_value_dict[i]
        # save pi_alpha
        if self.model.get_type() == "mmvt":
            pi_fig, ax = plt.subplots()
            plt.errorbar(
                anchor_values, self.pi_alpha.flatten()[0:len(anchor_indices)],
                yerr=self.pi_alpha_error, ecolor="k", capsize=2)
            plt.xticks(anchor_values, anchor_values, rotation=90)
            plt.ylabel("\u03C0_{\u03B1}")
            plt.xlabel("anchor value")
            plt.yscale("log", nonpositive="mask")
            plt.tight_layout()
            pi_fig.savefig(os.path.join(image_directory, "pi_alpha.png"))
            
        # save p_i
        pi_fig, ax = plt.subplots()
        plt.errorbar(np.round(milestone_values, 3), self.p_i, 
                     yerr=self.p_i_error, ecolor="k", capsize=2)
        plt.xticks(np.round(milestone_values, 3), np.round(milestone_values, 3),
                   rotation=90)
        plt.ylabel("p_i")
        plt.xlabel("milestones")
        plt.yscale("log", nonpositive="mask")
        plt.tight_layout()
        pi_fig.savefig(os.path.join(image_directory, "p_i.png"))
        # save free energy milestone profile
        pi_fig, ax = plt.subplots()
        plt.errorbar(np.round(milestone_values, 3), self.free_energy_profile, 
                 yerr=self.free_energy_profile_err, ecolor="k", capsize=2)
        plt.xticks(np.round(milestone_values, 3), np.round(milestone_values, 3),
                   rotation=90)
        plt.ylabel("\u0394G(milestone) (kcal/mol)")
        plt.xlabel("milestones")
        plt.tight_layout()
        pi_fig.savefig(os.path.join(
            image_directory, "free_energy_profile_milestones.png"))
        
        if self.free_energy_anchors is not None:
            # save free energy anchor profile
            pi_fig, ax = plt.subplots()
            plt.errorbar(anchor_values, self.free_energy_anchors, 
                     yerr=self.free_energy_anchors_err, ecolor="k", capsize=2)
            plt.xticks(anchor_values, anchor_values, rotation=90)
            plt.ylabel("\u0394G(anchor) (kcal/mol)")
            plt.xlabel("anchor")
            plt.tight_layout()
            pi_fig.savefig(os.path.join(
                image_directory, "free_energy_profile_anchor.png"))
        return
        
def analyze(model, force_warning=False, num_error_samples=1000, 
            stride_error_samples=None, skip_error_samples=None,
            skip_checks=False, min_time=0.0, max_time=None):
    """Perform all the analysis steps at once."""
    curdir = os.getcwd()
    analysis = Analysis(model, force_warning, num_error_samples, 
                        stride_error_samples, skip_error_samples)
    analysis.extract_data(min_time=min_time, max_time=max_time)
    if not skip_checks:
        analysis.check_extraction()
    analysis.fill_out_data_samples()
    analysis.process_data_samples()
    os.chdir(curdir)
    return analysis

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument(
        "model_file", metavar="MODEL_FILE", type=str, 
        help="name of model file for SEEKR2 calculation. This would be the "\
        "XML file generated in the prepare stage.")
    argparser.add_argument(
        "-f", "--force_warning", dest="force_warning", default=False, 
        help="By default, missing statistics for any anchors will generate "\
        "fatal errors. This option will instead raise a warning and attempt "\
        "the calculation anyway.", 
        action="store_true")
    argparser.add_argument(
        "-n", "--num_error_samples", dest="num_error_samples", 
        default=100, type=int, help="Specify the number of error samples" \
        " to generate for estimating error/uncertainty of computed "\
        "values. Default: 100")
    argparser.add_argument(
        "-S", "--stride_error_samples", dest="stride_error_samples", 
        default=None, type=int, help="Specify the number of strides between" \
        "saved error samples. An argument of None automatically assigns the "\
        "quantity at the number of milestones in the model squared. Default: "\
        "None")
    argparser.add_argument(
        "-K", "--skip_error_samples", dest="skip_error_samples", 
        default=None, type=int, help="Specify the number of error samples. "\
        "An argument of None automatically assigns the quantity at ten times "\
        "the number of milestones in the model squared. Default: None.")
    argparser.add_argument(
        "-d", "--image_directory", dest="image_directory", 
        default=None, type=str,
        help="Define the directory where all plots and images will be saved. "\
            "By default, graphics will be saved to the "\
            "'%s' directory in the model's anchor root directory."\
            % common_analyze.DEFAULT_IMAGE_DIR)
    argparser.add_argument(
        "-s", "--skip_checks", dest="skip_checks", default=False, 
        help="By default, post-simulation checks will be run before the "\
        "analysis is started, and if the checks fail, the analysis "\
        "will not proceed. This argument bypasses those "\
        "checks and allows the analysis to proceed anyways.",
        action="store_true")
    argparser.add_argument(
        "-t", "--minimum_time", dest="minimum_time", default=0.0, type=float,
        help="A user may wish to skip an amount of simulation time for each "\
        "anchor before counting the transitions for milestoning analysis. "\
        "Enter the time (in ps) to skip a portion of the production "\
        "simulations when performing analysis.")
    argparser.add_argument(
        "-T", "--maximum_time", dest="maximum_time", default=None, type=float,
        help="A user may wish to stop the analysis of simulation time for "\
        "each anchor at a particular time. Enter the time (in ps) at which to "\
        "end the analysis at a given anchor if the simulation time exceeds it.")
    
    args = argparser.parse_args() # parse the args into a dictionary
    args = vars(args)
    model_file = args["model_file"]
    force_warning = args["force_warning"]
    num_error_samples = args["num_error_samples"]
    stride_error_samples = args["stride_error_samples"]
    skip_error_samples = args["skip_error_samples"]
    image_directory = args["image_directory"]
    skip_checks = args["skip_checks"]
    min_time = args["minimum_time"]
    max_time = args["maximum_time"]
    model = base.load_model(model_file)
    image_directory = common_analyze.make_image_directory(
        model, image_directory)
    if not skip_checks:
        check.check_post_simulation_all(model, long_check=True)
    
    analysis = analyze(model, force_warning=force_warning, 
            num_error_samples=num_error_samples, 
            stride_error_samples=stride_error_samples,
            skip_error_samples=skip_error_samples, 
            skip_checks=skip_checks, min_time=min_time, max_time=max_time)
    analysis.print_results()
        
    print("All plots being saved to:", image_directory)
    analysis.save_plots(image_directory)