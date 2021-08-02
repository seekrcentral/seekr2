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

import seekr2.modules.common_base as base
import seekr2.modules.common_analyze as common_analyze
import seekr2.modules.mmvt_analyze as mmvt_analyze
import seekr2.modules.elber_analyze as elber_analyze
import seekr2.modules.check as check

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
        
    N_ij : defaultdict
        The N_ij quantity in MMVT theory - It represents an estimate of
        the relative number of transitions starting from milestone i 
        and ending in milestone j, regardless of the anchor the system
        is in. In this defaultdict object, the keys are tuples of ints
        (i,j), which represent the milestone indices of the source and
        destination milestones, respectively, and the values are floats
        representing the relative 'counts' of transitions from i to j.
    
    R_i : defaultdict
        The R_i quantity in MMVT theory - it represents an estimate of
        the relative time spent after encountering milestone i, 
        regardless of the anchor the system is in. In this defaultdict
        object, the keys are ints i representing the milestone index
        and the values are floats representing time.
        
    T : float
        The T quantity in MMVT theory - it represents the "total" time
        spent in simulation. It is merely used as a normalization and
        to cancel units in the theory.
        
    Q : numpy.array
        The rate matrix Q in MMVT theory - most of the important
        quantities in MMVT are computed using Q. Each element in row
        i and column j represents the rate of transition between i and
        j.
        
    K : numpy.array
        The transition matrix K, which appears in MMVT theory as well
        as classical milestoning theory. K represents probabilities of
        transitions (whereas Q represents rates of transitions). Each
        element in row i and column j represent the probability of 
        transition between i and j.
        
    p_i : numpy.array
        The probability vector p_i, which appears in MMVT theory as 
        well as classical milestoning theory. p_i represents 
        the probabilites of finding the system within the vicinity of
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
        system. The keys are the indeces of the end state milestones
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
    
    def __init__(self, model, force_warning=False, num_error_samples=0):
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
        self.MFPTs = {}
        self.MFPTs_error = {}
        self.k_off = None
        self.k_off_error = None
        self.k_ons = {}
        self.k_ons_error = {}
        self.force_warning = force_warning
        self.num_error_samples = num_error_samples
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
        
        anchors_missing_statistics = []
        for i, anchor in enumerate(self.model.anchors):
            if anchor.bulkstate:
                continue
            anchor_stats = self.anchor_stats_list[i]
            existing_alias_ids = []
            existing_alias_transitions = []
            for key in anchor_stats.k_alpha_beta:
                existing_alias_ids.append(key)
            for key in anchor_stats.N_i_j_alpha:
                existing_alias_transitions.append(key)
            
            for milestone in anchor.milestones:
                found_problem = False
                if milestone.alias_index not in existing_alias_ids:
                    anchors_missing_statistics.append(anchor.index)
                    break
                
                for milestone2 in anchor.milestones:
                    if milestone.alias_index == milestone2.alias_index:
                        continue
                    if (milestone.alias_index, milestone2.alias_index) \
                            not in existing_alias_transitions:
                        anchors_missing_statistics.append(anchor.index)
                        found_problem = True
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
    
    def extract_data(self, max_step_list=None, silence_errors=True):
        """
        Extract the data from simulations used in this analysis.
        """
        
        # If possible, avoid expensive I/O
        files_already_read = False
        if len(self.anchor_stats_list) > 0:
            files_already_read = True
            
        if self.model.openmm_settings is not None:
            timestep = self.model.openmm_settings.langevin_integrator.timestep
        elif self.model.namd_settings is not None:
            timestep = self.model.namd_settings.langevin_integrator.timestep
        else:
            raise Exception("No OpenMM or NAMD simulation settings in model.")
        
        for alpha, anchor in enumerate(self.model.anchors):
            if anchor.bulkstate:
                continue
            if max_step_list is not None:
                max_time = max_step_list[alpha] * timestep
            else:
                max_time = None
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
                if not silence_errors:
                    assert len(output_file_list) > 0, \
                        "Files not found: %s" % output_file_glob
                if self.model.openmm_settings is not None:
                    anchor_stats.read_output_file_list(
                        "openmm", output_file_list, max_time, anchor, timestep)
                elif self.model.namd_settings is not None:
                    anchor_stats.read_output_file_list(
                        "namd", output_file_list, max_time, anchor, timestep)
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
        R_i_alpha_average = []
        R_i_alpha_std_dev = []
        R_i_alpha_count = []
        T_alpha_total = []
        T_alpha_average = []
        T_alpha_std_dev = []
        T_alpha_count = []
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
                if id_alias in anchor_N_alpha_beta:
                    N_alpha_beta[(alpha, beta)] = anchor_N_alpha_beta[id_alias]
                    k_alpha_beta[(alpha, beta)] = anchor_k_alpha_beta[id_alias]
                else:
                    N_alpha_beta[(alpha, beta)] = 0
                    k_alpha_beta[(alpha, beta)] = 0.0
                
            anchor_N_i_j_alpha = self.anchor_stats_list[alpha].N_i_j_alpha
            N_i_j_alpha_element = defaultdict(int)
            for key in anchor_N_i_j_alpha:
                (alias_id_i, alias_id_j) = key
                id_i = anchor1.id_from_alias(alias_id_i)
                id_j = anchor1.id_from_alias(alias_id_j)
                new_key = (id_i, id_j)
                N_i_j_alpha_element[new_key] = anchor_N_i_j_alpha[key]
            N_i_j_alpha.append(N_i_j_alpha_element)
            
            anchor_R_i_alpha = self.anchor_stats_list[alpha].R_i_alpha_total
            anchor_R_i_alpha_std = self.anchor_stats_list[alpha].R_i_alpha_std_dev
            anchor_R_i_alpha_list = self.anchor_stats_list[alpha].R_i_alpha_list
            R_i_alpha_element = defaultdict(float)
            R_i_alpha_count_element = defaultdict(int)
            R_i_alpha_std_element = defaultdict(float)
            for key in anchor_R_i_alpha:
                alias_id_i = key
                id_i = anchor1.id_from_alias(alias_id_i)
                R_i_alpha_element[id_i] = anchor_R_i_alpha[key]
                R_i_alpha_std_element[id_i] = anchor_R_i_alpha_std[key]
                R_i_alpha_count_element[id_i] = len(anchor_R_i_alpha_list[key])
                
            R_i_alpha_total.append(R_i_alpha_element)
            R_i_alpha_std_dev.append(R_i_alpha_std_element)
            R_i_alpha_count.append(R_i_alpha_count_element)            
            anchor_T_alpha = self.anchor_stats_list[alpha].T_alpha_total
            anchor_T_alpha_std = self.anchor_stats_list[alpha].T_alpha_std_dev
            anchor_T_alpha_list = self.anchor_stats_list[alpha].T_alpha_list
            T_alpha_total.append(anchor_T_alpha)
            T_alpha_std_dev.append(anchor_T_alpha_std)
            T_alpha_count.append(len(anchor_T_alpha_list))
            
        self.main_data_sample = mmvt_analyze.MMVT_data_sample(
            self.model, N_alpha_beta, k_alpha_beta, N_i_j_alpha, 
            R_i_alpha_total, T_alpha_total)
        
        for i in range(self.num_error_samples):
            sampled_k_alpha_beta, sampled_N_i_j_alpha, \
            sampled_R_i_alpha_total, sampled_T_alpha_total \
                = self.resample_k_N_R_T(
                    N_alpha_beta, N_i_j_alpha, R_i_alpha_total,
                    R_i_alpha_average, R_i_alpha_std_dev, R_i_alpha_count,
                    T_alpha_total, T_alpha_average, T_alpha_std_dev, 
                    T_alpha_count)
            data_sample = mmvt_analyze.MMVT_data_sample(
                self.model, N_alpha_beta, sampled_k_alpha_beta, 
                sampled_N_i_j_alpha, sampled_R_i_alpha_total, 
                sampled_T_alpha_total)
            self.data_sample_list.append(data_sample)
        return

    def process_data_samples_mmvt(self, pre_equilibrium_approx=False):
        """
        Since the global, system-side statistics have been gathered, 
        compute the thermodynamic and kinetic quantities and their
        uncertainties. Applies to systems using MMVT milestoning.
        """
        self.main_data_sample.calculate_pi_alpha()
        self.main_data_sample.fill_out_data_quantities()
        if self.model.k_on_info is not None:
            self.main_data_sample.parse_browndye_results()
        self.main_data_sample.compute_rate_matrix()
        self.main_data_sample.calculate_thermodynamics()
        self.main_data_sample.calculate_kinetics(pre_equilibrium_approx)
        # do data_sample_list here
                     
        k_offs = []
        p_i_list = []
        pi_alpha_list = []
        free_energy_profile_list = []
        MFPTs_list = defaultdict(list)
        k_ons_list = defaultdict(list)
        for i in range(self.num_error_samples):
            data_sample = self.data_sample_list[i]
            data_sample.calculate_pi_alpha()
            data_sample.fill_out_data_quantities()
            if self.model.k_on_info is not None:
                data_sample.parse_browndye_results(bd_sample_from_normal=True)
            data_sample.compute_rate_matrix()
            data_sample.calculate_thermodynamics()
            data_sample.calculate_kinetics(pre_equilibrium_approx)
            k_offs.append(data_sample.k_off)
            p_i_list.append(data_sample.p_i)
            pi_alpha_list.append(data_sample.pi_alpha)
            free_energy_profile_list.append(data_sample.free_energy_profile)
            for key in data_sample.MFPTs:
                MFPTs_list[key].append(data_sample.MFPTs[key])
            for key in data_sample.k_ons:
                k_ons_list[key].append(data_sample.k_ons[key])
        
        pi_alpha_error = np.zeros(self.main_data_sample.pi_alpha.shape[0])
        p_i_error = np.zeros(self.main_data_sample.p_i.shape)
        free_energy_profile_err = np.zeros(
            self.main_data_sample.free_energy_profile.shape)
        k_off_error = None
        MFPTs_error = {}
        k_ons_error = {}
        if len(k_offs) > 0:
            k_off_error = np.std(k_offs)
            for i in range(pi_alpha_error.shape[0]):
                pi_alpha_val_list = []
                for j in range(len(pi_alpha_list)):
                    pi_alpha_val_list.append(pi_alpha_list[j][i])
                pi_alpha_error[i] = np.std(pi_alpha_val_list)
            
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
        
            for key in self.main_data_sample.MFPTs:
                MFPTs_error[key] = np.std(MFPTs_list[key])
            
            for key in self.main_data_sample.k_ons:
                k_ons_error[key] = np.std(k_ons_list[key])
        
        self.pi_alpha = self.main_data_sample.pi_alpha
        self.pi_alpha_error = pi_alpha_error
        self.p_i = self.main_data_sample.p_i
        self.p_i_error = p_i_error
        self.free_energy_profile = self.main_data_sample.free_energy_profile
        self.free_energy_profile_err = free_energy_profile_err
        self.MFPTs = self.main_data_sample.MFPTs
        self.MFPTs_error = MFPTs_error
        self.k_off = self.main_data_sample.k_off
        self.k_off_error = k_off_error
        self.k_ons = self.main_data_sample.k_ons
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
        R_i_average = []
        R_i_std_dev = []
        R_i_count = []
        bulkstate = None
        for i, anchor1 in enumerate(self.model.anchors):
            if anchor1.bulkstate:
                bulkstate = i
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
            anchor_R_i_std = self.anchor_stats_list[i].R_i_std_dev
            anchor_R_i_list = self.anchor_stats_list[i].R_i_list
            R_i_element = defaultdict(float)
            R_i_count_element = defaultdict(int)
            R_i_std_element = defaultdict(float)
            
            
            R_i_element[i] = anchor_R_i
            R_i_std_element[i] = anchor_R_i_std
            R_i_count_element[i] = len(anchor_R_i_list)
                
            R_i_total.append(R_i_element)
            R_i_std_dev.append(R_i_std_element)
            R_i_count.append(R_i_count_element)            
        
        self.main_data_sample = elber_analyze.Elber_data_sample(
            self.model, N_i_j_list, R_i_total)
        self.main_data_sample.fill_out_data_quantities()
        error_sample = elber_analyze.Elber_data_sample(
            self.model, N_i_j_list, R_i_total)
        error_sample.fill_out_data_quantities()
        self.data_sample_list.append(error_sample)
        return
        
        
    def process_data_samples_elber(self, pre_equilibrium_approx=False):
        """
        Since the global, system-side statistics have been gathered, 
        compute the thermodynamic and kinetic quantities and their
        uncertainties. Applies to systems using Elber milestoning.
        """
        if self.model.k_on_info is not None:
            self.main_data_sample.parse_browndye_results()
        self.main_data_sample.compute_rate_matrix()
        #self.main_data_sample.Q = common_analyze.minor2d(
        #    self.main_data_sample.Q, bulkstate, bulkstate)
        #self.main_data_sample.K = common_analyze.minor2d(
        #    self.main_data_sample.K, bulkstate, bulkstate)
        self.main_data_sample.calculate_thermodynamics()
        self.main_data_sample.calculate_kinetics(pre_equilibrium_approx)
        """
        error_sample = self.data_sample_list[0]
        error_sample.compute_rate_matrix()
        #self.main_data_sample.Q = common_analyze.minor2d(
        #    self.main_data_sample.Q, bulkstate, bulkstate)
        #self.main_data_sample.K = common_analyze.minor2d(
        #    self.main_data_sample.K, bulkstate, bulkstate)
        error_sample.calculate_thermodynamics()
        error_sample.calculate_kinetics(pre_equilibrium_approx)
        if self.num_error_samples > 0:
            p_i_error, free_energy_profile_err, MFPTs_error, k_off_error, \
                k_ons_error = error_sample.monte_carlo_milestoning_error(
                    num=self.num_error_samples,
                    pre_equilibrium_approx=pre_equilibrium_approx)
        """
        self.p_i = self.main_data_sample.p_i
        self.free_energy_profile = self.main_data_sample.free_energy_profile
        self.MFPTs = self.main_data_sample.MFPTs
        self.k_off = self.main_data_sample.k_off
        self.k_ons = self.main_data_sample.k_ons
        if self.num_error_samples > 0:
            p_i_error, free_energy_profile_err, MFPTs_error, k_off_error, \
                k_ons_error = self.main_data_sample.monte_carlo_milestoning_error(
                    num=self.num_error_samples,
                    pre_equilibrium_approx=pre_equilibrium_approx)
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
    
    def process_data_samples(self, pre_equilibrium_approx=False):
        """
        Based on the type of milestoning, use the data samples to 
        compute thermo and kinetics quantities and their uncertainties.
        """
        if self.model.get_type() == "mmvt":
            self.process_data_samples_mmvt(pre_equilibrium_approx)
        elif self.model.get_type() == "elber":
            self.process_data_samples_elber(pre_equilibrium_approx)
        return
    
    def resample_k_N_R_T(self, N_alpha_beta, N_i_j_alpha, R_i_alpha_total,
                         R_i_alpha_average, R_i_alpha_std_dev, R_i_alpha_count,
                         T_alpha_total, T_alpha_average, T_alpha_std_dev, 
                         T_alpha_count):
        """
        Create data samples from a distribution for computing the
        uncertainties of the thermo and kinetics.
        """
        sampled_k_alpha_beta = {}
        sampled_T_alpha_total = []
        sampled_R_i_alpha_total = []
        for alpha, anchor in enumerate(self.model.anchors):
            if anchor.bulkstate:
                continue
            element_R_i_alpha_total = {}
            for key in R_i_alpha_total[alpha]:
                n_R_i = R_i_alpha_count[alpha][key]
                if n_R_i != 0:
                    R_i_total_std_dev = R_i_alpha_std_dev[alpha][key] * np.sqrt(n_R_i)
                    R_fluctuation = np.random.normal(scale=R_i_total_std_dev)
                else:
                    R_fluctuation = 0.0
                element_R_i_alpha_total[key] = abs(R_i_alpha_total[alpha][key] + R_fluctuation)
            
            n_T = T_alpha_count[alpha]
            if n_T != 0:
                T_total_std_dev = T_alpha_std_dev[alpha] * np.sqrt(n_T)
                T_fluctuation = np.random.normal(scale=T_total_std_dev)
            else:
                T_fluctuation = 0.0
            element_T_alpha_total = abs(T_alpha_total[alpha] + T_fluctuation)
            # This way preserves T = sum of R_i
            sampled_T_alpha_total.append(np.sum(list(element_R_i_alpha_total.values())))
            # In contrast, this way samples T and R_i independently
            #sampled_T_alpha_total.append(element_T_alpha_total)
            sampled_R_i_alpha_total.append(element_R_i_alpha_total)
            
            for beta, anchor2 in enumerate(self.model.anchors):
                key = (alpha, beta)
                sampled_k_alpha_beta[key] = N_alpha_beta[key] / sampled_T_alpha_total[alpha]
            
        sampled_N_alpha_beta = N_alpha_beta
        sampled_N_i_j_alpha = N_i_j_alpha
        
        return sampled_k_alpha_beta, sampled_N_i_j_alpha, \
            sampled_R_i_alpha_total, sampled_T_alpha_total
    
    def print_results(self):
        """Print all results of the analysis calculation."""
        print("Printing results from MMVT SEEKR calculation")
        print("k_off (1/s):", common_analyze.pretty_string_value_error(
            self.k_off, self.k_off_error))
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
                diss_constant_err = diss_constant * common_analyze.quadriture(
                    k_on_err/k_on, self.k_off_error/self.k_off)
                delta_G_err = diss_constant_err*common_analyze.GAS_CONSTANT\
                    *self.model.temperature/diss_constant
            else:
                print("  k_on (1/s * 1/M) to state", key, ":", 
                      common_analyze.pretty_string_value_error(k_on, None))
                diss_constant_err = None
            
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
                print("  MFPT from state", state1, "to state", state2, ":",
                      common_analyze.pretty_string_value_error(
                          float(self.MFPTs[key]*1.0e-12), 
                          float(self.MFPTs_error[key]*1.0e-12)))
            else:
                print("  MFPT from state", state1, "to state", state2, ":",
                      common_analyze.pretty_string_value_error(
                          float(self.MFPTs[key]*1.0e-12), None))
        return
                
    def save_plots(self, image_directory):
        """
        Save a potentially useful series of plots of some quantities
        obtained during the analysis.
        
        TODO: interact with model, because the way these plots are saved
        depends on the structure of the CVs.
        """
        
        anchor_indices = np.zeros(len(self.model.anchors), dtype=np.int8)
        for i, anchor in enumerate(self.model.anchors):
            anchor_indices[i] = anchor.index
        milestone_indices = np.zeros(self.p_i.shape[0], dtype=np.int8)
        for i in range(self.p_i.shape[0]):
            milestone_indices[i] = i
        # save pi_alpha
        if self.model.get_type() == "mmvt":
            pi_fig, ax = plt.subplots()
            plt.errorbar(anchor_indices, self.pi_alpha, yerr=self.pi_alpha_error, 
                         ecolor="k", capsize=2)
            #ax.plot(anchor_indices, self.pi_alpha, linestyle='-', 
            #        marker="o", markersize = 1)
            plt.ylabel("\u03C0_\u03B1")
            plt.xlabel("anchors")
            pi_fig.savefig(os.path.join(image_directory, "pi_alpha.png"))
            
        # save p_i
        pi_fig, ax = plt.subplots()
        plt.errorbar(milestone_indices, self.p_i, yerr=self.p_i_error, 
                     ecolor="k", capsize=2)
        plt.ylabel("p_i")
        plt.xlabel("milestones")
        pi_fig.savefig(os.path.join(image_directory, "p_i.png"))
        # save free energy profile
        pi_fig, ax = plt.subplots()
        plt.errorbar(milestone_indices, self.free_energy_profile, 
                 yerr=self.free_energy_profile_err, ecolor="k", capsize=2)
        plt.ylabel("\u0394G(milestone) (kcal/mol)")
        plt.xlabel("milestones")
        pi_fig.savefig(os.path.join(image_directory, "free_energy_profile.png"))
        return

def analyze(model, force_warning=False, num_error_samples=1000, 
            pre_equilibrium_approx=False, skip_checks=False):
    """Perform all the analysis steps at once."""
    analysis = Analysis(model, force_warning, num_error_samples)
    analysis.extract_data()
    if not skip_checks:
        analysis.check_extraction()
    analysis.fill_out_data_samples()
    analysis.process_data_samples(pre_equilibrium_approx)
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
        default=1000, type=int, help="Specify the number of error samples" \
        " to generate for estimating error/uncertainty of computed "\
        "values. Default: 1000")
    argparser.add_argument(
        "-p", "--pre_equilibrium_approx", dest="pre_equilibrium_approx", 
        default=False, help="This option uses the pre-equilibrium "\
        "approximation when computing system kinetics. This setting may "\
        "be desirable for very long-timescale kinetic processes, which might "\
        "cause the poor matrix conditioning in the milestoning rate matrix, "\
        "causing the typical SEEKR2 analysis approach to fail.", 
        action="store_true")
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
    
    args = argparser.parse_args() # parse the args into a dictionary
    args = vars(args)
    xmlfile = args["model_file"]
    force_warning = args["force_warning"]
    num_error_samples = args["num_error_samples"]
    pre_equilibrium_approx = args["pre_equilibrium_approx"]
    image_directory = args["image_directory"]
    skip_checks = args["skip_checks"]
    
    model = base.Model()
    model.deserialize(xmlfile)
    if model.anchor_rootdir == ".":
        model_dir = os.path.dirname(xmlfile)
        model.anchor_rootdir = os.path.abspath(model_dir)
    
    if image_directory is None:
        image_directory = os.path.join(model.anchor_rootdir, 
                                       common_analyze.DEFAULT_IMAGE_DIR)
    if not os.path.exists(image_directory):
        os.mkdir(image_directory)
        
    if not skip_checks:
        check.check_post_simulation_all(model, long_check=True)
    
    analysis = analyze(model, force_warning=force_warning, 
            num_error_samples=num_error_samples, 
            pre_equilibrium_approx=pre_equilibrium_approx, 
            skip_checks=skip_checks)
    
    analysis.print_results()
    
    print("All plots being saved to:", image_directory)
    analysis.save_plots(image_directory)