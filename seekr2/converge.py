"""
converge.py

Functions and objects for extracting and predicting convergence of 
SEEKR2 simulation outputs and results.
"""

import os
import argparse

import matplotlib.pyplot as plt

import seekr2.modules.common_base as base
import seekr2.modules.common_analyze as common_analyze
import seekr2.modules.common_converge as common_converge

AVG_TRANSITION_TIME_MINIMUM = 0.1 # ps

def converge(model, k_on_state=None, image_directory=None, verbose=False,
             long_converge=True):
    """
    Perform all convergence steps: a generic analysis of convergence 
    of quantities such as N_ij, R_i, k_off, and k_on.
    """
    curdir = os.getcwd()
    k_on_conv, k_off_conv, N_alpha_beta_conv, T_alpha_conv, \
        k_alpha_beta_conv, pi_alpha_conv, N_ij_alpha_conv, R_i_alpha_conv, \
        N_ij_conv, R_i_conv, max_step_list, timestep_in_ns, data_sample_list, \
        times_dict \
        = common_converge.check_milestone_convergence(
            model, k_on_state=k_on_state, verbose=verbose, 
            long_converge=long_converge)
    
    if image_directory is None or image_directory == "" or not long_converge:
        return data_sample_list, times_dict
    
    k_off_fig, ax = common_converge.plot_scalar_conv(
        k_off_conv, max_step_list, title="$k_{off}$ Convergence", 
        label="k_{off} (s^{-1})", timestep_in_ns=timestep_in_ns)
    if k_off_fig is not None:
        k_off_fig.savefig(os.path.join(image_directory, "k_off_convergence.png"))
        plt.close(k_off_fig)
        
    k_on_fig, ax = common_converge.plot_scalar_conv(
        k_on_conv, max_step_list, title="$k_{on}$ Convergence", 
        label="k_{on} (s^{-1} M^{-1})", timestep_in_ns=timestep_in_ns)
    if k_on_fig is not None:
        k_on_fig.savefig(os.path.join(image_directory, "k_on_convergence.png"))
        plt.close(k_on_fig)
    
    N_alpha_beta_fig_list, ax, N_alpha_beta_title_list, N_alpha_beta_name_list \
        = common_converge.plot_dict_conv(
        N_alpha_beta_conv, max_step_list, label_base="N", unit="",
        timestep_in_ns=timestep_in_ns, y_axis_logarithmic=False, 
        title_suffix="", draw_double=True)
    for i, name in enumerate(N_alpha_beta_name_list):
        N_alpha_beta_fig = N_alpha_beta_fig_list[i]
        if N_alpha_beta_fig is not None:
            N_alpha_beta_filename = os.path.join(
                image_directory, common_analyze.N_ALPHA_BETA_DIR, name+".png")
            N_alpha_beta_fig.savefig(N_alpha_beta_filename)
            plt.close(N_alpha_beta_fig)
    
    T_alpha_fig_list, ax, T_alpha_title_list, T_alpha_name_list \
        = common_converge.plot_dict_conv(
        T_alpha_conv, max_step_list, label_base="T", unit="ps",
        timestep_in_ns=timestep_in_ns, y_axis_logarithmic=False, 
        title_suffix="", draw_double=True)
    
    for i, name in enumerate(T_alpha_name_list):
        T_alpha_fig = T_alpha_fig_list[i]
        if T_alpha_fig is not None:
            T_alpha_filename = os.path.join(
                image_directory, common_analyze.T_ALPHA_DIR, name+".png")
            T_alpha_fig.savefig(T_alpha_filename)
            plt.close(T_alpha_fig)
    
    k_alpha_beta_fig_list, ax, k_alpha_beta_title_list, k_alpha_beta_name_list \
        = common_converge.plot_dict_conv(
        k_alpha_beta_conv, max_step_list, label_base="k", unit="",
        timestep_in_ns=timestep_in_ns, title_suffix="")
    for i, name in enumerate(k_alpha_beta_name_list):
        k_alpha_beta_fig = k_alpha_beta_fig_list[i]
        if k_alpha_beta_fig is not None:
            k_alpha_beta_filename = os.path.join(
                image_directory, common_analyze.K_ALPHA_BETA_DIR, name+".png")
            k_alpha_beta_fig.savefig(k_alpha_beta_filename)
            plt.close(k_alpha_beta_fig)
    
    pi_alpha_fig_list, ax, pi_alpha_title_list, pi_alpha_name_list \
        = common_converge.plot_dict_conv(
        pi_alpha_conv, max_step_list, label_base="\pi", unit="",
        timestep_in_ns=timestep_in_ns, title_suffix="", name_base="pi",
        draw_double=False)
    for i, name in enumerate(pi_alpha_name_list):
        pi_alpha_fig = pi_alpha_fig_list[i]
        if pi_alpha_fig is not None:
            pi_alpha_filename = os.path.join(
                image_directory, common_analyze.PI_ALPHA_DIR, name+".png")
            pi_alpha_fig.savefig(pi_alpha_filename)
            plt.close(pi_alpha_fig)
            
    N_ij_alpha_fig_list, ax, N_ij_alpha_title_list, N_ij_alpha_name_list \
        = common_converge.plot_dict_conv(
        N_ij_alpha_conv, max_step_list, label_base="N", unit="",
        timestep_in_ns=timestep_in_ns, y_axis_logarithmic=False, 
        title_suffix="", draw_double=True)
    for i, name in enumerate(N_ij_alpha_name_list):
        N_ij_alpha_fig = N_ij_alpha_fig_list[i]
        if N_ij_alpha_fig is not None:
            N_ij_alpha_filename = os.path.join(
                image_directory, common_analyze.N_IJ_ALPHA_DIR, name+".png")
            N_ij_alpha_fig.savefig(N_ij_alpha_filename)
            plt.close(N_ij_alpha_fig)
    
    R_i_alpha_fig_list, ax, R_i_alpha_title_list, R_i_alpha_name_list \
        = common_converge.plot_dict_conv(
        R_i_alpha_conv, max_step_list, label_base="R", unit="ps",
        timestep_in_ns=timestep_in_ns, title_suffix="", 
        y_axis_logarithmic=False, draw_double=True)
    for i, name in enumerate(R_i_alpha_name_list):
        R_i_alpha_fig = R_i_alpha_fig_list[i]
        if R_i_alpha_fig is not None:
            R_i_alpha_filename = os.path.join(
                image_directory, common_analyze.R_I_ALPHA_DIR, name+".png")
            R_i_alpha_fig.savefig(R_i_alpha_filename)
            plt.close(R_i_alpha_fig)
    
    N_ij_fig_list, ax, N_ij_title_list, N_ij_name_list \
        = common_converge.plot_dict_conv(
        N_ij_conv, max_step_list, label_base="N", unit="",
        timestep_in_ns=timestep_in_ns, y_axis_logarithmic=False, 
        draw_double=True)
    for i, name in enumerate(N_ij_name_list):
        N_ij_fig = N_ij_fig_list[i]
        if N_ij_fig is not None:
            N_ij_filename = os.path.join(
                image_directory, common_analyze.N_IJ_DIR, name+".png")
            N_ij_fig.savefig(N_ij_filename)
            plt.close(N_ij_fig)
    
    R_i_fig_list, ax, R_i_title_list, R_i_name_list \
        = common_converge.plot_dict_conv(
        R_i_conv, max_step_list, label_base="R", unit="ps",
        timestep_in_ns=timestep_in_ns, y_axis_logarithmic=False, 
        title_suffix="", draw_double=True)
    for i, name in enumerate(R_i_name_list):
        R_i_fig = R_i_fig_list[i]
        if R_i_fig is not None:
            R_i_filename = os.path.join(
                image_directory, common_analyze.R_I_DIR, name+".png")
            R_i_fig.savefig(R_i_filename)
            plt.close(R_i_fig)
    
    print("All plots have been saved to:", image_directory)
    os.chdir(curdir)
    return data_sample_list, times_dict

def print_convergence_results(model, convergence_results, cutoff, 
                              transition_results, transition_time_results,
                              minimum_anchor_transitions, 
                              bd_transition_counts={}, times_dict={}):
    """
    Print the results of a convergence test.
    """
    print("Molecular dynamics results:")
    for alpha, anchor in enumerate(model.anchors):
        if anchor.bulkstate:
            continue
        if convergence_results[alpha] < cutoff \
                and min(transition_results[alpha].values()) \
                > minimum_anchor_transitions:
            is_converged = True
        else:
            is_converged = False
            
        transition_detail = transition_results[alpha]
        transition_avg_times = transition_time_results[alpha]
        transition_string = ""
        time_string = ""
        for key in transition_detail:
            transition_string += " {}->{} : {},".format(
                key[0], key[1], transition_detail[key])
        
        warnstr = ""
        for key in transition_avg_times:
            if transition_avg_times[key] < AVG_TRANSITION_TIME_MINIMUM \
                    and len(transition_detail) > 1:
                warnstr = "\n    WARNING: transition time(s) for this "\
                    "milestone are below the recommended minimum "\
                    "({} ps) ".format(AVG_TRANSITION_TIME_MINIMUM)
            time_string += " {} : {:.3f},".format(
                key, transition_avg_times[key])
        if anchor.index in times_dict:
            anchor_time = times_dict[anchor.index]
        else:
            anchor_time = 0.0
        anchor_string = " - Anchor {}: ".format(anchor.index) \
            +"\n     Milestone transitions:{} ".format(transition_string) \
            +"\n     Milestone avg. transition time (ps):{}"\
                .format(time_string) + warnstr \
            +"\n     Time simulated in anchor (ps): {:.3f}"\
                .format(anchor_time) \
            +"\n     Convergence value: " \
            +"{:.4e}. ".format(convergence_results[alpha]) \
            +"\n     Converged? {}".format(is_converged)
        print(anchor_string)
        
    if "b_surface" in bd_transition_counts:
        print("Brownian dynamics results:")
        print(" - b-surface reaction counts:")
        for bd_milestone_id in bd_transition_counts["b_surface"]:
            if bd_milestone_id in ["escaped", "stuck", "total"]:
                print("     {}:".format(bd_milestone_id), 
                      bd_transition_counts["b_surface"][bd_milestone_id])
            else:
                print("     to milestone {}:".format(bd_milestone_id), 
                    bd_transition_counts["b_surface"][bd_milestone_id])
        
    else:
        print("No BD results found.")
    
    return

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument(
        "model_file", metavar="MODEL_FILE", type=str, 
        help="name of model file for SEEKR2 calculation. This would be the "\
        "XML file generated in the prepare stage.")
    argparser.add_argument(
        "-s", "--k_on_state", dest="k_on_state", default=None, type=int,
        help="Define the bound state that will be used to compute k-on. If "\
        "left blank, and k-on statistics exist, then the first end state will "\
        "be chosen by default.")
    argparser.add_argument(
        "-d", "--image_directory", dest="image_directory", 
        default=None, type=str,
        help="Define the directory where all plots and images will be saved. "\
            "By default, graphics will be saved to the "\
            "'%s' directory in the model's anchor root directory."\
            % common_analyze.DEFAULT_IMAGE_DIR)
    argparser.add_argument(
        "-c", "--cutoff", dest="cutoff", default=0.1, type=float, 
        help="The minimum convergence that must be achieved before concluding "\
        "that the calculations have converged for a given anchor. "\
        "Default: 0.1.")
    argparser.add_argument(
        "-m", "--minimum_anchor_transitions", dest="minimum_anchor_transitions",
        default=100, type=int, help="Enter a minimum number of transitions "\
        "that must be observed per milestone in a given anchor as a criteria "\
        "for the simulations. Default: 100")
    argparser.add_argument(
        "-l", "--long_converge", dest="long_converge", default=False, 
        help="Whether to run a full, long convergence analysis. If set to "\
        "True, 100 intervals will be sampled for convergence values and "\
        "plots will be generated. Default: False.", action="store_true")
    
    args = argparser.parse_args() # parse the args into a dictionary
    args = vars(args)
    model_file = args["model_file"]
    k_on_state = args["k_on_state"]
    image_directory = args["image_directory"]
    cutoff = args["cutoff"]
    minimum_anchor_transitions = args["minimum_anchor_transitions"]
    long_converge = args["long_converge"]
    
    model = base.load_model(model_file)
    
    image_directory = common_analyze.make_image_directory(
        model, image_directory)
    
    if not model.using_bd():
        assert k_on_state is None, "--k_on_state cannot be defined for "\
            "this model. The model was not initialized to compute k-on "\
            "quantities."
        print("This model has no information for computing k-on - therefore "\
              "k-on convergence will be skipped.")
    else:
        end_milestones = []
        for anchor in model.anchors:
            if anchor.endstate:
                for milestone_id in anchor.get_ids():
                    if model.get_type() == "mmvt":
                        end_milestones.append(milestone_id)
                    else:
                        if anchor.milestones[milestone_id].is_source_milestone:
                            end_milestones.append(milestone_id)
                    continue
        assert len(end_milestones) > 0, "No end-state milestones for this "\
            "model: k-on convergence cannot be computed."
        print("All available options for --k_on_state include:", end_milestones)
        if k_on_state is None:
            k_on_state = end_milestones[0]
            print("No BD milestone has been provided for k_on_state. "\
                  "The BD milestone %d has been chosen by default." \
                  % k_on_state)
        else:
            assert k_on_state in end_milestones, "The provided "\
                "BD milestone of %d for k_on_state is not available." \
                % k_on_state
    
    data_sample_list, times_dict = converge(model, k_on_state, image_directory, 
                                verbose=True, long_converge=long_converge)
    main_data_sample = data_sample_list[-1]
    rmsd_convergence_results = common_converge.calc_RMSD_conv_amount(
        model, data_sample_list)
    transition_minima, transition_prob_results, transition_time_results \
        = common_converge.calc_transition_steps(
        model, main_data_sample)
    bd_transition_counts = main_data_sample.bd_transition_counts
    print_convergence_results(model, rmsd_convergence_results, cutoff, 
                              transition_prob_results, transition_time_results,
                              minimum_anchor_transitions,
                              bd_transition_counts, times_dict)