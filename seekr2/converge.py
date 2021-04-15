
"""
Functions and objects for extracting and predicting convergence of 
SEEKR simulation outputs and results.
"""

import os
import argparse

import seekr2.modules.common_base as base
import seekr2.modules.common_analyze as common_analyze
import seekr2.modules.common_converge as common_converge

N_IJ_DIR = "N_ij/"
R_I_DIR = "R_i/"
MIN_PLOT_NORM = 1e-8
WINDOW_SIZE = 30

def converge(model, k_on_state, image_directory):
    """
    Perform a generic analysis of convergence of emergenct quantities
    such as N_ij, R_i, k_off, and k_on.
    """
    k_on_conv, k_off_conv, N_ij_conv, R_i_conv, conv_intervals, \
        timestep_in_ns, data_sample_list \
        = common_converge.check_milestone_convergence(
            mymodel, k_on_state=k_on_state, 
            pre_equilibrium_approx=pre_equilibrium_approx)
    
    k_off_fig, ax = common_converge.plot_scalar_conv(
        k_off_conv, conv_intervals, label="k_{off} (s^{-1})", 
        timestep_in_ns=timestep_in_ns)
    k_on_fig, ax = common_converge.plot_scalar_conv(
        k_on_conv, conv_intervals, label="k_{on} (s^{-1} M^{-1})", 
        timestep_in_ns=timestep_in_ns)
    N_ij_fig_list, ax, N_ij_title_list, N_ij_name_list \
        = common_converge.plot_dict_conv(
        N_ij_conv, conv_intervals, label_base="N", 
        timestep_in_ns=timestep_in_ns)
    R_i_fig_list, ax, R_i_title_list, R_i_name_list \
        = common_converge.plot_dict_conv(
        R_i_conv, conv_intervals, label_base="R", 
        timestep_in_ns=timestep_in_ns)
    
    k_off_fig.savefig(os.path.join(image_directory, "k_off_convergence.png"))
    k_on_fig.savefig(os.path.join(image_directory, "k_on_convergence.png"))
    for i, (name, title) in enumerate(zip(N_ij_name_list, N_ij_title_list)):
        N_ij_fig = N_ij_fig_list[i]
        N_ij_filename = os.path.join(image_directory, N_IJ_DIR, name+".png")
        N_ij_fig.savefig(N_ij_filename)
    for i, (name, title) in enumerate(zip(R_i_name_list, R_i_title_list)):
        R_i_fig = R_i_fig_list[i]
        R_i_filename = os.path.join(image_directory, R_I_DIR, name+".png")
        R_i_fig.savefig(R_i_filename)
    
    print("All plots have been saved to:", image_directory)
    
    return data_sample_list

def print_convergence_results(model, convergence_results, cutoff):
    """
    Print the results of a convergence test.
    """
    for alpha, anchor in enumerate(model.anchors):
        if anchor.bulkstate:
            continue
        if convergence_results[alpha] < cutoff:
            is_converged = True
        else:
            is_converged = False
        anchor_string = "Anchor {}: ".format(anchor.index) \
            +"Convergence value: " \
            +"{}. ".format(convergence_results[alpha]) \
            +"Converged? {}".format(is_converged)
        print(anchor_string)

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument(
        "input_file", metavar="INPUT_FILE", type=str, 
        help="name of input file for OpenMMVT calculation. This would be the "\
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
        "-c", "--cutoff", dest="cutoff", default=0.01, type=float, 
        help="")
    argparser.add_argument(
        "-p", "--pre_equilibrium_approx", dest="pre_equilibrium_approx", 
        default=False, help="Optionally use the pre-equilibrium approximation"\
        "when computing system kinetics. This setting may be desirable for "\
        "very long-timescale kinetic processes, which would cause the typical"\
        "SEEKR2 analysis approach to fail.", action="store_true")
    
    args = argparser.parse_args() # parse the args into a dictionary
    args = vars(args)
    input_file = args["input_file"]
    k_on_state = args["k_on_state"]
    image_directory = args["image_directory"]
    cutoff = args["cutoff"]
    pre_equilibrium_approx = args["pre_equilibrium_approx"]
    
    mymodel = base.Model()
    mymodel.deserialize(input_file)
    if mymodel.anchor_rootdir == ".":
        model_dir = os.path.dirname(input_file)
        mymodel.anchor_rootdir = os.path.abspath(model_dir)
    
    if image_directory is None:
        image_directory = os.path.join(mymodel.anchor_rootdir, 
                                       common_analyze.DEFAULT_IMAGE_DIR)
    
    if not os.path.exists(image_directory):
        os.mkdir(image_directory)
    
    if not os.path.exists(os.path.join(image_directory, N_IJ_DIR)):
        os.mkdir(os.path.join(image_directory, N_IJ_DIR))
        
    if not os.path.exists(os.path.join(image_directory, R_I_DIR)):
        os.mkdir(os.path.join(image_directory, R_I_DIR))
    
    if mymodel.k_on_info is None:
        assert k_on_state is None, "--k_on_state cannot be defined for "\
            "this model. The model was not initialized to compute k-on "\
            "quantities."
        print("This model has no information for computing k-on - therefore "\
              "k-on convergence will be skipped.")
    else:
        end_milestones = []
        for anchor in mymodel.anchors:
            if anchor.endstate:
                for milestone_id in anchor.get_ids():
                    end_milestones.append(milestone_id)
        assert len(end_milestones) > 0, "No end-state milestones for this "\
            "model: k-on convergence cannot be computed."
        print("All available options for --k_on_state include:", end_milestones)
        if k_on_state is None:
            k_on_state = end_milestones[0]
            print("No milestone has been provided for k_on_state. "\
                  "The milestone %d has been chosen by default." % k_on_state)
        else:
            assert k_on_state in end_milestones, "The provided "\
                "milestone of %d for k_on_state is not available." % k_on_state
    
    data_sample_list = converge(mymodel, k_on_state, image_directory)
        
    rmsd_convergence_results = common_converge.calc_mmvt_RMSD_conv_amount(
        mymodel, data_sample_list)
    print_convergence_results(mymodel, rmsd_convergence_results, cutoff)