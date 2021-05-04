"""
Run entire SEEKR2 calculations to test for problems in the 
pipeline.
"""
import os
import sys
import time
import tempfile

import seekr2.tests.create_model_input as create_model_input
import seekr2.prepare as prepare
import seekr2.modules.check as check
import seekr2.run as run
import seekr2.converge as converge
import seekr2.modules.common_converge as common_converge
import seekr2.analyze as analyze

def run_short_ci(model_input):
    start_dir = os.getcwd()
    model, xml_path = prepare.generate_openmmvt_model_and_filetree(
        model_input, force_overwrite=False)
    model_dir = os.path.dirname(xml_path)
    model.anchor_rootdir = os.path.abspath(model_dir)
    check.check_pre_simulation_all(model)
    run.run(model, "any", min_b_surface_simulation_length=1000,
        min_bd_milestone_simulation_length=100, 
        max_b_surface_trajs_to_extract=10)
    data_sample_list = converge.converge(model, k_on_state=0)
    rmsd_convergence_results = common_converge.calc_RMSD_conv_amount(
        model, data_sample_list)
    transition_minima, transition_details \
        = common_converge.calc_transition_steps(
        model, data_sample_list[-1])
    converge.print_convergence_results(
        model, rmsd_convergence_results, cutoff=0.1, 
        transition_results=transition_details, 
        minimum_anchor_transitions=10, 
        bd_transition_counts=data_sample_list[-1].bd_transition_counts)
    analysis = analyze.analyze(model)
    analysis.print_results()
    os.chdir(start_dir)
    return

def run_generic_hostguest_ci():
    with tempfile.TemporaryDirectory() as temp_dir:
        host_guest_model_input \
            = create_model_input.create_host_guest_mmvt_model_input(temp_dir)
        run_short_ci(host_guest_model_input)
        
    return

def run_generic_namd_hostguest_ci():
    with tempfile.TemporaryDirectory() as temp_dir:
        host_guest_model_input \
            = create_model_input.create_host_guest_mmvt_model_input(temp_dir)
        host_guest_model_input.md_program = "namd"
        for input_anchor in host_guest_model_input.cv_inputs[0].input_anchors:
            input_anchor.starting_amber_params.prmtop_filename \
                = "../data/hostguest_files/hostguest_for_NAMD.parm7"
        run_short_ci(host_guest_model_input)
        
    return

if __name__ == "__main__":
    if len(sys.argv) <= 1:
        argument = "short"
    else:
        argument = sys.argv[1]
    
    starttime = time.time()
    if argument == "short":
        run_generic_hostguest_ci()
    elif argument == "namd":
        run_generic_namd_hostguest_ci()
    elif argument == "long":
        run_generic_hostguest_ci()
        run_generic_namd_hostguest_ci()
    
    print("Time elapsed: {:.3f}".format(time.time() - starttime))
    print("Continuous Integration Tests Passed Successfully.")