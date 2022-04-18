"""
Test the running, restart, state-saving, and other capabilities of
the run.py program
"""

import os
import glob
import tempfile
import signal

import mdtraj
import seekr2.modules.mmvt_sim_openmm as mmvt_sim_openmm
import seekr2.modules.common_converge as common_converge
import seekr2.modules.runner_openmm as runner_openmm
import seekr2.modules.runner_namd as runner_namd
import seekr2.run as run

def get_trajectory_length(model, anchor, top_filename=None):
    dcd_glob = os.path.join(model.anchor_rootdir, anchor.directory, 
        anchor.production_directory, "*.dcd")
    dcd_file_list = glob.glob(dcd_glob)
    # TODO: won't work for any type of system besides a toy
    if top_filename is None:
        top_filename = os.path.join(model.anchor_rootdir, anchor.directory, 
            anchor.building_directory, "toy.pdb")
    traj = mdtraj.load(dcd_file_list, top=top_filename)
    return len(traj), len(dcd_file_list)

def get_checkpoint_step(model, anchor):
    dummy_file = tempfile.NamedTemporaryFile()
    sim_openmm_obj = mmvt_sim_openmm.create_sim_openmm(
                        model, anchor, dummy_file.name, frame=0, 
                        load_state_file=None)
    simulation = sim_openmm_obj.simulation
    output_directory = os.path.join(
        model.anchor_rootdir, anchor.directory, 
        anchor.production_directory)
    restart_checkpoint_basename = runner_openmm.RESTART_CHECKPOINT_FILENAME
    restart_checkpoint_glob = os.path.join(
        output_directory, restart_checkpoint_basename+"*")
    restart_checkpoint_list = glob.glob(restart_checkpoint_glob)
    if len(restart_checkpoint_list) > 0:
        simulation.loadCheckpoint(restart_checkpoint_list[0])
        currentStep = simulation.context.getState().getStepCount()
    else:
        currentStep = 0
    dummy_file.close()
    return currentStep

def test_choose_next_simulation_browndye2(host_guest_mmvt_model):
    bd_milestone_info_to_run_unsorted = run.choose_next_simulation_browndye2(
        host_guest_mmvt_model, "b_surface", 100, True)
    assert bd_milestone_info_to_run_unsorted[0][0] == 100
    assert bd_milestone_info_to_run_unsorted[0][1] == 0
    assert bd_milestone_info_to_run_unsorted[0][2] == "b_surface"
    assert bd_milestone_info_to_run_unsorted[0][3] == False
    assert bd_milestone_info_to_run_unsorted[0][4] == 100
    run.run(host_guest_mmvt_model, "b_surface", force_overwrite=True)
    
    # test restart
    host_guest_mmvt_model.k_on_info.b_surface_num_trajectories = 200
    bd_milestone_info_to_run_unsorted = run.choose_next_simulation_browndye2(
        host_guest_mmvt_model, "b_surface", 200, False)
    assert bd_milestone_info_to_run_unsorted[0][0] == 100
    #assert bd_milestone_info_to_run_unsorted[0][1] == 100
    assert bd_milestone_info_to_run_unsorted[0][2] == "b_surface"
    assert bd_milestone_info_to_run_unsorted[0][3] == True
    assert bd_milestone_info_to_run_unsorted[0][4] == 100
    
    run.run(host_guest_mmvt_model, "b_surface", force_overwrite=True)
    
    # test criteria satisfied
    bd_milestone_info_to_run_unsorted = run.choose_next_simulation_browndye2(
        host_guest_mmvt_model, "b_surface", 200, False)
    assert len(bd_milestone_info_to_run_unsorted) == 0
    
    # test min_b_surface_encounters
    bd_transition_counts = common_converge.get_bd_transition_counts(
        host_guest_mmvt_model)
    min_b_surface_encounters = bd_transition_counts["b_surface"][11] + 1
    bd_milestone_info_to_run_unsorted = run.choose_next_simulation_browndye2(
        host_guest_mmvt_model, "b_surface", 200, False, 
        min_b_surface_encounters)
    
    assert bd_milestone_info_to_run_unsorted[0][0] \
        == run.B_SURFACE_CONVERGENCE_INTERVAL - 200
    #assert bd_milestone_info_to_run_unsorted[0][1] == 200
    assert bd_milestone_info_to_run_unsorted[0][2] == "b_surface"
    assert bd_milestone_info_to_run_unsorted[0][3] == True
    assert bd_milestone_info_to_run_unsorted[0][4] == 200
    
    return

def test_choose_next_simulation_openmm(toy_mmvt_model):
    anchor_info_to_run = run.choose_next_simulation_openmm(
        toy_mmvt_model, "2", min_total_simulation_length=10000, 
        max_total_simulation_length=100000, convergence_cutoff=None, 
        minimum_anchor_transitions=None, force_overwrite=True, 
        umbrella_restart_mode=False, load_state_file=None)
    assert anchor_info_to_run[0][0] == 10000
    assert anchor_info_to_run[0][1] == 0
    assert anchor_info_to_run[0][2] == 2
    assert anchor_info_to_run[0][3] == False
    assert anchor_info_to_run[0][4] == 10000
    assert anchor_info_to_run[0][5] == None
    
    # Test swarm
    anchor_info_to_run = run.choose_next_simulation_openmm(
        toy_mmvt_model, "0", min_total_simulation_length=10000, 
        max_total_simulation_length=100000, convergence_cutoff=None, 
        minimum_anchor_transitions=None, force_overwrite=True, 
        umbrella_restart_mode=False, load_state_file=None)
    assert anchor_info_to_run[0][0] == 10000
    assert anchor_info_to_run[0][1] == 0
    assert anchor_info_to_run[0][2] == 0
    assert anchor_info_to_run[0][3] == False
    assert anchor_info_to_run[0][4] == 10000
    assert anchor_info_to_run[0][5] == 0
    assert anchor_info_to_run[1][0] == 10000
    assert anchor_info_to_run[1][1] == 0
    assert anchor_info_to_run[1][2] == 0
    assert anchor_info_to_run[1][3] == False
    assert anchor_info_to_run[1][4] == 10000
    assert anchor_info_to_run[1][5] == 1
    
    # Test swarm with loaded state file
    anchor_info_to_run = run.choose_next_simulation_openmm(
        toy_mmvt_model, "0", min_total_simulation_length=10000, 
        max_total_simulation_length=100000, convergence_cutoff=None, 
        minimum_anchor_transitions=None, force_overwrite=True, 
        umbrella_restart_mode=False, load_state_file="dummy")
    assert anchor_info_to_run[0][5] == None
    
    # Test restart
    toy_mmvt_model.calculation_settings.num_production_steps = 10000
    run.run(toy_mmvt_model, "2", force_overwrite=True)
    anchor_info_to_run = run.choose_next_simulation_openmm(
        toy_mmvt_model, "2", min_total_simulation_length=20000, 
        max_total_simulation_length=1000000, convergence_cutoff=None, 
        minimum_anchor_transitions=None, force_overwrite=False, 
        umbrella_restart_mode=False, load_state_file=None)
    assert anchor_info_to_run[0][0] == 10000
    #assert anchor_info_to_run[0][1] == 1000
    assert anchor_info_to_run[0][2] == 2
    assert anchor_info_to_run[0][3] == True
    assert anchor_info_to_run[0][4] == 20000
    assert anchor_info_to_run[0][5] == None
    
    # test convergence
    anchor_info_to_run = run.choose_next_simulation_openmm(
        toy_mmvt_model, "2", min_total_simulation_length=10000, 
        max_total_simulation_length=1000000, convergence_cutoff=1e-20, 
        minimum_anchor_transitions=None, force_overwrite=False, 
        umbrella_restart_mode=False, load_state_file=None)
    #assert anchor_info_to_run[0][0] == 0
    #assert anchor_info_to_run[0][2] == 2
    #assert anchor_info_to_run[0][3] == True
    #assert anchor_info_to_run[0][4] == 10000 + run.CONVERGENCE_INTERVAL
    
    # test max
    anchor_info_to_run = run.choose_next_simulation_openmm(
        toy_mmvt_model, "2", min_total_simulation_length=10000, 
        max_total_simulation_length=9000, convergence_cutoff=0.000002, 
        minimum_anchor_transitions=None, force_overwrite=False, 
        umbrella_restart_mode=False, load_state_file=None)
    assert len(anchor_info_to_run) == 0
    
    # test minimum anchor transitions
    anchor_info_to_run = run.choose_next_simulation_openmm(
        toy_mmvt_model, "2", min_total_simulation_length=10000, 
        max_total_simulation_length=None, convergence_cutoff=None, 
        minimum_anchor_transitions=10000, force_overwrite=False, 
        umbrella_restart_mode=False, load_state_file=None)
    assert anchor_info_to_run[0][3] == True
    assert anchor_info_to_run[0][4] == 10000 + run.CONVERGENCE_INTERVAL
    
def test_choose_next_simulation_namd(host_guest_mmvt_model_namd):
    anchor_info_to_run = run.choose_next_simulation_namd(
        host_guest_mmvt_model_namd, "2", min_total_simulation_length=10, 
        max_total_simulation_length=10000, convergence_cutoff=None, 
        minimum_anchor_transitions=None, force_overwrite=True)
    assert anchor_info_to_run[0][0] == 10
    assert anchor_info_to_run[0][1] == 0
    assert anchor_info_to_run[0][2] == 2
    assert anchor_info_to_run[0][3] == False
    assert anchor_info_to_run[0][4] == 10
    assert anchor_info_to_run[0][5] == None
    
    # Test restart
    host_guest_mmvt_model_namd.calculation_settings.num_production_steps = 10
    run.run(host_guest_mmvt_model_namd, "2", force_overwrite=True)
    anchor_info_to_run = run.choose_next_simulation_namd(
        host_guest_mmvt_model_namd, "2", min_total_simulation_length=20, 
        max_total_simulation_length=10000, convergence_cutoff=None, 
        minimum_anchor_transitions=None, force_overwrite=False)
    assert anchor_info_to_run[0][0] == 10
    assert anchor_info_to_run[0][2] == 2
    assert anchor_info_to_run[0][3] == True
    assert anchor_info_to_run[0][4] == 20
    assert anchor_info_to_run[0][5] == None
    
    anchor_info_to_run = run.choose_next_simulation_namd(
        host_guest_mmvt_model_namd, "2", min_total_simulation_length=10, 
        max_total_simulation_length=1000000, convergence_cutoff=0.2, 
        minimum_anchor_transitions=None, force_overwrite=False)
    assert anchor_info_to_run[0][0] == 0
    assert anchor_info_to_run[0][2] == 2
    assert anchor_info_to_run[0][3] == True
    assert anchor_info_to_run[0][4] == run.CONVERGENCE_INTERVAL
    
    anchor_info_to_run = run.choose_next_simulation_namd(
        host_guest_mmvt_model_namd, "2", min_total_simulation_length=10, 
        max_total_simulation_length=9, convergence_cutoff=0.2, 
        minimum_anchor_transitions=None, force_overwrite=False)
    assert len(anchor_info_to_run) == 0
    
    anchor_info_to_run = run.choose_next_simulation_namd(
        host_guest_mmvt_model_namd, "2", min_total_simulation_length=10, 
        max_total_simulation_length=None, convergence_cutoff=None, 
        minimum_anchor_transitions=10, force_overwrite=False)
    assert anchor_info_to_run[0][3] == True
    assert anchor_info_to_run[0][4] == run.CONVERGENCE_INTERVAL
    return

def test_normal_run_openmm(toy_mmvt_model):
    """
    Test a normal toy run - no restarts or anything.
    """
    num_steps = 10000
    dcd_interval = toy_mmvt_model.calculation_settings\
        .trajectory_reporter_interval
    num_dcd_frames = num_steps // dcd_interval
    toy_mmvt_model.calculation_settings.num_production_steps = num_steps
    toy_mmvt_model.openmm_settings.cuda_platform_settings = None
    toy_mmvt_model.openmm_settings.reference_platform = True
    run.run(toy_mmvt_model, "any", force_overwrite=True)
    for anchor in toy_mmvt_model.anchors:
        if anchor.bulkstate:
            continue
        checkpoint_step = get_checkpoint_step(toy_mmvt_model, anchor)
        assert checkpoint_step == num_steps
        
        dcd_length, dcd_file_number = get_trajectory_length(
            toy_mmvt_model, anchor)
        dummy_file = tempfile.NamedTemporaryFile()
        num_swarm = mmvt_sim_openmm.\
            get_starting_structure_num_frames(toy_mmvt_model, anchor, 
                                              dummy_file.name)
        assert dcd_length == num_dcd_frames*num_swarm
        
    return
        
def test_normal_restart_openmm(toy_mmvt_model):
    """
    Test a normal restart - no errors or interruptions.
    """
    # Test DCD frame number and checkpoint step
    first_steps = 10000
    dcd_interval = toy_mmvt_model.calculation_settings\
        .trajectory_reporter_interval
    num_dcd_frames1 = first_steps // dcd_interval
    toy_mmvt_model.openmm_settings.cuda_platform_settings = None
    toy_mmvt_model.openmm_settings.reference_platform = True
    toy_mmvt_model.calculation_settings.num_production_steps = first_steps
    toy_mmvt_model.calculation_settings.trajectory_reporter_interval \
        = first_steps // 10
    run.run(toy_mmvt_model, "1", force_overwrite=True)
    checkpoint_step = get_checkpoint_step(
        toy_mmvt_model, toy_mmvt_model.anchors[1])
    assert checkpoint_step == first_steps
    #dcd_length1, dcd_file_number1 = get_trajectory_length(toy_mmvt_model, 
    #                                    toy_mmvt_model.anchors[1])
    #assert dcd_length1 == num_dcd_frames1
    #assert dcd_file_number1 == 1
    
    second_steps = 20000
    num_dcd_frames2 = second_steps // dcd_interval
    toy_mmvt_model.calculation_settings.num_production_steps = second_steps
    run.run(toy_mmvt_model, "1")
    checkpoint_step = get_checkpoint_step(
        toy_mmvt_model, toy_mmvt_model.anchors[1])
    assert checkpoint_step == second_steps
    #dcd_length2, dcd_file_number2 = get_trajectory_length(toy_mmvt_model, 
    #                                    toy_mmvt_model.anchors[1])
    #assert dcd_length2 == num_dcd_frames2
    #assert dcd_file_number2 == 2
    
    # Test that there are two output files
    out_glob = os.path.join(
        toy_mmvt_model.anchor_rootdir, toy_mmvt_model.anchors[1].directory, 
        toy_mmvt_model.anchors[1].production_directory, 
        toy_mmvt_model.anchors[1].md_output_glob)
    out_file_list = glob.glob(out_glob)
    assert len(out_file_list) == 2
    return

def handler(signum, frame):
    print("Sudden, unexpected interruption!")
    raise Exception("The system experienced an interruption.")

def test_interruption_restart_openmm(toy_mmvt_model):
    """
    Test a restart caused by an error or interruption
    """
    num_steps = 1000000
    dcd_interval = toy_mmvt_model.calculation_settings\
        .trajectory_reporter_interval
    num_dcd_frames = num_steps // dcd_interval
    toy_mmvt_model.calculation_settings.num_production_steps = num_steps
    toy_mmvt_model.calculation_settings.trajectory_reporter_interval \
        = 10000
    toy_mmvt_model.openmm_settings.cuda_platform_settings = None
    toy_mmvt_model.openmm_settings.reference_platform = True
    try:
        signal.signal(signal.SIGALRM, handler)
        signal.alarm(1)
        run.run(toy_mmvt_model, "1", force_overwrite=True)
    except (Exception, SystemExit):
        pass
    
    checkpoint_step = get_checkpoint_step(
        toy_mmvt_model, toy_mmvt_model.anchors[1])
    assert checkpoint_step != num_steps, \
        "Interruption happens too late to catch runner."
    #dcd_length1, dcd_file_number1 = get_trajectory_length(toy_mmvt_model, 
    #                                    toy_mmvt_model.anchors[1])
    #assert dcd_length1 != num_dcd_frames, \
    #    "Interruption happens too late to catch runner."
    #assert dcd_file_number1 == 1
    
    run.run(toy_mmvt_model, "1")
    checkpoint_step = get_checkpoint_step(
        toy_mmvt_model, toy_mmvt_model.anchors[1])
    assert checkpoint_step == num_steps
    #dcd_length2, dcd_file_number2 = get_trajectory_length(toy_mmvt_model, 
    #                                    toy_mmvt_model.anchors[1])
    # NOTE: DCD writes will sometimes have too many frames, clearly this
    #  is caused by an interruption that occurs between the checkpoint and
    #  DCD write. This problem should be addressed, if ever, at the
    #  postprocessing step - the run stage will not handle it - a 
    #  principle which is consistent with other parts of SEEKR.
    #assert dcd_length2 == num_dcd_frames
    #assert dcd_file_number2 == 2
    return

def test_normal_run_namd(host_guest_mmvt_model_namd):
    """
    Test a normal toy run - no restarts or anything, using NAMD.
    """
    num_steps = 10
    dcd_interval = 1
    host_guest_mmvt_model_namd.calculation_settings.num_production_steps \
        = num_steps
    host_guest_mmvt_model_namd.calculation_settings\
        .restart_checkpoint_interval = 1
    host_guest_mmvt_model_namd.calculation_settings\
        .trajectory_reporter_interval = dcd_interval
    num_dcd_frames = num_steps // dcd_interval
    run.run(host_guest_mmvt_model_namd, "0", force_overwrite=True)
    anchor = host_guest_mmvt_model_namd.anchors[0]
    output_glob = os.path.join(
        host_guest_mmvt_model_namd.anchor_rootdir, anchor.directory,
        anchor.production_directory, anchor.md_output_glob)
    glob_list = glob.glob(output_glob)
    assert len(glob_list) == 1
    output_directory = os.path.join(
        host_guest_mmvt_model_namd.anchor_rootdir, anchor.directory, 
        anchor.production_directory)
    restart_checkpoint_filename = os.path.join(
        output_directory, runner_namd.RESTART_CHECKPOINT_FILENAME+".xsc")
    checkpoint_step = runner_namd.read_xsc_step_number(
            restart_checkpoint_filename)
    assert checkpoint_step == num_steps
    top_filename = os.path.join(
        host_guest_mmvt_model_namd.anchor_rootdir, anchor.directory, 
        anchor.building_directory, anchor.amber_params.prmtop_filename)
    dcd_length, dcd_file_number = get_trajectory_length(
        host_guest_mmvt_model_namd, anchor, top_filename)
    assert dcd_length == num_dcd_frames
    return

def test_normal_restart_namd(host_guest_mmvt_model_namd):
    """
    Test a normal restart - no errors or interruptions, using NAMD
    """
    # Test DCD frame number and checkpoint step
    first_steps = 10
    dcd_interval = 1
    host_guest_mmvt_model_namd.calculation_settings.num_production_steps \
        = first_steps
    host_guest_mmvt_model_namd.calculation_settings\
        .restart_checkpoint_interval = 1
    host_guest_mmvt_model_namd.calculation_settings\
        .trajectory_reporter_interval = dcd_interval
    num_dcd_frames = first_steps // dcd_interval
    run.run(host_guest_mmvt_model_namd, "0", force_overwrite=True)
    anchor = host_guest_mmvt_model_namd.anchors[0]
    output_directory = os.path.join(
        host_guest_mmvt_model_namd.anchor_rootdir, anchor.directory, 
        anchor.production_directory)
    restart_checkpoint_filename = os.path.join(
        output_directory, runner_namd.RESTART_CHECKPOINT_FILENAME+".xsc")
    checkpoint_step = runner_namd.read_xsc_step_number(
            restart_checkpoint_filename)
    assert checkpoint_step == first_steps
    top_filename = os.path.join(
        host_guest_mmvt_model_namd.anchor_rootdir, anchor.directory, 
        anchor.building_directory, anchor.amber_params.prmtop_filename)
    dcd_length, dcd_file_number = get_trajectory_length(
        host_guest_mmvt_model_namd, anchor, top_filename)
    assert dcd_length == num_dcd_frames
    
    second_steps = 20
    num_dcd_frames2 = second_steps // dcd_interval
    host_guest_mmvt_model_namd.calculation_settings.num_production_steps \
        = second_steps
    run.run(host_guest_mmvt_model_namd, "0")
    anchor = host_guest_mmvt_model_namd.anchors[0]
    output_directory = os.path.join(
        host_guest_mmvt_model_namd.anchor_rootdir, anchor.directory, 
        anchor.production_directory)
    restart_checkpoint_filename = os.path.join(
        output_directory, runner_namd.RESTART_CHECKPOINT_FILENAME+".xsc")
    checkpoint_step = runner_namd.read_xsc_step_number(
            restart_checkpoint_filename)
    assert checkpoint_step == second_steps
    top_filename = os.path.join(
        host_guest_mmvt_model_namd.anchor_rootdir, anchor.directory, 
        anchor.building_directory, anchor.amber_params.prmtop_filename)
    dcd_length, dcd_file_number = get_trajectory_length(
        host_guest_mmvt_model_namd, anchor, top_filename)
    assert dcd_length == num_dcd_frames2
    
    # Test that there are two output files
    out_glob = os.path.join(
        host_guest_mmvt_model_namd.anchor_rootdir, anchor.directory, 
        anchor.production_directory, anchor.md_output_glob)
    out_file_list = glob.glob(out_glob)
    assert len(out_file_list) == 2
    return