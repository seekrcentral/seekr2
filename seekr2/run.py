"""
run.py

Run any and all simulations required for the SEEKR2 calculation.
"""

import os
import argparse
import tempfile
import math

from parmed import unit

import seekr2.modules.common_base as base
import seekr2.converge as converge
import seekr2.modules.common_converge as common_converge

# run 20 ns between convergence checks
CONVERGENCE_INTERVAL = 100000
B_SURFACE_CONVERGENCE_INTERVAL = 10000
BD_MILESTONE_CONVERGENCE_INTERVAL = 100
MAX_ITER = 99999

def choose_next_simulation_browndye2(
        model, instruction, min_b_surface_simulation_length, 
        min_bd_milestone_simulation_length, max_b_surface_trajs_to_extract,
        force_overwrite, min_b_surface_encounters=None, 
        min_bd_milestone_encounters=None):
    """
    Examine the model and all Browndye2 simulations that have run so
    far (if any), then construct a list of which BD milestones to
    run, in order.
    """
    if instruction == "any_md":
        return []
    
    import seekr2.modules.common_sim_browndye2 as sim_browndye2
    import seekr2.modules.runner_browndye2 as runner_browndye2
    
    if min_b_surface_simulation_length is None:
        min_b_surface_simulation_length \
            = model.k_on_info.b_surface_num_trajectories
        
    # if overwrite is True, cleanse all BD directories for new runs
    b_surface_directory = os.path.join(
        model.anchor_rootdir, model.k_on_info.b_surface_directory)
    if force_overwrite in ["any", "any_bd", "b_surface"]:
        runner_browndye2.cleanse_bd_outputs(b_surface_directory, 
                                            check_mode=False)
        for bd_milestone in model.k_on_info.bd_milestones:
            bd_milestone_directory = os.path.join(
                model.anchor_rootdir, bd_milestone.directory)
            
    bd_milestone_info_to_run_unsorted = []
    data_sample_list = converge.converge(model)
    data_sample = data_sample_list[-1]
    bd_transition_counts = data_sample.bd_transition_counts
    if instruction in ["any", "any_bd", "b_surface"]:
        # then try to run the b-surface simulations
        if "b_surface" not in bd_transition_counts:
            bd_milestone_info_to_run_unsorted.append(
                [min_b_surface_simulation_length, 0, "b_surface", False, 
                 min_b_surface_simulation_length, 
                 max_b_surface_trajs_to_extract])
        else:
            # see how many simulations need to be run in b_surface
            total_b_surface_sims = sum(
                bd_transition_counts["b_surface"].values())
            total_b_surface_encounters = min(
                bd_transition_counts["b_surface"].values())
            steps_to_go_minimum = min_b_surface_simulation_length \
                - total_b_surface_sims
            if steps_to_go_minimum > 0:
                bd_milestone_info_to_run_unsorted.append(
                    [steps_to_go_minimum, total_b_surface_encounters, \
                     "b_surface", True, steps_to_go_minimum, \
                     max_b_surface_trajs_to_extract])
            
            elif min_b_surface_encounters is not None:
                if total_b_surface_encounters < min_b_surface_encounters:
                    total_bd_simulation_length \
                        = (total_b_surface_sims \
                        // B_SURFACE_CONVERGENCE_INTERVAL + 1) \
                        * B_SURFACE_CONVERGENCE_INTERVAL
                    bd_milestone_info_to_run_unsorted.append(
                        [total_bd_simulation_length-total_b_surface_sims, \
                         total_b_surface_encounters, "b_surface", True, \
                         total_bd_simulation_length, \
                         max_b_surface_trajs_to_extract])
    
    for i, bd_milestone in enumerate(model.k_on_info.bd_milestones):
        if instruction == "b_surface": break
        if instruction not in ["any", "any_bd"]:
            if instruction.startswith("b"):
                try:
                    integer_instruction = int(instruction[1:])
                except ValueError:
                    return []
                    
                if bd_milestone.index != integer_instruction:
                    continue
            
            else:
                return []
        if min_bd_milestone_simulation_length is None:
            min_bd_milestone_simulation_length = bd_milestone.num_trajectories
            
        if max_b_surface_trajs_to_extract is None:
            max_b_surface_trajs_to_extract \
                = bd_milestone.max_b_surface_trajs_to_extract
        
        if bd_milestone.index not in bd_transition_counts:
            bd_milestone_info_to_run_unsorted.append(
                [min_bd_milestone_simulation_length, 0, bd_milestone.index, 
                 False, min_bd_milestone_simulation_length, 
                 max_b_surface_trajs_to_extract])
        else:
            # see how many simulations need to be run in b_surface
            total_b_surface_sims = sum(
                bd_transition_counts[bd_milestone.index].values())
            total_b_surface_encounters = min(
                bd_transition_counts[bd_milestone.index].values())
            steps_to_go_minimum = min_bd_milestone_simulation_length \
                * max_b_surface_trajs_to_extract - total_b_surface_sims
            if steps_to_go_minimum > 0:
                bd_milestone_info_to_run_unsorted.append(
                    [steps_to_go_minimum, total_b_surface_encounters, \
                     bd_milestone.index, True, \
                     min_bd_milestone_simulation_length, \
                     max_b_surface_trajs_to_extract])
            
            elif min_bd_milestone_encounters is not None:
                if total_b_surface_encounters < min_bd_milestone_encounters:
                    total_bd_simulation_length \
                        = (total_b_surface_sims \
                        // BD_MILESTONE_CONVERGENCE_INTERVAL + 1) \
                        * BD_MILESTONE_CONVERGENCE_INTERVAL
                    bd_milestone_info_to_run_unsorted.append(
                        [0, total_b_surface_encounters, \
                         bd_milestone.index, True, total_bd_simulation_length, \
                         max_b_surface_trajs_to_extract])
            
    return bd_milestone_info_to_run_unsorted

def choose_next_simulation_openmm(
        model, instruction, min_total_simulation_length, 
        max_total_simulation_length, convergence_cutoff, 
        minimum_anchor_transitions, force_overwrite, umbrella_restart_mode):
    """
    Examine the model and all MD simulations that have run so far.
    Using this information, as well as the specified criteria (minimum
    number of steps, minimum convergence, maximum number of steps,
    etc.), construct a list of anchors to run, in order.
    """
    import seekr2.modules.runner_openmm as runner_openmm
    import seekr2.modules.mmvt_sim_openmm as mmvt_sim_openmm
    import seekr2.modules.elber_sim_openmm as elber_sim_openmm
        
    if instruction == "any_bd":
        return []
    anchor_info_to_run_unsorted = []
    for alpha, anchor in enumerate(model.anchors):
        if anchor.bulkstate:
            continue
        if instruction not in ["any", "any_md"]:
            try:
                integer_instruction = int(instruction)
            except ValueError:
                return []
                
            if alpha != integer_instruction:
                continue
            
        if min_total_simulation_length is None:
            if model.get_type() == "mmvt":
                min_total_simulation_length \
                    = model.calculation_settings.num_production_steps
            elif model.get_type() == "elber":
                min_total_simulation_length \
                    = model.calculation_settings.num_umbrella_stage_steps
                
        output_directory = os.path.join(
            model.anchor_rootdir, anchor.directory, 
            anchor.production_directory)
        restart_checkpoint_filename = os.path.join(
            output_directory, runner_openmm.RESTART_CHECKPOINT_FILENAME)
        if os.path.exists(restart_checkpoint_filename) and not force_overwrite\
                and not umbrella_restart_mode:
            dummy_file = tempfile.NamedTemporaryFile()
            if model.get_type() == "mmvt":
                sim_openmm_obj = mmvt_sim_openmm.create_sim_openmm(
                    model, anchor, dummy_file.name)
                simulation = sim_openmm_obj.simulation
            elif model.get_type() == "elber":
                sim_openmm_obj = elber_sim_openmm.create_sim_openmm(
                    model, anchor, dummy_file.name)
                simulation = sim_openmm_obj.umbrella_simulation
            
            simulation.loadCheckpoint(restart_checkpoint_filename)
            currentStep = int(math.ceil(simulation.context.getState().getTime()\
                              .value_in_unit(unit.picoseconds) \
                              / sim_openmm_obj.timestep))
            dummy_file.close()
            restart = True
        else:
            currentStep = 0
            restart = False
        
        if max_total_simulation_length is not None:
            if currentStep > max_total_simulation_length:
                # this simulation has run long enough
                continue
        
        # make sure that all simulations reach minimum number of steps
        # before sorting by convergence
        steps_to_go_to_minimum = min_total_simulation_length - currentStep
        if steps_to_go_to_minimum <= 0:
            data_sample_list = converge.converge(model)
            steps_to_go_to_minimum = 0
            transition_minima, dummy \
                = common_converge.calc_transition_steps(
                model, data_sample_list[-1])
            
            num_transitions = transition_minima[alpha]
            if convergence_cutoff is not None:
                rmsd_convergence_results \
                    = common_converge.calc_RMSD_conv_amount(
                        model, data_sample_list)
                convergence = rmsd_convergence_results[alpha]
                if convergence < float(convergence_cutoff):
                    continue
                else:
                    print("anchor", alpha, "has not reached the point of "\
                          "convergence:", convergence, "of", convergence_cutoff)
                    total_simulation_length \
                        = (currentStep // CONVERGENCE_INTERVAL + 1) \
                        * CONVERGENCE_INTERVAL
                    anchor_info = [steps_to_go_to_minimum, num_transitions, 
                                   alpha, restart, total_simulation_length]
                    anchor_info_to_run_unsorted.append(anchor_info)
                    continue
                    
            if minimum_anchor_transitions is not None:
                minimum_anchor_transitions = int(minimum_anchor_transitions)
                if num_transitions >= minimum_anchor_transitions:
                    continue
                else:
                    print("anchor", alpha, "has not had the minimum number of "\
                          "transitions:", num_transitions, "of", 
                          minimum_anchor_transitions)
                    total_simulation_length = (
                        currentStep // CONVERGENCE_INTERVAL + 1) \
                        * (CONVERGENCE_INTERVAL)
                    anchor_info = [steps_to_go_to_minimum, num_transitions, 
                                   alpha, restart, total_simulation_length]
                    anchor_info_to_run_unsorted.append(anchor_info)
                    
        else:
            print("anchor", alpha, "has not run the minimum number of steps",
                  currentStep, "of", min_total_simulation_length)
            total_simulation_length = min_total_simulation_length
            num_transitions = 0
            anchor_info = [steps_to_go_to_minimum, num_transitions, 
                           alpha, restart, total_simulation_length]
            anchor_info_to_run_unsorted.append(anchor_info)
        
    # sort anchors first by how many steps still to run, then by how many
    # transitions have been completed.    
    anchor_info_to_run = sorted(
        anchor_info_to_run_unsorted, 
        key=lambda item: (min_total_simulation_length-item[0], item[1]))
    return anchor_info_to_run

def choose_next_simulation_namd(
        model, instruction, min_total_simulation_length, 
        max_total_simulation_length, convergence_cutoff, 
        minimum_anchor_transitions, force_overwrite):
    """
    Examine the model and all MD simulations that have run so far.
    Using this information, as well as the specified criteria (minimum
    number of steps, minimum convergence, maximum number of steps,
    etc.), construct a list of anchors to run, in order.
    """
    import seekr2.modules.runner_namd as runner_namd
    import seekr2.modules.mmvt_sim_namd as mmvt_sim_namd
        
    if instruction == "any_bd":
        return []
    anchor_info_to_run_unsorted = []
    for alpha, anchor in enumerate(model.anchors):
        if anchor.bulkstate:
            continue
        if instruction not in ["any", "any_md"]:
            try:
                integer_instruction = int(instruction)
            except ValueError:
                #print("Invalid argument for INSTRUCTION provided: "\
                #      "'{}'. ".format(instruction)\
                #      +"Allowed arguments: 'any', 'any_md', 'any_bd'.")
                #exit()
                return []
                
            if alpha != integer_instruction:
                continue
        if min_total_simulation_length is None:
            min_total_simulation_length \
                = model.calculation_settings.num_production_steps
                
        output_directory = os.path.join(
            model.anchor_rootdir, anchor.directory, 
            anchor.production_directory)
        restart_checkpoint_filename = os.path.join(
            output_directory, runner_namd.RESTART_CHECKPOINT_FILENAME+".xsc")
        if os.path.exists(restart_checkpoint_filename) and not force_overwrite:
            currentStep = runner_namd.read_xsc_step_number(
                restart_checkpoint_filename)
            restart = True
        else:
            currentStep = 0
            restart = False
        
        if max_total_simulation_length is not None:
            if currentStep > max_total_simulation_length:
                # this simulation has run long enough
                continue
        
        # make sure that all simulations reach minimum number of steps
        # before sorting by convergence
        steps_to_go_to_minimum = min_total_simulation_length - currentStep
        if steps_to_go_to_minimum <= 0:
            data_sample_list = converge.converge(model)
            steps_to_go_to_minimum = 0
            transition_results, dummy \
                = common_converge.calc_transition_steps(
                model, data_sample_list[-1])
            
            num_transitions = transition_results[alpha]
            if convergence_cutoff is not None:
                if model.get_type() == "mmvt":
                    rmsd_convergence_results \
                        = common_converge.calc_RMSD_conv_amount(
                            model, data_sample_list)
                elif model.get_type() == "elber":
                    raise Exception("Elber auto-run not available.")
                convergence = rmsd_convergence_results[alpha]
                if convergence < float(convergence_cutoff):
                    continue
                else:
                    print("anchor", alpha, "has not reached the point of "\
                          "convergence:", convergence, "of", convergence_cutoff)
                    total_simulation_length \
                        = (currentStep // CONVERGENCE_INTERVAL + 1) \
                        * CONVERGENCE_INTERVAL
                    anchor_info = [steps_to_go_to_minimum, num_transitions, 
                                   alpha, restart, total_simulation_length]
                    anchor_info_to_run_unsorted.append(anchor_info)
                    continue
                    
            if minimum_anchor_transitions is not None:
                minimum_anchor_transitions = int(minimum_anchor_transitions)
                if num_transitions >= minimum_anchor_transitions:
                    continue
                else:
                    print("anchor", alpha, "has not had the minimum number of "\
                          "transitions:", num_transitions, "of", 
                          minimum_anchor_transitions)
                    total_simulation_length = (
                        currentStep // CONVERGENCE_INTERVAL + 1) \
                        * (CONVERGENCE_INTERVAL)
                    anchor_info = [steps_to_go_to_minimum, num_transitions, 
                                   alpha, restart, total_simulation_length]
                    anchor_info_to_run_unsorted.append(anchor_info)
                    
        else:
            print("anchor", alpha, "has not run the minimum number of steps",
                  currentStep, "of", min_total_simulation_length)
            total_simulation_length = min_total_simulation_length
            num_transitions = 0
            anchor_info = [steps_to_go_to_minimum, num_transitions, 
                           alpha, restart, total_simulation_length]
            anchor_info_to_run_unsorted.append(anchor_info)
        
    # sort anchors first by how many steps still to run, then by how many
    # transitions have been completed.    
    anchor_info_to_run = sorted(
        anchor_info_to_run_unsorted, 
        key=lambda item: (min_total_simulation_length-item[0], item[1]))
    return anchor_info_to_run

def run_browndye2(model, bd_milestone_index, restart, n_trajectories, 
                  force_overwrite=False, max_b_surface_trajs_to_extract=1000):
    """Run a Browndye2 simulation."""
    import seekr2.modules.runner_browndye2 as runner_browndye2
    
    if bd_milestone_index == "b_surface":
        bd_milestone_directory = os.path.join(
            model.anchor_rootdir, model.k_on_info.b_surface_directory)
        bd_directory_list = [bd_milestone_directory]
    
    else:
        assert bd_milestone_index >= 0, "only positive indices allowed."
        try:
            bd_milestone = model.k_on_info.bd_milestones[bd_milestone_index]
        except IndexError:
            print("Invalid bd_milestone index provided.")
            exit()
        bd_milestone_directory = os.path.join(
            model.anchor_rootdir, bd_milestone.directory)
        #assert not b_surface, "Extraction may not be performed on the "\
        #    "b-surface."
        lig_pqr_filenames, rec_pqr_filenames \
            = runner_browndye2.extract_bd_surface(
            model, bd_milestone, max_b_surface_trajs_to_extract, 
            force_overwrite, restart)
        bd_directory_list = runner_browndye2.make_fhpd_directories(
            model, bd_milestone, lig_pqr_filenames, rec_pqr_filenames, 
            force_overwrite)
        
    for bd_directory in bd_directory_list:
        runner_browndye2.run_bd_top(
            model.browndye_settings.browndye_bin_dir, bd_directory, restart, 
            force_overwrite)
        n_trajectories_per_output = runner_browndye2.DEFAULT_N_TRAJ_PER_OUT
        if n_trajectories is not None:
            if n_trajectories_per_output > n_trajectories:
                n_trajectories_per_output = n_trajectories
        runner_browndye2.modify_variables(
            bd_directory, model.k_on_info.bd_output_glob, n_trajectories, 
            restart=restart, n_trajectories_per_output\
            =n_trajectories_per_output)
        runner_browndye2.run_nam_simulation(
            model.browndye_settings.browndye_bin_dir, bd_directory, 
            model.k_on_info.bd_output_glob)
    
    if bd_milestone_index != "b_surface":
        runner_browndye2.combine_fhpd_results(
            model, bd_milestone, bd_directory_list)
    return

def run_openmm(model, anchor_index, restart, total_simulation_length, 
               cuda_device_index=False, force_overwrite=False,
               save_state_file=False, num_rev_launches=1, 
               umbrella_restart_mode=False):
    """Run an OpenMM simulation."""
    import seekr2.modules.runner_openmm as runner_openmm
    import seekr2.modules.mmvt_sim_openmm as mmvt_sim_openmm
    import seekr2.modules.elber_sim_openmm as elber_sim_openmm
        
    assert anchor_index >= 0, "only positive indices allowed."
    try:
        myanchor = model.anchors[anchor_index]
    except IndexError:
        print("Invalid anchor index provided.")
        exit()
    
    if cuda_device_index is not None:
        assert model.openmm_settings.cuda_platform_settings is not None
        model.openmm_settings.cuda_platform_settings.cuda_device_index = \
            cuda_device_index
            
    if total_simulation_length is not None:
        if model.get_type() == "mmvt":
            old_num_production_steps = \
                model.calculation_settings.num_production_steps
            model.calculation_settings.num_production_steps = \
                total_simulation_length
        elif model.get_type() == "elber":
            old_num_production_steps = \
                model.calculation_settings.num_umbrella_stage_steps
            model.calculation_settings.num_umbrella_stage_steps = \
                total_simulation_length
            model.calculation_settings.num_rev_launches = num_rev_launches
    
    runner = runner_openmm.Runner_openmm(model, myanchor)
    default_output_file, state_file_prefix, restart_index = runner.prepare(
        restart, save_state_file, force_overwrite, umbrella_restart_mode)
    
    if model.get_type() == "mmvt":
        sim_openmm_obj = mmvt_sim_openmm.create_sim_openmm(
            model, myanchor, default_output_file, state_file_prefix)
    elif model.get_type() == "elber":
        sim_openmm_obj = elber_sim_openmm.create_sim_openmm(
            model, myanchor, default_output_file, state_file_prefix)
    
    runner.run(sim_openmm_obj, restart, None, restart_index)
    if model.get_type() == "mmvt":
        model.calculation_settings.num_production_steps \
            = old_num_production_steps
    elif model.get_type() == "elber":
        model.calculation_settings.num_umbrella_stage_steps \
            = old_num_production_steps

def run_namd(model, anchor_index, restart, total_simulation_length, 
                    cuda_device_index=False, force_overwrite=False,
                    save_state_file=False, namd_command="namd2", namd_arguments=""):
    """Run a NAMD simulation."""
    import seekr2.modules.runner_namd as runner_namd
    import seekr2.modules.mmvt_sim_namd as mmvt_sim_namd
        
    assert anchor_index >= 0, "only positive indices allowed."
    try:
        myanchor = model.anchors[anchor_index]
    except IndexError:
        print("Invalid anchor index provided.")
        exit()
    
    if cuda_device_index is not None:
        assert model.openmm_settings.cuda_platform_settings is not None
        model.openmm_settings.cuda_platform_settings.cuda_device_index = \
            cuda_device_index
            
    if total_simulation_length is not None:
        old_num_production_steps = \
            model.calculation_settings.num_production_steps
        model.calculation_settings.num_production_steps = \
            total_simulation_length
    
    runner = runner_namd.Runner_namd(model, myanchor, namd_command, 
                                     namd_arguments)
    default_output_file, output_basename, state_file_prefix, restart_index \
        = runner.prepare(restart, save_state_file, force_overwrite)
    
    if model.get_type() == "mmvt":
        sim_namd_obj = mmvt_sim_namd.create_sim_namd(
            model, myanchor, output_basename)
    else:
        raise Exception("Method not implemented for NAMD:", model.get_type())
    
    sim_namd_obj.seekr_namd_settings.save_state = save_state_file
    sim_namd_obj.seekr_namd_settings.save_one_state_for_all_boundaries\
         = runner.save_one_state_for_all_boundaries
    sim_namd_obj.seekr_namd_settings.check_state_interval\
         = runner.check_state_interval
    
    runner.run(sim_namd_obj, default_output_file, restart, None, restart_index)
    model.calculation_settings.num_production_steps = old_num_production_steps

def run(model, instruction, min_total_simulation_length=None, 
        max_total_simulation_length=None, convergence_cutoff=None, 
        directory=None, minimum_anchor_transitions=None, 
        cuda_device_index=None, force_overwrite=False, save_state_file=False,
        namd_command="namd2", namd_arguments="", min_b_surface_simulation_length=None,
        min_bd_milestone_simulation_length=None, 
        max_b_surface_trajs_to_extract=None, min_b_surface_encounters=None, 
        min_bd_milestone_encounters=None, num_rev_launches=1, 
        umbrella_restart_mode=False):
    """
    Run all simulations, both MD and BD, as specified by the user 
    inputs.
    """
    md_complete = False
    bd_complete = False
    
    # Only cleanse BD files if BD is being run in this instance
    if (instruction in ["any", "any_bd"]) or instruction.startswith("b"):
        bd_force_overwrite = force_overwrite
    else:
        bd_force_overwrite = False
        
    if directory is not None:
        model.anchor_rootdir = os.path.abspath(directory)
    elif model.anchor_rootdir == ".":
        model_dir = os.path.dirname(input_file)
        model.anchor_rootdir = os.path.abspath(model_dir)
        
    assert os.path.exists(model.anchor_rootdir), "An incorrect anchor "\
        "root directory was provided."
    
    counter = 0
    while not md_complete:
        if model.openmm_settings is not None:
            anchor_info_to_run = choose_next_simulation_openmm(
                model, instruction, min_total_simulation_length, 
                max_total_simulation_length, convergence_cutoff, 
                minimum_anchor_transitions, force_overwrite, 
                umbrella_restart_mode)
        elif model.namd_settings is not None:
            anchor_info_to_run = choose_next_simulation_namd(
                model, instruction, min_total_simulation_length, 
                max_total_simulation_length, convergence_cutoff, 
                minimum_anchor_transitions, force_overwrite)
        else:
            raise Exception("No valid MD engine settings provided.")
        
        for anchor_info in anchor_info_to_run:
            steps_to_go_to_minimum, num_transitions, anchor_index, \
            restart, total_simulation_length = anchor_info
            if (force_overwrite or umbrella_restart_mode) and restart:
                restart = False
            print("running anchor_index:", anchor_index, "restart:", 
                  restart, "total_simulation_length:", 
                  total_simulation_length, "num_transitions:", 
                  num_transitions)
            if model.openmm_settings is not None:
                # import OpenMM libraries only if it will be used
                run_openmm(model, anchor_index, restart, 
                           total_simulation_length, 
                           cuda_device_index=cuda_device_index, 
                           force_overwrite=force_overwrite, 
                           save_state_file=save_state_file,
                           num_rev_launches=num_rev_launches,
                           umbrella_restart_mode=umbrella_restart_mode)
            elif model.namd_settings is not None:
                run_namd(model, anchor_index, restart, 
                         total_simulation_length, 
                         cuda_device_index=cuda_device_index, 
                         force_overwrite=force_overwrite, 
                         save_state_file=save_state_file, 
                         namd_command=namd_command,
                         namd_arguments=namd_arguments)
            else:
                raise Exception("No OpenMM or NAMD simulation settings in "\
                                "model.")
                
        if len(anchor_info_to_run) > 0:
            md_complete = False
        else:
            md_complete = True
        counter += 1
        if counter > MAX_ITER:
            raise Exception("MD while loop appears to be stuck.")
        force_overwrite = False
        if umbrella_restart_mode:
            md_complete = True
        umbrella_restart_mode = False
    
    counter = 0
    while not bd_complete:
        if model.k_on_info is None:
            break
        
        #if max_b_surface_trajs_to_extract is None:
        #    #max_b_surface_trajs_to_extract = 1e99
        #    max_b_surface_trajs_to_extract = model.
            
        bd_milestone_info_to_run = choose_next_simulation_browndye2(
            model, instruction, min_b_surface_simulation_length, 
            min_bd_milestone_simulation_length, 
            max_b_surface_trajs_to_extract, bd_force_overwrite,
            min_b_surface_encounters, min_bd_milestone_encounters)
        
        for bd_milestone_info in bd_milestone_info_to_run:
            steps_to_go_to_minimum, num_transitions, bd_milestone_index, \
                restart, total_num_trajs, max_b_surface_trajs_to_extract \
                = bd_milestone_info
            if bd_force_overwrite and restart:
                restart = False
            print("running BD:", bd_milestone_index, "restart:", 
                  restart, "trajectories to run:", steps_to_go_to_minimum, 
                  "transitions so far:", num_transitions)
            run_browndye2(
                model, bd_milestone_index, restart, steps_to_go_to_minimum, 
                force_overwrite=bd_force_overwrite, 
                max_b_surface_trajs_to_extract=max_b_surface_trajs_to_extract)
        
        if len(anchor_info_to_run) > 0:
            bd_complete = False
        else:
            bd_complete = True
        counter += 1
        if counter > MAX_ITER:
            raise Exception("BD while loop appears to be stuck.")
        
    return

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument(
        "instruction", metavar="INSTRUCTION",
        help="The index of the anchor whose simulation to start or restart. "\
        "The arguments 'any_md', 'any_bd', and 'any' are also allowed. The "\
        "argument 'any_md' will run any unfinished MD anchors. The argument "\
        "'any_bd' will run any unfinished BD calculations. The argument 'any' "\
        "will run either MD or BD calculations that still need to finish.")
    argparser.add_argument(
        "input_file", metavar="INPUT_FILE", type=str, 
        help="The name of the input file for SEEKR2 calculation. This would "\
        "be the XML file generated in the prepare stage.")
    argparser.add_argument("-c", "--cuda_device_index", 
                           dest="cuda_device_index", default=None,
                           help="modify which cuda_device_index to run the "\
                           "simulation on. For example, the number 0 or 1 "\
                           "would suffice. To run on multiple GPU indices, "
                           "simply enter comma separated indices. Example: "\
                           "'0,1'. If a value is not supplied, the value in "\
                           "the INPUT_FILE will be used by default.", type=str)
    argparser.add_argument("-t", "--minimum_total_simulation_length", 
                           dest="min_total_simulation_length", default=None,
                           help="Enter a minimum MD simulation length (in "\
                           "time steps) to run the simulation if a different "\
                           "number of steps are desired than what is in the "\
                           "INPUT_FILE. A longer simulation may be run if "
                           "other specified criteria are not met, such as the "\
                           "--convergence_cutoff, for example.", type=int)
    argparser.add_argument("-T", "--maximum_total_simulation_length", 
                           dest="max_total_simulation_length", default=None,
                           help="Enter the maximum number of MD simulation "\
                           "timesteps. Simulations will not exceed this "\
                           "length, even if other criteria, such as "\
                           "--convergence_cutoff are not met.", type=int)
    argparser.add_argument("-C", "--convergence_cutoff", 
                           dest="convergence_cutoff", default=None,
                           help="Enter a target convergence for each " \
                           "anchor. See the documentation of converge.py " \
                           "for more detail about how convergence is " \
                           "calculated.")
    argparser.add_argument("-m", "--minimum_anchor_transitions", 
                           dest="minimum_anchor_transitions", default=None,
                           help="Enter a minimum number of transitions that " \
                           "must be observed per milestone in a given anchor "\
                           "as a criteria for the simulations.")
    argparser.add_argument("-d", "--directory", dest="directory",
                           default=None, help="Provide the location of the "\
                           "directory which contains the anchors. If this "\
                           "argument is omitted, then the directory within "\
                           "the anchor_rootdir setting of the INPUT_FILE will "\
                           "be used.", type=str)
    argparser.add_argument("-f", "--force_overwrite", dest="force_overwrite",
                           default=False, help="Toggle whether to overwrite "\
                           "existing files within an anchor. If this option "\
                           "is enabled then the "\
                           "anchor's production directory will be emptied of "\
                           "all output, trajectory, and backup files for the "\
                           "new simulation.", action="store_true")
    argparser.add_argument("-s", "--save_state_file", dest="save_state_file",
                           default=False, help="Toggle whether to save a "\
                           "state file whenever a bounce occurs.", 
                           action="store_true")
    argparser.add_argument("-n", "--namd_command", dest="namd_command",
                           default="namd2", help="Define a different NAMD "\
                           "command. This allows users to define a path to "\
                           "NAMD or to use a 'charmrun namd2' command."\
                           " By default, 'namd2' is used. Only valid for NAMD "\
                           "runs.", type=str)
    argparser.add_argument("-a", "--namd_arguments", dest="namd_arguments",
                           default="", help="Additional arguments for NAMD can"\
                           " be entered here. Note: this will require "\
                           "quotation marks around the arguments. Example: "\
                           "-a '+p8 +devices 0,2'. Only valid for NAMD "\
                           "runs.", type=str)
    argparser.add_argument("-b", "--minimum_b_surface_trajectories", 
                           dest="minimum_b_surface_trajectories", default=None,
                           help="Enter a minimum number of BD trajectories to "\
                           "run for the b-surface if a different number of "\
                           "trajectories are desired than what is indicated "\
                           "in the <b_surface_num_trajectories> tag in the "\
                           "INPUT_FILE. A longer simulation may be run if "\
                           "other specified criteria are not met, such as the "\
                           "--min_b_surface_encounters, for example.",
                           type=int)
    argparser.add_argument("-B", "--minimum_bd_milestone_trajectories", 
                           dest="minimum_bd_milestone_trajectories", 
                           default=None,
                           help="Enter a minimum number of BD trajectories to "\
                           "run for a BD milestone if a different number of "\
                           "trajectories are desired than what is indicated "\
                           "in the <num_trajectories> tag in the "\
                           "<bd_milestone> tag of the "\
                           "INPUT_FILE. A longer simulation may be run if "\
                           "other specified criteria are not met, such as the "\
                           "--min_bd_milestone_encounters, for example.",
                           type=int)
    argparser.add_argument("-e", "--max_b_surface_trajs_to_extract", 
                           dest="max_b_surface_trajs_to_extract", default=None,
                           help="Enter a maximum number of b-surface "\
                           "encounter complexes to extract to construct the "\
                           "first hitting point distribution (FHPD) to run "\
                           "the simulations for any BD milestones that might "\
                           "be run. By default, all b-surface encounter "\
                           "complexes will be extracted.", type=int)
    argparser.add_argument("-y", "--min_b_surface_encounters", 
                           dest="min_b_surface_encounters", default=None,
                           help="Enter a minimum number of encounters that " \
                           "must be observed for the b-surface "\
                           "as a criteria to finish running simulations.",
                           type=int)
    argparser.add_argument("-z", "--min_bd_milestone_encounters", 
                           dest="min_bd_milestone_encounters", default=None,
                           help="Enter a minimum number of encounters that " \
                           "must be observed for the BD milestone "\
                           "as a criteria to finish running simulations.",
                           type=int)
    argparser.add_argument("-l", "--num_rev_launches", dest="num_rev_launches",
                          default=1, help="In Elber milestoning, this "\
                          "parameter defines how many reversals to launch "\
                          "for each equilibrium configuration generated by "\
                          "the umbrella stage. For each launch, the positions "\
                          "will be identical, but the velocities will be "\
                          "resampled from a Maxwell-Boltzmann distribution.",
                          type=int)
    argparser.add_argument("-u", "--umbrella_restart_mode", 
                           dest="umbrella_restart_mode", default=False,
                          help="In Elber milestoning, this option allows one"\
                          "to use the umbrella simulations that already exist "\
                          "in the anchor, and just re-run the reversals and "\
                          "forwards simulations.", action="store_true")
    
    args = argparser.parse_args()
    args = vars(args)
    instruction = args["instruction"]
    input_file = args["input_file"]
    cuda_device_index = args["cuda_device_index"]
    min_total_simulation_length = args["min_total_simulation_length"]
    max_total_simulation_length = args["max_total_simulation_length"]
    convergence_cutoff = args["convergence_cutoff"]
    minimum_anchor_transitions = args["minimum_anchor_transitions"]
    directory = args["directory"]
    force_overwrite = args["force_overwrite"]
    save_state_file = args["save_state_file"]
    namd_command = args["namd_command"]
    namd_arguments = args["namd_arguments"]
    minimum_b_surface_trajectories = args["minimum_b_surface_trajectories"]
    minimum_bd_milestone_trajectories = args[
        "minimum_bd_milestone_trajectories"]
    max_b_surface_trajs_to_extract = args["max_b_surface_trajs_to_extract"]
    min_b_surface_encounters = args["min_b_surface_encounters"]
    min_bd_milestone_encounters = args["min_bd_milestone_encounters"]
    num_rev_launches = args["num_rev_launches"]
    umbrella_restart_mode = args["umbrella_restart_mode"]
    
    assert os.path.exists(input_file), "A nonexistent input file was provided."
    model = base.Model()
    model.deserialize(input_file)
    
    run(model, instruction, min_total_simulation_length, 
        max_total_simulation_length, convergence_cutoff, 
        directory, minimum_anchor_transitions, cuda_device_index, 
        force_overwrite, save_state_file, namd_command, namd_arguments,
        minimum_b_surface_trajectories, minimum_bd_milestone_trajectories,
        max_b_surface_trajs_to_extract, min_b_surface_encounters, 
        min_bd_milestone_encounters, num_rev_launches, umbrella_restart_mode)
        