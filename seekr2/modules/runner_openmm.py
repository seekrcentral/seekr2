"""
runner_openmm.py

A command-line tool for running SEEKR calculations using the OpenMM
engine either locally or on a computing resource.
"""

import os
import re
import sys
import argparse
import glob
import shutil
import time
from collections import defaultdict

import numpy as np
from parmed import unit

import seekr2.modules.common_base as base
import seekr2.modules.mmvt_cvs.mmvt_cv_base as mmvt_cv_base
import seekr2.modules.elber_cvs.elber_cv_base as elber_cv_base
import seekr2.modules.mmvt_sim_openmm as mmvt_sim_openmm
import seekr2.modules.elber_sim_openmm as elber_sim_openmm
import seekr2.modules.check as check

RESTART_CHECKPOINT_FILENAME = "backup.checkpoint"
SAVE_STATE_DIRECTORY = "states/"
SAVE_STATE_PREFIX = "openmm"

MAX_REVERSE_ITER = 1000000 # 2000 ns
MAX_FORWARD_ITER = 1000000 # 2000 ns
REV_STAGE_STEPS_PER_BLOCK = 100
FWD_STAGE_STEPS_PER_BLOCK = 100

assert "_" not in SAVE_STATE_PREFIX, "State files will not be properly "\
    "handled if they contain an underscore character: '_'"

def min_number_of_states_per_boundary(myglob, anchor):
    """
    Given a glob that should match all state files from bounces,
    determine the minimum number of bounces against any of the
    milestones.
    """
    all_milestone_aliases = []
    milestone_state_count = {}
    for milestone in anchor.milestones:
        all_milestone_aliases.append(str(milestone.alias_index))
        milestone_state_count[str(milestone.alias_index)] = 0
        
    orig_all_milestoning_aliases = all_milestone_aliases[:]
    all_state_files = glob.glob(myglob)
    for state_file in all_state_files:
        state_alias_index = os.path.basename(state_file).split("_")[2]
        assert state_alias_index in orig_all_milestoning_aliases, \
            "Unknown boundary alias index found: "+state_alias_index\
            +" orig_all_milestoning_aliases: "+str(orig_all_milestoning_aliases)
        milestone_state_count[str(state_alias_index)] += 1
    
    min_state_count = min(milestone_state_count.values())
    return min_state_count

def all_boundaries_have_state(myglob, anchor):
    """
    Given a glob that should match all state files from bounces, 
    check whether bounces have been made against all milestones
    in this anchor.
    """
    min_state_count = min_number_of_states_per_boundary(myglob, anchor)
    if min_state_count == 0:
        return False
    else:
        return True
    
def elber_anchor_has_umbrella_files(model, anchor):
    """
    Determine whether umbrella files already exist, which can save time
    by skipping the Elber umbrella stage.
    """
    umbrella_glob = elber_cv_base.ELBER_UMBRELLA_BASENAME+"*.dcd"
    anchor_prod_directory = os.path.join(
        model.anchor_rootdir, anchor.directory, anchor.production_directory)
    umbrella_files_glob = os.path.join(anchor_prod_directory, umbrella_glob)
    umbrella_files_list = glob.glob(umbrella_files_glob)
    if len(umbrella_files_list) > 0:
        return True
    else:
        return False

def search_for_state_to_load(model, anchor):
    """
    Find an existing state in an adjacent anchor that could be used to
    seed this anchor as a set of starting positions.
    """
    for milestone in anchor.milestones:
        if milestone.neighbor_anchor_index is None:
            continue
        neighbor_anchor = model.anchors[milestone.neighbor_anchor_index]
        for neighbor_milestone in neighbor_anchor.milestones:
            if neighbor_milestone.neighbor_anchor_index == anchor.index:
                target_alias_index = str(neighbor_milestone.alias_index)
                glob_prefix = SAVE_STATE_PREFIX + "_*_" + target_alias_index
                state_glob = os.path.join(
                    model.anchor_rootdir, neighbor_anchor.directory, 
                    neighbor_anchor.production_directory, SAVE_STATE_DIRECTORY,
                    glob_prefix)
                state_glob_files = glob.glob(state_glob)
                if len(state_glob_files) > 0:
                    state_file = state_glob_files[0]
                    assert os.path.exists(state_file)
                    return state_file
    return None

def get_data_file_length(data_file_name):
    """
    Open a file, and read the number of lines in it.
    
    Parameters
    -----------
    data_file_name : str
        A path to the file name to read the length of.
    
    Returns
    -------
    file_length : int
        The number of lines within the file.
    """
    
    if not os.path.exists(data_file_name):
        return 0
    with open(data_file_name, "r") as f:
        file_length = len(f.readlines())
    return file_length 

def read_reversal_data_file_last(data_file_name):
    """
    Open the data file, read the final transition, and return whether the
    latest reversal was successful.
    
    Parameters
    -----------
    data_file_name : str
        A path to the file name to read the transitions from
    
    Returns
    -------
    success : bool
        Whether the latest trajectory was successful
    """
    with open(data_file_name, "r") as f:
        line = f.readlines()[-1].strip().split(",")
    if line[0] != "2": # TODO: HACKY!
        return True
    else:
        return False

def cleanse_anchor_outputs(model, anchor, skip_umbrella_files=False):
    """
    Delete all simulation outputs for an existing anchor to make way
    for new outputs.
    """
    output_directory = os.path.join(
        model.anchor_rootdir, anchor.directory, anchor.production_directory)
    output_files_glob = os.path.join(output_directory, anchor.md_output_glob)
    output_restarts_list = glob.glob(output_files_glob)
    for output_file in output_restarts_list:
        os.remove(output_file)
    mmvt_dcd_glob = os.path.join(output_directory, "mmvt*.dcd")
    for dcd_file in glob.glob(mmvt_dcd_glob):
        os.remove(dcd_file)
    if not skip_umbrella_files:
        umbrella_dcd_glob = os.path.join(output_directory, "umbrella*.dcd")
        for dcd_file in glob.glob(umbrella_dcd_glob):
            os.remove(dcd_file)
    rev_dcd_glob = os.path.join(output_directory, "reverse*.dcd")
    for dcd_file in glob.glob(rev_dcd_glob):
        os.remove(dcd_file)
    fwd_dcd_glob = os.path.join(output_directory, "forward*.dcd")
    for dcd_file in glob.glob(fwd_dcd_glob):
        os.remove(dcd_file)
    backup_file = os.path.join(output_directory, 
                           RESTART_CHECKPOINT_FILENAME)
    if os.path.exists(backup_file):
        os.remove(backup_file)
    backup_swarm_glob = os.path.join(output_directory, 
                           RESTART_CHECKPOINT_FILENAME+".swarm_*")
    for backup_file in glob.glob(backup_swarm_glob):
        if os.path.exists(backup_file):
            os.remove(backup_file)
    states_dir = os.path.join(output_directory, 
                              SAVE_STATE_DIRECTORY)
    if os.path.exists(states_dir):
        shutil.rmtree(states_dir)
    
    elber_rev_glob = os.path.join(
        output_directory, elber_cv_base.ELBER_REV_GLOB)
    for elber_rev_file in glob.glob(elber_rev_glob):
        os.remove(elber_rev_file)
        
    elber_fwd_glob = os.path.join(
        output_directory, elber_cv_base.ELBER_FWD_GLOB)
    for elber_fwd_file in glob.glob(elber_fwd_glob):
        os.remove(elber_fwd_file)
    return

def get_last_bounce(data_file_name):
    """
    Find the index of the last bounce and return it.
    """
    if not os.path.exists(data_file_name):
        return None
    with open(data_file_name, 'r') as data_file:
        lines = data_file.readlines()
        if len(lines) == 0:
            return None
        line = lines[-1]
    if line.startswith("#") or line.startswith("CHECKPOINT"):
        return None
    else:
        count = int(line.split(",")[1])+1
    return count

def saveCheckpoint(sim_openmm, checkpoint_filename):
    """
    Write a checkpoint file and save a line in the output to indicate
    that a checkpoint occurred.
    """
    sim_openmm.simulation.saveCheckpoint(checkpoint_filename)
    sim_time = sim_openmm.simulation.context.getTime().value_in_unit(
        unit.picoseconds)
    with open(sim_openmm.output_filename, "a") as f:
        f.write("CHECKPOINT,{:.3f}\n".format(sim_time))
    return

def make_lock_string():
    """
    Make a string that can be put into the LOCK file to detect race
    conditions.
    """
    pid = os.getpid()
    timestamp = time.time()
    lock_string = "{},{:.3f}".format(pid, timestamp)
    return lock_string

def create_lock(model, anchor, verbose=False):
    """
    Create a lock file to prevent running if a SEEKR anchor is currently being
    run by another process.
    """
    lock_filename = os.path.join(
        model.anchor_rootdir, anchor.directory, anchor.production_directory,
        "LOCK")
    if os.path.exists(lock_filename) and verbose:
        print("Warning: a lock file was detected in anchor {}. "\
                .format(anchor.index)+"This may happen because the a previous "\
                "process failed suddenly, or because a process is currently "\
                "running on this anchor.")
    lock_string = make_lock_string()
    with open(lock_filename, "w") as f:
        f.write(lock_string)
    return lock_filename, lock_string

def update_lock(lock_filename, lock_string, anchor):
    """
    Check the lock file for changes, and update the lock file.
    """
    if not os.path.exists(lock_filename):
        raise Exception("The lock file was deleted for anchor {}. "\
                        .format(anchor.index))
    
    with open(lock_filename, "r") as f:
        old_lock_string = f.readline()
    if old_lock_string != lock_string:
        with open(lock_filename, "w") as f:
            # kill the other process
            f.write("ABORT")
        raise Exception("A lock file modification was detected for anchor {}. "\
                        .format(anchor.index)+"Another SEEKR process was "\
                        "running on the same anchor at the same time. Both "\
                        "processes will be killed. You must force this "\
                        "anchor to restart from the beginning to obtain "\
                        "correct results.")
    # otherwise all is well
    lock_string = make_lock_string()
    with open(lock_filename, "w") as f:
        f.write(lock_string)
    return lock_string

class Runner_openmm():
    """
    Prepare and run SEEKR2 simulation using the OpenMM simulation 
    engine.
    
    Attributes:
    -----------
    model : Model()
        The model contains all the settings, parameters, directories,
        and file names used to perform a SEEKR2 simulation.
        
    sim_openmm : Sim_openmm()
        The sim_openmm object contains all the OpenMM simulation
        objects, including the Simulation(), System(), Integrator(),
        etc.
        
    anchor : Anchor()
        The anchor is the spatial object for this simulation.
    
    output_directory : str
        A path to the directory where output will be written.
    
    state_prefix : str or None
        If str, then this will be the prefix for the MMVT openMM plugin
        to generate state files.
        
    save_one_state_for_all_boundaries : bool, Default True
        If True, then states will be automatically saved from all
        bounces until every boundary in this runner's anchor has
        been encountered. This parameter is merely for the 
        convenience of starting simulations on adjacent anchors
        if starting structures do not exist. A slight (or possibly
        substantial) speed gain may be obtained if set to False but
        starting structures must exist for all adjacent anchors.
    """
    
    def __init__(self, model, anchor):
        self.model = model
        self.sim_openmm = None
        self.anchor = anchor
        self.output_directory = os.path.join(
                self.model.anchor_rootdir, self.anchor.directory, 
                self.anchor.production_directory)
        if model.get_type() == "mmvt":
            self.glob = mmvt_cv_base.OPENMMVT_GLOB
            self.basename = mmvt_cv_base.OPENMMVT_BASENAME
            self.extension = mmvt_cv_base.OPENMMVT_EXTENSION
            
        elif model.get_type() == "elber":
            self.glob = elber_cv_base.ELBER_REV_GLOB
            self.basename = elber_cv_base.OPENMM_ELBER_BASENAME
            self.extension = elber_cv_base.OPENMM_ELBER_EXTENSION
            
        self.state_prefix = None
        self.save_one_state_for_all_boundaries=False
        self.restart_checkpoint_filename = None
        self.restart_checkpoint_interval = None
        self.start_chunk = None
        self.end_chunk = None
        #self.num_production_steps = None
        self.steps_per_chunk = None
        self.start_bounce_counter = 0
        self.save_all_states = False
        self.umbrellas_already_exist_mode = False
        self.swarm_string = ""
    
    def prepare(self, restart=False, save_state_file=False, 
                save_state_boundaries=False, force_overwrite=False, 
                umbrella_restart_mode=False, swarm_index=None):
        """
        This function gets run before the sim_openmm object is created
        so that the proper paths can be found, etc.
        """
        settings = self.model.openmm_settings
        assert settings is not None, "This model was not prepared for OpenMM."
        restart_index = 1
        if swarm_index is None:
            self.swarm_string = ""
        else:
            self.swarm_string = ".swarm_{}".format(swarm_index)
            self.glob = "%s%s*.%s" % (mmvt_cv_base.OPENMMVT_BASENAME, 
                                      self.swarm_string, 
                                      mmvt_cv_base.OPENMMVT_EXTENSION)
        
        output_files_glob = os.path.join(
            self.output_directory, self.glob)
        output_restarts_list = glob.glob(output_files_glob)
        if restart:
            output_restarts_list = base.order_files_numerically(
                output_restarts_list)
            assert len(output_restarts_list) > 0, \
                "No simulation has yet been run: cannot use restart mode."
            if self.model.get_type() == "mmvt":
                self.start_bounce_counter = get_last_bounce(
                    output_restarts_list[-1])
            if self.start_bounce_counter is None:
                self.start_bounce_counter = 0
            restart_index = len(output_restarts_list) + 1
            default_output_filename = os.path.join(
                self.output_directory, 
                "%s%s.restart%d.%s" % (self.basename, self.swarm_string, 
                                       restart_index, self.extension))
        else:
            if len(output_restarts_list) > 0:
                if not force_overwrite and not umbrella_restart_mode:
                    print("This anchor already has existing output files "\
                          "and the entered command would overwrite them. "\
                          "If you desire to overwrite the existing files, "\
                          "then use the --force_overwrite (-f) option, and "\
                          "all outputs will be deleted and replace by a new "\
                          "run.")
                    raise Exception("Cannot overwrite existing outputs.")
                elif force_overwrite:
                    cleanse_anchor_outputs(self.model, self.anchor)
                else:
                    cleanse_anchor_outputs(self.model, self.anchor,
                                           skip_umbrella_files=True)
            
            # check if umbrellas exist
            if self.model.get_type() == "elber":
                anchor_has_umbrella_files = elber_anchor_has_umbrella_files(
                    self.model, self.anchor)
                assert not force_overwrite or not umbrella_restart_mode, \
                    "The options force_overwrite and umbrella_restart_mode "\
                    "may not both be activated at the same time."
                if umbrella_restart_mode:
                    assert anchor_has_umbrella_files, "Cannot use umbrella "\
                        "restart mode if umbrella files don't exist for "\
                        "anchor {}.".format(self.anchor.index)
                        
                if anchor_has_umbrella_files and (not force_overwrite \
                        or umbrella_restart_mode):
                    self.umbrellas_already_exist_mode = True
                    
            default_output_filename = os.path.join(
                self.output_directory, 
                "%s%s.restart%d.%s" % (self.basename, self.swarm_string, 1, 
                               self.extension))
        
        state_dir = os.path.join(self.output_directory, 
                                 SAVE_STATE_DIRECTORY)
        self.save_one_state_for_all_boundaries = save_state_boundaries
        if self.save_one_state_for_all_boundaries:
            if not os.path.exists(state_dir):
                os.mkdir(state_dir)
        self.state_prefix = os.path.join(state_dir, SAVE_STATE_PREFIX)
        if save_state_file:
            state_prefix = self.state_prefix
            self.save_all_states = True
            if not os.path.exists(state_dir):
                os.mkdir(state_dir)
        else:
            state_prefix = None
            self.save_all_states = False
        
        restart_checkpoint_basename \
            = RESTART_CHECKPOINT_FILENAME + self.swarm_string
        self.restart_checkpoint_filename = os.path.join(
            self.output_directory, restart_checkpoint_basename)
        
        return default_output_filename, state_prefix, restart_index
    
    def run(self, sim_openmm_obj, restart=False, load_state_file=None, 
            restart_index=1):
        """Run the SEEKR simulation."""
        self.sim_openmm = sim_openmm_obj
        calc_settings = self.model.calculation_settings
        
        if self.model.get_type() == "mmvt":
            simulation = self.sim_openmm.simulation
            self.restart_checkpoint_interval \
                = calc_settings.restart_checkpoint_interval
        elif self.model.get_type() == "elber":
            simulation = self.sim_openmm.umbrella_simulation
            self.restart_checkpoint_interval = calc_settings.fwd_rev_interval

        if restart:
            simulation.loadCheckpoint(self.restart_checkpoint_filename)
            currentStep = simulation.context.getState().getStepCount()
            self.start_chunk = int(
                currentStep // self.restart_checkpoint_interval)
            
            simulation.currentStep = currentStep
            print("restarting from saved checkpoint:", 
                  self.restart_checkpoint_filename, "at step:", currentStep)
            # see how many restart files have already been created
            traj_filename = os.path.join(
                self.output_directory, 
                "%s.restart%d%s.%s" % (self.basename, restart_index, 
                                       self.swarm_string, "dcd"))
        else:
            self.start_chunk = 0
            traj_filename = os.path.join(
                self.output_directory, 
                "%s%d%s.%s" % (self.basename, 1, self.swarm_string, "dcd"))
            
            if load_state_file is not None:
                simulation.loadState(load_state_file)
                simulation.context.setTime(0.0)
                simulation.context.setStepCount(0)
                simulation.currentStep = 0
            elif self.sim_openmm.try_to_load_state:
                load_state_file = search_for_state_to_load(
                    self.model, self.anchor)
                if load_state_file is not None:
                    simulation.loadState(load_state_file)
                    simulation.context.setTime(0.0)
                    simulation.context.setStepCount(0)
                    simulation.currentStep = 0
                else:
                    print("No states found to load for anchor: %d" \
                            % self.anchor.index)
        
        starttime = time.time()
        if self.model.get_type() == "mmvt":
            self.run_mmvt(traj_filename)
        elif self.model.get_type() == "elber":
            self.run_elber(traj_filename)
        total_time = time.time() - starttime
        ns_per_day = (self.end_chunk - self.start_chunk) \
            * self.steps_per_chunk  * self.sim_openmm.timestep * 1e-3 \
            * 86400 / total_time
        print("Benchmark (ns/day):", ns_per_day)
        return
    
    def run_mmvt(self, traj_filename):
        """Run the SEEKR2 MMVT calculation."""
        calc_settings = self.model.calculation_settings
        simulation = self.sim_openmm.simulation
        trajectory_reporter_interval \
            = calc_settings.trajectory_reporter_interval
        energy_reporter_interval = calc_settings.energy_reporter_interval
        traj_reporter = self.sim_openmm.traj_reporter
        if trajectory_reporter_interval is not None:
            simulation.reporters.append(traj_reporter(
                traj_filename, trajectory_reporter_interval, 
                #enforcePeriodicBox=False))
                # Turning this on
                enforcePeriodicBox=None))
            if self.restart_checkpoint_interval is not None:
                assert trajectory_reporter_interval >= \
                    self.restart_checkpoint_interval
            if trajectory_reporter_interval \
                    > calc_settings.num_production_steps:
                trajectory_reporter_interval \
                    = calc_settings.num_production_steps
            
        if energy_reporter_interval is not None:
            simulation.reporters.append(
                self.sim_openmm.energy_reporter(
                    sys.stdout, energy_reporter_interval, step=True, 
                    potentialEnergy=True, temperature=True, volume=True))
            if energy_reporter_interval > calc_settings.num_production_steps:
                energy_reporter_interval = calc_settings.num_production_steps
            
        if self.restart_checkpoint_interval is not None:
            if self.restart_checkpoint_interval \
                    > calc_settings.num_production_steps:
                self.restart_checkpoint_interval \
                    = calc_settings.num_production_steps
            self.end_chunk = calc_settings.num_production_steps \
                // self.restart_checkpoint_interval
            self.steps_per_chunk = self.restart_checkpoint_interval
        else:
            self.end_chunk = self.start_chunk + 1
            self.steps_per_chunk = calc_settings.num_production_steps
        
        lock_filename, lock_string = create_lock(self.model, self.anchor)
        for chunk in range(self.start_chunk, self.end_chunk):
            if self.save_one_state_for_all_boundaries:
                if all_boundaries_have_state(self.state_prefix+"*", 
                                             self.anchor):
                    if self.save_all_states:
                        self.sim_openmm.integrator.setSaveStateFileName(
                            self.state_prefix)
                    else:
                        self.sim_openmm.integrator.setSaveStateFileName("")
                        self.save_one_state_for_all_boundaries = False
                    bounce_counter = get_last_bounce(
                        self.sim_openmm.output_filename)
                    if bounce_counter is None:
                        bounce_counter = self.start_bounce_counter
                    self.sim_openmm.integrator.setBounceCounter(bounce_counter)
                    self.sim_openmm.simulation.context.reinitialize(
                        preserveState=True)
                    
                else:
                    self.sim_openmm.integrator.setSaveStateFileName(
                        self.state_prefix)
                    bounce_counter = get_last_bounce(
                        self.sim_openmm.output_filename)
                    if bounce_counter is None:
                        bounce_counter = self.start_bounce_counter
                    self.sim_openmm.integrator.setBounceCounter(bounce_counter)
                    self.sim_openmm.simulation.context.reinitialize(
                        preserveState=True)
            
            if chunk > 0:    
                saveCheckpoint(self.sim_openmm, self.restart_checkpoint_filename)
            
            self.sim_openmm.simulation.step(self.steps_per_chunk)
            lock_string = update_lock(lock_filename, lock_string, self.anchor)
        
        saveCheckpoint(self.sim_openmm, self.restart_checkpoint_filename)
        os.remove(lock_filename)
        
        return
    
    def run_elber(self, traj_filename):
        """Run the SEEKR2 Elber calculation."""
        openmm_settings = self.model.openmm_settings
        calc_settings = self.model.calculation_settings
        umbrella_simulation = self.sim_openmm.umbrella_simulation
        umbrella_trajectory_reporter_interval \
            = calc_settings.umbrella_trajectory_reporter_interval
        umbrella_energy_reporter_interval \
            = calc_settings.umbrella_energy_reporter_interval
        umbrella_traj_reporter = self.sim_openmm.umbrella_traj_reporter
        directory = os.path.dirname(traj_filename)
        basename = os.path.basename(traj_filename)
        suffix = re.sub(elber_cv_base.OPENMM_ELBER_BASENAME, "", basename)
        umbrella_basename = elber_cv_base.ELBER_UMBRELLA_BASENAME+suffix
        umbrella_traj_filename = os.path.join(directory, umbrella_basename)
        if umbrella_trajectory_reporter_interval is not None \
                and not self.umbrellas_already_exist_mode:
            umbrella_simulation.reporters.append(umbrella_traj_reporter(
                umbrella_traj_filename, umbrella_trajectory_reporter_interval, 
                #enforcePeriodicBox=False))
                enforcePeriodicBox=None))
            if calc_settings.num_umbrella_stage_steps \
                    < umbrella_trajectory_reporter_interval:
                umbrella_trajectory_reporter_interval \
                    = calc_settings.num_umbrella_stage_steps
            #assert umbrella_trajectory_reporter_interval <= calc_settings.num_umbrella_stage_steps, \
            #    "The umbrella trajectory reporter interval must be less or equal to "\
            #    "the length of the simulation."
        
        if umbrella_energy_reporter_interval is not None:
            umbrella_simulation.reporters.append(
                self.sim_openmm.umbrella_energy_reporter(
                    sys.stdout, umbrella_energy_reporter_interval, step=True, 
                    potentialEnergy=True, temperature=True, volume=True))
            if calc_settings.num_umbrella_stage_steps \
                    < umbrella_energy_reporter_interval:
                umbrella_energy_reporter_interval \
                    = calc_settings.num_umbrella_stage_steps
            #assert umbrella_energy_reporter_interval \
            #        <= calc_settings.num_umbrella_stage_steps, \
            #    "The umbrella energy reporter interval must be less "\
            #    "than or equal to the length of the simulation."
        
        rev_simulation = self.sim_openmm.rev_simulation
        rev_trajectory_reporter_interval \
            = calc_settings.rev_trajectory_reporter_interval
        rev_energy_reporter_interval \
            = calc_settings.rev_energy_reporter_interval
        rev_traj_reporter = self.sim_openmm.rev_traj_reporter
                
        fwd_simulation = self.sim_openmm.fwd_simulation
        fwd_trajectory_reporter_interval \
            = calc_settings.fwd_trajectory_reporter_interval
        fwd_energy_reporter_interval \
            = calc_settings.fwd_energy_reporter_interval
        fwd_traj_reporter = self.sim_openmm.fwd_traj_reporter
        
        assert calc_settings.fwd_rev_interval \
            < calc_settings.num_umbrella_stage_steps, \
            "The interval between reversal/forward trajectories must be less "\
            "than the length of the umbrella simulation."
        assert calc_settings.num_umbrella_stage_steps \
            % calc_settings.fwd_rev_interval == 0, \
            "The number of umbrella stage steps should be a multiple of the "\
            "interval between reversal/forward trajectories."
        
        self.end_chunk = calc_settings.num_umbrella_stage_steps \
            // calc_settings.fwd_rev_interval
        self.steps_per_chunk = calc_settings.fwd_rev_interval
        rev_data_file_name = self.sim_openmm.rev_output_filename
        if os.path.exists(rev_data_file_name):
            os.remove(rev_data_file_name)
        fwd_data_file_name = self.sim_openmm.fwd_output_filename
        if os.path.exists(fwd_data_file_name):
            os.remove(fwd_data_file_name)
        num_errors = 0
        
        # If the user provides an existing umbrella trajectory, don't re-
        # simulation, just load the frames and skip to rev and fwd.
        if self.umbrellas_already_exist_mode:
            print("Umbrella trajectories found already existing in "\
                  "anchor {}. Skipping umbrella simulations.".format(
                      self.anchor.index))
            # load DCD using mdtraj and use it as the umbrella state
            umbrella_traj = check.load_structure_with_mdtraj(
                self.model, self.anchor, mode="elber_umbrella")
            assert umbrella_traj.n_frames == self.end_chunk, \
                "Umbrella trajectories found have total length of "\
                "{} frames ".format(umbrella_traj.n_frames) \
                +"while this calculation is configured to run "\
                "{} frames.".format(self.end_chunk)
            assert openmm_settings.barostat is None, "Cannot use existing" \
                "umbrella trajectories in constant pressure mode."
        for chunk in range(self.start_chunk, self.end_chunk):
            if self.umbrellas_already_exist_mode:
                # extract frame from trajectory and load the coordinates
                print("loading umbrella frame {}.".format(chunk))
                umbrella_simulation.context.setPositions(
                    umbrella_traj.openmm_positions(chunk))
                umbrella_state = \
                    umbrella_simulation.context.getState(
                        getPositions=True, getVelocities=True)
            else:
                umbrella_simulation.saveCheckpoint(
                    self.restart_checkpoint_filename)
                umbrella_simulation.step(self.steps_per_chunk)
                umbrella_state = \
                    umbrella_simulation.context.getState(
                        getPositions=True, getVelocities=True)
                
            # for each velocity restart
            for launch_id in range(calc_settings.num_rev_launches):
                crossing_counter = chunk * calc_settings.num_rev_launches \
                    + launch_id
                counter_str = "_%d" % crossing_counter
                if self.save_one_state_for_all_boundaries:
                    if all_boundaries_have_state(self.state_prefix+"*", 
                                                 self.anchor):
                        if self.save_all_states:
                            self.sim_openmm.rev_integrator.setSaveStateFileName(
                                self.state_prefix+counter_str+"r")
                            self.sim_openmm.fwd_integrator.setSaveStateFileName(
                                self.state_prefix+counter_str+"f")
                        else:
                            self.sim_openmm.rev_integrator\
                                .setSaveStateFileName("")
                            self.sim_openmm.fwd_integrator\
                                .setSaveStateFileName("")
                            self.save_one_state_for_all_boundaries = False
                        self.sim_openmm.rev_simulation.context.reinitialize(
                            preserveState=True)
                        self.sim_openmm.fwd_simulation.context.reinitialize(
                            preserveState=True)
                        
                    else:
                        self.sim_openmm.rev_integrator.setSaveStateFileName(
                            self.state_prefix+counter_str+"r")
                        self.sim_openmm.fwd_integrator.setSaveStateFileName(
                            self.state_prefix+counter_str+"f"+str(launch_id))
                        self.sim_openmm.rev_simulation.context.reinitialize(
                            preserveState=True)
                        self.sim_openmm.fwd_simulation.context.reinitialize(
                            preserveState=True)
                rev_simulation.context.setPositions(
                    umbrella_state.getPositions())
                rev_simulation.context.setVelocitiesToTemperature(
                    self.model.openmm_settings.initial_temperature \
                    * unit.kelvin)
                rev_simulation.context.setPeriodicBoxVectors(
                    *umbrella_state.getPeriodicBoxVectors())
                self.sim_openmm.rev_integrator.setCrossingCounter(
                    crossing_counter)
                if rev_trajectory_reporter_interval is not None:
                    rev_traj_filename = os.path.join(
                        self.output_directory, 
                        "reverse_%d.dcd" % crossing_counter)
                    print("rev_traj_filename", rev_traj_filename)
                    rev_simulation.reporters = [rev_traj_reporter(
                        rev_traj_filename, rev_trajectory_reporter_interval, 
                        #enforcePeriodicBox=False)]
                        enforcePeriodicBox=None)]
                if rev_energy_reporter_interval is not None:
                    rev_simulation.reporters.append(
                        self.sim_openmm.rev_energy_reporter(
                            sys.stdout, rev_energy_reporter_interval, step=True, 
                            potentialEnergy=True, temperature=True, 
                            volume=True))
                rev_simulation.context.reinitialize(preserveState=True)
                rev_start_state = rev_simulation.context.getState(
                    getPositions=True, getVelocities=True)
                rev_data_file_length = get_data_file_length(rev_data_file_name)
                rev_block_counter = 0
                had_error = False
                while get_data_file_length(rev_data_file_name) \
                        == rev_data_file_length:
                    rev_data_file_length = get_data_file_length(
                        rev_data_file_name)
                    try:
                        rev_simulation.step(REV_STAGE_STEPS_PER_BLOCK)
                    except Exception: 
                        # if there was a NAN error
                        print("Error encountered. Continuing with the next "\
                              "umbrella frame.")
                        num_errors += 1
                        had_error = True
                        # don't want to log this as a success
                        break 
                    rev_block_counter += 1
                    if rev_block_counter > MAX_REVERSE_ITER:
                        print("maximum iterations exceeded.")
                        break
                
                rev_simulation.context.setTime(0.0)
                rev_simulation.context.setStepCount(0)
                rev_simulation.currentStep = 0
                if had_error == True:
                    # move on to the next frame
                    break 
                
                if read_reversal_data_file_last(rev_data_file_name):
                    fwd_simulation.context.setPositions(
                        rev_start_state.getPositions())
                    fwd_simulation.context.setVelocities(
                        -1.0*rev_start_state.getVelocities())
                    fwd_simulation.context.setPeriodicBoxVectors(
                        *rev_start_state.getPeriodicBoxVectors())
                    self.sim_openmm.fwd_integrator.setCrossingCounter(
                        crossing_counter)
                    fwd_simulation.context.reinitialize(preserveState=True)
                    # TODO: add reporter update here
                    fwd_data_file_length = get_data_file_length(
                        fwd_data_file_name)
                    if fwd_trajectory_reporter_interval is not None:
                        fwd_traj_filename = os.path.join(
                            self.output_directory, 
                            "forward_%d.dcd" % crossing_counter)
                        fwd_simulation.reporters = [fwd_traj_reporter(
                            fwd_traj_filename, fwd_trajectory_reporter_interval, 
                            #enforcePeriodicBox=False)]
                            enforcePeriodicBox=None)]
                    if fwd_energy_reporter_interval is not None:
                        fwd_simulation.reporters.append(
                            self.sim_openmm.fwd_energy_reporter(
                                sys.stdout, fwd_energy_reporter_interval, 
                                step=True, potentialEnergy=True, 
                                temperature=True, volume=True))
                    fwd_block_counter = 0
                    had_error = False
                    while get_data_file_length(fwd_data_file_name) \
                            == fwd_data_file_length:
                        fwd_data_file_length = get_data_file_length(
                            fwd_data_file_name)
                        try:
                            fwd_simulation.step(FWD_STAGE_STEPS_PER_BLOCK)
                        except Exception: 
                            # if there was a NAN error
                            print("Error encountered. Continuing with the "\
                                  "next umbrella frame.")
                            num_errors += 1
                            had_error = True
                            break
                        fwd_block_counter += 1
                        if fwd_block_counter > MAX_FORWARD_ITER:
                            print("maximum iterations exceeded.")
                            break
                    
                    fwd_simulation.context.setTime(0.0)
                    fwd_simulation.context.setStepCount(0)
                    fwd_simulation.currentStep = 0
                                        
        umbrella_simulation.saveCheckpoint(self.restart_checkpoint_filename)
        return
    
if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument(
        "index", metavar="ANCHOR_INDEX", type=int, 
        help="The index of the anchor whose simulation to start or restart.")
    argparser.add_argument(
        "input_file", metavar="INPUT_FILE", type=str, 
        help="name of input file for OpenMMVT calculation. This would be the "\
        "XML file generated in the prepare stage.")
    argparser.add_argument("-r", "--restart", dest="restart", default=False,
                           help="Restart simulation from backup checkpoint in "\
                           "input file. Overrides '-l' argument to load a "\
                           "state.", action="store_true")
    argparser.add_argument("-c", "--cuda_device_index", 
                           dest="cuda_device_index", default=None,
                           help="modify which cuda_device_index to run the "\
                           "simulation on. For example, the number 0 or 1 "\
                           "would suffice. To run on multiple GPU indices, "
                           "simply enter comma separated indices. Example: "\
                           "'0,1'. If a value is not supplied, the value in "\
                           "the INPUT_FILE will be used by default.", type=str)
    argparser.add_argument("-t", "--total_simulation_length", 
                           dest="total_simulation_length", default=None,
                           help="Enter a different simulation length (in "\
                           "time steps) to run the simulation if a different "\
                           "number of steps are desired than what is in the "\
                           "INPUT_FILE.", type=int)
    argparser.add_argument("-o", "--output_file", dest="output_file",
                           default=None,
                           help="Enter a file name different from the default "\
                           "if the MMVT bounces should be recorded to a "\
                           "different location.")
    argparser.add_argument("-s", "--save_state_file", dest="save_state_file",
                           default=False, help="Toggle whether to save a "\
                           "state file whenever a bounce occurs.", 
                           action="store_true")
    argparser.add_argument("-l", "--load_state_file", dest="load_state_file",
                           default=None, help="Load a state file from a "\
                           "provided file path, presumably from another "\
                           "anchor as initial conditions for this calculation.",
                           type=str)
    argparser.add_argument("-d", "--directory", dest="directory",
                           default=None, help="Provide the location of the "\
                           "directory which contains the anchors. If this "\
                           "argument is omitted, then the directory within "\
                           "the anchor_rootdir setting of the INPUT_FILE will "\
                           "be used.", type=str)
    argparser.add_argument("-f", "--force_overwrite", dest="force_overwrite",
                           default=False, help="Toggle whether to overwrite "\
                           "existing files within an anchor. If this option "\
                           "is enabled and --restart is False, then the "\
                           "anchor's production directory will be emptied of "\
                           "all output, trajectory, and backup files for the "\
                           "new simulation.", action="store_true")
    argparser.add_argument("-L", "--num_rev_launches", dest="num_rev_launches",
                          default=1, help="In Elber milestoning, this "\
                          "parameter defines how many reversals to launch "\
                          "for each equilibrium configuration generated by "\
                          "the umbrella stage. For each launch, the positions "\
                          "will be identical, but the velocities will be "\
                          "resampled from a Maxwell-Boltzmann distribution.",
                          type=int)
    argparser.add_argument("-u", "--umbrella_restart_mode", 
                           dest="num_rev_launches", default=False,
                          help="In Elber milestoning, this option allows one"\
                          "to use the umbrella simulations that already exist "\
                          "in the anchor, and just re-run the reversals and "\
                          "forwards simulations.", action="store_true")
        
    args = argparser.parse_args()
    args = vars(args)
    anchor_index = int(args["index"])
    input_file = args["input_file"]
    restart = args["restart"]
    cuda_device_index = args["cuda_device_index"]
    total_simulation_length = args["total_simulation_length"]
    output_file = args["output_file"]
    save_state_file = args["save_state_file"]
    load_state_file = args["load_state_file"]
    directory = args["directory"]
    force_overwrite = args["force_overwrite"]
    num_rev_launches = args["num_rev_launches"]
    umbrella_restart_mode = args["umbrella_restart_mode"]
    
    assert os.path.exists(input_file), "A nonexistent input file was provided."
    model = base.Model()
    model.deserialize(input_file)
    
    if directory is not None:
        model.anchor_rootdir = os.path.abspath(directory)
    elif model.anchor_rootdir == ".":
        model_dir = os.path.dirname(input_file)
        model.anchor_rootdir = os.path.abspath(model_dir)
        
    assert os.path.exists(model.anchor_rootdir), "An incorrect anchor "\
        "root directory was provided."
    
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
        model.calculation_settings.num_production_steps = \
            total_simulation_length
    
    runner = Runner_openmm(model, myanchor)
    default_output_file, state_file_prefix, restart_index = runner.prepare(
        restart, save_state_file, force_overwrite, umbrella_restart_mode)
    if output_file is None:
        output_file = default_output_file
    
    if model.get_type() == "mmvt":
        sim_openmm_obj = mmvt_sim_openmm.create_sim_openmm(
            model, myanchor, output_file, state_file_prefix, frame=0, 
            load_state_file=load_state_file)
    elif model.get_type() == "elber":
        sim_openmm_factory = elber_sim_openmm.create_sim_openmm(
            model, myanchor, output_file, state_file_prefix)
        model.calculation_settings.num_rev_launches = num_rev_launches
    
    runner.run(sim_openmm_obj, restart, load_state_file, restart_index)