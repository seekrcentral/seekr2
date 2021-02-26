"""
runner_openmm.py

A command-line tool for running SEEKR calculations using the OpenMM
engine either locally or on a computing resource.
"""

import os
import sys
import argparse
import glob
import shutil
import time

from simtk import unit

import openmmvt.base as base
import openmmvt.sim_openmm as sim_openmm

RESTART_CHECKPOINT_FILENAME = "backup.checkpoint"
SAVE_STATE_DIRECTORY = "states/"
SAVE_STATE_PREFIX = "openmmvt"

class Runner_openmm():
    """
    Prepare and run an MMVT simulation using the OpenMM simulation 
    engine.
    
    Attributes:
    -----------
    model : Model()
        The model contains all the settings, parameters, directories,
        and file names used to perform an MMVT simulation.
        
    sim_openmm : Sim_openmm()
        The sim_openmm object contains all the OpenMM simulation
        objects, including the Simulation(), System(), Integrator(),
        etc.
        
    anchor : Anchor()
        The anchor is the Voronoi Cell within which the simulation
        will be run.
    """
    def __init__(self, model, anchor):
        self.model = model
        self.sim_openmm = None
        self.anchor = anchor
    
    def prepare(self, restart=False, save_state_file=False, 
                force_overwrite=False):
        """
        This function gets run before the sim_openmm object is created
        so that the proper paths can be found, etc.
        """
        settings = self.model.openmm_settings
        
        output_directory = os.path.join(
            self.model.anchor_rootdir, self.anchor.directory, 
            self.anchor.production_directory)
        restart_index = 1
        if restart:
            mmvt_output_files_glob = os.path.join(
                output_directory, base.OPENMMVT_GLOB)
            mmvt_output_restarts_list = glob.glob(mmvt_output_files_glob)
            restart_index = len(mmvt_output_restarts_list) + 1
            default_output_filename = os.path.join(
                output_directory, 
                "%s.restart%d.%s" % (base.OPENMMVT_BASENAME, restart_index, 
                                     base.OPENMMVT_EXTENSION))
        else:
            mmvt_output_files_glob = os.path.join(
                output_directory, base.OPENMMVT_GLOB)
            mmvt_output_restarts_list = glob.glob(mmvt_output_files_glob)
            if len(mmvt_output_restarts_list) > 0:
                if not force_overwrite:
                    print("This anchor already has existing output files "\
                          "and the entered command would overwrite them. "\
                          "If you desire to overwrite the existing files, "\
                          "then use the --force_overwrite (-f) option, and "\
                          "all outputs will be deleted and replace by a new "\
                          "run.")
                    raise Exception("Cannot overwrite existing MMVT outputs.")
                else:
                    for mmvt_output_file in mmvt_output_restarts_list:
                        os.remove(mmvt_output_file)
                    dcd_glob = os.path.join(
                        output_directory, '%s*.dcd' % base.OPENMMVT_BASENAME)
                    for dcd_file in glob.glob(dcd_glob):
                        os.remove(dcd_file)
                    backup_file = os.path.join(output_directory, 
                                           RESTART_CHECKPOINT_FILENAME)
                    if os.path.exists(backup_file):
                        os.remove(backup_file)
                    states_dir = os.path.join(output_directory, 
                                              SAVE_STATE_DIRECTORY)
                    if os.path.exists(states_dir):
                        shutil.rmtree(states_dir)
                        
            default_output_filename = os.path.join(
                output_directory, 
                "%s%d.%s" % (base.OPENMMVT_BASENAME, 1, 
                             base.OPENMMVT_EXTENSION))
            
        state_prefix = None 
        if save_state_file:
            state_dir = os.path.join(output_directory, SAVE_STATE_DIRECTORY)
            if not os.path.exists(state_dir):
                os.mkdir(state_dir)
            state_prefix = os.path.join(state_dir, SAVE_STATE_PREFIX)
        
        return default_output_filename, state_prefix, restart_index
    
    def run(self, sim_openmm_obj, restart=False, load_state_file=None, 
            restart_index=1):
        """
        Run the MMVT OpenMM simulation.
        """
        self.sim_openmm = sim_openmm_obj
        settings = self.model.openmm_settings
        output_directory = os.path.join(
            self.model.anchor_rootdir, self.anchor.directory, 
            self.anchor.production_directory)
        restart_checkpoint_frequency = settings.restart_checkpoint_frequency
        restart_checkpoint_filename = os.path.join(
            output_directory, RESTART_CHECKPOINT_FILENAME)
        trajectory_reporter_frequency = settings.trajectory_reporter_frequency
        energy_reporter_frequency = settings.energy_reporter_frequency
        
        system = self.sim_openmm.system
        simulation = self.sim_openmm.simulation
        integrator = self.sim_openmm.integrator
        traj_reporter = self.sim_openmm.traj_reporter
        
        if restart:
            simulation.loadCheckpoint(restart_checkpoint_filename)
            currentStep = int(simulation.context.getState().getTime()\
                              .value_in_unit(unit.picoseconds) \
                              // self.sim_openmm.timestep)
            start_chunk = int(currentStep // restart_checkpoint_frequency)
            simulation.currentStep = currentStep
            print("restarting from saved checkpoint:", 
                  restart_checkpoint_filename, "at step:", currentStep)
            # see how many restart files have already been created
            mmvt_output_files_glob = os.path.join(
                output_directory, base.OPENMMVT_GLOB)
            mmvt_output_restarts_list = glob.glob(mmvt_output_files_glob)
            
            traj_filename = os.path.join(
                output_directory, 
                '%s.restart%d.%s' % (base.OPENMMVT_BASENAME, restart_index,
                                         "dcd"))
        else:
            start_chunk = 0
            traj_filename = os.path.join(
                output_directory, 
                '%s%d.%s' % (base.OPENMMVT_BASENAME, 1,
                                         "dcd"))
            
            if load_state_file is not None:
                simulation.loadState(load_state_file)
                simulation.context.setTime(0.0)
        
        if restart_checkpoint_frequency is not None:
            end_chunk = settings.total_simulation_length // \
                restart_checkpoint_frequency
            steps_per_chunk = restart_checkpoint_frequency
        else:
            end_chunk = start_chunk + 1
            steps_per_chunk = settings.total_simulation_length
            
        if trajectory_reporter_frequency is not None:
            simulation.reporters.append(traj_reporter(
                traj_filename, trajectory_reporter_frequency))
            assert trajectory_reporter_frequency >= restart_checkpoint_frequency
        
        if energy_reporter_frequency is not None:
            simulation.reporters.append(
                self.sim_openmm.energy_reporter(
                    sys.stdout, energy_reporter_frequency, step=True, 
                    potentialEnergy=True, temperature=True, volume=True))
        
        starttime = time.time()
        for chunk in range(start_chunk, end_chunk):
            simulation.saveCheckpoint(restart_checkpoint_filename)
            simulation.step(steps_per_chunk)
        total_time = time.time() - starttime
        ns_per_day = settings.total_simulation_length \
            * self.sim_openmm.timestep * 1e-3 * 86400 / total_time
        print("Benchmark (ns/day):", ns_per_day)
        simulation.saveCheckpoint(restart_checkpoint_filename)
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
    
    assert os.path.exists(input_file), "A nonexistent input file was provided."
    mymodel = base.Model()
    mymodel.deserialize(input_file)
    
    if directory is not None:
        mymodel.anchor_rootdir = os.path.abspath(directory)
    elif mymodel.anchor_rootdir == ".":
        model_dir = os.path.dirname(input_file)
        mymodel.anchor_rootdir = os.path.abspath(model_dir)
        
    assert os.path.exists(mymodel.anchor_rootdir), "An incorrect anchor "\
        "root directory was provided."
    
    assert anchor_index >= 0, "only positive indices allowed."
    try:
        myanchor = mymodel.anchors[anchor_index]
    except IndexError:
        print("Invalid anchor index provided.")
        exit()
    
    if cuda_device_index is not None:
        assert mymodel.openmm_settings.cuda_platform_settings is not None
        mymodel.openmm_settings.cuda_platform_settings.cuda_device_index = \
            cuda_device_index
            
    if total_simulation_length is not None:
        mymodel.openmm_settings.total_simulation_length = \
            total_simulation_length
    
    runner = Runner_openmm(mymodel, myanchor)
    default_output_file, state_file_prefix, restart_index = runner.prepare(
        restart, save_state_file, force_overwrite)
    if output_file is None:
        output_file = default_output_file
    
    sim_openmm_factory = sim_openmm.Sim_openmm_factory()
    sim_openmm_obj = sim_openmm_factory.create_sim_openmm(
        mymodel, myanchor, output_file, state_file_prefix)
        
    runner.run(sim_openmm_obj, restart, load_state_file, restart_index)