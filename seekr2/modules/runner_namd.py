"""
runner_namd.py

A command-line tool for running SEEKR calculations using the NAMD
engine either locally or on a computing resource.
"""

import os
import argparse
import glob
import shutil

import seekr2.modules.common_base as base
import seekr2.modules.mmvt_base as mmvt_base
import seekr2.modules.elber_base as elber_base
import seekr2.modules.common_sim_namd as common_sim_namd
import seekr2.modules.mmvt_sim_namd as mmvt_sim_namd
import seekr2.modules.elber_sim_namd as elber_sim_namd

RESTART_CHECKPOINT_FILENAME = "backup_checkpoint"
SAVE_STATE_DIRECTORY = "states/"
SAVE_STATE_PREFIX = "namdmmvt"

def search_for_state_to_load(model, anchor):
    """
    Find an existing state in an adjacent anchor that could be used to
    seed this anchor as a set of starting positions.
    """
    for milestone in anchor.milestones:
        neighbor_anchor = model.anchors[milestone.neighbor_anchor_index]
        for neighbor_milestone in neighbor_anchor.milestones:
            if neighbor_milestone.neighbor_anchor_index == anchor.index:
                target_alias_index = str(neighbor_milestone.alias_index)
                glob_prefix = SAVE_STATE_PREFIX + "_" + target_alias_index \
                    + "_*"
                state_glob = os.path.join(
                    model.anchor_rootdir, neighbor_anchor.directory, 
                    neighbor_anchor.production_directory, SAVE_STATE_DIRECTORY,
                    glob_prefix)
                state_glob_files = glob.glob(state_glob)
                if len(state_glob_files) > 0:
                    state_file = state_glob_files[0]
                    assert os.path.exists(state_file)
                    state_file_no_ext = os.path.splitext(state_file)[0]
                    return state_file_no_ext
    return None

def read_xsc_step_number(xsc_filename):
    """
    Read a NAMD .xsc file to extract the latest step to use as input
    for restarts and such.
    """
    last_step = 0
    with open(xsc_filename, "r") as f:
        for line in f.readlines():
            if line.startswith("#"):
                continue
            last_step = int(line.strip().split(" ")[0])
    
    return last_step

def cleanse_anchor_outputs(model, anchor):
    """
    Delete all simulation outputs for an existing anchor to make way
    for new outputs.
    """
    if model.get_type() == "mmvt":
        basename = mmvt_base.NAMDMMVT_BASENAME
        extension = mmvt_base.NAMDMMVT_EXTENSION
    # Elber not supported in NAMD
    #elif model.get_type() == "elber":
    #    basename = elber_base.OPENMM_ELBER_BASENAME
    #    extension = elber_base.OPENMM_ELBER_EXTENSION
    output_directory = os.path.join(
        model.anchor_rootdir, anchor.directory, anchor.production_directory)
    output_files_glob = os.path.join(output_directory, anchor.md_output_glob)
    output_restarts_list = glob.glob(output_files_glob)
    for mmvt_output_file in output_restarts_list:
        os.remove(mmvt_output_file)
    dcd_glob = os.path.join(output_directory, "*.dcd*")
    for dcd_file in glob.glob(dcd_glob):
        os.remove(dcd_file)
    colvars_glob = os.path.join(
        output_directory, "%s*.colvars.*" % basename)
    for colvars_file in glob.glob(colvars_glob):
        os.remove(colvars_file)
    coor_glob = os.path.join(
        output_directory, "%s*.coor*" % basename)
    for coor_file in glob.glob(coor_glob):
        os.remove(coor_file)
    vel_glob = os.path.join(
        output_directory, "%s*.vel*" % basename)
    for vel_file in glob.glob(vel_glob):
        os.remove(vel_file)
    xsc_glob = os.path.join(
        output_directory, "%s*.xsc*" % basename)
    for xsc_file in glob.glob(xsc_glob):
        os.remove(xsc_file)
    input_glob = os.path.join(
        output_directory, 
        common_sim_namd.NAMD_INPUT_FILENAME.format("*"))
    for input_file in glob.glob(input_glob):
        os.remove(input_file)
    checkpoint_glob = os.path.join(
        output_directory, RESTART_CHECKPOINT_FILENAME+"*")
    for checkpoint_file in glob.glob(checkpoint_glob):
        os.remove(checkpoint_file)
    states_dir = os.path.join(
        output_directory, SAVE_STATE_DIRECTORY)
    if os.path.exists(states_dir):
        shutil.rmtree(states_dir)
    return

class Runner_namd():
    """
    Prepare and run a SEEKR2 simulation using the NAMD simulation 
    engine.
    
    Attributes:
    -----------
    model : Model()
        The model contains all the settings, parameters, directories,
        and file names used to perform a SEEKR2 simulation.
        
    sim_namd : Sim_namd()
        The sim_namd object contains all the NAMD simulation
        settings.
        
    anchor : Anchor()
        The anchor is the Voronoi Cell (MMVT) or milestoning (Elber) 
        within which the simulation will be run.
    """
    def __init__(self, model, anchor, namd_command, namd_arguments):
        self.model = model
        self.sim_namd = None
        self.anchor = anchor
        self.namd_command = namd_command
        self.namd_arguments = namd_arguments
        self.seekr2_output_files_glob = ""
        self.basename = ""
        self.extension = ""
        self.namd_input_filename = ""
        
        if model.get_type() == "mmvt":
            self.glob = mmvt_base.NAMDMMVT_GLOB
            self.basename = mmvt_base.NAMDMMVT_BASENAME
            self.extension = mmvt_base.NAMDMMVT_EXTENSION
            self.output_directory = os.path.join(
                self.model.anchor_rootdir, self.anchor.directory, 
                self.anchor.production_directory)
            self.header = mmvt_sim_namd.SEEKR_MMVT_PROD_HEADER
        elif model.get_type() == "elber":
            self.glob = elber_base.OPENMM_ELBER_GLOB
            self.basename = elber_base.OPENMM_ELBER_BASENAME
            self.extension = elber_base.OPENMM_ELBER_EXTENSION
            #assert stage in ["temp_equil", "umbrella", "fwd_rev"], \
            #    "stage not allowed for Elber milestoning: {1}.".format(stage) \
            #    + "Stage must be one of ['temp_equil', 'umbrella', 'fwd_rev']"
            #self.output_directory = os.path.join(
            #    self.model.anchor_rootdir, self.anchor.directory, 
            #    stage)
            self.output_directory = os.path.join(
                self.model.anchor_rootdir, self.anchor.directory)
            self.header = elber_sim_namd.SEEKR_ELBER_PROD_HEADER
        self.save_one_state_for_all_boundaries = True
        self.check_state_interval = 1000
    
    def prepare(self, restart=False, save_state=False, 
                force_overwrite=False):
        """
        This function gets run before the sim_namd object is created
        so that the proper paths can be found, etc.
        """
        settings = self.model.namd_settings
        assert settings is not None, "This model was not prepared for NAMD."
        
        #output_directory = os.path.join(
        #    self.model.anchor_rootdir, self.anchor.directory, 
        #    self.anchor.production_directory)
        restart_index = 1
        seekr2_output_files_glob = os.path.join(
            self.output_directory, self.anchor.md_output_glob)
        self.seekr2_output_files_glob = seekr2_output_files_glob
        mmvt_output_restarts_list = glob.glob(seekr2_output_files_glob)
        if restart:
            restart_index = len(mmvt_output_restarts_list) + 1
            default_output_filename = "%s.restart%d.%s" \
                    % (self.basename, restart_index, self.extension)
            output_basename = "%s%d" % (mmvt_base.NAMDMMVT_BASENAME, 
                                        restart_index)
        else:
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
                    """
                    for mmvt_output_file in mmvt_output_restarts_list:
                        os.remove(mmvt_output_file)
                    dcd_glob = os.path.join(
                        self.output_directory, "%s*.dcd*" % self.basename)
                    for dcd_file in glob.glob(dcd_glob):
                        os.remove(dcd_file)
                    colvars_glob = os.path.join(
                        self.output_directory, "%s*.colvars.*" % self.basename)
                    for colvars_file in glob.glob(colvars_glob):
                        os.remove(colvars_file)
                    coor_glob = os.path.join(
                        self.output_directory, "%s*.coor*" % self.basename)
                    for coor_file in glob.glob(coor_glob):
                        os.remove(coor_file)
                    vel_glob = os.path.join(
                        self.output_directory, "%s*.vel*" % self.basename)
                    for vel_file in glob.glob(vel_glob):
                        os.remove(vel_file)
                    xsc_glob = os.path.join(
                        self.output_directory, "%s*.xsc*" % self.basename)
                    for xsc_file in glob.glob(xsc_glob):
                        os.remove(xsc_file)
                    input_glob = os.path.join(
                        self.output_directory, 
                        common_sim_namd.NAMD_INPUT_FILENAME.format("*"))
                    for input_file in glob.glob(input_glob):
                        os.remove(input_file)
                    checkpoint_glob = os.path.join(
                        self.output_directory, RESTART_CHECKPOINT_FILENAME+"*")
                    for checkpoint_file in glob.glob(checkpoint_glob):
                        os.remove(checkpoint_file)
                    states_dir = os.path.join(
                        self.output_directory, SAVE_STATE_DIRECTORY)
                    if os.path.exists(states_dir):
                        shutil.rmtree(states_dir)
                    """
                        
            default_output_filename = "%s%d.%s" % (
                self.basename, 1, self.extension)
            output_basename = "%s%d" % (self.basename, 1)
            
        state_dir = os.path.join(self.output_directory, 
                                 SAVE_STATE_DIRECTORY)
        if self.save_one_state_for_all_boundaries:
            if not os.path.exists(state_dir):
                os.mkdir(state_dir)
        self.state_prefix = os.path.join(state_dir, SAVE_STATE_PREFIX)
        if save_state:
            state_prefix = self.state_prefix
            self.save_one_state_for_all_boundaries=False
            if not os.path.exists(state_dir):
                os.mkdir(state_dir)
        else:
            state_prefix = None
            
        return default_output_filename, output_basename, state_prefix, \
            restart_index
    
    def run(self, sim_namd_obj, output_file, restart=False, 
            load_state_file=None, restart_index=1):
        """
        Run the MMVT NAMD simulation.
        """
        self.sim_namd = sim_namd_obj
        settings = self.model.namd_settings
        calc_settings = self.model.calculation_settings
        output_directory = os.path.join(
            self.model.anchor_rootdir, self.anchor.directory, 
            self.anchor.production_directory)
        
        self.sim_namd.namd_root.output_files.restartfreq = \
            calc_settings.restart_checkpoint_interval
        self.sim_namd.namd_root.output_files.restartname = \
            RESTART_CHECKPOINT_FILENAME
        self.sim_namd.header = self.header
        
        if restart:
            self.sim_namd.namd_root.input_files.coordinates = ""
            self.sim_namd.namd_root.input_files.binCoordinates = \
                self.sim_namd.namd_root.output_files.restartname + ".coor"
            self.sim_namd.namd_root.input_files.binVelocities = \
                self.sim_namd.namd_root.output_files.restartname + ".vel"
            self.sim_namd.namd_root.input_files.extendedSystem = \
                self.sim_namd.namd_root.output_files.restartname + ".xsc"
            self.sim_namd.namd_root.simulation_parameters.temperature = None
            xsc_file_abs_path = os.path.join(
                output_directory, 
                self.sim_namd.namd_root.input_files.extendedSystem)
            firsttimestep = read_xsc_step_number(xsc_file_abs_path)
            self.sim_namd.namd_root.simulation_parameters.firsttimestep = \
                firsttimestep
            
            print("restarting from saved checkpoint:", 
                  self.sim_namd.namd_root.output_files.restartname, 
                  "at step number:", firsttimestep)
            
        else:
            if load_state_file is None \
                    and self.sim_namd.namd_root.input_files.try_to_load_state:
                load_state_file = search_for_state_to_load(
                    self.model, self.anchor)
                print("loading state file:", load_state_file)
                if load_state_file is None:
                    print("No states found to load for anchor: %d" \
                            % self.anchor.index)
            
            if load_state_file is not None:
                self.sim_namd.namd_root.input_files.coordinates = ""
                self.sim_namd.namd_root.input_files.binCoordinates = \
                    load_state_file + ".coor"
                self.sim_namd.namd_root.input_files.binVelocities = \
                    load_state_file + ".vel"
                self.sim_namd.namd_root.input_files.extendedSystem = \
                    load_state_file + ".xsc"
                self.sim_namd.namd_root.simulation_parameters.temperature = None
                xsc_file_abs_path = os.path.join(
                    output_directory, 
                    self.sim_namd.namd_root.input_files.extendedSystem)
                self.sim_namd.namd_root.simulation_parameters.firsttimestep = 0
                            
        curdir = os.getcwd()
        print("moving to directory:", output_directory)
        os.chdir(output_directory)
        namd_filename = sim_namd_obj.write_namd_script(self.model, self.anchor,
                                                       restart_index)
        sim_namd_obj.write_colvar_script(
            self.model, self.anchor)
        if self.namd_arguments:
            self.namd_arguments = " " + self.namd_arguments
        run_command = self.namd_command + self.namd_arguments + " " \
            + namd_filename + " +stdout " + output_file + ".%d"
        print("Running command:", run_command)
        os.system(run_command)
        results_glob = glob.glob(output_file + ".*")
        assert len(results_glob) > 0, "Problem occurred running "\
            "namd: results file(s) was not generated."
        os.chdir(curdir)
        return

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument(
        "index", metavar="ANCHOR_INDEX", type=int, 
        help="The index of the anchor whose simulation to start or restart.")
    argparser.add_argument(
        "input_file", metavar="INPUT_FILE", type=str, 
        help="name of input file for NAMD MMVT calculation. This would be the "\
        "XML file generated in the prepare stage.")
    argparser.add_argument("-r", "--restart", dest="restart", default=False,
                           help="Restart simulation from backup checkpoint in "\
                           "input file. Overrides '-l' argument to load a "\
                           "state.", action="store_true")
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
                           "set of .coor, .vel, and .xsc state files "\
                           "whenever a bounce occurs.", 
                           action="store_true")
    argparser.add_argument("-l", "--load_state_file", dest="load_state_file",
                           default=None, help="Load state files from a "\
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
    argparser.add_argument("-n", "--namd_command", dest="namd_command",
                           default="namd2", help="Define a different NAMD "\
                           "command. This allows users to define a path to "\
                           "NAMD or to use a 'charmrun namd2' command."\
                           " By default, 'namd2' is used.",
                           type=str)
    argparser.add_argument("-a", "--namd_arguments", dest="namd_arguments",
                           default="", help="Additional arguments for NAMD can"\
                           " be entered here. Note: this will require "\
                           "quotation marks around the arguments. Example: "\
                           "-a '+p8 +devices 0,2'.", type=str)
        
    args = argparser.parse_args()
    args = vars(args)
    anchor_index = int(args["index"])
    input_file = args["input_file"]
    restart = args["restart"]
    total_simulation_length = args["total_simulation_length"]
    output_file = args["output_file"]
    save_state_file = args["save_state_file"]
    load_state_file = args["load_state_file"]
    directory = args["directory"]
    force_overwrite = args["force_overwrite"]
    namd_command = args["namd_command"]
    namd_arguments = args["namd_arguments"]
    
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
    
    if total_simulation_length is not None:
        mymodel.namd_settings.total_simulation_length = \
            total_simulation_length
    runner = Runner_namd(mymodel, myanchor, namd_command, namd_arguments)
    default_output_file, output_basename, state_file_prefix, restart_index \
        = runner.prepare(restart, save_state_file, force_overwrite)
    if output_file is None:
        output_file = default_output_file
    
    if mymodel.get_type() == "mmvt":
        sim_namd_factory = mmvt_sim_namd.MMVT_sim_namd_factory()
    elif mymodel.get_type() == "elber":
        pass
    else:
        raise Exception("Calculation type not supported: {1}".format(
            mymodel.get_type()))
        
    sim_namd_obj = sim_namd_factory.create_sim_namd(
        mymodel, myanchor, output_basename)
    sim_namd_obj.seekr_namd_settings.save_state = save_state_file
    sim_namd_obj.seekr_namd_settings.save_one_state_for_all_boundaries\
         = runner.save_one_state_for_all_boundaries
    sim_namd_obj.seekr_namd_settings.check_state_interval\
         = runner.check_state_interval
    
    runner.run(sim_namd_obj, output_file, restart, load_state_file, 
               restart_index)