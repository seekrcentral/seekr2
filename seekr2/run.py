"""
run.py

Run any and all simulations required for the SEEKR2 calculation.
"""

import os
import argparse

import seekr2.common.base as base

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument(
        "index", metavar="ANCHOR_INDEX", type=int, 
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
    argparser.add_argument("-t", "--total_simulation_length", 
                           dest="total_simulation_length", default=None,
                           help="Enter a different simulation length (in "\
                           "time steps) to run the simulation if a different "\
                           "number of steps are desired than what is in the "\
                           "INPUT_FILE.", type=int)
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
        
    args = argparser.parse_args()
    args = vars(args)
    anchor_index = int(args["index"])
    input_file = args["input_file"]
    restart = args["restart"]
    cuda_device_index = args["cuda_device_index"]
    total_simulation_length = args["total_simulation_length"]
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
        
    if mymodel.get_type() == "mmvt":
        if mymodel.openmm_settings is not None:
            pass
        elif mymodel.namd_settings is not None:
            pass
        else:
            raise Exception("No OpenMM or NAMD simulation settings in model.")
    elif mymodel.get_type() == "elber":
        if mymodel.openmm_settings is not None:
            pass
        elif mymodel.namd_settings is not None:
            pass
        else:
            raise Exception("No OpenMM or NAMD simulation settings in model.")
        