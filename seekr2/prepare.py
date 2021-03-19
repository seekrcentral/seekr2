"""
Accept arguments and settings for creating a seekr2 calculation, 
including its model.xml file as well as its entire filetree.
"""

import os
import argparse

import seekr2.modules.common_base as base
import seekr2.modules.common_prepare as common_prepare
import seekr2.modules.filetree as filetree
from seekr2.modules.common_base import Amber_params, Forcefield_params 
from seekr2.modules.common_cv import *
from seekr2.modules.common_prepare import Browndye_settings_input, Ion, \
    MMVT_input_settings, Elber_input_settings

def generate_openmmvt_model_and_filetree(model_input):
    """
    Using the Model_input from the user, prepare the Model
    object and the filetree. Then prepare all building files
    for each anchor and serialize the Model to XML.
    """
    model_factory = common_prepare.Model_factory()
    model = model_factory.create_model(model_input)
    common_prepare.prepare_model_cvs_and_anchors(model, model_input)
    root_directory = os.path.expanduser(model_input.root_directory)
    filetree.generate_filetree(model, root_directory)
    filetree.copy_building_files(model, model_input, root_directory)
    common_prepare.generate_bd_files(model, root_directory)
    xml_path = os.path.join(root_directory, "model.xml")
    model.serialize(xml_path)
    return model

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument(
        "model_input_file", metavar="MODEL_INPUT_FILE", type=str, 
        help="The name of input XML file for a SEEKR2 calculation.")
    """ # TODO: remove
    argparser.add_argument("-r", "--run", dest="run", default=False,
        help="Run the BrownDye bd_top and nam_simulation programs for the "\
        "BD_MILESTONE indicated. By default, this mode is activated if the "\
        "extract argument is False, otherwise, if extract is True, this "\
        "argument must be specified in order to do both.",
        action="store_true")
    """
        
    args = argparser.parse_args() # parse the args into a dictionary
    args = vars(args)
    model_input_filename = args["model_input_file"]
    #run = args["run"]
    model_input = common_prepare.Model_input()
    model_input.deserialize(model_input_filename, user_input = True)
    generate_openmmvt_model_and_filetree(model_input)