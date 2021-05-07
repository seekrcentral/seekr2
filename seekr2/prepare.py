"""
prepare.py

Accept arguments and settings for creating a seekr2 calculation, 
including its model.xml file as well as its entire filetree.
"""

import os
import argparse

import seekr2.modules.common_base as base
import seekr2.modules.common_prepare as common_prepare
import seekr2.modules.filetree as filetree
import seekr2.modules.check as check
# Don't remove the following library imports - needed by deserializer
from seekr2.modules.common_prepare import Browndye_settings_input, \
    MMVT_input_settings, Elber_input_settings
from seekr2.modules.common_base import Ion, Amber_params, Forcefield_params
from seekr2.modules.common_cv import Spherical_cv_anchor, Spherical_cv_input

def generate_seekr2_model_and_filetree(model_input, force_overwrite):
    """
    Using the Model_input from the user, prepare the Model
    object and the filetree. Then prepare all building files
    for each anchor and serialize the Model to XML.
    """
    model = common_prepare.model_factory(model_input)
    common_prepare.prepare_model_cvs_and_anchors(model, model_input)
    root_directory = os.path.expanduser(model_input.root_directory)
    xml_path = os.path.join(root_directory, "model.xml")
    if os.path.exists(xml_path):
        # then a model file already exists at this location: update
        # the anchor directories.
        old_model = base.Model()
        old_model.deserialize(xml_path)
        common_prepare.modify_model(old_model, model, root_directory,
                                    force_overwrite)
    
    filetree.generate_filetree(model, root_directory)
    filetree.copy_building_files(model, model_input, root_directory)
    common_prepare.generate_bd_files(model, root_directory)
    model.serialize(xml_path)
    return model, xml_path

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument(
        "model_input_file", metavar="MODEL_INPUT_FILE", type=str, 
        help="The name of input XML file for a SEEKR2 calculation.")
    
    argparser.add_argument(
        "-f", "--force_overwrite", dest="force_overwrite", default=False,
        help="Toggle whether to overwrite existing simulation output files "\
        "within any anchor that might have existed in an old model that would "\
        "be overwritten by generating this new model. If not toggled, this "\
        "program will throw an exception instead of performing any such "\
        "overwrite.", action="store_true")
    argparser.add_argument(
        "-s", "--skip_checks", dest="skip_checks", default=False, 
        help="By default, pre-simulation checks will be run after the "\
        "preparation is complete, and if the checks fail, the SEEKR2 "\
        "model will not be saved. This argument bypasses those "\
        "checks and allows the model to be generated anyways.", 
        action="store_true")
        
    args = argparser.parse_args() # parse the args into a dictionary
    args = vars(args)
    model_input_filename = args["model_input_file"]
    force_overwrite = args["force_overwrite"]
    skip_checks = args["skip_checks"]
    model_input = common_prepare.Model_input()
    model_input.deserialize(model_input_filename, user_input = True)
    model, xml_path = generate_seekr2_model_and_filetree(
        model_input, force_overwrite)
    if model.anchor_rootdir == ".":
        model_dir = os.path.dirname(xml_path)
        model.anchor_rootdir = os.path.abspath(model_dir)
    if not skip_checks:
        check.check_pre_simulation_all(model)