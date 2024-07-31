"""
filetree.py

Generate the necessary directories and files to execute SEEKR2 
simulations.
"""

import os
import shutil
from shutil import copyfile

import parmed

import seekr2.modules.common_base as base

class Filetree():
    """
    Define a file tree: a framework of directories to be populated with
    files.
    
    Attributes:
    -----------
    tree : dict
        a nested dictionary object defining the file structure to create.
        example: {'file1':{}, 'file2':{'file3':{}}}
    
    """
    
    def __init__(self, tree):
        self.tree = tree
        return

    def make_tree(self, rootdir, branch={}):
        """
        Construct the file tree given a root directory, populating
        it with all subsequent branches and leaves in the file tree.
        
        Parameters:
        -----------
        rootdir : str
            a string indicating the path in which to place the root of
            the file tree 
        
        branch: optional nested dictionary structure
            of the tree. (see Filetree.__doc__)
        
        Returns:
        --------
        None
        """
        assert os.path.isdir(rootdir), "rootdir argument must be a "\
            "real directory"
        if not branch: branch = self.tree
        for subbranch in list(branch.keys()):
            # first create each subbranch
            subbranchpath = os.path.join(rootdir,subbranch)
            if not os.path.exists(subbranchpath):
                os.mkdir(subbranchpath)
            if not branch[subbranch]: 
                # if its an empty list, then we have a leaf
                continue
            else: 
                # then we can descend further, recursively
                self.make_tree(subbranchpath, branch=branch[subbranch])
        return

def generate_filetree_root(model, rootdir, empty_rootdir=False):
    """
    Create (or recreate) the lowest directory of a filetree within 
    a model.
    
    Parameters:
    -----------
    model : Model()
        The MMVT Model to generate an entire filetree for
    
    rootdir : str
        A string path to the Model's root directory, where all the
        anchors and other directories will be located.
    
    empty_rootdir : bool, default False
        Whether to delete an existing Model root directory, if it
        exists.
    
    """
    if empty_rootdir and os.path.isdir(rootdir):
        print("Deleting all subdirectories/files in rootdir:", rootdir)
        shutil.rmtree(rootdir)
        
    if not os.path.isdir(rootdir):
        os.mkdir(rootdir)
    
    return

def generate_filetree_by_anchor(anchor, rootdir):
    """
    Create (or recreate) the filetree within a model, including the
    anchor directories and subdirectories.
    
    Parameters:
    -----------
    anchor : Anchor()
        The anchor to generate the filetree for
    
    rootdir : str
        A string path to the Model's root directory, where all the
        anchors and other directories will be located.
    
    """
    if not anchor.md:
        return
    anchor_dict = {}
    anchor_dict[anchor.production_directory] = {}
    anchor_dict[anchor.building_directory] = {}    
    anchor_filetree = Filetree({anchor.directory:anchor_dict})
    anchor_filetree.make_tree(rootdir)
    return

def generate_filetree_bd(model, rootdir):
    """
    Create (or recreate) the filetree used in BD calculations.
    
    Parameters:
    -----------
    model : Model()
        The MMVT Model to generate an entire filetree for
    
    rootdir : str
        A string path to the Model's root directory, where all the
        anchors and other directories will be located.
    
   
    """
    if model.using_bd():
        b_surface_dict = {}
        b_surface_filetree = Filetree(
            {model.k_on_info.b_surface_directory:b_surface_dict})
        b_surface_filetree.make_tree(rootdir)
    
    return

def copy_building_files_by_anchor(anchor, input_anchor, rootdir):
    """
    For each of the anchors and other directories, copy the necessary
    initial simulation (building) files for the simulations into those
    directories.
    
    Parameters:
    -----------
    anchor : Anchor
        the Anchor to copy the building files into.
        
    input_anchor : Input_anchor
        Input_model's Input_anchor() objects which contain the input 
        files to copy.
        
    rootdir : str
        A path to the model's root directory.
    
    """
    if not anchor.md:
        return
    
    if anchor.__class__.__name__ in ["MMVT_toy_anchor", "Elber_toy_anchor"]:
        anchor.starting_positions = input_anchor.starting_positions
        return
    
    anchor_building_dir = os.path.join(rootdir, anchor.directory, 
                                    anchor.building_directory)
    assert os.path.exists(anchor_building_dir)
    
    amber = input_anchor.starting_amber_params
    new_prmtop_filename = None
    if amber is not None:
        anchor.amber_params = base.Amber_params()
        if amber.prmtop_filename is not None and \
                amber.prmtop_filename != "":
            amber.prmtop_filename = os.path.expanduser(
                amber.prmtop_filename)
            assert os.path.exists(amber.prmtop_filename), \
                "Provided file does not exist: {}".format(
                    amber.prmtop_filename)
            prmtop_filename = os.path.basename(amber.prmtop_filename)
            new_prmtop_filename = os.path.join(anchor_building_dir, 
                                               prmtop_filename)
            copyfile(amber.prmtop_filename, new_prmtop_filename)
            anchor.amber_params.prmtop_filename = prmtop_filename
            
        if amber.pdb_coordinates_filename is not None and \
                amber.pdb_coordinates_filename != "":
            pdb_filename = os.path.basename(amber.pdb_coordinates_filename)
            new_pdb_filename = os.path.join(anchor_building_dir, 
                                            pdb_filename)
            copyfile(os.path.expanduser(amber.pdb_coordinates_filename), 
                     new_pdb_filename)
            anchor.amber_params.pdb_coordinates_filename = pdb_filename
            anchor.amber_params.box_vectors = amber.box_vectors
            if anchor.amber_params.box_vectors is None:
                anchor.amber_params.box_vectors = base.Box_vectors()
                box_vectors = base.get_box_vectors_from_pdb(new_pdb_filename)
                anchor.amber_params.box_vectors.from_quantity(
                    box_vectors)
        
    forcefield = input_anchor.starting_forcefield_params
        
    if forcefield is not None:
        assert amber is None, "Parameters may not be included for both "\
            "Amber and Forcefield inputs."
        anchor.forcefield_params = base.Forcefield_params()
        if forcefield.built_in_forcefield_filenames is not None and \
                len(forcefield.built_in_forcefield_filenames) > 0:
            for filename in forcefield.built_in_forcefield_filenames:
                anchor.forcefield_params.built_in_forcefield_filenames.\
                    append(os.path.expanduser(filename))
                    
        if forcefield.custom_forcefield_filenames is not None and \
                len(forcefield.custom_forcefield_filenames) > 0:
            for filename in forcefield.custom_forcefield_filenames:
                ff_filename = os.path.basename(filename)
                new_ff_filename = os.path.join(anchor_building_dir, 
                                               ff_filename)
                copyfile(filename, new_ff_filename)
                anchor.forcefield_params.custom_forcefield_filenames.\
                    append(os.path.expanduser(ff_filename))
        
        if forcefield.pdb_coordinates_filename is not None and \
                forcefield.pdb_coordinates_filename != "":
            forcefield.pdb_coordinates_filename = os.path.expanduser(
                forcefield.pdb_coordinates_filename)
            assert os.path.exists(forcefield.pdb_coordinates_filename), \
                "Provided file does not exist: {}".format(
                    forcefield.pdb_coordinates_filename)
            pdb_filename = os.path.basename(forcefield.pdb_coordinates_filename)
            new_pdb_filename = os.path.join(anchor_building_dir, 
                                            pdb_filename)
            copyfile(forcefield.pdb_coordinates_filename, new_pdb_filename)
            anchor.forcefield_params.pdb_coordinates_filename = pdb_filename
            anchor.forcefield_params.box_vectors = forcefield.box_vectors
            if anchor.forcefield_params.box_vectors is None:
                anchor.forcefield_params.box_vectors = base.Box_vectors()
                box_vectors = base.get_box_vectors_from_pdb(new_pdb_filename)
                anchor.forcefield_params.box_vectors.from_quantity(
                    box_vectors)
                
        if forcefield.system_filename is not None and \
                forcefield.system_filename != "":
            forcefield.system_filename = os.path.expanduser(
                forcefield.system_filename)
            assert os.path.exists(forcefield.system_filename), \
                "Provided file does not exist: {}".format(
                    forcefield.system_filename)
            system_filename = os.path.basename(forcefield.system_filename)
            new_system_filename = os.path.join(anchor_building_dir, 
                                               system_filename)
            copyfile(forcefield.system_filename, new_system_filename)
            anchor.forcefield_params.system_filename = system_filename
                
    charmm = input_anchor.starting_charmm_params
    new_psf_filename = None
    if charmm is not None:
        anchor.charmm_params = base.Charmm_params()
        if charmm.psf_filename is not None and \
                charmm.psf_filename != "":
            charmm.psf_filename = os.path.expanduser(
                charmm.psf_filename)
            assert os.path.exists(charmm.psf_filename), \
                "Provided file does not exist: {}".format(
                    charmm.psf_filename)
            psf_filename = os.path.basename(charmm.psf_filename)
            new_psf_filename = os.path.join(anchor_building_dir, 
                                               psf_filename)
            copyfile(charmm.psf_filename, new_psf_filename)
            anchor.charmm_params.psf_filename = psf_filename
            
        if charmm.charmm_ff_files is not None \
               and len(charmm.charmm_ff_files) > 0:
            anchor.charmm_params.charmm_ff_files = []
            for input_charmm_ff_filename in charmm.charmm_ff_files:
                input_charmm_ff_filename = os.path.expanduser(
                    input_charmm_ff_filename)
                assert os.path.exists(input_charmm_ff_filename), \
                    "Provided file does not exist: {}".format(
                        input_charmm_ff_filename)
                charmm_ff_filename = os.path.basename(input_charmm_ff_filename)
                new_ff_filename = os.path.join(anchor_building_dir, 
                                                   charmm_ff_filename)
                copyfile(input_charmm_ff_filename, new_ff_filename)
                anchor.charmm_params.charmm_ff_files.append(charmm_ff_filename)
            
        if charmm.pdb_coordinates_filename is not None and \
                charmm.pdb_coordinates_filename != "":
            pdb_filename = os.path.basename(charmm.pdb_coordinates_filename)
            new_pdb_filename = os.path.join(anchor_building_dir, 
                                            pdb_filename)
            copyfile(os.path.expanduser(charmm.pdb_coordinates_filename), 
                     new_pdb_filename)
            anchor.charmm_params.pdb_coordinates_filename = pdb_filename
            anchor.charmm_params.box_vectors = charmm.box_vectors
            if anchor.charmm_params.box_vectors is None:
                anchor.charmm_params.box_vectors = base.Box_vectors()
                box_vectors = base.get_box_vectors_from_pdb(new_pdb_filename)
                anchor.charmm_params.box_vectors.from_quantity(
                    box_vectors)
    
    return
    
def copy_bd_files(model, input_model, rootdir):
    """
    Copy the necessary files for the BD simulations into those 
    directories.
    
    Parameters:
    -----------
    model : Model
        the Model to prepare files for.
        
    input_model : Input_Model
        The Input_model's object which contains the input files
        information to copy over.
        
    rootdir : str
        A path to the model's root directory.
    
    """
    if model.using_bd():
        bd_settings = model.browndye_settings
        bd_input_settings = input_model.browndye_settings_input
        k_on_info = model.k_on_info
        b_surface_dir = os.path.join(rootdir, k_on_info.b_surface_directory)
        
        ligand_pqr_filename = os.path.basename(bd_settings.ligand_pqr_filename)
        ligand_pqr_dest_filename = os.path.join(
            b_surface_dir, ligand_pqr_filename)
        copyfile(os.path.expanduser(bd_input_settings.ligand_pqr_filename), 
                 ligand_pqr_dest_filename)
        
        receptor_pqr_filename = os.path.basename(
            bd_settings.receptor_pqr_filename)
        receptor_pqr_dest_filename = os.path.join(
            b_surface_dir, receptor_pqr_filename)
        copyfile(os.path.expanduser(bd_input_settings.receptor_pqr_filename), 
                 receptor_pqr_dest_filename)
    
    return