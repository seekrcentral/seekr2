"""
filetree.py

Generate the necessary directories and files to execute SEEKR2 
simulations.
"""

import os
import shutil
from shutil import copyfile
import numpy as np

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

def generate_filetree(model, rootdir, empty_rootdir=False):
    """
    Create (or recreate) the filetree within a model, including the
    anchor directories and subdirectories.
    
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
        
    for anchor in model.anchors:
        if not anchor.md:
            continue
        anchor_dict = {}
        anchor_dict[anchor.production_directory] = {}
        anchor_dict[anchor.building_directory] = {}
        
        anchor_filetree = Filetree({anchor.name:anchor_dict})
        anchor_filetree.make_tree(rootdir)
        
    if model.k_on_info is not None:
        
        b_surface_dict = {}
        b_surface_filetree = Filetree(
            {model.k_on_info.b_surface_directory:b_surface_dict})
        b_surface_filetree.make_tree(rootdir)
        for bd_milestone in model.k_on_info.bd_milestones:
            bd_milestone_dict = {}
            bd_milestone_dict[bd_milestone.extracted_directory] = {}
            bd_milestone_dict[bd_milestone.fhpd_directory] = {}
            
            bd_milestone_filetree = Filetree(
                {bd_milestone.directory:bd_milestone_dict})
            bd_milestone_filetree.make_tree(rootdir)
    
    return

def copy_building_files(model, input_model, rootdir):
    """
    For each of the anchors and other directories, copy the necessary
    files for the simulations into those directories.
    
    Parameters:
    -----------
    model : list
        the Model to prepare files for
        
    input_anchors : list
        An equivalent list of the Input_model's Input_anchor() objects.
        Which contain the input files to copy.
        
    rootdir : str
        A path to the model's root directory.
    
    """
    input_anchors = input_model.cv_inputs[0].input_anchors
    for anchor, input_anchor in zip(model.anchors, input_anchors):
        if not anchor.md:
            continue
        anchor_building_dir = os.path.join(rootdir, anchor.directory, 
                                        anchor.building_directory)
        assert os.path.exists(anchor_building_dir)
        
        try: # TODO: fix simple XML parser so this isn't necessary
            amber = input_anchor.starting_amber_params
        except AttributeError:
            amber = None
        
        if amber is not None:
            anchor.amber_params = base.Amber_params()
            #assert amber.prmtop_filename is not None, \
            #    "Amber PRMTOP file must be provided for anchor {}".format(
            #        anchor.index)
            #assert amber.prmtop_filename != "", \
            #    "Amber PRMTOP file must be provided for anchor {}".format(
            #        anchor.index)
            if amber.prmtop_filename is not None and \
                    amber.prmtop_filename != "":
                assert os.path.exists(amber.prmtop_filename)
                prmtop_filename = os.path.basename(amber.prmtop_filename)
                new_prmtop_filename = os.path.join(anchor_building_dir, 
                                                   prmtop_filename)
                copyfile(amber.prmtop_filename, new_prmtop_filename)
                anchor.amber_params.prmtop_filename = prmtop_filename
                
            """# TODO: remove
            if amber.inpcrd_filename is not None and \
                    amber.inpcrd_filename != "":
                assert os.path.exists(amber.inpcrd_filename)
                inpcrd_filename = os.path.basename(amber.inpcrd_filename)
                new_inpcrd_filename = os.path.join(anchor_building_dir, 
                                                   inpcrd_filename)
                copyfile(amber.inpcrd_filename, new_inpcrd_filename)
                anchor.amber_params.inpcrd_filename = inpcrd_filename
                if anchor.amber_params.box_vectors is None:
                    assert new_prmtop_filename is not None
                    anchor.amber_params.box_vectors = base.Box_vectors()
                    inpcrd_structure = parmed.load_file(new_prmtop_filename, 
                                                xyz=new_inpcrd_filename)
                    anchor.amber_params.box_vectors.from_quantity(
                        inpcrd_structure.box_vectors)
            """
            
            #assert amber.pdb_coordinates_filename is not None and \
            #    amber.pdb_coordinates_filename != "", \
            #    "PDB file must be provided for anchor {}".format(
            #        anchor.index)
            #amber.pdb_coordinates_filename = os.path.expanduser(
            #    amber.pdb_coordinates_filename)
            #assert os.path.exists(amber.pdb_coordinates_filename), \
            #    "Provided file does not exist: {}".format(
            #        amber.pdb_coordinates_filename)
            if amber.pdb_coordinates_filename is not None and \
                    amber.pdb_coordinates_filename != "":
                pdb_filename = os.path.basename(amber.pdb_coordinates_filename)
                new_pdb_filename = os.path.join(anchor_building_dir, 
                                                pdb_filename)
                copyfile(amber.pdb_coordinates_filename, new_pdb_filename)
                anchor.amber_params.pdb_coordinates_filename = pdb_filename
                if amber.box_vectors is not None:
                    flattened_box_vectors = np.array([
                        val for row in amber.box_vectors for val in row
                    ]) 
                    box_vectors = base.Box_vectors()
                    box_vectors.ax, box_vectors.ay, box_vectors.az, \
                    box_vectors.bx, box_vectors.by, box_vectors.bz, \
                    box_vectors.cx, box_vectors.cy, box_vectors.cz = \
                    flattened_box_vectors
                    amber.box_vectors = box_vectors
                anchor.amber_params.box_vectors = amber.box_vectors
                if anchor.amber_params.box_vectors is None:
                    anchor.amber_params.box_vectors = base.Box_vectors()
                    pdb_structure = parmed.load_file(new_pdb_filename)
                    anchor.amber_params.box_vectors.from_quantity(
                        pdb_structure.box_vectors)
            
        try: # TODO: fix simple XML parser so this isn't necessary
            forcefield = input_anchor.starting_forcefield_params
        except AttributeError:
            forcefield = None
            
        if forcefield is not None:
            assert amber is None, "Parameters may not be included for both "\
                "Amber and Forcefield inputs."
            anchor.forcefield_params = base.Forcefield_params()
            if forcefield.built_in_forcefield_filenames is not None and \
                    len(forcefield.built_in_forcefield_filenames) > 0:
                for filename in forcefield.built_in_forcefield_filenames:
                    anchor.forcefield_params.built_in_forcefield_filenames.\
                        append(filename)
                        
            if forcefield.custom_forcefield_filenames is not None and \
                    len(forcefield.custom_forcefield_filenames) > 0:
                for filename in forcefield.custom_forcefield_filenames:
                    ff_filename = os.path.basename(filename)
                    new_ff_filename = os.path.join(anchor_building_dir, 
                                                   ff_filename)
                    copyfile(filename, new_ff_filename)
                    anchor.forcefield_params.custom_forcefield_filenames.\
                        append(ff_filename)
            
            if forcefield.pdb_filename is not None and \
                    forcefield.pdb_filename != "":
                assert os.path.exists(forcefield.pdb_filename)
                pdb_filename = os.path.basename(forcefield.pdb_filename)
                new_pdb_filename = os.path.join(anchor_building_dir, 
                                                pdb_filename)
                copyfile(forcefield.pdb_filename, new_pdb_filename)
                anchor.forcefield_params.pdb_filename = pdb_filename
                anchor.forcefield_params.box_vectors = forcefield.box_vectors
                if anchor.forcefield_params.box_vectors is None:
                    pdb_structure = parmed.load_file(new_pdb_filename)
                    anchor.forcefield_params.box_vectors.from_quantity(
                        pdb_structure.box_vectors)
                
    if model.k_on_info is not None:
        bd_settings = model.browndye_settings
        bd_input_settings = input_model.browndye_settings_input
        k_on_info = model.k_on_info
        b_surface_dir = os.path.join(rootdir, k_on_info.b_surface_directory)
        
        ligand_pqr_filename = os.path.basename(bd_settings.ligand_pqr_filename)
        ligand_pqr_dest_filename = os.path.join(
            b_surface_dir, ligand_pqr_filename)
        copyfile(bd_input_settings.ligand_pqr_filename, 
                 ligand_pqr_dest_filename)
        
        receptor_pqr_filename = os.path.basename(
            bd_settings.receptor_pqr_filename)
        receptor_pqr_dest_filename = os.path.join(
            b_surface_dir, receptor_pqr_filename)
        copyfile(bd_input_settings.receptor_pqr_filename, 
                 receptor_pqr_dest_filename)
    return
