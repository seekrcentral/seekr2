"""
check.py

Check SEEKR2 calculations for common problems both before and after
simulations.

check_pre_simulation_all : 
    These are checks that should be done before running any
    simulations (run stage). They are automatically run by prepare.py.
 - Check for bubbles in starting structures, or unusually low densities
   (high densities will probably fail on their own). Check that box
   boundaries aren't significantly larger than starting structure 
   boxes.
 - Check to ensure that the same (or similar) salt concentrations exist
   between the MD and BD stages.
 - Check that the system exists within the expected Voronoi cell and
   give detailed information of what is wrong if it is not in there.
 - Check that atom selections exist on the same molecule, and that
   atoms from water molecules are not included accidentally.
 - Check that BD milestone's atom selections are the same between
   MD and BD.

check_post_simulation_all : 
    These are checks that should be performed after any or all
    simulations are performed, and before the analysis stage.
    They are automatically run by analyze.py
 - Check that all simulation trajectories keep the system where
   expected: 
    - Elber umbrella stage stays close to milestone.
    - Elber reversal/forward stage ends close to the milestone it was
      recorded to have crossed.
    - MMVT trajectories stay within the Voronoi cell, and they end
      close to the boundary recorded by the crossing.
  - Re-check simulation density? ***
  - Check that BD simulations end on the correct milestone
  
This module maybe also be run from the command line. Example:

$ python check.py /path/to/model.xml
"""

import os
import sys
import argparse
import collections
import glob
import warnings
warnings.filterwarnings("ignore")
import tempfile
#import bubblebuster

import numpy as np
import parmed
import mdtraj

import seekr2.modules.common_base as base
import seekr2.modules.elber_base as elber_base
import seekr2.modules.mmvt_base as mmvt_base
import seekr2.modules.mmvt_sim_openmm as mmvt_sim_openmm
import seekr2.modules.elber_sim_openmm as elber_sim_openmm
import seekr2.modules.common_converge as common_converge

# The charges of common ions
ION_CHARGE_DICT = {"li":1.0, "na":1.0, "k":1.0, "rb":1.0, "cs":1.0, "fr":1.0,
                   "be":2.0, "mg":2.0, "ca":2.0, "sr":2.0, "ba":2.0, "ra":2.0,
                   "f":-1.0, "cl":-1.0, "br":-1.0, "i":-1.0, "at":-1.0,
                   "he":0.0, "ne":0.0, "ar":0.0, "kr":0.0, "xe":0.0, "rn":0.0}

AVOGADROS_NUMBER = 6.022e23
SAVE_STATE_DIRECTORY = "states/"
MAX_STRUCTURES_TO_CHECK = 100
RECURSION_LIMIT = 100000
MAX_BD_STUCK = 5

def load_structure_with_parmed(model, anchor):
    """
    Given the simulation inputs, load an anchor's structure for one of
    the checks and return the parmed structure object.
    """
    
    if anchor.amber_params is not None:
        building_directory = os.path.join(
            model.anchor_rootdir, anchor.directory, 
            anchor.building_directory)
        if anchor.amber_params.prmtop_filename is not None:
            prmtop_filename = os.path.join(
                building_directory, anchor.amber_params.prmtop_filename)
        else:
            return None
        
        if anchor.amber_params.pdb_coordinates_filename is not None \
                and anchor.amber_params.pdb_coordinates_filename != "":
            pdb_filename = os.path.join(
                building_directory, 
                anchor.amber_params.pdb_coordinates_filename)
            structure = parmed.load_file(pdb_filename)
            if anchor.amber_params.box_vectors is not None:
                structure.box_vectors = anchor.amber_params.box_vectors.to_quantity()
        else:
            return None

        return structure
        
    elif anchor.forcefield_params is not None:
        forcefield_filenames = []
        for forcefield_filename in \
                anchor.forcefield_params.built_in_forcefield_filenames:
            forcefield_filenames.append(forcefield_filename)
        for forcefield_filename in \
                anchor.forcefield_params.custom_forcefield_filenames:
            forcefield_filenames.append(os.path.join(
                building_directory, forcefield_filename))
        parameter_set = parmed.openmm.parameters.OpenMMParameterSet(
            forcefield_filenames)
        pdb_filename = os.path.join(building_directory, 
                               anchor.forcefield_params.pdb_filename)
        structure = parmed.load_file(parameter_set)
        pdb_structure = parmed.load_file(pdb_filename)
        structure.coordinates = pdb_structure.coordinates
        if anchor.forcefield_params.box_vectors is not None:
                structure.box_vectors \
                    = anchor.forcefield_params.box_vectors.to_quantity()
        return structure
    
    elif anchor.charmm_params is not None:
        building_directory = os.path.join(
            model.anchor_rootdir, anchor.directory, 
            anchor.building_directory)
        if anchor.charmm_params.psf_filename is not None:
            psf_filename = os.path.join(
                building_directory, anchor.charmm_params.psf_filename)
        else:
            return None
        
        if anchor.charmm_params.pdb_coordinates_filename is not None \
                and anchor.charmm_params.pdb_coordinates_filename != "":
            pdb_filename = os.path.join(
                building_directory, 
                anchor.charmm_params.pdb_coordinates_filename)
            structure = parmed.load_file(pdb_filename)
            if anchor.charmm_params.box_vectors is not None:
                structure.box_vectors = anchor.charmm_params\
                    .box_vectors.to_quantity()
        else:
            return None

        return structure
        
    else:
        return None

def load_structure_with_mdtraj(model, anchor, mode="pdb", coords_filename=None):
    """
    Given the simulation inputs, load an anchor's structure for one of
    the checks and return the mdtraj Trajectory() object.
    """
    
    building_directory = os.path.join(
        model.anchor_rootdir, anchor.directory, anchor.building_directory)
    prod_directory = os.path.join(
        model.anchor_rootdir, anchor.directory, anchor.production_directory)
    if mode == "pdb":
        pass
    elif mode == "elber_umbrella":
        umbrella_dcd_path = os.path.join(prod_directory)
        umbrella_basename = elber_base.ELBER_UMBRELLA_BASENAME+"*.dcd"
        umbrella_traj_glob = os.path.join(umbrella_dcd_path, umbrella_basename)
        umbrella_traj_filenames = glob.glob(umbrella_traj_glob)
        if len(umbrella_traj_filenames) == 0:
            return None
        
        # search for and remove any empty trajectory files
        indices_to_pop = []
        for i, umbrella_traj_filename in enumerate(umbrella_traj_filenames):
            if os.path.getsize(umbrella_traj_filename) == 0:
                indices_to_pop.append(i)
        
        for i in indices_to_pop[::-1]:
            umbrella_traj_filenames.pop(i)
            
        assert len(umbrella_traj_filenames) > 0, "Only empty umbrella " \
        "trajectories were found. You can force SEEKR to skip these checks "\
        "by using the --skip_checks (-s) argument"
    
    elif mode == "state_xml":
        assert coords_filename is not None
        
    elif mode == "mmvt_traj":
        mmvt_traj_basename = mmvt_base.OPENMMVT_BASENAME+"*.dcd"
        mmvt_traj_glob = os.path.join(prod_directory, mmvt_traj_basename)
        mmvt_traj_filenames = glob.glob(mmvt_traj_glob)
        if len(mmvt_traj_filenames) == 0:
            return None
        
        # search for and remove any empty trajectory files
        indices_to_pop = []
        for i, mmvt_traj_filename in enumerate(mmvt_traj_filenames):
            if os.path.getsize(mmvt_traj_filename) == 0:
                indices_to_pop.append(i)
        
        for i in indices_to_pop[::-1]:
            mmvt_traj_filenames.pop(i)
    else:
        raise Exception("Check mode not implemented: {}".format(mode))
        
    if model.using_toy():
        pdb_filename = os.path.join(building_directory, "toy.pdb")
        if mode == "pdb":
            traj = mdtraj.load(pdb_filename)
            
        elif mode == "elber_umbrella":
            assert len(umbrella_traj_filenames) > 0, "Only empty umbrella " \
                "trajectories were found. You can force SEEKR to skip these "\
                "checks by using the --skip_checks (-s) argument. Anchor: {}"\
                .format(anchor.index)
            traj = mdtraj.load(umbrella_traj_filenames, top=pdb_filename)
            
        elif mode == "state_xml":
            traj = mdtraj.load_xml(coords_filename, top=pdb_filename)
        
        elif mode == "mmvt_traj":
            if not len(mmvt_traj_filenames) > 0:
                warnings.warn("Empty mmvt trajectories were found in "\
                      "anchor: {}.".format(anchor.index))
                return None
                
            traj = mdtraj.load(mmvt_traj_filenames, top=pdb_filename)
            
        return traj
    
    elif anchor.amber_params is not None:
        
        if anchor.amber_params.prmtop_filename is not None:
            prmtop_filename = os.path.join(
                building_directory, anchor.amber_params.prmtop_filename)
        else:
            return None
        if mode == "pdb":
            if anchor.amber_params.pdb_coordinates_filename is not None \
                    and anchor.amber_params.pdb_coordinates_filename != "":
                pdb_filename = os.path.join(
                    building_directory, 
                    anchor.amber_params.pdb_coordinates_filename)
                traj = mdtraj.load(pdb_filename) #, top=prmtop_filename)
            else:
                # anchor has no structure files
                return None
        
        elif mode == "elber_umbrella":
            assert len(umbrella_traj_filenames) > 0, "Only empty umbrella " \
                "trajectories were found. You can force SEEKR to skip these "\
                "checks by using the --skip_checks (-s) argument. Anchor: {}"\
                .format(anchor.index)
            traj = mdtraj.load(umbrella_traj_filenames, top=prmtop_filename)
        
        elif mode == "state_xml":
            traj = mdtraj.load_xml(coords_filename, top=prmtop_filename)
        
        elif mode == "mmvt_traj":
            if not len(mmvt_traj_filenames) > 0:
                warnings.warn("Empty mmvt trajectories were found in "\
                      "anchor: {}.".format(anchor.index))
                return None
                
            traj = mdtraj.load(mmvt_traj_filenames, top=prmtop_filename)
            
        return traj
        
    elif anchor.forcefield_params is not None:
        forcefield_filenames = []
        for forcefield_filename in \
                anchor.forcefield_params.built_in_forcefield_filenames:
            forcefield_filenames.append(forcefield_filename)
        for forcefield_filename in \
                anchor.forcefield_params.custom_forcefield_filenames:
            forcefield_filenames.append(os.path.join(
                building_directory, forcefield_filename))
        parameter_set = parmed.openmm.parameters.OpenMMParameterSet(
            forcefield_filenames)
        pdb_filename = os.path.join(building_directory, 
                               anchor.forcefield_params.pdb_filename)
        if mode == "pdb":
            traj = mdtraj.load(pdb_filename)
        elif mode == "elber_umbrella":
            assert len(umbrella_traj_filenames) > 0, "Only empty umbrella " \
                "trajectories were found. You can force SEEKR to skip these "\
                "checks by using the --skip_checks (-s) argument"
            traj = mdtraj.load(umbrella_traj_filenames, top=pdb_filename)
        elif mode == "state_xml":
            traj = mdtraj.load_xml(coords_filename, top=prmtop_filename)
        elif mode == "mmvt_traj":
            assert len(mmvt_traj_filenames) > 0, "Only empty mmvt " \
                "trajectories were found. You can force SEEKR to skip these "\
                "checks by using the --skip_checks (-s) argument"
            traj = mdtraj.load(mmvt_traj_filenames, top=pdb_filename)
        
        return traj
    
    elif anchor.charmm_params is not None:
        if anchor.charmm_params.psf_filename is not None:
            psf_filename = os.path.join(
                building_directory, anchor.charmm_params.psf_filename)
        else:
            return None
        if mode == "pdb":
            if anchor.charmm_params.pdb_coordinates_filename is not None \
                    and anchor.charmm_params.pdb_coordinates_filename != "":
                pdb_filename = os.path.join(
                    building_directory, 
                    anchor.charmm_params.pdb_coordinates_filename)
                traj = mdtraj.load(pdb_filename) #, top=prmtop_filename)
            else:
                # anchor has no structure files
                return None
        
        elif mode == "elber_umbrella":
            assert len(umbrella_traj_filenames) > 0, "Only empty umbrella " \
                "trajectories were found. You can force SEEKR to skip these "\
                "checks by using the --skip_checks (-s) argument. Anchor: {}"\
                .format(anchor.index)
            traj = mdtraj.load(umbrella_traj_filenames, top=psf_filename)
        
        elif mode == "state_xml":
            traj = mdtraj.load_xml(coords_filename, top=psf_filename)
        
        elif mode == "mmvt_traj":
            if not len(mmvt_traj_filenames) > 0:
                warnings.warn("Empty mmvt trajectories were found in "\
                      "anchor: {}.".format(anchor.index))
                return None
                
            traj = mdtraj.load(mmvt_traj_filenames, top=psf_filename)
            
        return traj
        
    else:
        return None

def is_ion(atom):
    """If a lone atom has no bonds, assume it's an ion."""
    
    if len(atom.bond_partners) == 0:
        return True
    else:
        return False

def check_pre_sim_bubbles(model):
    """
    Checks starting pdb structures for water box bubbles.
    
    Parameters
    ----------
    model : Model
        The SEEKR2 model object containing all calculation information.
    
    Returns
    -------
    bool
        False if a bubble is detected and True otherwise.
    """
    pass
    """
    for anchor in model.anchors:
        building_directory = os.path.join(
            model.anchor_rootdir, 
            anchor.directory, 
            anchor.building_directory
        )
        if anchor.amber_params is not None:
            if (anchor.amber_params.pdb_coordinates_filename is not None \
                and anchor.amber_params.pdb_coordinates_filename != ""
            ):
                pdb_filename = os.path.join(
                    building_directory,
                    anchor.amber_params.pdb_coordinates_filename
                )
                bvecs = anchor.amber_params.box_vectors
                if bvecs is not None:
                    bvecs = bvecs.to_quantity().value_in_unit(
                        parmed.unit.nanometers)
                    
                    pdb_periodic_box_properties = \
                        bubblebuster.periodic_box_properties(
                            pdb_filename,
                            mesh=0.7,
                            box_vectors=np.array(bvecs, dtype=np.float32)
                        )
                    if pdb_periodic_box_properties.has_bubble:
                        warnstr = "CHECK FAILURE: Water box bubble detected "\
                        "in starting structure: {}".format(pdb_filename) 
                        print(warnstr)
                    else:
                        pdb_periodic_box_properties = \
                            bubblebuster.periodic_box_properties(
                                pdb_filename,
                                mesh=0.7,
                                box_vectors=np.array(bvecs, dtype=np.float32)
                            )
                        if pdb_periodic_box_properties.has_bubble:
                            warnstr = "CHECK FAILURE: Water box bubble detected "\
                            "in starting structure: {}.".format(pdb_filename) 
                            print(warnstr)
                            return False

                else:
                    pdb_periodic_box_properties = \
                        bubblebuster.periodic_box_properties(
                            pdb_filename,
                            mesh=0.7
                        )
                    if pdb_periodic_box_properties.has_bubble:
                        warnstr = "CHECK FAILURE: Water box bubble detected "\
                        "in one or more starting structures." 
                        print(warnstr)
                        return False
    return True
    """

def check_pre_sim_MD_and_BD_salt_concentration(model):
    """
    Users might inadvertently define different salt concentrations 
    between the MD and BD stages. 
    Examine BD settings and count the numbers of ions in MD starting
    structures to ensure that ion concentrations are relatively 
    consistent to avoid situations where different salt concentrations
    exist between the various anchors and scales.
    """
    
    RELATIVE_TOLERANCE = 0.1
    ABSOLUTE_TOLERANCE = 0.01
    if model.using_bd():
        bd_ionic_strength = 0.0
        for ion in model.k_on_info.ions:
            bd_ionic_strength += 0.5 * ion.conc * ion.charge**2
    else:
        # No BD to check
        return True
    
    for anchor in model.anchors:
        md_ionic_strength = 0.0
        structure = load_structure_with_parmed(model, anchor)
        if structure is None:
            continue
        box_6_vector = structure.get_box()
        assert box_6_vector is not None, "Unable to compute box volume for "\
            "structures in anchor {}".format(anchor.index)
        box_vectors = base.Box_vectors()
        box_vectors.from_6_vector(box_6_vector[0])
        box_volume = box_vectors.get_volume() # volume in nm^3
        particle_concentration = 1.0e24 / (AVOGADROS_NUMBER * box_volume)
        for index, atom in enumerate(structure.atoms):
            if is_ion(atom):
                found_ion = False
                for key in ION_CHARGE_DICT:
                    if atom.name.lower().startswith(key):
                        charge = ION_CHARGE_DICT[key]
                        md_ionic_strength += 0.5 * particle_concentration \
                            * charge**2
                        found_ion = True
                    
                if not found_ion:
                    print("check_pre_sim_MD_and_BD_salt_concentration - Found "\
                          "unbonded atom with anchor {}, ".format(anchor.index)\
                          +"index: {} and name: {}. ".format(index, atom.name)\
                          +"Charge is uncertain. Assuming its contribution to "\
                          "solvent ionic strength is zero.")
        if not np.isclose(md_ionic_strength, bd_ionic_strength, 
                          rtol=RELATIVE_TOLERANCE, atol=ABSOLUTE_TOLERANCE):
            print("""CHECK FAILURE: BD simulation has significantly different
    ionic strength of {:.3f} M*e^2 than MD simulation ionic strength of 
    {:.3f} M*e^2 for anchor {}. Please check the ion concentrations in the 
    BD simulation settings, and also count the number of ions 
    in the MD simulations. Counterions added to neutralize the protein in an
    MD simulation may also trigger this failure. If that is the case, simply 
    skip this check.""".format(bd_ionic_strength, 
                                 md_ionic_strength,
                                 anchor.index))
            return False
        
    return True

def check_systems_within_Voronoi_cells(model):
    """
    Users might provide the wrong starting structure for a given 
    anchor, and the structure may actually belong in a different one.
    When the model defines Voronoi cells (such as in MMVT), check that
    the starting structures lie within the expected Voronoi cells, and
    suggest corrections if the check fails. Otherwise, the SEEKR2
    backend would fail with a non-helpful error message.
    """
    returning_result = True
    if model.get_type() != "mmvt":
        # only apply to MMVT systems
        return True
    
    curdir = os.getcwd()
    os.chdir(model.anchor_rootdir)
    
    for anchor in model.anchors:
        if model.using_toy():
            if anchor.starting_positions is None: continue
            if len(anchor.starting_positions) == 0: continue
            tmp_path = tempfile.NamedTemporaryFile()
            output_file = tmp_path.name
            model.openmm_settings.cuda_platform_settings = None
            model.openmm_settings.reference_platform = True
            if model.get_type() == "mmvt":
                my_sim_openmm = mmvt_sim_openmm.create_sim_openmm(
                    model, anchor, output_file)
                context = my_sim_openmm.simulation.context
            else:
                my_sim_openmm = elber_sim_openmm.create_sim_openmm(
                    model, anchor, output_file)
                context = my_sim_openmm.umbrella_simulation.context
            
            
        else:
            traj = load_structure_with_mdtraj(model, anchor)
            if traj is None:
                continue
            
        for milestone in anchor.milestones:
            cv = model.collective_variables[milestone.cv_index]
            if model.using_toy():
                result = cv.check_openmm_context_within_boundary(
                    context, milestone.variables, verbose=True)
                
            else:
                result = cv.check_mdtraj_within_boundary(
                    traj, milestone.variables, verbose=True)
                
            if result == False:
                correct_anchor = None
                for anchor2 in model.anchors:
                    within_milestones = True
                    for milestone in anchor2.milestones:
                        if anchor.__class__.__name__ \
                                in ["MMVT_toy_anchor", "Elber_toy_anchor"]:
                            result2 = cv.check_openmm_context_within_boundary(
                                context, milestone.variables, verbose=True)
                            
                        else:
                            result2 = cv.check_mdtraj_within_boundary(
                                traj, milestone.variables)
                        if not result2:
                            within_milestones = False
                            break
                    if within_milestones:
                        correct_anchor = anchor2
                        break
                
                warnstr = """CHECK FAILURE: The provided initial starting 
    structure for anchor {} does not lie within the 
    anchor boundaries. The simulation would fail. Please 
    carefully check this anchor's starting structure, as 
    well as the collective variable's (CV) atom selections, 
    and anchor/milestone variables.""".format(anchor.index)
                print(warnstr)
                if correct_anchor is not None:
                    print("It looks like the failed structure might belong in "\
                          "anchor {}.".format(correct_anchor.index))
                returning_result = False
    
    os.chdir(curdir)
    return returning_result

def recurse_atoms(atom, _visited_indices=set()):
    """
    Recursively visit all atoms within a molecule for the purposes
    of determining the molecules (all sets of atoms connected by bonds)
    in the system.
    """
    
    _visited_indices.add(atom.idx)
    for bonded_atom in atom.bond_partners:
        if not bonded_atom.idx in _visited_indices:
            branch_indices = recurse_atoms(bonded_atom, _visited_indices)
            _visited_indices.update(branch_indices)
    return _visited_indices

def check_atom_selections_on_same_molecule(model):
    """
    The user might accidentally define atom selections that span
    multiple molecules. Check this possibility by finding all molecules
    in the system and ensure that atom selections only exist on one
    molecule.
    """
    sys.setrecursionlimit(RECURSION_LIMIT)
    warnstr1 = """CHECK FAILURE: the atom selection for collective variable 
    (CV) number {} is split over multiple molecules. Atom index {}, 
    which has the name {} and serial id {}, was the first atom to 
    be detected in a different molecule than the others. Please
    check the structure for anchor {} to ensure that atom indexing 
    is correct. Keep in mind: SEEKR2 atom indexing starts at 0 (but 
    PDB files often start with a different atom serial index, such
    as 1 or another number)."""
    warnstr2 = """CHECK FAILURE: the atom index {} for collective variable
    (CV) number {} does not exist in the structure for anchor {}."""
    for anchor in model.anchors:
        structure = load_structure_with_parmed(model, anchor)
        if structure is None:
            continue
        molecules = []
        traversed_atoms = set()
        for index, atom in enumerate(structure.atoms):
            if index in traversed_atoms:
                continue
            molecule = recurse_atoms(atom, set())
            traversed_atoms.update(molecule)
            molecules.append(molecule)
        
        for cv in model.collective_variables:
            cv_in_anchor = False
            for milestone in anchor.milestones:
                if milestone.cv_index == cv.index:
                    cv_in_anchor = True
                    
            if not cv_in_anchor:
                continue
            
            if isinstance(cv, mmvt_base.MMVT_closest_pair_CV) \
                    or isinstance(cv, mmvt_base.MMVT_count_contacts_CV):
                # Closest pair and count contacts CVs use nonbonded pairs.
                continue
            
            atom_groups = cv.get_atom_groups()
            for atom_group in atom_groups:
                in_molecule = None
                for atom_index in atom_group:
                    index_found = False
                    for mol_id, molecule in enumerate(molecules):
                        if atom_index in molecule :
                            if in_molecule is None:
                                in_molecule = mol_id
                                index_found = True
                            elif in_molecule == mol_id:
                                index_found = True
                            else:
                                print(warnstr1.format(
                                    cv.index, atom_index, 
                                    structure.atoms[atom_index].name, 
                                    structure.atoms[atom_index].number, 
                                    anchor.index))
                                return False
                    if not index_found:
                        print(warnstr2.format(atom_index, cv.index, 
                                              anchor.index))
                        return False
    return True

def check_atom_selections_MD_BD(model):
    """
    Users might accidentally select atoms that are not equivalent
    between the MD and BD stages. Detect whether this is the case.
    """
    
    warnstr1 = """CHECK FAILURE: The atom selection as defined for BD
    milestone {} includes atom index {} that does not exist 
    in the file {}. Please check the atom indices in Browndye
    inputs and also check the atom indices in your input PQR
    files."""
    warnstr2 = """CHECK FAILURE: Different atoms are selected between
    the MD and BD stages. The BD selection contains {} while the
    MD selection contains {}. Please check the input molecular
    structure files, as well as the input atom indices (keep in
    mind that SEEKR2 requires atom indexing to start from 0 in 
    each molecular structure input file, but PDB and PQR files
    might start their numbering at 1 or another number."""
    if model.using_bd():
        b_surface_dir = os.path.join(
            model.anchor_rootdir, model.k_on_info.b_surface_directory)
        rec_pqr_path = os.path.join(
            b_surface_dir, model.browndye_settings.receptor_pqr_filename)
        lig_pqr_path = os.path.join(
            b_surface_dir, model.browndye_settings.ligand_pqr_filename)
        rec_pqr_structure = parmed.load_file(rec_pqr_path)
        lig_pqr_structure = parmed.load_file(lig_pqr_path)
        for bd_index, bd_milestone in enumerate(model.k_on_info.bd_milestones):
            receptor_indices = bd_milestone.receptor_indices
            for receptor_index in receptor_indices:
                if receptor_index < 0 \
                        or receptor_index > len(rec_pqr_structure.atoms):
                    print(warnstr1.format(bd_index, receptor_index, 
                                          rec_pqr_path))
                    return False
            ligand_indices = bd_milestone.ligand_indices
            for ligand_index in ligand_indices:
                if ligand_index < 0 \
                        or ligand_index > len(lig_pqr_structure.atoms):
                    print(warnstr1.format(bd_index, ligand_index, 
                                          lig_pqr_path))
                    return False
    else:
        # No BD to check
        return True
    
    for anchor in model.anchors:
        md_structure = load_structure_with_parmed(model, anchor)
        if md_structure is None:
            continue
        md_bd_cv_pairs = []
        for md_cv in model.collective_variables:
            for bd_index, bd_milestone in enumerate(
                    model.k_on_info.bd_milestones):
                if bd_milestone.outer_milestone.cv_index == md_cv.index:
                    md_bd_cv_pairs.append([md_cv.index, bd_index])
        
        for md_bd_pair in md_bd_cv_pairs:
            md_cv = model.collective_variables[md_bd_pair[0]]
            bd_milestone = model.k_on_info.bd_milestones[md_bd_pair[1]]
            md_atom_groups = sorted(md_cv.get_atom_groups(), 
                                    key=lambda L: (len(L), L))
            bd_atom_groups = sorted([bd_milestone.receptor_indices, 
                                    bd_milestone.ligand_indices],
                                    key=lambda L: (len(L), L))
            # check atom selection lengths
            for md_atom_group, bd_atom_group in zip(
                    md_atom_groups, bd_atom_groups):
                if len(md_atom_group) != len(bd_atom_group):
                    bd_error_str = "{} atom(s)".format(len(bd_atom_group))
                    md_error_str = "{} atom(s)".format(len(md_atom_group))
                    print(warnstr2.format(bd_error_str, md_error_str))
                    return False
                md_atomic_name_dict = collections.defaultdict(int)
                bd_atomic_name_dict = collections.defaultdict(int)
                for md_atom_idx in md_atom_group:
                    md_atom = md_structure.atoms[md_atom_idx]
                    md_atomic_name_dict[md_atom.name] += 1
                rec_pqr_indices = bd_milestone.receptor_indices
                lig_pqr_indices = bd_milestone.ligand_indices
                if bd_atom_group == rec_pqr_indices:
                    for bd_atom_idx in rec_pqr_indices:
                        bd_atom = rec_pqr_structure.atoms[bd_atom_idx]
                        bd_atomic_name_dict[bd_atom.name] += 1
                elif bd_atom_group == lig_pqr_indices:
                    for bd_atom_idx in lig_pqr_indices:
                        bd_atom = lig_pqr_structure.atoms[bd_atom_idx]
                        bd_atomic_name_dict[bd_atom.name] += 1
                else:
                    raise Exception("BD/MD mismatch in atom groups")
                
                if md_atomic_name_dict != bd_atomic_name_dict:
                    bd_err_str = ""
                    for key in bd_atomic_name_dict:
                        this_str = "{} atoms with atomic number {}".format(
                            bd_atomic_name_dict[key], key)
                        bd_err_str += this_str
                    md_err_str = ""
                    for key in md_atomic_name_dict:
                        this_str = "{} atoms with atomic number {}".format(
                            md_atomic_name_dict[key], key)
                        md_err_str += this_str
                    print(warnstr2.format(bd_err_str, md_err_str))
                    return False
    return True

def check_pqr_residues(model):
    """
    Browndye will create lumped charges (test charges) out of entire
    residues. This check will suggest that users split their PQR
    atoms into individual residues to increase accuracy.
    """
    if model.using_bd():
        b_surface_dir = os.path.join(
            model.anchor_rootdir, model.k_on_info.b_surface_directory)
        rec_pqr_path = os.path.join(
            b_surface_dir, model.browndye_settings.receptor_pqr_filename)
        lig_pqr_path = os.path.join(
            b_surface_dir, model.browndye_settings.ligand_pqr_filename)
        rec_pqr_structure = parmed.load_file(rec_pqr_path)
        lig_pqr_structure = parmed.load_file(lig_pqr_path)
        rec_residues = []
        for residue in rec_pqr_structure.residues:
            if residue.name != "GHO":
                rec_residues.append(residue)
        lig_residues = []
        for residue in lig_pqr_structure.residues:
            if residue.name != "GHO":
                lig_residues.append(residue)
        
    else:
        # No BD to check
        return True
    
    warnstr1 = """CHECK FAILURE: The PQR file {} is comprised of a single
    residue, named {}. By default, Browndye2 will lump all the charges in
    a residue into a single "test charge". This feature can be overridden by
    numbering every atom in your PQR file(s) as a different residue. If you
    don't mind this, and wish to proceed, then run this script again, skipping
    all checks."""
    assert rec_residues != 0, "Empty PQR file: {}".format(
        rec_pqr_path)
    if len(rec_residues) == 1:
        print(warnstr1.format(rec_pqr_path, rec_residues[0].name))
        return False
        
    assert lig_residues != 0, "Empty PQR file: {}".format(
        lig_pqr_path)
    if len(lig_residues) == 1:
        print(warnstr1.format(lig_pqr_path, lig_residues[0].name))
        return False
    
    return True

def check_for_one_bulk_anchor(model):
    """
    In order for the model to work properly, one or more bulk anchors is
    required if BD is enabled.
    """
    if not model.using_bd():
        return True
    
    num_bulk_anchors = 0
    for anchor in model.anchors:
        if anchor.bulkstate:
            num_bulk_anchors += 1
            
    if num_bulk_anchors >= 1:
        return True
    else:
        warnstr = """CHECK FAILURE: {} bulk anchors were found. There needs
        to be one or more bulk anchors per model if BD is enabled."""
        print(warnstr.format(num_bulk_anchors))
        return False

def check_mutual_neighbor_anchors(model):
    """
    Make sure that every anchor's neighbor recognizes the anchor itself
    as a neighbor
    """
    anchors_missing_neighbors = set()
    double_booked_neighbors = set()
    for anchor in model.anchors:
        for milestone in anchor.milestones:
            anchor2 = model.anchors[milestone.neighbor_anchor_index]
            # now check this anchor to ensure that it has me as a neighbor
            has_me_as_neighbor = False
            for milestone2 in anchor2.milestones:
                if anchor.index == milestone2.neighbor_anchor_index:
                    if has_me_as_neighbor:
                        double_booked_neighbors.add(anchor.index)
                    has_me_as_neighbor = True
            
            if not has_me_as_neighbor:
                anchors_missing_neighbors.add(anchor.index)
    
    success = True
    if len(anchors_missing_neighbors) > 0:
        warnstr = """CHECK FAILURE: Anchors {} were found to have mismatched
        neighbor anchors. This is likely a bug. Please contact the 
        developers."""
        print(warnstr.format(anchors_missing_neighbors))
        success = False
        
    if len(double_booked_neighbors) > 0:
        warnstr = """CHECK FAILURE: Anchors {} were found to have double-
        booked neighbor anchors. This is likely a bug. Please contact the 
        developers."""
        print(warnstr.format(double_booked_neighbors))
        success = False
        
    return success

def check_pre_simulation_all(model):
    """
    After the completion of the prepare stage, check inputs for some
    of the most common problems and mistakes a user is likely to 
    encounter. If a check fails, raise an Exception.
    
    Parameters:
    -----------
    model : Model()
        The SEEKR2 model object containing all calculation information.
        
    Returns:
    None
    """
    curdir = os.getcwd()
    check_passed_list = []
    #check_passed_list.append(check_pre_sim_bubbles(model))
    if not model.using_toy():
        # Skipping MD/BD salt conc. check because the best results seem to 
        # come from using no salt in BD.
        #check_passed_list.append(check_pre_sim_MD_and_BD_salt_concentration(model))
        check_passed_list.append(check_atom_selections_on_same_molecule(model))
        check_passed_list.append(check_atom_selections_MD_BD(model))
        check_passed_list.append(check_pqr_residues(model))
    else:
        # Checks for toy systems
        pass
    
    check_passed_list.append(check_systems_within_Voronoi_cells(model))
    check_passed_list.append(check_for_one_bulk_anchor(model))
    check_passed_list.append(check_mutual_neighbor_anchors(model))
    
    no_failures = True
    for check_passed in check_passed_list:
        if not check_passed:
            no_failures = False
            
    if no_failures:
        print("All pre-simulation checks passed.")
        return
    else:
        check_fail_str = "One or more fatal pre-simulation checks failed. It "\
        "is highly recommended that you address and correct each of these "\
        "problems. However, you can force SEEKR to skip these checks by using "\
        "the --skip_checks (-s) argument on prepare.py. If these checks were "\
        "performed using prepare.py, the model was saved, and though it is "\
        "not recommended, you may choose to proceed with the calculation "\
        "regardless of these failed checks."
        print(check_fail_str)
        raise Exception("The SEEKR2 calculation should not proceed due "\
                        "to failed pre-simulation checks.")
    os.chdir(curdir)
    return

def check_elber_umbrella_stage(model):
    """
    If an umbrella force is improperly constructed or because of a bug,
    a system may deviate too far from the milestone during an umbrella
    stage. Detect this situation to alert the user.
    """
    
    warnstr = """CHECK FAILURE: Elber umbrella stage trajectory for anchor {} 
    deviates significantly from the location of the central milestone.
    This could be caused by a bad umbrella force constant (too large or
    too small) umbrella_force_constant value.
    If may also happen if the model.xml file was improperly modified,
    or perhaps due to a bug."""
    if model.get_type() != "elber":
        # this test would be irrelevant
        return True
    if model.openmm_settings is not None:
        if model.openmm_settings.hydrogenMass is not None:
            return True
    for anchor in model.anchors:
        traj = load_structure_with_mdtraj(model, anchor, mode="elber_umbrella")
        if traj is None:
            continue
        for milestone in anchor.milestones:
            if milestone.alias_index == 2:
                cv = model.collective_variables[milestone.cv_index]
                result = cv.check_mdtraj_close_to_boundary(
                    traj, milestone.variables, verbose=True)
                if result == False:
                    print(warnstr.format(anchor.index))
                    return False
    return True

def check_xml_boundary_states(model):
    """
    SEEKR2 calculations generate OpenMM XML state files when 
    boundaries are encountered. Ensure that these boundary encounters
    are properly accounted and that they exist where expected.
    """
    warnstr = """CHECK FAILURE: Saved boundary states in anchor {} were not
    saved close to the expected milestone(s). File name: {}. 
    This could be caused if the model.xml file was improperly 
    modified, or perhaps due to a bug."""
    for anchor in model.anchors:
        states_dir = os.path.join(
            model.anchor_rootdir, anchor.directory, 
            anchor.production_directory, SAVE_STATE_DIRECTORY)
        for milestone in anchor.milestones:
            if model.get_type() == "elber" and milestone.alias_index == 2:
                continue
            state_file_glob = os.path.join(
                states_dir, "*_*_%d" % milestone.alias_index)
            state_files = glob.glob(state_file_glob)
            if len(state_files) == 0:
                continue
            for i, state_file in enumerate(state_files):
                traj = load_structure_with_mdtraj(
                    model, anchor, mode="state_xml", coords_filename=state_file)
                cv = model.collective_variables[milestone.cv_index]
                result = cv.check_mdtraj_close_to_boundary(
                    traj, milestone.variables, verbose=True)
                save_pdb_filename = state_file+".pdb"
                # Do not save these PDBs
                # This messes up the automatic state loader in the runners
                #traj.save_pdb(save_pdb_filename)
                if result == False:
                    save_pdb_filename = os.path.join(states_dir, "problem.pdb")
                    traj.save_pdb(save_pdb_filename)
                    print(warnstr.format(anchor.index, state_file))
                    print("Saved problematic structure to:", save_pdb_filename)
                    return False
                
                if i > MAX_STRUCTURES_TO_CHECK:
                    break
    
    return True

def check_mmvt_in_Voronoi_cell(model):
    """
    For some reason, a MMVT system may drift out of a Voronoi cell.
    This would most likely indicate a bug, so detect this situation.
    """
    
    if model.get_type() != "mmvt":
        # irrelevant test
        return True
    if model.openmm_settings is not None:
        if model.openmm_settings.hydrogenMass is not None:
            return True
    for anchor in model.anchors:
        traj = load_structure_with_mdtraj(model, anchor, mode="mmvt_traj")
        if traj is None:
            continue
        for milestone in anchor.milestones:
            cv = model.collective_variables[milestone.cv_index]
            curdir = os.getcwd()
            os.chdir(model.anchor_rootdir)
            result = cv.check_mdtraj_within_boundary(traj, milestone.variables, 
                                                     verbose=True)
            os.chdir(curdir)
            if result == False:
                warnstr = """CHECK FAILURE: The MMVT trajectory(ies)
    for anchor {} do not lie within the 
    anchor boundaries. This could be caused by
    incorrect modification of files or the model.xml
    file, or possibly a bug.""".format(anchor.index)
                print(warnstr)
                return False
    
    return True

def find_parmed_structure_com(structure, indices):
    """
    For a parmed structure, find the center of mass (COM) of a set
    of atoms defined by "indices".
    """
    
    total_mass = 0.0
    com_vector = np.zeros(3)
    for index in indices:
        atom = structure.atoms[index]
        total_mass += atom.mass
        com_vector += atom.mass * structure.coordinates[index]
    
    assert total_mass > 0.0
    return com_vector / total_mass

# TODO: remove? This check no longer used since BD milestones have changed
def check_bd_simulation_end_state(model):
    """
    The BD stage simulation program, Browndye2, can save encounter
    complexes when the system encounters a milestone. Check these
    encounter complex structures to ensure that they exist close to
    the expected milestone.
    """
    
    ATOL = 0.1
    if not model.using_bd():
        return True
    
    for bd_index, bd_milestone in enumerate(model.k_on_info.bd_milestones):
        outer_radius = bd_milestone.outer_milestone.variables["radius"]
        fhpd_directory = os.path.join(
            model.anchor_rootdir, bd_milestone.directory, 
            bd_milestone.fhpd_directory)
        directory_glob = os.path.join(fhpd_directory, "lig*")
        directories = glob.glob(directory_glob)
        for directory in directories:
            lig_pqr_filename = os.path.join(directory, "ligand.pqr")
            rec_pqr_filename = os.path.join(directory, "receptor.pqr")
            lig_struct = parmed.load_file(lig_pqr_filename)
            rec_struct = parmed.load_file(rec_pqr_filename)
            lig_com = find_parmed_structure_com(
                lig_struct, bd_milestone.ligand_indices)
            rec_com = find_parmed_structure_com(
                rec_struct, bd_milestone.receptor_indices)
            distance = np.linalg.norm(rec_com - lig_com) / 10.0 # A_to_nm
            if not np.isclose(distance, outer_radius, atol=ATOL):
                warnstr = """CHECK FAILURE: The BD milestone number {}
    has saved FHPD structures from the b-surface simulations
    that were significantly far away from the outermost 
    milestone of this binding site. Observed distance: {:.3f}. 
    Expected distance: {:.3f}.""".format(bd_milestone.index, distance, outer_radius)
                print(warnstr)
                return False
    return True

def check_output_files_for_stuck_anchors(model):
    """
    Systems started outside the milestone boundaries may bounce every
    time step. There are checks before the simulation to ensure that this 
    doesn't happen, but restarts and states have been observed to cause 
    this phenomenon unexpectedly. This check makes sure that no output files
    after the simulation are stuck behind a milestone.
    """
    MAX_SEQUENTIAL_BOUNCES = 1000
    anchors_stuck = []
    timestep = model.get_timestep() + 1e-6
    for anchor in model.anchors:
        output_file_glob = os.path.join(
            model.anchor_rootdir, anchor.directory, 
            anchor.production_directory, anchor.md_output_glob)
        
        output_file_list = glob.glob(output_file_glob)
        
        for output_file_name in output_file_list:
            if anchor.index in anchors_stuck:
                break
            sequential_bounces = 0
            last_time = 0.0
            with open(output_file_name, "r") as f:
                for line in f.readlines():
                    if line.startswith("#"):
                        continue
                    if line.startswith("CHECKPOINT"):
                        cur_time = float(line.strip().split(",")[1])
                    else:
                        cur_time = float(line.strip().split(",")[2])
                    if cur_time - last_time <= timestep:
                        sequential_bounces += 1
                    else:
                        sequential_bounces = 0
                        
                    if sequential_bounces > MAX_SEQUENTIAL_BOUNCES:
                        print("max sequential bounces exceeded for anchor "\
                              "{} at time {}".format(anchor.index, cur_time))
                        anchors_stuck.append(anchor.index)
                        break
                    last_time = cur_time
                    
    if len(anchors_stuck) > 0:
        warnstr = """CHECK FAILURE: Anchors {} had simulations that were
        found to be stuck behind a milestone, and never stopped bouncing.
        Any analysis performed on this model will be incorrect. This problem
        can possibly be caused by an incorrect starting state file."""
        print(warnstr.format(anchors_stuck))
        return False
    else:
        return True
            
def check_for_stuck_BD(model):
    """
    Examine BD output to ensure that too many BD trajectories weren't
    stuck.
    """
    if model.using_bd():
        bd_transition_counts = common_converge.get_bd_transition_counts(model)
        if "stuck" in bd_transition_counts:
            if bd_transition_counts["stuck"] > MAX_BD_STUCK:
                warnstr = """CHECK FAILURE: A large number of BD
                trajectories stuck. {} stuck BD trajectories found."""
                print(warnstr.format(bd_transition_counts["stuck"]))
                return False
            else:
                return True
    return True

def check_post_simulation_all(model, long_check=False):
    """
    After the completion of the run stage, check simulation files
    for some of the most common problems and mistakes a user is likely 
    to encounter. If a check fails, raise an Exception.
    
    Parameters:
    -----------
    model : Model()
        The SEEKR2 model object containing all calculation information.
    
    long_check : bool, Default False
        Whether to conduct a detailed check of post-simulation outputs.
        If set to True, the checks are likely to take a significantly
        longer amount of time, but is more likely to detect problems.
    
    Returns:
    None
    """
    curdir = os.getcwd()
    check_passed_list = []
    if model.using_toy():
        pass
    else:
        check_passed_list.append(check_elber_umbrella_stage(model))
        check_passed_list.append(check_mmvt_in_Voronoi_cell(model))
        # TODO: remove?
        #check_passed_list.append(check_bd_simulation_end_state(model))
        if long_check:
            check_passed_list.append(check_xml_boundary_states(model))
    
    check_passed_list.append(check_output_files_for_stuck_anchors(model))
    check_passed_list.append(check_for_stuck_BD(model))
    
    no_failures = True
    for check_passed in check_passed_list:
        if not check_passed:
            no_failures = False
            
    if no_failures:
        print("All post-simulation checks passed.")
        return
    else:
        check_fail_str = "One or more fatal post-simulation checks failed. It "\
        "is highly recommended that you address and correct each of these "\
        "problems. However, you can force SEEKR to skip these checks by using "\
        "the --skip_checks (-s) argument on analyze.py."
        print(check_fail_str)
        raise Exception("The SEEKR2 calculation can not proceed due "\
                        "to failed post-simulation checks.")
    os.chdir(curdir)
    return

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(
        description="Check SEEKR2 calculations for common problems both "\
        "before and after simulations.")
    
    argparser.add_argument(
        "model_file", metavar="MODEL_FILE", type=str, 
        help="name of model file for OpenMMVT calculation. This would be the "\
        "XML file generated in the prepare stage.")
    
    args = argparser.parse_args() # parse the args into a dictionary
    args = vars(args)
    xmlfile = args["model_file"]
    model = base.Model()
    model.deserialize(xmlfile)
    if model.anchor_rootdir == ".":
        model_dir = os.path.dirname(xmlfile)
        model.anchor_rootdir = os.path.abspath(model_dir)
    check_pre_simulation_all(model)
    check_post_simulation_all(model, long_check=True)
    