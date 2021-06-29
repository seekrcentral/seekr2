"""
runner_browndye.py

A command-line tool for running SEEKR calculations using the BrownDye
engine either locally or on a computing resource.
"""
import os
import subprocess
import re
import argparse
import glob
import shutil
import xml.etree.ElementTree as ET
from xml.dom import minidom
from collections import defaultdict

import parmed
import mdtraj

import seekr2.modules.common_base as base
import seekr2.modules.common_sim_browndye2 as sim_browndye2

REACTION_FILENAME = "rxns.xml"
DEFAULT_N_TRAJ_PER_OUT = 1000

def make_empty_pqrxml(directory, filename="empty.pqrxml"):
    """
    Create an empty pqrxml which is used in the Extraction phase to 
    generate new ligand encounter complex PQRs.
    
    Returns
    -------
    empty_pqrxml_path: str
        Generated path of empty pqrxml file.
    """
    
    empty_pqrxml_path = os.path.join(directory, filename)
    with open(empty_pqrxml_path, "w") as f:
        empty_pqrxml_string = "<roottag> \n</roottag>"
        f.write(empty_pqrxml_string)
    return empty_pqrxml_path

def cleanse_bd_outputs(directory, check_mode=True):
    """
    Check for simulation outputs in an existing anchor. Optionally
    delete them to make way for new outputs.
    """
    files_will_be_removed = False
    output_files_glob = os.path.join(directory, base.BROWNDYE_OUTPUT)
    output_files_list = glob.glob(output_files_glob)
    for output_file in output_files_list:
        files_will_be_removed = True
        if not check_mode:
            os.remove(output_file)
    
    simulation_file_glob = os.path.join(directory, "*_simulation.xml")
    for simulation_file in glob.glob(simulation_file_glob):
        files_will_be_removed = True
        if not check_mode:
            os.remove(simulation_file)
    
    traj_file_glob = os.path.join(directory, "traj*.xml")
    for traj_file in glob.glob(traj_file_glob):
        files_will_be_removed = True
        if not check_mode:
            os.remove(traj_file)
    return files_will_be_removed

def make_browndye_input_xml(model, rootdir, receptor_xml_filename, 
                            ligand_xml_filename, num_bd_steps, 
                            bd_directory=None, make_apbs_mode=True):
    """
    Creates Browndye input file with all input parameters set.
    If bd_directory is None, then b_surface is assumed.
    
    Parameters
    ----------
    model : Model()
        The SEEKR2 model.
        
    rootdir : str
        The model's root directory.
           
    receptor_xml_filename : str
        receptor PQRXML structure.
    
    ligand_xml_filename : str
        ligand PQRXML structure.
    
    num_bd_steps : int
        The number of BD steps to run
    
    bd_directory : Directory
        Browndye directory
    
    make_apbs_mode : bool
        Whether to use Browndye2's automatic APBS grid generator.
    
    Returns
    -------
    debye_length : float
        The Debye length of the system based on calculation by APBS.
        
    reaction_filename : str
        A path to the reaction file.
    """
    
    root = sim_browndye2.Root()
    if bd_directory is None:
        # B-surface
        bd_directory = os.path.join(
            rootdir, model.k_on_info.b_surface_directory)
        root.system.start_at_site = "False"
    else:
        # BD Milestone
        bd_directory = os.path.join(rootdir, bd_directory)
        root.system.start_at_site = "True"
        
    root.n_trajectories = num_bd_steps
    if root.n_trajectories_per_output > root.n_trajectories:
        root.n_trajectories_per_output = root.n_trajectories
    root.n_threads = model.browndye_settings.n_threads
    root.system.reaction_file = REACTION_FILENAME
    reaction_filename = root.system.reaction_file
    if make_apbs_mode:
        root.system.solvent.debye_length = -1.0
        input_xml_filename = os.path.join(bd_directory, "apbs_input.xml")
    else:
        root.system.solvent.debye_length = model.browndye_settings.debye_length
        assert root.system.solvent.debye_length > 0.0, "The Debye length must "\
            "be set if make_apbs_mode=False"
        input_xml_filename = os.path.join(bd_directory, "input.xml")
        
    root.system.solvent.kT = model.temperature / 298.0
    root.system.time_step_tolerances.minimum_core_dt = 0.2
    root.system.time_step_tolerances.minimum_core_reaction_dt = 0.05
    
    root.system.solvent.ions = []
    for model_ion in model.k_on_info.ions:
        bd_ion = sim_browndye2.Ion()
        bd_ion.radius = model_ion.radius
        bd_ion.charge = model_ion.charge
        bd_ion.conc = model_ion.conc
        root.system.solvent.ions.append(bd_ion)
        
    receptor_group = sim_browndye2.Group()
    receptor_group.name = sim_browndye2.BROWNDYE_RECEPTOR
    receptor_core = sim_browndye2.Core()
    receptor_core.name = sim_browndye2.BROWNDYE_RECEPTOR
    receptor_pqrxml_filename = os.path.basename(receptor_xml_filename)
    receptor_core.atoms = receptor_pqrxml_filename
    receptor_core.grid_spacing = model.browndye_settings.apbs_grid_spacing
    receptor_group.core_list.append(receptor_core)
    ligand_group = sim_browndye2.Group()
    ligand_group.name = sim_browndye2.BROWNDYE_LIGAND
    ligand_core = sim_browndye2.Core()
    ligand_core.name = sim_browndye2.BROWNDYE_LIGAND
    ligand_pqrxml_filename = os.path.basename(ligand_xml_filename)
    ligand_core.atoms = ligand_pqrxml_filename
    ligand_core.grid_spacing = model.browndye_settings.apbs_grid_spacing
    ligand_group.core_list.append(ligand_core)
    root.system.group_list.append(receptor_group)
    root.system.group_list.append(ligand_group)
    if make_apbs_mode:
        root.write(input_xml_filename)
        debye_length = sim_browndye2.make_and_run_apbs(
            root, input_xml_filename, 
            browndye2_bin=model.browndye_settings.browndye_bin_dir,
            new_input_xml_base=sim_browndye2.BROWNDYE_INPUT_FILENAME)
    else:
        root.write(input_xml_filename)
        debye_length = model.browndye_settings.debye_length
                
    return debye_length, reaction_filename
    
def make_browndye_reaction_xml(model, abs_reaction_path, bd_milestone=None):
    """
    Creates Browndye2 reaction file.
    """
    rxnroot = sim_browndye2.Reaction_root()
    if bd_milestone is None:
        rxnroot.first_state = "b_surface"
        bd_milestone_name = "b_surface"
    else:
        rxnroot.first_state = str(
            bd_milestone.outer_milestone.index)
        bd_milestone_name = bd_milestone.name
    
    ghost_indices_rec = model.browndye_settings.ghost_indices_rec
    ghost_indices_lig = model.browndye_settings.ghost_indices_lig
    
    for i, bd_milestone2 in enumerate(model.k_on_info.bd_milestones):
        ghost_index_rec = ghost_indices_rec[i] # comes from the model?
        ghost_index_lig = ghost_indices_lig[i]
        rxn = sim_browndye2.Reaction()
        pair = sim_browndye2.Pair()
        if bd_milestone_name == bd_milestone2.name:
            rxn.name = str(bd_milestone2.inner_milestone.index)
            rxn.state_after = str(bd_milestone2.inner_milestone.index)
            pair.distance = bd_milestone2.inner_milestone.variables['radius'] \
                * 10.0
        else:
            rxn.name = str(bd_milestone2.outer_milestone.index)
            rxn.state_after = str(bd_milestone2.outer_milestone.index)
            pair.distance = bd_milestone2.outer_milestone.variables['radius'] \
                * 10.0
        rxn.state_before = rxnroot.first_state
        rxn.molecule0_group = sim_browndye2.BROWNDYE_RECEPTOR
        rxn.molecule0_core = sim_browndye2.BROWNDYE_RECEPTOR
        rxn.molecule1_group = sim_browndye2.BROWNDYE_LIGAND
        rxn.molecule1_core = sim_browndye2.BROWNDYE_LIGAND
        rxn.n_needed = 1
        pair.atom1_index = ghost_index_rec
        pair.atom2_index = ghost_index_lig
        
        rxn.pair_list.append(pair)
        rxnroot.reaction_list.append(rxn)
    rxnroot.write(abs_reaction_path)
    return

def run_bd_top(browndye_bin_dir, bd_directory, restart=False,
               force_overwrite=False):
    """
    Prepare Browndye2 input files for simulation.
    """
    curdir = os.getcwd()
    print("moving to directory:", bd_directory)
    os.chdir(bd_directory)
    simulation_filename = sim_browndye2.BROWNDYE_RECEPTOR + "_" \
        + sim_browndye2.BROWNDYE_LIGAND + "_simulation.xml"
    if os.path.exists(simulation_filename) and not restart:
        if force_overwrite:
            print("force_overwrite set to True: existing files will be "\
                  "overwritten.")
            os.remove(simulation_filename)
        else:
            print("This anchor already has existing output files and the "\
                  "entered command would overwrite them. If you desire to "\
                  "overwrite the existing files, then use the "\
                  "--force_overwrite (-f) option, and all outputs will be "\
                  "deleted and replaced by a new run.")
            raise Exception("Cannot overwrite existing Browndye outputs.")
        
    bd_command = os.path.join(browndye_bin_dir, "bd_top")
    command = bd_command + " " + sim_browndye2.BROWNDYE_INPUT_FILENAME
    assert os.path.exists(sim_browndye2.BROWNDYE_INPUT_FILENAME), \
        "Necessary file doesn't exist: %s" % \
        sim_browndye2.BROWNDYE_INPUT_FILENAME
    print("running command:", command)
    os.system(command)
    assert os.path.exists(simulation_filename), "Problem occurred running "\
        "bd_top: simulation file %s was not generated." % simulation_filename
    os.chdir(curdir)
    return

def modify_variables(bd_milestone_directory, bd_output_glob, 
                     n_trajectories, n_threads=None, seed=None, 
                     output_file=None, restart=False, 
                     n_trajectories_per_output=None):
    """Modify several variables at runtime based on user inputs."""
    simulation_filename_base = sim_browndye2.BROWNDYE_RECEPTOR + "_" \
        + sim_browndye2.BROWNDYE_LIGAND + "_simulation.xml"
    simulation_filename = os.path.join(bd_milestone_directory, 
                                       simulation_filename_base)
    sim_file_old_lines = []
    sim_file_new_lines = []
    with open(simulation_filename, 'r') as f:
        for line in f.readlines():
            sim_file_old_lines.append(line)
            
    results_file_glob = os.path.join(bd_milestone_directory, bd_output_glob)
    results_file_list = glob.glob(results_file_glob)
    if restart:
        assert len(results_file_list) > 0, \
            "No simulation has yet been run: cannot use restart mode."
        restart_index = len(results_file_list) + 1
        if seed is None:
            seed = restart
        else:
            seed += restart
    else:
        restart_index = 1
        cleanse_bd_outputs(bd_milestone_directory, check_mode=False)
        
    trajectory_file = "traj%d_" % restart_index
    
    if output_file is None:
        output_file = "results%d.xml" % restart_index
    
    for line in sim_file_old_lines:
        new_line = line
        if n_trajectories is not None:
            new_line = re.sub(r"(?is)<n_trajectories>.+</n_trajectories>", 
                   "<n_trajectories> %d </n_trajectories>" % n_trajectories, 
                   new_line)
        
            if n_trajectories_per_output is None:
                n_trajectories_per_output = n_trajectories
            
            new_line = re.sub(r"(?is)<n_trajectories_per_output>.+"\
                              "</n_trajectories_per_output>", 
                   "<n_trajectories_per_output> %d "\
                   "</n_trajectories_per_output>" % n_trajectories_per_output, 
                   new_line)
        
        if n_threads is not None:
            new_line = re.sub(r"(?is)<n_threads>.+</n_threads>", 
                   "<n_threads> %d </n_threads>" % n_threads, new_line)
            
        if seed is not None:
            new_line = re.sub(r"(?is)<seed>.+</seed>", 
                   "<seed> %d </seed>" % seed, new_line)
            
        if output_file is not None:
            new_line = re.sub(r"(?is)<output>.+</output>", 
                   "<output> %s </output>" % output_file, new_line)
            
        new_line = re.sub(r"(?is)<trajectory_file>.+</trajectory_file>", 
            "<trajectory_file> %s </trajectory_file>" % trajectory_file, 
            new_line)
            
        sim_file_new_lines.append(new_line)
    
    with open(simulation_filename, 'w') as f:
        for line in sim_file_new_lines:
            f.write(line)
            
    return

def run_nam_simulation(browndye_bin_dir, bd_directory, bd_output_glob):
    """Run the Browndye2 simulation."""
    curdir = os.getcwd()
    print("moving to directory:", bd_directory)
    os.chdir(bd_directory)
    bd_command = os.path.join(browndye_bin_dir, "nam_simulation")
    simulation_filename = sim_browndye2.BROWNDYE_RECEPTOR + "_" \
        + sim_browndye2.BROWNDYE_LIGAND + "_simulation.xml"
    command = bd_command + " " + simulation_filename
    print("running command:", command)
    os.system(command)
    results_glob = glob.glob(bd_output_glob)
    assert len(results_glob) > 0, "Problem occurred running "\
        "nam_simulation: results file was not generated."
    os.chdir(curdir)
    return

def make_proc_file_last_frame(input_filename, output_filename, 
                              pqrxml_path_1, pqrxml_path_2):
    """
    Extract the last frame from a process_trajectories output XML
    and write a new XML containing only the last frame.
    """
    input_tree = ET.parse(input_filename)
    output_trajectory = ET.Element("trajectory")
    output_trajectory.text = "\n  "
    input_trajectory = input_tree.getroot()
    last_state = None
    for item in input_trajectory:
        if item.tag == "state":
            last_state = item
        elif item.tag == "atom_files":
            new_atom_files = ET.SubElement(output_trajectory, "atom_files")
            new_atom_files.text = " %s %s " % (pqrxml_path_1, pqrxml_path_2)
            new_atom_files.tail = "\n  "
        else:
            output_trajectory.append(item)
    assert last_state is not None
    output_trajectory.append(last_state)
    xmlstr = ET.tostring(output_trajectory).decode("utf-8")
    with open(output_filename, "w") as f:
        f.write(xmlstr)
        
    return

def extract_bd_surface(model, bd_milestone, max_b_surface_trajs_to_extract,
                       force_overwrite=False, restart=False):
    """
    Given the processed trajectories output XML files, extract the
    encounter complex PQR files from the b-surface simulation(s).
    """
    lig_pqr_filenames = []
    rec_pqr_filenames = []
    b_surface_dir = os.path.join(model.anchor_rootdir, 
                                 model.k_on_info.b_surface_directory)
    b_surface_ligand_pqr = model.browndye_settings.ligand_pqr_filename
    b_surface_ligand_pqrxml = os.path.splitext(
        b_surface_ligand_pqr)[0] + ".xml"
    b_surface_ligand_pqrxml_full_path = os.path.join(b_surface_dir, 
                                                     b_surface_ligand_pqrxml)
    assert os.path.exists(b_surface_ligand_pqrxml_full_path), "PQRXML file %s "\
        "not found for b-surface." % b_surface_ligand_pqrxml_full_path
    b_surface_receptor_pqr = model.browndye_settings.receptor_pqr_filename
    b_surface_receptor_pqrxml = os.path.splitext(
        b_surface_receptor_pqr)[0] + ".xml"
    b_surface_receptor_pqrxml_full_path = os.path.join(
        b_surface_dir, b_surface_receptor_pqrxml)
    assert os.path.exists(b_surface_receptor_pqrxml_full_path), "PQRXML file "\
        "%s not found for b-surface." % b_surface_receptor_pqrxml_full_path
    bd_milestone_directory = os.path.join(model.anchor_rootdir, 
                                          bd_milestone.directory)
    extract_directory = os.path.join(bd_milestone_directory, 
                                     bd_milestone.extracted_directory)
    assert os.path.exists(extract_directory)
    sitename = str(bd_milestone.outer_milestone.index)
    empty_pqrxml_path = make_empty_pqrxml(extract_directory)
    process_trajectories = os.path.join(
        model.browndye_settings.browndye_bin_dir, "process_trajectories")
    vtf_trajectory = os.path.join(
        model.browndye_settings.browndye_bin_dir, "vtf_trajectory")
    
    quitting = False
    counter = 0
    
    results_file_glob = os.path.join(b_surface_dir, base.BROWNDYE_OUTPUT)
    results_file_list = glob.glob(results_file_glob)
    assert len(results_file_list) > 0, \
        "No b-surface simulation has yet been run: cannot extract to FHPD."
    num_restarts = len(results_file_list)
    
    for i in range(1,num_restarts+1):
        for j in range(model.browndye_settings.n_threads):
            if quitting: break
            print("extracting trajectories from traj number:", i, "thread:", j)
            output_filename = os.path.join(extract_directory, 
                                          "rxn_output%d_%d.txt" % (i,j))
            if os.path.exists(output_filename):
                if force_overwrite:
                    print("force_overwrite set to True: existing files will "\
                          "be overwritten.")
                elif restart:
                    print("restarting extraction.")
                else:
                    print("This folder already has existing output files and "\
                          "the entered command would overwrite them. If you "\
                          "desire to overwrite the existing files, then use "\
                          "the --force_overwrite (-f) option, and all "\
                          "outputs will be deleted and replaced by a new run.")
                    raise Exception("Cannot overwrite existing Browndye "\
                                    "outputs.")
            traj_filename = os.path.join(b_surface_dir, "traj%d_%d.xml" % (i,j))
            trajindex_filename = os.path.join(
                b_surface_dir, "traj%d_%d.index.xml" % (i,j))
            #assert os.path.exists(traj_filename), "trajectory output file "\
            #    "%s not found for b-surface. Are you sure you ran b-surface "\
            #    "simulations?" % traj_filename
            #assert os.path.exists(trajindex_filename), "trajectory output "\
            #    "file %s not found for b-surface. Are you sure you ran "\
            #    "b-surface simulations?" % traj_filename
            if not os.path.exists(traj_filename) \
                    or not os.path.exists(trajindex_filename):
                continue
            
            command = "echo 'Browndye Trajectory number'; "\
                +process_trajectories+" -traj %s -index %s -srxn %s > %s" \
                % (traj_filename, trajindex_filename, sitename, 
                   output_filename)
            print("running command:", command)
            std_out = subprocess.check_output(command, stderr=subprocess.STDOUT,
                                               shell=True)
            assert os.path.exists(output_filename) and \
                    os.stat(output_filename).st_size > 0.0, "Problem running "\
                "process_trajectories: reaction list file not generated."
            number_list = []
            subtraj_list = []
            with open(output_filename, "r") as f:
                for line in f.readlines():
                    if re.search("<number>",line):
                        number_list.append(int(line.strip().split()[1]))
                    elif re.search("<subtrajectory>",line):
                        subtraj_list.append(int(line.strip().split()[1]))
            #print("number_list:", number_list)
            #print("subtraj_list:", subtraj_list)
            if len(number_list) == 0 or len(subtraj_list) == 0:
                continue
            # sort both lists concurrently so that structures are in order
            number_list, subtraj_list = zip(*sorted(zip(
                number_list, subtraj_list)))
            for k, rxn_number in enumerate(number_list):
                if counter > max_b_surface_trajs_to_extract:
                    quitting = True
                    break
                rxn_subtraj = subtraj_list[k]
                proc_traj_basename = os.path.join(
                    extract_directory, "proc_traj%d_%d_%d" % (i, j, k))
                xml_traj_filename = proc_traj_basename + ".xml"
                command = process_trajectories+" -traj %s -index %s -n %d "\
                    "-sn %d -nstride 1 > %s" % (
                        traj_filename, trajindex_filename, rxn_number, 
                        rxn_subtraj, xml_traj_filename)
                print("running command:", command)
                std_out = subprocess.check_output(
                    command, stderr=subprocess.STDOUT, shell=True)
                assert os.path.exists(xml_traj_filename) and \
                    os.stat(xml_traj_filename).st_size > 0.0, \
                    "Problem running process_trajectories: trajectory XML "\
                    "file not generated."
                
                last_frame_name = proc_traj_basename + "_last.xml"
                make_proc_file_last_frame(xml_traj_filename, last_frame_name,
                                          empty_pqrxml_path, 
                                          b_surface_ligand_pqrxml_full_path)
                
                # write the last frame as a pqr file
                pqr_filename = os.path.join(extract_directory, 
                                            "lig%d_%d_%d.pqr" % (i,j,k))
                lig_pqr_filenames.append(pqr_filename)
                command = vtf_trajectory+" -traj %s -pqr > %s"\
                    % (last_frame_name, pqr_filename)
                print("running command:", command)
                std_out = subprocess.check_output(
                    command, stderr=subprocess.STDOUT, shell=True)
                assert os.path.exists(pqr_filename) and \
                    os.stat(pqr_filename).st_size > 0.0, "Problem running "\
                    "vtf_trajectory: ligand%d_%d_%d PQR file not generated." \
                    % (i, j, k)
                    
                make_proc_file_last_frame(xml_traj_filename, last_frame_name,
                                      b_surface_receptor_pqrxml_full_path, 
                                      empty_pqrxml_path)
                
                pqr_rec_filename = os.path.join(
                    extract_directory, "receptor%d_%d_%d.pqr" % (i,j,k))
                rec_pqr_filenames.append(pqr_rec_filename)
                command = vtf_trajectory+" -traj %s -pqr > "\
                    "%s" % (last_frame_name, pqr_rec_filename)
                print("running command:", command)
                std_out = subprocess.check_output(
                    command, stderr=subprocess.STDOUT, shell=True)
                assert os.path.exists(pqr_filename) and \
                    os.stat(pqr_filename).st_size > 0.0, "Problem running "\
                    "vtf_trajectory: receptor%d_%d_%d PQR file not generated." \
                    % (i,j,k)
                
                os.remove(xml_traj_filename)
                counter += 1
    
    assert len(lig_pqr_filenames) > 0, "No trajectories found in b_surface "\
            "simulations. Consider using larger outermost milestone or "\
            "simulating more b-surface trajectories."
    assert len(lig_pqr_filenames) == len(rec_pqr_filenames)
    return lig_pqr_filenames, rec_pqr_filenames

def make_big_fhpd_trajectory(directory, lig_pqr_filenames, rec_pqr_filenames):
    """
    Combine all extracted PQR files from the b-surface simulation 
    encounter complexes, and make one big PDB trajectory representing
    the first hitting point distribution (FHPD).
    """
    fhpd_traj_filename = os.path.join(directory, "fhpd.pdb")
    temporary_pdb_filenames = []
    for i, (lig_pqr_filename, rec_pqr_filename) in enumerate(zip(
            lig_pqr_filenames, rec_pqr_filenames)):
        lig_frame = parmed.load_file(lig_pqr_filename)
        rec_frame = parmed.load_file(rec_pqr_filename)
        combined_frame = rec_frame + lig_frame
        temp_pdb_filename = os.path.join(directory, "fhpd_TEMP%d.pdb" % i)
        combined_frame.save(temp_pdb_filename, overwrite=True)
        temporary_pdb_filenames.append(temp_pdb_filename)
        
    traj = mdtraj.load(temporary_pdb_filenames)
    traj.save_pdb(fhpd_traj_filename)
    for temp_pdb_filename in temporary_pdb_filenames:
        os.remove(temp_pdb_filename)
        
    return
    

def make_fhpd_directories(model, bd_milestone, lig_pqr_filenames, 
                          rec_pqr_filenames, force_overwrite=False):
    """
    For each of the extracted first hitting point distribution (FHPD)
    PQR files, generate a directory to run its own Browndye2 
    simulation(s).
    """
    b_surface_dir = os.path.join(model.anchor_rootdir, 
                                 model.k_on_info.b_surface_directory)
    bd_milestone_directory = os.path.join(model.anchor_rootdir, 
                                          bd_milestone.directory)
    extract_directory = os.path.join(bd_milestone_directory, 
                                     bd_milestone.extracted_directory)
    fhpd_directory = os.path.join(bd_milestone_directory, 
                                  bd_milestone.fhpd_directory)
    make_big_fhpd_trajectory(fhpd_directory, lig_pqr_filenames, 
                             rec_pqr_filenames)
    sitename = str(bd_milestone.outer_milestone.index)
    
    if force_overwrite:
        existing_fhpd_directories_glob = os.path.join(fhpd_directory, "lig*")
        for fhpd_lig_directory in glob.glob(existing_fhpd_directories_glob):
            shutil.rmtree(fhpd_lig_directory)
    
    directories = []
    for i, (lig_orig_pqr_filename, rec_orig_pqr_filename) in enumerate(zip(
            lig_pqr_filenames, rec_pqr_filenames)):
        ligand_basename = os.path.splitext(
            os.path.basename(lig_orig_pqr_filename))[0]
        directory_name = os.path.join(fhpd_directory, ligand_basename)
        if not os.path.exists(directory_name):
            os.mkdir(directory_name)
        lig_pqr_filename = os.path.join(directory_name, "ligand.pqr")
        shutil.copyfile(lig_orig_pqr_filename, lig_pqr_filename)
        ligand_xml_filename = os.path.join(directory_name, "ligand.xml")
        rec_pqr_filename = os.path.join(directory_name, "receptor.pqr")
        shutil.copyfile(rec_orig_pqr_filename, rec_pqr_filename)
        receptor_xml_filename = os.path.join(
            directory_name, "receptor.xml")
        
        abs_reaction_path = os.path.join(directory_name, 
                                         REACTION_FILENAME)
        make_browndye_reaction_xml(model, abs_reaction_path, 
                                   bd_milestone=bd_milestone)
        
        sim_browndye2.make_pqrxml(lig_pqr_filename, 
                    browndye2_bin=model.browndye_settings.browndye_bin_dir,
                    output_xml_filename=ligand_xml_filename)
        sim_browndye2.make_pqrxml(rec_pqr_filename, 
                    browndye2_bin=model.browndye_settings.browndye_bin_dir,
                    output_xml_filename=receptor_xml_filename)
        
        if model.browndye_settings.recompute_ligand_electrostatics:
            debye_length, reaction_filename = make_browndye_input_xml(
                model, model.anchor_rootdir, receptor_xml_filename, 
                ligand_xml_filename, bd_milestone.num_trajectories, 
                bd_directory=directory_name, make_apbs_mode=True)
        else:
            raise Exception("recompute_ligand_electrostatics=False is not yet"\
                            "implemented.")
            debye_length, reaction_filename = make_browndye_input_xml(
                model, model.anchor_rootdir, receptor_xml_filename, 
                ligand_xml_filename, bd_milestone.num_trajectories, 
                bd_directory=directory_name, make_apbs_mode=False)
        
        directories.append(directory_name)
        
    return directories

def combine_fhpd_results(model, bd_milestone, fhpd_directories):
    """
    For each of the simulations performed within a BD milestone,
    combine their outputs and results into a single results file
    for analysis.
    """
    reaction_dict = defaultdict(float)
    number_escaped = 0
    number_stuck = 0
    number_total = 0
    number_total_check = 0
    results_filename_list = []
    for fhpd_directory in fhpd_directories:
        results_glob = os.path.join(fhpd_directory, 
                                    bd_milestone.bd_output_glob)
        results_filename_list += glob.glob(results_glob)
    
    assert len(results_filename_list) > 0, "No BD output files found."
    for results_filename in results_filename_list:
        tree = ET.parse(results_filename)
        root = tree.getroot()
        reactions_XML = root.find("reactions")
        number_total += int(reactions_XML.find("n_trajectories").text.strip())
        number_stuck += int(reactions_XML.find("stuck").text.strip())
        number_escaped += int(reactions_XML.find("escaped").text.strip())
        for completed_XML in reactions_XML.iter("completed"):
            name = completed_XML.find("name").text.strip()
            n = int(completed_XML.find("n").text.strip())
            reaction_dict[name] += n
            number_total_check += n
        
    assert number_total == number_total_check + number_stuck + number_escaped
    for completed_XML in reactions_XML.iter("completed"):
        reactions_XML.remove(completed_XML)
    
    reactions_XML.find("n_trajectories").text = str(number_total)
    reactions_XML.find("stuck").text = str(number_stuck)
    reactions_XML.find("escaped").text = str(number_escaped)
    for key in reaction_dict:
        completed_XML = ET.SubElement(reactions_XML, "completed")
        completed_XML.text = "\n      "
        completed_XML.tail = "\n  "
        name_XML = ET.SubElement(completed_XML, "name")
        name_XML.text = key
        name_XML.tail = "\n      "
        n_XML = ET.SubElement(completed_XML, "n")
        n_XML.text = str(int(reaction_dict[key]))
        n_XML.tail = "\n    "
    
    xmlstr = ET.tostring(root).decode("utf-8")
    bd_milestone_directory = os.path.join(model.anchor_rootdir, 
                                          bd_milestone.directory)
    dest_filename = os.path.join(bd_milestone_directory, "results.xml")
    with open(dest_filename, 'w') as f:
        f.write(xmlstr)
        
    return

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument(
        "bd_milestone", metavar="BD_MILESTONE", type=str, 
        help="Which milestone to run BrownDye for. Arguments may be the "\
        "string 'b_surface' or a numerical index of a given BD_milestone.")
    argparser.add_argument(
        "input_file", metavar="INPUT_FILE", type=str, 
        help="name of input file for OpenMMVT calculation. This would be the "\
        "XML file generated in the prepare stage.")
    argparser.add_argument(
        "-r", "--restart", dest="restart", default=False,
        help="Restart simulations to obtain more results.", action="store_true")
    argparser.add_argument("-n", "--n_trajectories", dest="n_trajectories", 
        default=None, help="Enter a different number of trajectories to run "\
        "the simulation if a different number of trajectories are desired "\
        "than what is in the input.xml file.", type=int)
    argparser.add_argument("-t", "--n_threads", dest="n_threads", 
        default=None, help="Enter a different number of threads to run "\
        "the simulation if a different number of threads are desired "\
        "than what is in the input.xml file.", type=int)
    argparser.add_argument("-o", "--output_file", dest="output_file",
        default=None, help="Enter a new output file name different from the "\
        "output file specified in the input.xml file.", type=str)
    argparser.add_argument("-s", "--seed", dest="seed", default=None, 
        help="Enter a new random number seed if different from the "\
        "seed specified in the input.xml file.", type=int)
    argparser.add_argument("-d", "--directory", dest="directory",
        default=None, help="Provide the location of the directory which "\
        "contains the anchors. If this argument is omitted, then the "\
        "directory within the anchor_rootdir setting of the INPUT_FILE will "\
        "be used.", type=str)
    argparser.add_argument("-m", "--max_b_surface_trajs", 
        dest="max_b_surface_trajs", help="In extraction mode, enter the "\
        "maximum number of successful encounter states to extract from the "\
        "b-surface simulations to start the bd_milestone simulations from.", 
        type=int, default=1000)
    argparser.add_argument("-f", "--force_overwrite", dest="force_overwrite",
        default=False, help="Toggle whether to overwrite existing files "\
        "within milestone. If this option is enabled, "\
        "then the anchor's production directory will be emptied of all "\
        "output, trajectory, and results files for the new simulation.", 
        action="store_true")
        
    args = argparser.parse_args() # parse the args into a dictionary
    args = vars(args)
    bd_milestone_index = args["bd_milestone"]
    input_file = args["input_file"]
    restart = args["restart"]
    n_trajectories = args["n_trajectories"]
    n_threads = args["n_threads"]
    output_file = args["output_file"]
    seed = args["seed"]
    directory = args["directory"]
    max_b_surface_trajs_to_extract = args["max_b_surface_trajs"]
    force_overwrite = args["force_overwrite"]
    
    b_surface = False
    if bd_milestone_index.lower() == "b_surface":
        bd_milestone_index = bd_milestone_index.lower()
        b_surface = True
    else:
        try:
            bd_milestone_index = int(bd_milestone_index)
        except ValueError:
            print("The bd_milestone argument must either be the string "\
                  "'b_surface' or an integer.")
            exit()
    
    assert os.path.exists(input_file), "A nonexistent input file was provided."
    model = base.Model()
    model.deserialize(input_file)
    
    if directory is not None:
        model.anchor_rootdir = os.path.abspath(directory)
    elif model.anchor_rootdir == ".":
        model_dir = os.path.dirname(input_file)
        model.anchor_rootdir = os.path.abspath(model_dir)
        
    if n_threads is not None:
        model.browndye_settings.n_threads = n_threads
        
    assert os.path.exists(model.anchor_rootdir), "An incorrect directory "\
        "was provided."
    
    if b_surface:
        bd_milestone_directory = os.path.join(
            model.anchor_rootdir, model.k_on_info.b_surface_directory)
        
    else:
        assert bd_milestone_index >= 0, "only positive indices allowed."
        try:
            bd_milestone = model.k_on_info.bd_milestones[bd_milestone_index]
        except IndexError:
            print("Invalid bd_milestone index provided.")
            exit()
        bd_milestone_directory = os.path.join(
            model.anchor_rootdir, bd_milestone.directory)
    
    if not b_surface:
        assert max_b_surface_trajs_to_extract is not None, "The argument "\
            "'--max_b_surface_trajs_to_extract' (-m) must be set for the "\
            "specified calculation."
        lig_pqr_filenames, rec_pqr_filenames = extract_bd_surface(
            model, bd_milestone, max_b_surface_trajs_to_extract, 
            force_overwrite)
        bd_directory_list = make_fhpd_directories(
            model, bd_milestone, lig_pqr_filenames, rec_pqr_filenames,
            force_overwrite)
        
    else:
        bd_directory_list = [bd_milestone_directory]
        
    for bd_directory in bd_directory_list:
        run_bd_top(model.browndye_settings.browndye_bin_dir, 
                   bd_directory, restart, force_overwrite)
        n_trajectories_per_output = DEFAULT_N_TRAJ_PER_OUT
        assert n_trajectories is not None, "The argument '--n_trajectories' "\
            "(-n) must be set for the specified calculation."
        if n_trajectories_per_output > n_trajectories:
            n_trajectories_per_output = n_trajectories
        modify_variables(bd_directory, model.k_on_info.bd_output_glob, 
                         n_trajectories, n_threads, seed, output_file, restart,
                         n_trajectories_per_output=n_trajectories_per_output)
        run_nam_simulation(model.browndye_settings.browndye_bin_dir, 
                           bd_directory, 
                           model.k_on_info.bd_output_glob)
    
    if not b_surface:
        combine_fhpd_results(model, bd_milestone, bd_directory_list)