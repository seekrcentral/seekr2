"""
make_test_model.py

makes a Model object for testing purposes in sim_openmm.py,
runner_openmm.py, and others.
"""

import os
from shutil import copyfile

import seekr2.modules.common_base as base
import seekr2.modules.mmvt_base as mmvt_base
import seekr2.modules.elber_base as elber_base
import seekr2.modules.filetree as filetree

def make_browndye_params_b_surface(model, b_surface_dir):
    this_dir = os.path.dirname(os.path.realpath(__file__))
    receptor_src = os.path.join(
        this_dir, "../data/hostguest_files/hostguest_receptor.pqr")
    ligand_src = os.path.join(
        this_dir, "../data/hostguest_files/hostguest_ligand.pqr")
    
    receptor_dest = os.path.join(b_surface_dir, "hostguest_receptor.pqr")
    ligand_dest = os.path.join(b_surface_dir, "hostguest_ligand.pqr")
    
    copyfile(receptor_src, receptor_dest)
    copyfile(ligand_src, ligand_dest)
    
    model.browndye_settings.receptor_pqr_filename = "hostguest_receptor.pqr"
    model.browndye_settings.ligand_pqr_filename = "hostguest_ligand.pqr"
    return

def make_amber_params(anchor, building_dir, engine="openmm"):
    anchor.amber_params = base.Amber_params()
    
    this_dir = os.path.dirname(os.path.realpath(__file__))
    
    if engine == "openmm":
        prmtop_src = os.path.join(
            this_dir, "../data/hostguest_files/hostguest.parm7")
    elif engine == "namd":
        prmtop_src = os.path.join(
            this_dir, "../data/hostguest_files/hostguest_for_NAMD.parm7")
    else:
        raise Exception("Engine not implemented: %s" % engine)
    inpcrd_src = os.path.join(
        this_dir, "../data/hostguest_files/hostguest.rst7")
    pdb_coord_src = os.path.join(
        this_dir, "../data/hostguest_files/hostguest_at1.9.pdb")
    
    prmtop_dest = os.path.join(building_dir, "hostguest.parm7")
    inpcrd_dest = os.path.join(building_dir, "hostguest.rst7")
    pdb_coord_dest = os.path.join(building_dir, "hostguest_at1.9.pdb")
    
    copyfile(prmtop_src, prmtop_dest)
    copyfile(inpcrd_src, inpcrd_dest)
    copyfile(pdb_coord_src, pdb_coord_dest)
    
    anchor.amber_params.prmtop_filename = "hostguest.parm7"
    anchor.amber_params.inpcrd_filename = "hostguest.rst7"
    anchor.amber_params.box_vectors = None
    anchor.amber_params.pdb_coordinates_filename = "hostguest_at1.9.pdb"
    return

def make_forcefield_params(anchor, building_dir):
    anchor.forcefield_params = base.Forcefield_params()
    
    this_dir = os.path.dirname(os.path.realpath(__file__))
    
    xml_src = os.path.join(this_dir, "../data/hostguest_files/hostguest.xml")
    pdb_coord_src = os.path.join(
        this_dir, "../data/hostguest_files/hostguest_for_xml.pdb")
    
    xml_dest = os.path.join(building_dir, "hostguest.xml")
    pdb_coord_dest = os.path.join(building_dir, "hostguest_for_xml.pdb")
    
    copyfile(xml_src, xml_dest)
    copyfile(pdb_coord_src, pdb_coord_dest)
    
    anchor.forcefield_params.built_in_forcefield_filenames = \
        ["amber14/tip3pfb.xml"]
    anchor.forcefield_params.custom_forcefield_filenames = ["hostguest.xml"]
    anchor.forcefield_params.pdb_filename = "hostguest_for_xml.pdb"
    anchor.forcefield_params.box_vectors = None
    return

def make_test_model(tmp_dir, num_anchors=3, milestone_type="spherical", 
                    make_dir=True, mode="amber", engine="openmm", no_BD=False,
                    calculation_type="mmvt"):
    
    if milestone_type == "spherical":
        #group1 = [2478, 2489, 2499, 2535, 2718, 2745, 2769, 2787, 2794, 2867, 
        #          2926]
        group1 = list(range(147))
        #group2 = [3221, 3222, 3223, 3224, 3225, 3226, 3227, 3228, 3229]
        group2 = list(range(147, 162))
        groups = [group1, group2]
        if calculation_type == "mmvt":
            cv1 = mmvt_base.MMVT_spherical_CV(index=0, groups=groups)
        elif calculation_type == "elber":
            cv1 = elber_base.Elber_spherical_CV(index=0, groups=groups)
        else:
            raise Exception(
                "Calculation type not available: {1}".format(calculation_type))
        cvs = [cv1]
        
    else:
        raise Exception("milestone type not implemented")
    
    mymodel = base.Model()
    mymodel.temperature = 277.8
    if calculation_type == "mmvt":
        mymodel.calculation_type="MMVT"
        mymodel.calculation_settings = mmvt_base.MMVT_settings()
    elif calculation_type == "elber":
        mymodel.calculation_type="Elber"
        mymodel.calculation_settings = elber_base.Elber_settings()
    
    mymodel.anchor_rootdir = tmp_dir
    mymodel.num_anchors = num_anchors
    mymodel.num_milestones = num_anchors-1
    
    if engine == "openmm":
        mymodel.openmm_settings = base.Openmm_settings()
        mymodel.openmm_settings.energy_reporter_frequency = 50
        mymodel.openmm_settings.trajectory_reporter_frequency = 50
        mymodel.openmm_settings.restart_checkpoint_frequency = 50
        mymodel.openmm_settings.total_simulation_length = 50
        mymodel.openmm_settings.cuda_platform_settings = None
        mymodel.openmm_settings.reference_platform = True
    
    if engine == "namd":
        mymodel.namd_settings = base.Namd_settings()
        mymodel.namd_settings.energy_reporter_frequency = 50
        mymodel.namd_settings.trajectory_reporter_frequency = 50
        mymodel.namd_settings.restart_checkpoint_frequency = 50
        mymodel.namd_settings.total_simulation_length = 50
    
    
    if not no_BD:
        my_k_on_info = base.K_on_info()
        my_k_on_info.source_milestones = [2]
        mymodel.k_on_info = my_k_on_info
    
    mymodel.collective_variables = cvs
    
    '''
    for index in range(num_anchors):
        milestone_list = []
        if index < num_anchors-1:
            milestone1 = base.Milestone()
            milestone1.index = index
            milestone1.neighbor_index = index+1
            milestone1.alias_id = 2
            milestone1.cv_index = 0
            milestone1.variables = {"radius": 1.3}
            milestone_list.append(milestone1)
            bulkstate = False
        else:
            bulkstate = True
        
        if index > 0:
            milestone2 = base.Milestone()
            milestone2.index = index-1
            milestone2.neighbor_index = index-1
            milestone2.alias_id = 1
            milestone2.cv_index = 0
            milestone2.variables = {"radius": 1.1}
            milestone_list.append(milestone2)
            endstate = False
        else:
            endstate = True
            
        if bulkstate:
            endstate = True
    '''
    for i in range(num_anchors):
        index = i
        milestone_list = []
        num_milestones = 0
        
        if i > 0:
            milestone1 = base.Milestone()
            milestone1.index = i-1
            milestone1.neighbor_index = i-1
            milestone1.alias_index = num_milestones+1
            milestone1.cv_index = 0
            milestone1.variables = {"k": -1.0, "radius": 0.05 + 0.1*i}
            milestone_list.append(milestone1)
            num_milestones += 1
            endstate = False
        else:
            endstate = True
        
        if i < num_anchors - 1:
            milestone2 = base.Milestone()
            milestone2.index = i
            milestone2.neighbor_index = i+1
            milestone2.alias_index = num_milestones + 1
            milestone2.cv_index = 0
            milestone2.variables = {"k": 1.0, "radius": 0.05 + 0.1*(i+1)}
            milestone_list.append(milestone2)
            bulkstate = False
        else:
            bulkstate = True
        
        name = "anchor_%d" % index
        directory = name
        if make_dir:
            full_anchor_dir = os.path.join(tmp_dir, directory)
            os.mkdir(full_anchor_dir)
            prod_dir = os.path.join(full_anchor_dir, "prod")
            os.mkdir(prod_dir)
            building_dir = os.path.join(full_anchor_dir, "building")
            os.mkdir(building_dir)
        md_glob = "mmvt*.out"
        md = True
        bd = False
        
        #anchor = base.Anchor(name, index, directory, md_glob, md, bd, 
        #         endstate, bulkstate, len(milestone_list), milestone_list)
        
        if calculation_type == "mmvt":
            #cv1 = mmvt_base.MMVT_spherical_CV(index=0, groups=groups)
            anchor = mmvt_base.MMVT_anchor()
        elif calculation_type == "elber":
            anchor = elber_base.Elber_anchor()
        else:
            raise Exception(
                "Calculation type not available: {1}".format(calculation_type))
        
        anchor.name = name
        anchor.index = index
        anchor.directory = directory
        anchor.md_mmvt_output_glob = md_glob
        anchor.md = md
        anchor.bd = bd
        anchor.endstate = endstate
        anchor.bulkstate = bulkstate
        #anchor.num_milestones = len(milestone_list)
        anchor.milestones = milestone_list
        
        if mode == "amber":
            make_amber_params(anchor, building_dir, engine=engine)
        elif mode == "forcefield":
            make_forcefield_params(anchor, building_dir)
        else:
            raise Exception("mode not implemented:"+mode)
        mymodel.anchors.append(anchor)
    
    if not no_BD:
        mymodel.browndye_settings = base.Browndye_settings()
        mymodel.browndye_settings.apbs_grid_spacing = 0.5
        mymodel.browndye_settings.n_threads = 1
        mymodel.browndye_settings.recompute_ligand_electrostatics = 1
        mymodel.browndye_settings.debye_length = 7.5612
        mymodel.browndye_settings.ghost_indices_rec = [147]
        mymodel.browndye_settings.ghost_indices_lig = [15]
        mymodel.k_on_info.b_surface_directory = "b_surface"
        mymodel.k_on_info.b_surface_num_steps = 100
        ion1 = base.Ion()
        ion1.radius = 1.2
        ion1.charge = -1.0
        ion1.conc = 0.15
        ion2 = base.Ion()
        ion2.radius = 0.9
        ion2.charge = 1.0
        ion2.conc = 0.15
        mymodel.k_on_info.ions = [ion1, ion2]
        if make_dir:
            b_surface_abs_path = os.path.join(tmp_dir, "b_surface")
            os.mkdir(b_surface_abs_path)
            
        make_browndye_params_b_surface(mymodel, b_surface_abs_path)
        
        mymodel.k_on_info.bd_milestones = []
        bd_milestone = base.BD_milestone()
        bd_milestone.index = 0
        bd_milestone.name = "bd_milestone_0"
        bd_milestone.directory = bd_milestone.name
        bd_milestone.outer_milestone = milestone2 
        #bd_milestone.outer_milestone = mymodel.anchors[-1].milestones[-1]
        assert "radius" in bd_milestone.outer_milestone.variables, \
            "A BD outer milestone must be spherical."
        
        
        bd_milestone.inner_milestone = milestone1
        
        bd_milestone.num_steps = 1
        bd_milestone.receptor_indices = list(range(147))
        bd_milestone.ligand_indices = list(range(15))
        
        mymodel.k_on_info.bd_milestones.append(bd_milestone)
        
        for bd_milestone in mymodel.k_on_info.bd_milestones:
            bd_milestone_dict = {}
            bd_milestone_dict[bd_milestone.extracted_directory] = {}
            bd_milestone_dict[bd_milestone.fhpd_directory] = {}
            
            bd_milestone_filetree = filetree.Filetree(
                {bd_milestone.directory:bd_milestone_dict})
            bd_milestone_filetree.make_tree(tmp_dir)
        
    return mymodel