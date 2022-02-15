"""
common_sim_openmm.py

Create the Sim_openmm object which contains the objects needed for an
openmm simulation based on the settings provided by a user. These
objects and routines are common to both Elber and MMVT milestoning.
"""

import os

import numpy as np
try:
    import openmm
except ImportError:
    import simtk.openmm as openmm
    
try:
    import openmm.app as openmm_app
except ImportError:
    import simtk.openmm.app as openmm_app
    
import parmed

class Common_sim_openmm:
    """
    Contain superclass information necessary to run a SEEKR2 
    calculation in OpenMM.
    
    platform : openmm.Platform()
        The platform to use for the simulation. The CUDA platform
        is recommended, but Reference is also available.
        
    properties : dict
        A dictionary of properties passed to the Simulation() object
        that apply to the platform.
        
    timestep : float
        The timestep (in picoseconds) which was passed to the 
        integrator, but is also useful to know in other places within 
        the program.
        
    output_filename : str
        The file that the plugin will write SEEKR output to, and which
        will be parsed by the analysis software.
        
    try_to_load_state : bool
        Whether to attempt to load a state from an adjacent anchor.
    
    """
    def __init__(self):
        self.platform = None
        self.properties = None
        self.timestep = None
        self.output_filename = ""
        self.try_to_load_state = False
        
    def add_force_objects(self):
        raise Exception("This base class cannot be used for creating a "\
                        "Sim_openmm object's force objects.")


def fill_generic_parameters(sim_openmm, model, anchor, output_filename):
    """
    Enter parameters to the sim_openmm object that are common to
    both MMVT and Elber calculations.
    """
    sim_openmm.output_filename = output_filename
    return

def write_toy_pdb_file(topology, positions, out_file_name):
    """
    
    """
    out_file = open(out_file_name, "w")
    openmm_app.PDBFile.writeHeader(topology, out_file)
    openmm_app.PDBFile.writeModel(topology, positions, out_file)
    openmm_app.PDBFile.writeFooter(topology, out_file)
    return

def make_toy_system_object(model):
    """
    
    """
    system = openmm.System()
    topology = openmm_app.Topology()
    num_particles = model.toy_settings.num_particles
    assert model.toy_settings.num_particles == len(model.toy_settings.masses)
    new_expression = "k*("+model.toy_settings.potential_energy_expression+")"
    force = openmm.CustomCentroidBondForce(num_particles, new_expression)
    force.addGlobalParameter("k", 1.0*openmm.unit.kilocalories_per_mole)
    groups_list = []
    for i in range(num_particles):
        mass = model.toy_settings.masses[i] * openmm.unit.amu
        system.addParticle(mass)
        atom_name = "X{}".format(i)
        if not atom_name in openmm_app.element.Element._elements_by_symbol:
            element = openmm_app.element.Element(0, atom_name, atom_name, mass)
        else:
            element = openmm_app.element.Element._elements_by_symbol[atom_name]
        chain = topology.addChain()
        residue = topology.addResidue("UNK", chain)
        topology.addAtom(atom_name, element, residue)
        groups_list.append(force.addGroup([i]))
    
    force.addBond(groups_list, [])
    system.addForce(force)
    return system, topology

def create_openmm_system(sim_openmm, model, anchor, frame=0):
    """
    Create an openmm System() object.
    """
    building_directory = os.path.join(
        model.anchor_rootdir, anchor.directory, anchor.building_directory)
    box_vectors = None
    
    positions_obj = None
    if anchor.__class__.__name__ == "MMVT_toy_anchor":
        positions = np.array(anchor.starting_positions) * openmm.unit.nanometers
    else:
        if anchor.amber_params is not None:
            prmtop_filename = os.path.join(
                building_directory, anchor.amber_params.prmtop_filename)
            prmtop = openmm_app.AmberPrmtopFile(prmtop_filename)
            #assert anchor.amber_params.pdb_coordinates_filename is not None
            if anchor.amber_params.pdb_coordinates_filename is None \
                    or anchor.amber_params.pdb_coordinates_filename == "":
                positions_obj = None
                positions = None
                sim_openmm.try_to_load_state = True
            else:
                pdb_coordinates_filename = os.path.join(
                    building_directory, 
                    anchor.amber_params.pdb_coordinates_filename)
                positions_obj = openmm_app.PDBFile(pdb_coordinates_filename)
                
            #assert anchor.amber_params.box_vectors is not None
            box_vectors = anchor.amber_params.box_vectors
            topology = prmtop.topology
        
        elif anchor.forcefield_params is not None:
            forcefield_filenames = []
            for forcefield_filename in \
                    anchor.forcefield_params.built_in_forcefield_filenames:
                forcefield_filenames.append(forcefield_filename)
            for forcefield_filename in \
                    anchor.forcefield_params.custom_forcefield_filenames:
                forcefield_filenames.append(os.path.join(
                    building_directory, forcefield_filename))
            pdb_filename = os.path.join(building_directory, 
                                   anchor.forcefield_params.pdb_filename)
            pdb = openmm_app.PDBFile(pdb_filename)
            forcefield = openmm_app.ForceField(
                *forcefield_filenames)
            box_vectors = anchor.forcefield_params.box_vectors
            
            topology = pdb.topology
            positions_obj = pdb
        
        elif anchor.charmm_params is not None:
            raise Exception("Charmm systems not yet implemented")
        
        else:
            raise Exception("No Amber or Charmm input settings detected.")
    
    #assert box_vectors is not None, "No source of box vectors provided."
    nonbonded_method = model.openmm_settings.nonbonded_method.lower()
    if nonbonded_method == "pme":
        nonbondedMethod = openmm_app.PME
        
    elif nonbonded_method == "nocutoff":
        nonbondedMethod = openmm_app.NoCutoff
        
    elif nonbonded_method == "cutoffnonperiodic":
        nonbondedMethod = openmm_app.CutoffNonPeriodic
        
    elif nonbonded_method == "cutoffperiodic":
        nonbondedMethod = openmm_app.CutoffPeriodic
        
    elif nonbonded_method == "ewald":
        nonbondedMethod = openmm_app.Ewald
    
    else:
        raise Exception("nonbonded method not found: %s", 
                        model.openmm_settings.nonbonded_method)
    
    if model.openmm_settings.constraints is None:
        constraints_str = None
    else:
        constraints_str = model.openmm_settings.constraints
        
    if constraints_str is None:
        constraints = None
    
    elif constraints_str.lower() == "none":
        constraints = None
    
    elif constraints_str.lower() == "hbonds":
        constraints = openmm_app.HBonds
        
    elif constraints_str.lower() == "allbonds":
        constraints = openmm_app.AllBonds
        
    elif constraints_str.lower() == "hangles":
        constraints = openmm_app.HAngles
        
    else:
        raise Exception("constraints not found: %s", 
                        model.openmm_settings.constraints)
    
    if model.openmm_settings.hydrogenMass is not None:
        hydrogenMass = model.openmm_settings.hydrogenMass*openmm.unit.amu
    else:
        hydrogenMass = None
        
    rigidWater = model.openmm_settings.rigidWater
    
    if anchor.__class__.__name__ == "MMVT_toy_anchor":
        system, topology = make_toy_system_object(model)
        out_file_name = os.path.join(model.anchor_rootdir, anchor.directory, 
                                     anchor.building_directory, "toy.pdb")
        write_toy_pdb_file(topology, positions, out_file_name)
    else:
        if anchor.amber_params is not None:
            system = prmtop.createSystem(
                nonbondedMethod=nonbondedMethod, 
                nonbondedCutoff=model.openmm_settings.nonbonded_cutoff, 
                constraints=constraints, hydrogenMass=hydrogenMass, 
                rigidWater=rigidWater)
        
        elif anchor.forcefield_params is not None:
            system = forcefield.createSystem(
                pdb.topology, nonbondedMethod=nonbondedMethod, 
                nonbondedCutoff=model.openmm_settings.nonbonded_cutoff, 
                constraints=constraints, hydrogenMass=hydrogenMass, 
                rigidWater=rigidWater)
        
        elif anchor.charmm_params is not None:
            raise Exception("Charmm input settings not yet implemented")
            
        else:
            print("Settings for Amber or Charmm simulations not found")
    
    if positions_obj is not None:
        positions = positions_obj.getPositions(frame=frame)
        assert frame >= 0, "Cannot have negative frame index"
        assert frame < positions_obj.getNumFrames(), \
            "Frame index {} out of range.".format(frame)
        
    return system, topology, positions, box_vectors

def add_barostat(sim_openmm, model):
    """
    Add a barostat to the sim_openmm object to maintain constant
    simulation pressure.
    """
    if model.openmm_settings.barostat:
        barostat = openmm.MonteCarloBarostat(
            model.openmm_settings.barostat.target_pressure, 
            model.openmm_settings.barostat.target_temperature,
            model.openmm_settings.barostat.frequency)
        sim_openmm.system.addForce(barostat)
    return

def add_platform(sim_openmm, model):
    """
    Specify which platform to use for this OpenMM simulation, whether
    'reference', 'CUDA', or otherwise.
    """
    if model.openmm_settings.reference_platform:
        sim_openmm.platform = \
            openmm.Platform.getPlatformByName('Reference')
        sim_openmm.properties = {}
    elif model.openmm_settings.cuda_platform_settings is not None:
        sim_openmm.platform = openmm.Platform.getPlatformByName('CUDA')
        sim_openmm.properties = model.openmm_settings.\
            cuda_platform_settings.make_properties_dict()
    else:
        raise Exception("Platform settings not found for either CUDA "\
                        "or Reference")
    return