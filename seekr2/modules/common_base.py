"""
common_base.py

Contain base classes, constants, and other objects used for both 
MMVT and Elber milestoning.
"""

import os
import re
import math

import numpy as np
import parmed
from parmed import unit

#import seekr2.libraries.serializer.serializer as serializer
from abserdes import Serializer

# A glob for BrownDye output files
BROWNDYE_OUTPUT = "results*.xml"

def strBool(bool_str):
    """
    Take the string "true" or "false" of any case and return a 
    boolean object.
    """
    
    if bool_str.lower() == "true":
        return True
    elif bool_str.lower() == "false":
        return False
    else:
        raise Exception(
            "argument for strBool must be string either 'True' or 'False'.")

def order_files_numerically(file_list, func=float, use_basename=False):
    """
    If there is a list of files, order them numerically, not
    alphabetically and return the sorted list of files. Note that
    only the base name is sorted, not any part of the directory.
    
    Parameters
    ----------
    file_list : list
        A list of strings that contain one or more numerical values. 
        The strings will be sorted by only the numerical values within
        them.
        
    Returns
    -------
    sorted_file_list : list
        A new list of strings sorted numerically
    """
    
    sorted_file_list = []
    numerical_dict = {}
    for i, file_name in enumerate(file_list):
        file_name2 = os.path.splitext(file_name)[0]
        if use_basename:
            file_name2 = os.path.basename(file_name2)
            
        numbers = re.findall(r"(-?\d+(?:\.\d+)?)", file_name2)
        numbers = tuple([func(j) for j in numbers])
        numerical_dict[numbers] = i
        
    numerical_list = sorted(numerical_dict.keys())
    for number in numerical_list:
        index = numerical_dict[number]
        sorted_file_list.append(file_list[index])
        
    return sorted_file_list

def get_openmm_center_of_mass_com(system, positions, group1):
    """
    Compute the center of mass of a group of atoms given an OpenMM
    system, the system's atomic positions, and the indices of the
    atoms of interest.
    
    Parameters:
    -----------
    system : System()
        The OpenMM System object which contains the atomic masses.
        
    positions : Quantity
        A n-by-3 array of positions of all atoms in the system. Can
        be the output of State.getPositions()
        
    group1 : list
        A list of integers representing atom indices of the group
        the center of mass is being computed for.
        
    Returns:
    --------
    com : Quantity
        A 1x3 array representing the x,y,z coordinates of the center of
        mass of the group of atoms.
    """
    total_mass = 0.0 * unit.daltons
    com = np.array([0, 0, 0]) * unit.nanometers * unit.daltons
    for index in group1:
        atom_mass = system.getParticleMass(index)
        total_mass += atom_mass
        com += atom_mass * positions[index]
        
    com /= total_mass
    return com

def get_openmm_rmsd(system, positions, ref_positions, group1):
    """
    Compute the RMSD between two groups of atoms given an OpenMM
    system, the system's atomic positions, a set of reference 
    positions, and the indices of the atoms of interest.
    
    Parameters:
    -----------
    system : System()
        The OpenMM System object which contains the atomic masses.
        
    positions : Quantity
        A n-by-3 array of positions of all atoms in the system. Can
        be the output of State.getPositions()
        
    ref_positions : Quantity
        A n-by-3 array of reference positions of all atoms in the 
        system.
        
    group1 : list
        A list of integers representing atom indices of the group
        the center of mass is being computed for.
        
    Returns:
    --------
    com : Quantity
        A 1x3 array representing the x,y,z coordinates of the center of
        mass of the group of atoms.
    """
    sq_dist_sum = np.array([0, 0, 0]) * unit.nanometers
    for index in group1:
        sq_dist_sum += np.linalg.norm(positions[index], ref_positions[index])**2
    
    return np.sqrt(sq_dist_sum / len(group1))

def get_box_vectors_from_pdb(pdb_filename):
    """
    Extract the box_vectors from the CRYST line in a pdb file.
    """
    pdb_structure = parmed.load_file(pdb_filename, skip_bonds=True)
    assert pdb_structure.box_vectors is not None, "No box vectors "\
    "found in {}. ".format(pdb_filename) \
    + "Box vectors for an anchor must be defined with a CRYST "\
    "line within the PDB file, or explicitly set in the model "\
    "input XML file."
    length1 = np.linalg.norm(pdb_structure.box_vectors[0]\
                             .value_in_unit(unit.nanometers))
    length2 = np.linalg.norm(pdb_structure.box_vectors[1]\
                             .value_in_unit(unit.nanometers))
    length3 = np.linalg.norm(pdb_structure.box_vectors[2]\
                             .value_in_unit(unit.nanometers))
    assert (length1 > 0.0) and (length2 > 0.0) and (length3 > 0.0), \
        "Box vector length(s) of zero detected. Please check box vectors in " \
        "input PDB."
    return pdb_structure.box_vectors

def convert_openmm_to_python_expr(old_function_str):
    """
    Convert the expression from a form used by OpenMM
    to a form used by Python.
    """
    new_function_str = re.sub(r"\^", "**", old_function_str)
    return new_function_str

class Box_vectors(Serializer):
    """
    A box vector object that contains the a, b, and c vectors in units
    of nanometers.
    """
    
    def __init__(self):
        self.ax = 1.0
        self.ay = 0.0
        self.az = 0.0
        self.bx = 0.0
        self.by = 1.0
        self.bz = 0.0
        self.cx = 0.0
        self.cy = 0.0
        self.cz = 1.0
        
    def to_quantity(self):
        """
        Convert this object to a Quantity object that could be used
        in OpenMM or parmed.
        """
        box_vector = unit.Quantity(
            [[self.ax, self.ay, self.az], 
             [self.bx, self.by, self.bz], 
             [self.cx, self.cy, self.cz]], unit=unit.nanometer)
        return box_vector
    
    def from_quantity(self, quantity):
        """
        Fill out the entries within this object from a Quantity
        object representing box vectors.
        """
        
        values = quantity.value_in_unit(unit.nanometer)
        self.ax = float(values[0][0])
        self.ay = float(values[0][1])
        self.az = float(values[0][2])
        self.bx = float(values[1][0])
        self.by = float(values[1][1]) # TEST THESE FUNCTIONS
        self.bz = float(values[1][2])
        self.cx = float(values[2][0])
        self.cy = float(values[2][1])
        self.cz = float(values[2][2])
        return
    
    def from_6_vector(self, box_6_vector):
        """
        Sometimes box vectors are represented using 6 numbers, of which
        the first three are vector lengths, and the next three are the
        angles between them. Compute the full 3x3 box_vectors and assign
        to this object.
        
        box_6_vector is assumed to be in Angstroms and degrees.
        """
        
        # convert box_6_vector from Angstroms to nm
        box_6_vector[0] *= 0.1
        box_6_vector[1] *= 0.1
        box_6_vector[2] *= 0.1
        
        self.ax = float(box_6_vector[0])
        self.ay = 0.0
        self.az = 0.0
        self.bx = float(box_6_vector[1] \
                        * math.cos(box_6_vector[5]*math.pi/180.0))
        self.by = float(box_6_vector[1] \
                        * math.sin(box_6_vector[5]*math.pi/180.0))
        self.bz = 0.0
        self.cx = (box_6_vector[2] * math.cos(box_6_vector[4]*math.pi/180.0))
        if box_6_vector[5] == 90.0:
            self.cy = 0.0
        else:
            self.cy = float(box_6_vector[2] * (
                math.cos(box_6_vector[3]*math.pi/180.0) \
                - math.cos(box_6_vector[4]*math.pi/180.0) \
                * math.cos(box_6_vector[5]*math.pi/180.0)) \
                / math.sin(box_6_vector[5]*math.pi/180.0))
        self.cz = float(math.sqrt(box_6_vector[2]**2 - self.cx**2 \
            - self.cy**2))
        return
    
    def to_6_vector(self):
        """
        Sometimes box vectors are represented using 6 numbers, of which
        the first three are vector lengths, and the next three are the
        angles between them. Use this object's full 3x3 box_vectors and 
        return the equivalent 6-vector.
        """
        a = np.array([self.ax, self.ay, self.az])
        b = np.array([self.bx, self.by, self.bz])
        c = np.array([self.cx, self.cy, self.cz])
        a_length = np.linalg.norm(a)
        b_length = np.linalg.norm(b)
        c_length = np.linalg.norm(c)
        alpha = (180.0/math.pi) * math.acos(np.dot(b, c) / (b_length * c_length))
        beta = (180.0/math.pi) * math.acos(np.dot(c, a) / (c_length * a_length))
        gamma = (180.0/math.pi) * math.acos(np.dot(a, b) / (a_length * b_length))
        six_vector = [a_length, b_length, c_length, alpha, beta, gamma]
        return six_vector
    
    def get_volume(self):
        """
        Compute the volume of the box defined by these box vectors.
        Units in nm^3
        """
        
        A = np.array([self.ax, self.ay, self.az])
        B = np.array([self.bx, self.by, self.bz])
        C = np.array([self.cx, self.cy, self.cz])
        volume = abs(np.dot(A, np.cross(B, C)))
        return volume
    
    def get_min_length(self):
        """
        Get the dimension of largest magnitude.
        """
        lengths = []
        lengths.append(np.linalg.norm(np.array([self.ax, self.ay, self.az]))) 
        lengths.append(np.linalg.norm(np.array([self.bx, self.by, self.bz])))
        lengths.append(np.linalg.norm(np.array([self.cx, self.cy, self.cz])))
        return min(lengths)
        
class Langevin_integrator_settings(Serializer):
    """
    Contains settings used in the initialization of a Langevin-
    type integrator in either OpenMM or NAMD.
    
    Attributes:
    -----------
    friction_coefficient : float, Default 1.0
        The friction coefficient (gamma) value (in units of inverse
        picoseconds) used by Langevin integrators to create drag and
        stochastic 'kicks' to the particles.
        
    target_temperature : float, Default 298.15
        The temperature value (in units of Kelvin) which the integrator
        will fluctuate close to.
        
    random_seed : int, Default 0
        Set a seed for the random number generator.
        
    timestep : float, Default 0.002
        The length of time taken by a simulation step (in units of 
        picoseconds).
        
    rigid_tolerance :  float, Default 1e-6
        The acceptable tolerance when placing atoms of a rigid bond
        length or angle.
        
    integrator_type : str, Default "langevin"
        Which type of Langevin integrator to use. Options include "langevin",
        "langevinMiddle", ...
    """
    
    def __init__(self):
        self.friction_coefficient = 1.0
        self.target_temperature = 298.15
        self.random_seed = 0
        self.timestep = 0.002
        self.rigid_tolerance = 1e-6
        self.integrator_type = "langevin"
        return
    
class Barostat_settings_openmm(Serializer):
    """
    Contains settings used in the initialization of a 
    MonteCarloBarostat for the control of pressure within the 
    OpenMM simulation.
    
    Attributes:
    -----------
    target_pressure : float, Default 1.0
        The pressure (in units of bar) that this barostat will attempt
        to maintain.
        
    target_temperature : float, Default 298.15
        The temperature (in units of Kelvin) at which this system is
        being maintained.
        
    frequency : int, Default 25
        the Monte Carlo frequency (in time steps) at which pressure 
        changes should be attempted. 
    """
    
    def __init__(self):
        self.target_pressure = 1.0
        self.target_temperature = 298.15
        self.frequency = 25
        self.membrane = False
        return
    
class Barostat_settings_namd(Serializer):
    """
    Contains settings used in the initialization of a 
    LangevinPiston barostat for the control of pressure within the 
    NAMD simulation.
    
    Attributes:
    -----------
    target_pressure : float, Default 1.0
        The pressure (in units of bar) that this barostat will attempt
        to maintain.
        
    target_temperature : float, Default 298.15
        The temperature (in units of Kelvin) at which this system is
        being maintained.
        
    oscillation_frequency : float, Default 0.2
        describes barostat oscillation time scale (in ps) for the 
        piston method.
        
    decay_timescale : float, Default 0.1
        Describes the barostat damping time scale (in ps) for the 
        Langevin piston method.
    """
    
    def __init__(self):
        self.target_pressure = 1.0
        self.target_temperature = 298.15
        self.oscillation_frequency = 0.2
        self.decay_timescale = 0.1
        return
    
class Cuda_platform_settings(Serializer):
    """
    Contains settings used in the CUDA platform of OpenMM simulations.
    
    Attributes:
    -----------
    cuda_device_index : str, Default "0"
        The indices of the GPU device(s) on which to perform the 
        simulation. Examples: "0", "1", "0,1", etc... Inputs
        are the same as they would be for the CudaDeviceIndex
        argument of an OpenMM simulation.
        
    cuda_precision : str, Default "mixed"
        The precision to use within the CUDA kernels. Options which
        OpenMM accepts include: "single", "mixed", and "double".
    """
    
    def __init__(self):
        self.cuda_device_index = "0"
        self.cuda_precision = "mixed"
    
    def make_properties_dict(self):
        """
        Return a platform properties dictionary which goes into the 
        arguments of an OpenMM's Simulation() object.
        """
        properties = {'CudaDeviceIndex': str(self.cuda_device_index), 
                      'CudaPrecision': str(self.cuda_precision)}
        return properties
    
class Openmm_settings(Serializer):
    """
    Contains settings used by OpenMM MD simulation program.
    
    Attributes:
    -----------
    nonbonded_method : str, Default "PME"
        The method to use to account for long-range nonbonded forces.
        This argument may be any of the options for the nonbondedMethod
        argument in an OpenMM System() object.
        
    nonbonded_cutoff : float, Default 0.9
        The distance after which VDW nonbonded forces are cut off in
        units of nanometers. This argument is supplied to the 
        nonbondedCutoff argument of an OpenMM System() object.
        
    constraints : str, Default "hbonds"
        The constraints method to use for bonds and angles within the
        system. This argument is supplied to the constraints argument
        of the OpenMM System() object.
        
    hydrogenMass : float or None, Default None
        This parameter may be used for hydrogen mass repartitioning.
        If a float is provided, then mass will be shifted from every
        hydrogen atom's bonded heavy atom to the hydrogen, such that
        the hydrogen atoms will have this value (in AMU) for their
        masses. If None, then the default hydrogen mass is used.
        
    rigidWater : bool, Default True
        If True, then water bonds and angles will be made rigid.
        
    langevin_integrator : Langevin_integrator_settings
        Contains settings for the Langevin integrator.
        
    barostat : Barostat_settings_openmm
        Contains settings for the Barostat.
        
    cuda_platform_settings : Cuda_platform_settings or None
        Contains settings for the Cuda platform. If None, then the
        Reference platform is used instead.
        
    reference_platform : bool, Default False
        If cuda_platform_settings is not None, then this should be set
        to True to indicate that the reference platform must be used.
        
    run_minimization : bool, Default False
        Whether to run minimizations before doing any MMVT
        calculations.
        
    initial_temperature : float, Default 298.15
        The temperature (in units of Kelvin) to use to initialize the 
        atomic velocities randomly .
    """
    def __init__(self):
        self.nonbonded_method = "PME"
        self.nonbonded_cutoff = 0.9
        self.constraints = "hbonds"
        self.hydrogenMass = None
        self.rigidWater = True
        self.langevin_integrator = Langevin_integrator_settings()
        self.barostat = None
        self.cuda_platform_settings = Cuda_platform_settings()
        self.reference_platform = False
        self.run_minimization = False
        self.initial_temperature = 298.15
        return
    
class Namd_settings(Serializer):
    """
    Contains settings used by NAMD simulation program.
    
    Attributes:
    -----------
    watermodel : str
        Which water model to use in the simulation. In NAMD, options
        include: 'tip3p', 'tip4p', and 'swm4'.
        
    PMEGridSpacing : float, Default 0.1
        The spacing (in Angstroms) of the PME grid in NAMD.
        
    nonbonded_cutoff : float
        The cutoff (in Angstroms) for the Van der Waals interactions.
    
    constraints : str, Default "hbonds"
        The constraints method to use for bonds within the system. Can
        be "none", "water", or "hbonds".
        
    langevin_integrator : Langevin_integrator_settings
        Contains settings for the Langevin integrator.
        
    barostat : Barostat_settings_namd
        Contains settings for the Barostat.
        
    run_minimization : bool, Default False
        Whether to run minimizations before doing any calculations.
        
    initial_temperature : float, Default 298.15
        The temperature (in units of Kelvin) to use to initialize 
        the atomic velocities randomly .
    
    eval_stride : int, default 10
        How frequently (in timesteps) to check for boundary crossings.
    """
    
    def __init__(self):
        self.watermodel = ""
        self.PMEGridSpacing = 0.1
        self.nonbonded_cutoff = 0.9
        self.constraints = "hbonds"
        self.langevin_integrator = Langevin_integrator_settings()
        self.barostat = None #Barostat_settings()
        self.run_minimization = False
        self.initial_temperature = 298.15
        self.eval_stride = 10
        return

class Toy_settings(Serializer):
    """
    Contains settings used by Toy simulation program.
    
    Attributes:
    -----------
        
    potential_energy_expression : str
        The expression that can be used in OpenMM to direct particle 
        motions.
        
    num_particles : int
        The number of particles in the system.
    
    masses : list
        A list of particle masses
    """
    
    def __init__(self):
        super(Toy_settings, self).__init__()
        self.potential_energy_expression = None
        self.num_particles = -1
        self.masses = []
        return

class Browndye_settings(Serializer):
    """
    Read and parse the outputs from the BrownDye program, which runs
    the BD stage of the SEEKR2 calculation
    
    Attributes:
    -----------
    browndye_bin_dir : str, Default ""
        A path to the BrownDye programs. If Browndye's bin/ directory 
        has been added to system $PATH, then this string can be empty.
        
    receptor_pqr_filename : str
        The path to the receptor molecule's PQR-format file.
        
    ligand_pqr_filename : str
        The path to the ligand molecule's PQR-format file.
        
    apbs_grid_spacing : float
        The resolution (in Angstroms) of the APBS (electrostatics) 
        grid.
        
    n_threads : int, Default 1
        The number of cores to use for the BrownDye calculation.
        
    recompute_ligand_electrostatics : bool, Default True
        Whether ligand electrostatics should be computed for the
        bd_milestone portions of the calculation.
        
    debye_length : float
        The Debye length (in Angstroms) as computed by APBS and
        used by BrownDye.
        
    ghost_indices_rec : list
        The index(es) of the atom(s) which are the Ghost atoms - 
        that is, they are used to identify the center of mass of
        important sites.
        
    ghost_indices_lig : list
        Similar to ghost_indices_rec, only these exist on the ligand.
    """
    
    def __init__(self):
        self.browndye_bin_dir = ""
        self.receptor_pqr_filename = ""
        self.ligand_pqr_filename = ""
        self.apbs_grid_spacing = -1.0
        self.n_threads = 1
        self.recompute_ligand_electrostatics = True
        self.debye_length = -1.0
        self.ghost_indices_rec = []
        self.ghost_indices_lig = []
        return
    
class Amber_params(Serializer):
    """
    Contains parameters for an amber simulation.
    
    Attributes:
    -----------
    prmtop_filename : str
        The AMBER parameter/topology file for the system
        
    box_vectors : None or Box_vectors
        The 3 vectors which describe the shape of the simulation box.
        If empty, the box vectors are taken from the CRYST line of the 
        PDB file.
        
    pdb_coordinates_filename : str
        The path to the PDB file from which to obtain the atomic
        coordinates. 
    """
    
    def __init__(self):
        self.prmtop_filename = ""
        self.box_vectors = None
        self.pdb_coordinates_filename = ""
        return

def same_amber_params(amber_params1, amber_params2):
    """
    Returns True if both amber_params objects are the same
    """
    if amber_params1 is None:
        if amber_params2 is None:
            return True
        else:
            return False
    if amber_params2 is None:
        return False
    
    if amber_params1.prmtop_filename != amber_params2.prmtop_filename:
        return False
    if amber_params1.box_vectors != amber_params2.box_vectors:
        return False
    if amber_params1.pdb_coordinates_filename \
            != amber_params2.pdb_coordinates_filename:
        return False
    return True

class Forcefield_params(Serializer):
    """
    Contains parameters for an OpenMM simulation starting from a set
    of XML forcefields.
    
    Attributes:
    -----------
    built_in_forcefield_filenames : list
        A list of forcefield XML files that come with OpenMM, 
        used to start the simulation. Examples include: 
        'amber14-all.xml', 'amber14/tip3pfb.xml', ...
        
    custom_forcefield_filenames : list
        A list of forcefield XML files that *did not* come with 
        OpenMM. Presumably, these have been created custom for the
        system.
        
    pdb_filename : str
        The PDB input coordinates and box vectors for the system.
    
    system_filename : str
        The System() serialized object that can be used to reconstruct a
        complete OpenMM simulation.
    
    box_vectors : None or Box_vectors
        The 3 vectors which describe the shape of the simulation box.
        If None, then the box vectors are taken from the 
        PDB file defined by the pdb_filename argument.
        
    """
    def __init__(self):
        self.built_in_forcefield_filenames = []
        self.custom_forcefield_filenames = []
        self.pdb_coordinates_filename = ""
        self.system_filename = ""
        self.box_vectors = None
        return

class Charmm_params(Serializer):
    """
    Contains parameters for a CHARMM simulation.
    
    Attributes:
    -----------
    psf_filename : str
        The CHARMM psf file for the system
        
    charmm_ff_files : list
        The CHARMM forcefield files for the system
        
    box_vectors : None or Box_vectors
        The 3 vectors which describe the shape of the simulation box.
        
    pdb_coordinates_filename : str
        The path to the PDB file from which to obtain the atomic
        coordinates.
    """
    
    def __init__(self):
        self.psf_filename = ""
        self.charmm_ff_files = []
        self.box_vectors = None
        self.pdb_coordinates_filename = ""
        return

def same_charmm_params(charmm_params1, charmm_params2):
    """
    Returns True if both charmm_params objects are the same
    """
    if charmm_params1 is None:
        if charmm_params2 is None:
            return True
        else:
            return False
    if charmm_params2 is None:
        return False
    
    if charmm_params1.psf_filename != charmm_params2.psf_filename:
        return False
    if charmm_params1.charmm_ff_files != charmm_params2.charmm_ff_files:
        return False
    if charmm_params1.box_vectors != charmm_params2.box_vectors:
        return False
    if charmm_params1.pdb_coordinates_filename \
            != charmm_params2.pdb_coordinates_filename:
        return False
    return True

class Ion(Serializer):
    """
    An ion for input to BrownDye and APBS calculations.
    
    Attributes:
    -----------
    radius : float
        The radius (in Angstroms)
    
    charge : float
        The charge (in proton charge) of the ion.
        
    conc : float
        The concentration (in moles/liter) of the ion in the solution.
    """
    
    def __init__(self):
        self.radius = -1.0
        self.charge = -1.0
        self.conc = -1.0
        return

class K_on_info(Serializer):
    """
    Information needed to compute K-on-related quantities.
    
    Attributes:
    -----------
    bd_milestones : list
        A list of BD_milestone() objects which act as reaction 
        surfaces in the Browndye simulations.
        
    b_surface_directory : str
        The directory where to find the results XML file output by
        Browndye which applies to trajectories started from the
        b-surface.
        
    b_surface_num_trajectories : int
        The number of trajectories that will be run from the b-surface.
        
    bd_output_glob : str
        A glob which can be used to select the output XML files
        produced by Browndye within the directory argument above.
        
    ions : list
        A list of Ion() objects which will be provided for APBS
        calculations.
    """
    
    def __init__(self):
        self.bd_milestones = []
        self.b_surface_directory = "b_surface"
        self.b_surface_num_trajectories = -1
        self.bd_output_glob = BROWNDYE_OUTPUT
        self.ions = []
        return

class BD_milestone(Serializer):
    """
    BD simulations must end on a spherical milestone, and one
    additional set of BD simulations must be run for each of these.
    
    index : int
        The index of this BD_milestone, since there can be many of
        them. These are numbered starting from 0, and have no bearing
        on the numbering of the MD milestones.
        
    directory : str
        The directory for this BD_milestone and all of its sub-files.
        
    bd_output_glob : str
        A glob for BrownDye output files containing transition info
        related to this milestone.
        
    name : str
        A unique name for this milestone
    
    outer_milestone : None or Milestone()-inheriting class
        Each BD_milestone has an outer one, where the simulations
        start, 
        
    inner_milestone : None or Milestone()-inheriting class
        Each BD_milestone also has an inner one, where some of the 
        simulations may end.
        
    num_trajectories : int
        The number of BD_simulations per member of the FHPD to launch.
        
    receptor_indices : list
        The atom indices (starting from zero) whose center of mass
        will define the center of the binding site.
    
    ligand_indices : list
        The atom indices (starting from zero) whose center of mass
        defines the center of the ligand (or an important part of the 
        ligand).
        
    extracted_directory : str
        The directory where extractions from the b-surface simulations
        will be performed.
        
    fhpd_directory : str
        The directory where the FHPD simulations will be run from for
        this BD_milestone.
    
    max_b_surface_trajs_to_extract : int
        After the b-surface simulation, members of the encounter 
        complexes will be extracted to construct the FHPD. Then, these
        will be run as their own independent BD simulations. This
        parameter defines the maximum number of structures to extract
        from the FHPD and run simulations for.
    """
    
    def __init__(self):
        self.index = -1
        self.directory = ""
        self.bd_output_glob = BROWNDYE_OUTPUT
        self.name = ""
        self.outer_milestone = None
        self.inner_milestone = None
        #self.num_trajectories = -1
        self.receptor_indices = []
        self.ligand_indices = []
        # TODO: remove
        #self.extracted_directory = "extracted_from_b_surface"
        #self.fhpd_directory = "first_hitting_point_distribution"
        #self.max_b_surface_trajs_to_extract = -1
        return
    
class Milestone(Serializer):
    """
    Milestones represent the boundaries within the simulation.
    A given anchor object may contain any number of milestone
    objects.
    
    Attributes:
    -----------
    index : int
        This is index of a milestone instance across the entire model.
        Two milestone instances in different anchor objects should
        share the same index if the milestones represent the same
        boundary between two adjacent anchors.
        
    neighbor_anchor_index : int
        The index of the anchor instance on the other side of this
        milestone from the anchor instance that contains this
        milestone object.
        
    alias_index : int
        This field can be used to represent a simpler numbering scheme
        within an anchor. This is convenient for the backend 
        plugins.
        
    cv_index : int
        The index of which collective variable in the model gives rise
        to this milestone.
        
    variables : dict
        This milestone has many variables associated with it. The keys
        of the variables dictionary are the variable names that appear
        within the CV's expression, while the values of the dictionary
        are specific to each milestone's location.
    """
    def __init__(self):
        self.index = -1
        self.neighbor_anchor_index = -1
        self.alias_index = -1
        self.cv_index = -1
        self.variables = {}
        self.is_source_milestone = False
        
    def get_CV(self, model):
        """
        Obtain and return the collective variable object which
        applies to this milestone.
        """
        
        return model.collective_variables[self.cv_index]

class Model(Serializer): 
    """
    The Model object contains all the parameters, settings, 
    directory names, and file names which are relevant to a
    SEEKR2 calculation.
    
    Attributes
    ----------
    temperature : float
        The temperature in Kelvin of the system in all simulations 
        and calculations
        
    calculation_type : str
        A string representing the type of calculation that is being
        done by this Model object. Options include "MMVT" or "Elber"
        
    calculation_settings : MMVT_settings() or Elber_settings()
        Settings general to any MD engine, including the number of
        timesteps to run, backup intervals, etc.
        
    anchor_rootdir : str
        The directory which contains all the anchor directories. If
        this string is equal to "." in the XML, then the directory 
        in which the XML file was found is used as the anchor root 
        directory.
    
    num_anchors : int
        The number of anchor points and, equivalently, the number of Voronoi
        cells
        
    num_milestones : int
        The total number of milestones (boundaries) that separate the
        Voronoi cells of this model.
        
    openmm_settings : Openmm_settings
        The Openmm_settings() object for this model. It contains all
        the settings that could be used within an OpenMM simulation.
        
    namd_settings : Namd_settings
        The Namd_settings() object for this model. It contains all the
        settings that could be used in a NAMD simulation.
        
    toy_settings : Toy_settings
        The Toy_settings() object for this model. It contains all
        the settings that could be used within a Toy OpenMM simulation.
        
    browndye_settings : Browndye_settings
        The Browndye_settings() object for this model. It contains all
        the settings that could be used within a Browndye simulation.
        
    k_on_info : K_on_info
        The K_on_info() object which represents the settings and 
        parameters needed for calculating the k-on for this model.
        
    collective_variables : list
        A list of Collective_variable() (CV) instances which contain 
        information for each of the CVs in this model.
        
    anchors : list
        A list of Anchor() object instances for this model.
    """
    
    def __init__(self):
        self.temperature = 0.0
        self.calculation_type = "Unknown"
        self.calculation_settings = None
        self.anchor_rootdir = "."
        self.num_anchors = -1
        self.num_milestones = -1
        self.openmm_settings = None
        self.namd_settings = None
        self.toy_settings = None
        self.browndye_settings = None
        self.k_on_info = None
        self.collective_variables = []
        self.anchors = []
        return
    
    def get_type(self):
        """
        Return the calculation_type of this model, and check for
        erronious entries to this field.
        """
        if self.calculation_type.lower() == "elber":
            return "elber"
        elif self.calculation_type.lower() == "mmvt":
            return "mmvt"
        else:
            error_msg = "Calculation type not available: "\
                + "{}. Available types are 'elber' and 'mmvt'.".format(
                    self.calculation_type)
            raise Exception(error_msg)
        
    def using_bd(self):
        """
        Return whether this model is using BD
        """
        if self.k_on_info is None:
            return False
        else:
            return True
    
    def using_toy(self):
        """
        Return whether this model is using BD
        """
        if self.toy_settings is None:
            return False
        else:
            return True
    
    def get_timestep(self):
        """
        Return the timestep used in this model
        """
        if self.openmm_settings is not None:
            if self.openmm_settings.langevin_integrator is not None:
                timestep = self.openmm_settings.langevin_integrator.timestep
            else:
                raise Exception("Settings not provided for available "\
                        "integrator type(s).")
        elif self.namd_settings is not None:
            if self.namd_settings.langevin_integrator:
                timestep = self.namd_settings.langevin_integrator.timestep
            else:
                raise Exception("Settings not provided for available "\
                        "integrator type(s).")
        else:
            raise Exception("No settings provided for available simulators.")
        return timestep
    
    def get_bulk_index(self):
        """
        Get the anchor index of the bulk state. If there is no bulk state,
        return None.
        """
        bulk_index = None
        for alpha, anchor in enumerate(self.anchors):
            if anchor.bulkstate:
                assert bulk_index is None, "Only one bulk state is allowed "\
                    "in model"
                bulk_index = alpha
        return bulk_index
    
def load_model(model_file, directory=None):
    """
    
    """
    assert os.path.exists(model_file), \
        "No such file or directory: {}.".format(model_file)
    model = Model()
    model.deserialize(model_file)
    if directory is not None:
        model.anchor_rootdir = os.path.abspath(directory)
    elif model.anchor_rootdir == ".":
        model_dir = os.path.dirname(model_file)
        model.anchor_rootdir = os.path.abspath(model_dir)
    return model

def save_model(model, model_file):
    """
    
    """
    dont_modify_warning = "WARNING: this file is automatically generated and "\
        "should not be modified by hand. Instead, modify the input XML and "\
        "re-run prepare.py."
    model.serialize(model_file, xml_header=dont_modify_warning)
    return

def parse_xml_list(variable):
    """
    Take 'variable', and if it's a list, then return as-is, but if it's
    an argument representing a command, like range(), then execute
    the command to produce a list.
    """
    if isinstance(variable, list):
        return variable
    elif isinstance(variable, str):
        # Match the range() command
        m = re.match(r"range\((\d+),?(\d+)?,?(\d+)?\)", variable)
        if m:
            groups = m.groups()
            if groups[2] is None:
                if groups[1] is None:
                    # one argument
                    return list(range(0,int(m[1])))
                else:
                    # two arguments
                    return list(range(int(m[1]),int(m[2])))
            else:
                return list(range(int(m[1]), int(m[2]), int(m[3])))
            
        else:
            raise Exception("Invalid XML input: {}. ".format(variable) \
                            +"Available command(s): range(a[,b][,c]).")
    elif variable is None:
        return None
    else:
        raise Exception("Invalid XML input: {}".format(variable))
    
def get_anchor_pdb_filename(anchor):
    """
    Obtain and return an anchor's starting structure filename.
    
    Parameters
    ----------
    anchor : Anchor()
        For a given Anchor object, return a PDB file name from the input 
        parameters, whether they are Amber, Charmm, etc.
    """
    if anchor.amber_params is not None:
        return anchor.amber_params.pdb_coordinates_filename
    
    if anchor.forcefield_params is not None:
        return anchor.forcefield_params.pdb_coordinates_filename
        
    if anchor.charmm_params is not None:
        return anchor.charmm_params.pdb_coordinates_filename
    
    return ""