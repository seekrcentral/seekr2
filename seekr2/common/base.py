"""
common/base.py

This module of seekr2 contains base classes, constants, and other
objects used for both MMVT and Elber milestoning.
"""

import re

from parmed import unit

import seekr2.libraries.serializer.serializer as serializer

BROWNDYE_OUTPUT = "results.xml"

def strBool(bool_str):
    """
    Take the string "true" or "false" of any case and returns a 
    boolean object.
    """
    if bool_str.lower() == "true":
        return True
    elif bool_str.lower() == "false":
        return False
    else:
        raise Exception(
            "argument for strBool must be string either 'True' or 'False'.")

def order_files_numerically(file_list):
    """
    If there is a list of files, order them numerically, not
    alphabetically and return the sorted list of files.
    
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
        numbers = re.findall(r"\d+", file_name)
        numbers = tuple([int(j) for j in numbers])
        numerical_dict[numbers] = i
        
    numerical_list = sorted(numerical_dict.keys())
    for number in numerical_list:
        index = numerical_dict[number]
        sorted_file_list.append(file_list[index])
        
    return sorted_file_list

class Box_vectors(serializer.Serializer):
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
        in OpenMM or another simtk program.
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
        print("values:", values)
        self.ax = values[0][0]
        self.ay = values[0][1]
        self.az = values[0][2]
        self.bx = values[1][0]
        self.by = values[1][1] # TEST THESE FUNCTIONS
        self.bz = values[1][2]
        self.cx = values[2][0]
        self.cy = values[2][1]
        self.cz = values[2][2]
        return

class Langevin_integrator_settings(serializer.Serializer):
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
    """
    def __init__(self):
        self.friction_coefficient = 1.0
        self.target_temperature = 298.15
        self.random_seed = 0
        self.timestep = 0.002
        self.rigid_tolerance = 1e-6
        return
    
class Barostat_settings_openmm(serializer.Serializer):
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
        changes should be attempted. This only affects OpenMM's 
        barostat.
    """
    def __init__(self):
        self.target_pressure = 1.0
        self.target_temperature = 298.15
        self.frequency = 25
        return
    
class Barostat_settings_namd(serializer.Serializer):
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
    
class Cuda_platform_settings(serializer.Serializer):
    """
    Contains settings used in the CUDA platform of OpenMM simulations.
    
    Attributes:
    -----------
    cuda_device_index : str, Default "0"
        The indices of the GPU device(s) on which to perform the 
        simulation. Example: "0", "1", "0,1", etc...
        
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
        properties = {'CudaDeviceIndex': self.cuda_device_index, 
                      'CudaPrecision': self.cuda_precision}
        return properties
    
class Openmm_settings(serializer.Serializer):
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
        the hydrogen atoms will have that value (in AMU) for their
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
        The temperature to use to initialize the atomic velocities
        randomly in units of Kelvin.
        
    energy_reporter_frequency : int or None, Default None
        The frequency to report system state information, including
        energy. If set to None, then this information is never printed.
        
    restart_checkpoint_frequency : int or None, Default None
        The frequency to write a restart checkpoint for the system to
        be easily restarted. None means do not write backups.
        
    trajectory_reporter_frequency : int or None, Default None
        The frequency to write frames of the trajectory to file.
        None means do not write to a trajectory file.
        
    total_simulation_length : int, Default 30000
        The total number of time steps to simulate.
    """
    def __init__(self):
        self.nonbonded_method = "PME"
        self.nonbonded_cutoff = 0.9
        self.constraints = "hbonds"
        self.hydrogenMass = None
        self.rigidWater = True
        self.langevin_integrator = Langevin_integrator_settings()
        self.barostat = None #Barostat_settings()
        self.cuda_platform_settings = Cuda_platform_settings()
        self.reference_platform = False
        self.run_minimization = False
        self.initial_temperature = 298.15
        self.energy_reporter_frequency = None
        self.restart_checkpoint_frequency = None
        self.trajectory_reporter_frequency = None
        return
    
class Namd_settings(serializer.Serializer):
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
        The temperature to use to initialize the atomic velocities
        randomly in units of Kelvin.
        
    energy_reporter_frequency : int or None, Default None
        The frequency to report system state information, including
        energy. If set to None, then this information is never printed.
        
    restart_checkpoint_frequency : int or None, Default None
        The frequency to write a restart checkpoint for the system to
        be easily restarted. None means do not write backups.
        
    trajectory_reporter_frequency : int or None, Default None
        The frequency to write frames of the trajectory to file.
        None means do not write to a trajectory file.
        
    total_simulation_length : int, Default 30000
        The total number of time steps to simulate.
        
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
        self.energy_reporter_frequency = None
        self.restart_checkpoint_frequency = None
        self.trajectory_reporter_frequency = None
        self.eval_stride = 10
        return
    
class Browndye_settings(serializer.Serializer):
    """
    Read and parse the outputs from the BrownDye program, which runs
    the BD stage of the MMVT SEEKR calculation
    
    Attributes:
    -----------
    browndye_bin_dir : str, Default ""
        A path to the BrownDye programs. If added to sytem $PATH, then
        this string can be empty.
        
    receptor_pqr_filename : str
        The path to the receptor molecule's filename.
        
    ligand_pqr_filename : str
        The path to the ligand molecule's filename.
        
    apbs_grid_spacing : float
        The resolution (in Angstroms) of the APBS (electrostatics) 
        grid.
        
    n_threads : int, Default 1
        The number of cores to use for the BrownDye calculation.
        
    recompute_ligand_electrostatics : boo, Default True
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
    
class Amber_params(serializer.Serializer):
    """
    Contains parameters for an amber simulation.
    
    Attributes:
    -----------
    prmtop_filename : str
        The AMBER parameter/topology file for the system
        
    inpcrd_filename : str
        The AMBER input coordinates file for the system
        
    box_vectors : None or Box_vectors
        The 3 vectors which describe the shape of the simulation box.
        If None, then the box vectors are taken from the 
        inpcrd file defined by the inpcrd_filename argument.
        
    pdb_coordinates_filename : str
        The path to the PDB file from which to obtain the atomic
        coordinates. If it's an empty string or None, then the
        coordinates are taken from the inpcrd file.
    
    """
    def __init__(self):
        self.prmtop_filename = ""
        self.inpcrd_filename = ""
        self.box_vectors = None
        self.pdb_coordinates_filename = ""
        return
    
class Forcefield_params(serializer.Serializer):
    """
    Contains parameters for an openmm simulation starting from a set
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
        
    box_vectors : None or Box_vectors
        The 3 vectors which describe the shape of the simulation box.
        If None, then the box vectors are taken from the 
        PDB file defined by the pdb_filename argument.
        
    """
    def __init__(self):
        self.built_in_forcefield_filenames = []
        self.custom_forcefield_filenames = []
        self.pdb_filename = ""
        self.box_vectors = None
        return

class Ion(serializer.Serializer):
    """
    An input XML file ion for input to BrownDye and APBS calculations.
    
    Attributes:
    -----------
    radius : float
        The radius (in nanometers)
    
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

class K_on_info(serializer.Serializer):
    """
    Information needed to compute K-on-related quantities.
    
    Attributes:
    -----------
    source_milestones : list
        A list of milestones which act as 'sources' for the ligand. 
        This will be a list of integers whose indices apply to the 
        model's milestones.
        
    b_surface_directory : str
        The directory where to find the results XML file output by
        Browndye which applies to trajectories started from the
        b-surface.
        
    b_surface_num_trajectories : int
        The number of trajectories that will be run from the b-surface.
        
    bd_mmvt_output_glob : str
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
        self.bd_mmvt_output_glob = "results.xml"
        self.ions = []
        return

class BD_milestone(serializer.Serializer):
    """
    BD simulations must end on a spherical milestone, and one
    additional set of BD simulations must be run for each of these.
    
    index : int
        The index of this BD_milestone, since there can be many of
        them. These are numbered starting from 0, and have no bearing
        on the numbering of the MD milestones.
        
    directory : str
        The directory for this BD_milestone and all of its sub-files.
        
    bd_mmvt_output_glob : str
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
        defines the center of the ligand (or an import part of the 
        ligand).
        
    extracted_directory : str
        The directory where extractions from the b-surface simulations
        will be performed.
        
    fhpd_directory : str
        The directory where the FHPD simulations will be run from for
        this BD_milestone.
    
    """
    def __init__(self):
        self.index = -1
        self.directory = ""
        self.bd_mmvt_output_glob = BROWNDYE_OUTPUT
        self.name = ""
        self.outer_milestone = None
        self.inner_milestone = None
        self.num_trajectories = -1
        self.receptor_indices = []
        self.ligand_indices = []
        self.extracted_directory = "extracted_from_b_surface"
        self.fhpd_directory = "first_hitting_point_distribution"
        return
    
class Milestone(serializer.Serializer):
    """
    Milestones represent the boundaries within the simulation.
    A given Anchor() object may contain any number of Milestone() 
    objects.
    
    Attributes:
    -----------
    index : int
        This is index of a Milestone() instance across the entire model.
        Two Milestone() instances in different Anchor() objects should
        share the same index if the milestones represent the same
        boundary between two adjacent Anchors.
        
    neighbor_anchor_index : int
        The index of the Anchor() instance on the other side of this
        milestone from the Anchor() instance that contains this
        Milestone() object.
        
    alias_index : int
        This field can be used to represent a simpler numbering scheme
        within an Anchor(). This is convenient for the low-level 
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
        
    def get_CV(self, model):
        """
        Obtain and return the collective variable object which
        applies to this milestone.
        """
        return model.collective_variables[self.cv_index]

class Model(serializer.Serializer): 
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
        self.browndye_settings = None
        self.k_on_info = None
        self.collective_variables = []
        self.anchors = []
        return
    
    def get_type(self):
        if self.calculation_type.lower() == "elber":
            return "elber"
        elif self.calculation_type.lower() == "mmvt":
            return "mmvt"
        else:
            error_msg = "Calculation type not available: "\
                + "{1}. Available types are 'elber' and 'mmvt'.".format(
                    self.calculation_type)
            raise Exception(error_msg)