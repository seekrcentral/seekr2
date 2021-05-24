"""
common_sim_namd.py

Base objects and routines for preparing NAMD simulations
to run.
"""

import os
import datetime

from parmed import unit

dt_date = datetime.datetime.now()

RESTART_CHECKPOINT_FILENAME = "backup_checkpoint"
NAMD_INPUT_FILENAME = "prod_{0}.namd"

SEEKR_ELBER_PROD_HEADER = "#" + "*"*64 + "\n" \
    + "#  Name SEEKR Elber Production\n" \
    + "#  Date " + dt_date.strftime("%Y-%b-%d") + "\n" \
    + "#  Generated with the SEEKR2 automatic NAMD input file generator\n" \
    + "#" + "*"*64 + "\n"

def add_string_buffer(string1, string2, gap=25, add_newline=True):
    """
    Add a buffer of a constant number of spaces between two strings.
    This is used to make prettier output for the NAMD input files.
    """
    if len(string1) >= gap:
        new_string = str(string1) + " " + str(string2)
    else:
        new_string = "{0: <{gap}}{1}".format(string1, string2, gap=25)
    
    if add_newline:
        return new_string + "\n"
    else:
        return new_string

class Amber_input_files():
    """
    Input files needed for a simulation using the Amber forcefield.
    
    Attributes:
    -----------
    parmfile : str
        The name of the AMBER prmtop or parm7 file to use as topology.
    ambercoor : str
        The name of the AMBER inpcrd or rst7 file to use, possibly, as
        coordinates.
    coordinates : str
        The name of a coordinates file, usually in PDB format.
    readExclusions : str or None, Default None
        Whether to read 1-4 exclusions from the prmtop file. It is
        recommended to put "yes". None gives yes by default. Setting
        to "yes" or None will cause exclusions and scaling to be read
        from the prmtop directly - the easiest option.
    scnb : str or None, Default None
        This is the VDW 1-4 scaling factor: the number that the 1-4 VDW
        interactions are divided by. None gives 2.0 by default.
    exclude : str or None, Default "scaled1-4
        This parameter describes to NAMD which pairs of bonded atoms
        should be excluded from non-bonded interactions. While the 
        conventional AMBER entry would be "scaled1-4", this should
        be done automatically if readExclusions is left at the default.
    _1_4_scaling : float
        This is the scaling factor for the 1-4 electrostatic 
        interactions. Unless exclude is set to "scaled1-4", this
        parameter has no effect.
    """
    def __init__(self):
        self.parmfile = ""
        self.ambercoor = ""
        self.coordinates = ""
        self.readExclusions = None
        self.scnb = None
        self.exclude = "scaled1-4"
        self._1_4scaling = 0.833333
        return
    
    def fill_out_from_model(self, anchor, in_building=True):
        """
        Use the SEEKR Model() object settings to fill out this object
        with relevant quantities.
        """
        assert anchor.amber_params is not None, "Amber settings not provided "\
            "for this anchor"
        if in_building:
            building_relative_dir = "../building"
            self.parmfile = os.path.join(building_relative_dir,
                                         anchor.amber_params.prmtop_filename)
            #self.ambercoor = os.path.join(building_relative_dir,
            #                              anchor.amber_params.inpcrd_filename)
            if anchor.amber_params.pdb_coordinates_filename:
                self.coordinates = os.path.join(
                    building_relative_dir, 
                    anchor.amber_params.pdb_coordinates_filename)
            
        else:
            self.parmfile = anchor.amber_params.prmtop_filename
            self.ambercoor = anchor.amber_params.inpcrd_filename
            self.coordinates = anchor.amber_params.pdb_coordinates_filename
        return
    
    def to_string(self):
        """
        Write these settings to a string to be used as direct input
        into a NAMD input script.
        """
        my_string = ""
        my_string += add_string_buffer("amber", "yes")
        my_string += add_string_buffer("parmfile", self.parmfile)
        if self.coordinates:
            my_string += add_string_buffer("coordinates", self.coordinates)
        else:
            my_string += add_string_buffer("ambercoor", self.ambercoor)
        
        if self.readExclusions is not None:
            my_string += add_string_buffer("readExclusions", self.readExclusions)
            
        my_string += add_string_buffer("exclude", self.exclude)
        if self.exclude == "scaled1-4":
            my_string += add_string_buffer("1-4scaling", 
                                           self._1_4scaling)
        if self.scnb is not None:
            my_string += add_string_buffer("scnb", str(self.scnb))
        
        
        return my_string

class Input_files():
    """
    All NAMD input settings related to input files.
    
    Attributes:
    -----------
    amber_input_files : Amber_input_files or None, Default None
        Settings for the Amber forcefield in NAMD.
    charmm_input_files : Charmm_input_files or None, Default None
        Settings for the Charmm forcefield in NAMD.
    binCoordinates : str, Default ""
        NAMD produces a binary output file that may be used as
        input for atomic coordinates.
    binVelocities : str, Default ""
        NAMD produces a binary output file that may be used as
        input for atomic velocities.
    extendedSystem : str, Default ""
        NAMD can produce and read an extended system file which
        contains information about box origin and vectors.
    try_to_load_state : bool, Default False
        Whether to try and load a saved state from an adjacent 
        anchor.
    """
    def __init__(self):
        self.amber_input_files = None
        self.charmm_input_files = None
        self.binCoordinates = ""
        self.binVelocities = ""
        self.extendedSystem = ""
        self.try_to_load_state = False
        return
    
    def fill_out_from_model(self, anchor, in_building=True):
        """
        Use the SEEKR Model() object settings to fill out this object
        with relevant quantities.
        """
        if anchor.amber_params is not None:
            self.amber_input_files = Amber_input_files()
            self.amber_input_files.fill_out_from_model(anchor, in_building)
            if not self.amber_input_files.coordinates:
                self.try_to_load_state = True
                
        elif anchor.charmm_params is not None:
            raise Exception("Charmm forcefields not yet available for NAMD "\
                            "SEEKR")
            #self.charmm_input_files = Charmm_input_files()
        else:
            raise Exception("No FF parameters provided for anchor in model.")
        
        self.binCoordinates = ""
        self.binVelocities = ""
        self.extendedSystem = ""
        return
    
    def to_string(self):
        """
        Write these settings to a string to be used as direct input
        into a NAMD input script.
        """
        my_string = "\n# Input Files\n"
        if self.amber_input_files is not None:
            my_string += self.amber_input_files.to_string()
        elif self.charmm_input_files is not None:
            raise Exception("CHARMM forcefield not yet implemented.")
        else:
            raise Exception("No input file forcefield selected.")
        
        if self.binCoordinates:
            my_string += add_string_buffer(
                "binCoordinates", self.binCoordinates)
            
        if self.binVelocities:
            my_string += add_string_buffer(
                "binVelocities", self.binVelocities)
        
        if self.extendedSystem:
            my_string += add_string_buffer(
                "extendedSystem", self.extendedSystem)
        return my_string
    
class Output_files():
    """
    NAMD settings related to file outputs.
    
    Attributes:
    -----------
    xstFile : str, Default ""
        The name of a file to write box origin and vectors. If empty
        then no file will be written.
    xstFreq : int
        If xstFile is set, then define the frequency to write to it.
    outputname : str
        A prefix for the binary .coor and .vel written at the end of a
        NAMD simulation.
    binaryoutput : str, Default ""
        Enable the use of binary output files. If left at the default
        empty string, this will be enabled automatically by NAMD.
    dcdfile : str
        The name of the file to write the trajectory to in DCD format.
    dcdfreq : int
        The frequency that the DCD file should be written to.
    restartname : str, Default ""
        The name of the file to write restart information to. The
        default empty string will cause NAMD to write to 
        $outputname.restart.
    restartfreq : int
        The frequency to write to the restart file.
    outputEnergies : int
        How frequently to write energy information to standard output
    outputTiming : int
        How frequently to compute and report timing benchmarks.
    """
    def __init__(self):
        self.xstFile = ""
        self.xstFreq = -1
        self.outputname = ""
        self.binaryoutput = ""
        self.dcdfile = ""
        self.dcdfreq = -1
        self.restartname = ""
        self.restartfreq = -1
        self.outputEnergies = -1
        self.outputTiming = -1
        return
    
    def fill_out_from_model(self, model, output_filename):
        """
        Use the SEEKR Model() object settings to fill out this object
        with relevant quantities.
        """
        self.outputname = output_filename
        self.dcdfile = output_filename + ".dcd"
        self.dcdfreq = model.calculation_settings.trajectory_reporter_interval
        if self.dcdfreq is None:
            self.dcdfreq = 0
        self.restartname = RESTART_CHECKPOINT_FILENAME
        self.restartfreq = model.calculation_settings\
            .restart_checkpoint_interval
        if self.restartfreq is None:
            self.restartfreq = 0
        self.outputEnergies \
            = model.calculation_settings.energy_reporter_interval
        if self.outputEnergies is None or self.outputEnergies == 0:
            self.outputEnergies = 1000000
        self.outputTiming = 10 * self.outputEnergies
        return
    
    def to_string(self):
        """
        Write these settings to a string to be used as direct input
        into a NAMD input script.
        """
        my_string = "\n# Output Files\n"
        if self.xstFile:
            my_string += add_string_buffer("xstFile", self.xstFile)
            assert self.xstFreq >= 0
            my_string += add_string_buffer("xstFreq", str(self.xstFreq))
        my_string += add_string_buffer("outputname", self.outputname)
        if self.binaryoutput:
            my_string += add_string_buffer("binaryoutput", self.binaryoutput)
        my_string += add_string_buffer("dcdfile", self.dcdfile)
        assert self.dcdfreq >= 0
        my_string += add_string_buffer("dcdfreq", str(self.dcdfreq))
        if self.restartname and self.restartfreq is not None:
            my_string += add_string_buffer("restartname", self.restartname)
            assert self.restartfreq >= 0
            my_string += add_string_buffer("restartfreq", str(self.restartfreq))
        assert self.outputEnergies > 0, "outputEnergies must be positive."
        my_string += add_string_buffer("outputEnergies", 
                                       str(self.outputEnergies))
        assert self.outputTiming >= 0
        my_string += add_string_buffer("outputTiming", 
                                       str(self.outputTiming))
        return my_string

class Simulation_parameters():
    """
    Physical and algorithmic settings for the NAMD simulation.
    
    Attributes:
    -----------
    temperature : float
        The temperature (in K) of the simulation. This will be used to 
        assign random velocities at the start.
    watermodel : str
        Which water model to use in the simulation. In NAMD, options
        include: 'tip3p', 'tip4p', and 'swm4'.
    cutoff : float
        The cutoff (in Angstroms) for the Van der Waals interactions.
    switching : str, Default "off"
        If set to 'off', then a truncated cutoff is performed. If 
        turned to 'on', then smoothing functions are applied.
    margin : float, Default 0.0
        This parameter could be used to slightly improve performance
        in NAMD.
    zeroMomentum : str, Default "on"
        If 'on', remove center of mass drift.
    ljCorrection : str, Default "on"
        If 'on', apply an analytic tail correction to reported Van
        der Waals energy and virial that is equal to the amount lost
        due to switching and cutoff of Lennard-Jones potential.
    rigitBonds : str
        Control if and how SHAKE is used to constrain bonds involving 
        hydrogen. Options include 'none', 'water', and 'all'. If left 
        empty, 'none' is used by default.
    rigidTolerance : float, Default 1e-8
        The SHAKE algorithm is assumed to have converged when all 
        constrained bonds differ from the nominal length by less than
        this amount.
    rigidIterations : int, Default 100
        The maximum number of iterations SHAKE will perform before 
        giving up on constraining bond lengths.
    useSettle : str, Default 'on'
        Whether to use SETTLE instead of SHAKE to constrain bonds.
    seed : int or None
        Sets the random seed.
    firsttimestep : int, Default 0
        The value of the first time step, modify to cause the
        simulation to start at a different timestep than zero.
    """
    def __init__(self):
        self.temperature = -1.0
        self.watermodel = ""
        self.cutoff = -1.0
        self.switching = "off"
        self.margin = 0.0
        self.zeroMomentum = "on"
        self.ljCorrection = "on"
        self.rigidBonds = ""
        self.rigidTolerance = 1e-8
        self.rigidIterations = 100
        self.useSettle = "on"
        self.seed = 0
        self.firsttimestep = 0
        return
    
    def fill_out_from_model(self, model):
        """
        Use the SEEKR Model() object settings to fill out this object
        with relevant quantities.
        """
        self.temperature = model.temperature
        self.watermodel = model.namd_settings.watermodel
        self.cutoff = model.namd_settings.nonbonded_cutoff
        constraints_str = model.namd_settings.constraints.lower()
        if constraints_str == "none" or constraints_str is None:
            self.rigidBonds = "none"
        elif constraints_str == "water":
            self.rigidBonds = "water"
        elif constraints_str == "hbonds":
            self.rigidBonds = "all"
        else:
            raise Exception("This option is not allowed for bond constraints "\
                            "in NAMD: %s" % constraints_str)
        self.rigidTolerance = model.namd_settings.langevin_integrator\
            .rigid_tolerance
        self.seed = model.namd_settings.langevin_integrator.random_seed
        return
    
    def to_string(self):
        """
        Write these settings to a string to be used as direct input
        into a NAMD input script.
        """
        my_string = "\n# Simulation Parameters\n"
        if self.temperature is not None:
            assert self.temperature >= 0.0
            my_string += add_string_buffer("temperature", str(self.temperature))
        if self.watermodel:
            my_string += add_string_buffer("watermodel", self.watermodel)
        assert self.cutoff >= 0.0
        cutoff_in_nm = unit.Quantity(self.cutoff, unit.nanometers)
        cutoff_in_A = cutoff_in_nm.value_in_unit(unit.angstroms)
        my_string += add_string_buffer("cutoff", str(cutoff_in_A))
        my_string += add_string_buffer("switching", self.switching)
        assert self.margin >= 0.0
        margin_in_nm = unit.Quantity(self.margin, unit.nanometers)
        margin_in_A = margin_in_nm.value_in_unit(unit.angstroms)
        my_string += add_string_buffer("margin", str(margin_in_A))
        my_string += add_string_buffer("zeroMomentum", self.zeroMomentum)
        my_string += add_string_buffer("ljCorrection", self.ljCorrection)
        if self.rigidBonds:
            my_string += add_string_buffer("rigidBonds", self.rigidBonds)
            assert self.rigidTolerance > 0.0
            my_string += add_string_buffer("rigidTolerance", 
                                           str(self.rigidTolerance))
            assert self.rigidIterations > 0
            my_string += add_string_buffer("rigidIterations", 
                                           str(self.rigidIterations))
            my_string += add_string_buffer("useSettle", self.useSettle)
        if self.seed is not None and self.seed != 0:
            my_string += add_string_buffer("seed", str(self.seed))
        my_string += add_string_buffer("firsttimestep", str(self.firsttimestep))
        return my_string

class Integrator_parameters():
    """
    Settings involving the integrator for a NAMD calculation.
    
    Parameters:
    -----------
    timestep : float
        The length of time (in fs) between timesteps.
    nonbondedFreq : int, Default 1
        The frequency to evaluate nonbonded forces
    fullElectFrequency : int, Default 1
        The frequency to evaluate electrostatic forces
    stepspercycle : int, Default 10
        The number of timesteps between atom reassignments to pair lists.
    langevinTemp : float
        The temperature (in K) to maintain by using the thermostat.
    langevinDamping : float, 
        The friction coefficient (in ps^-1) of the thermostat.
    langvinHydrogen : str, Default 'off'
        If set to 'off', then turn off Langevin dynamics for hydrogen 
        atoms.
    """
    def __init__(self):
        self.timestep = -1.0
        self.nonbondedFreq = 1
        self.fullElectFrequency = 1
        self.stepspercycle = 10
        self.langevinTemp = -1.0
        self.langevinDamping = -1.0
        self.langevinHydrogen = "off"
        return
    
    def fill_out_from_model(self, model):
        """
        Use the SEEKR Model() object settings to fill out this object
        with relevant quantities.
        """
        self.timestep = model.namd_settings.langevin_integrator.timestep
        self.langevinTemp = model.namd_settings.langevin_integrator\
            .target_temperature
        self.langevinDamping = model.namd_settings.langevin_integrator\
            .friction_coefficient
        return

    def to_string(self):
        """
        Write all settings in this object to a format which can be
        provided to NAMD in an input file.
        """
        my_string = "\n# Integrator Parameters\n"
        assert self.timestep > 0.0
        timestep_in_ps = unit.Quantity(self.timestep, unit.picoseconds)
        timestep_in_fs = timestep_in_ps.value_in_unit(unit.femtoseconds)
        my_string += add_string_buffer("timestep", str(timestep_in_fs))
        assert self.nonbondedFreq > 0
        my_string += add_string_buffer("nonbondedFreq", str(self.nonbondedFreq))
        assert self.fullElectFrequency > 0
        my_string += add_string_buffer("fullElectFrequency", 
                                       str(self.fullElectFrequency))
        assert self.stepspercycle > 0
        my_string += add_string_buffer("stepspercycle", 
                                       str(self.stepspercycle))
        if self.langevinTemp is not None:
            my_string += add_string_buffer("langevin", "yes")
            assert self.langevinTemp >= 0.0
            my_string += add_string_buffer("langevinTemp", 
                                           str(self.langevinTemp))
            assert self.langevinDamping >= 0.0
            my_string += add_string_buffer("langevinDamping", 
                                           str(self.langevinDamping))
            my_string += add_string_buffer("langevinHydrogen", 
                                           self.langevinHydrogen)
        return my_string

class Pressure_control():
    """
    Settings involving the pressure control for a NAMD calculation.
    
    Parameters:
    -----------
    useGroupPressure : str, Default "no"
        Whether group of atom virials are used to calculate pressure.
        In NAMD, "no" fluctuates less, and is required for rigid bonds
        so it is used by Default.
    useFlexibleCell : str, Default "no"
        Whether to allow the three orthogonal dimensions of the 
        periodic cell to fluctuate independently.
    useConstantArea : str, Default "no"
        When enabled, NAMD keeps the dimensions of the unit cell in the
        X-Y plane constant while allowing fluctuations along the
        z-axis. This is useful for membrane simulations.
    langevinPistonTarget : float
        The target pressure (in bar) to maintain by the barostat.
    langevinPistonPeriod : float, Default 0.2
        Describes barostat oscillation time scale (in ps) for the 
        piston method.
    langevinPistonDecay : float, Default 0.1
        Describes the barostat damping time scale (in ps) for the 
        Langevin piston method.
    langevinPistonTemp : float
        Specifies barostat noise temperature (in K) for the Langevin
        piston. This should be the same temperature as the thermostat.
    """
    def __init__(self):
        self.useGroupPressure = "no"
        self.useFlexibleCell = "no"
        self.useConstantArea = "no"
        self.langevinPistonTarget = -1.0
        self.langevinPistonPeriod = 0.2
        self.langevinPistonDecay = 0.1
        self.langevinPistonTemp = -1.0
        return
    
    def fill_out_from_model(self, model):
        """
        Use the SEEKR Model() object settings to fill out this object
        with relevant quantities.
        """
        self.langevinPistonTarget = model.namd_settings.barostat.target_pressure
        self.langevinPistonPeriod = model.namd_settings.barostat\
            .oscillation_frequency
        self.langevinPistonDecay = model.namd_settings.barostat.decay_timescale
        self.langevinPistonTemp = model.namd_settings.barostat\
            .target_temperature
        return
        
    def to_string(self):
        """
        Write all settings in this object to a format which can be
        provided to NAMD in an input file.
        """
        my_string = "\n# Pressure Control\n"
        my_string += add_string_buffer("useGroupPressure", 
                                       self.useGroupPressure)
        my_string += add_string_buffer("useFlexibleCell", self.useFlexibleCell)
        my_string += add_string_buffer("useConstantArea", self.useConstantArea)
        if self.langevinPistonTarget is not None:
            my_string += add_string_buffer("langevinPiston", "yes")
            assert self.langevinPistonTarget >= 0.0
            my_string += add_string_buffer("langevinPistonTarget", 
                                           str(self.langevinPistonTarget))
            assert self.langevinPistonPeriod >= 0.0
            langevinPistonPeriod_in_ps = unit.Quantity(
                self.langevinPistonPeriod, unit.picoseconds)
            langevinPistonPeriod_in_fs = \
                langevinPistonPeriod_in_ps.value_in_unit(unit.femtoseconds)
            my_string += add_string_buffer("langevinPistonPeriod", 
                                           str(langevinPistonPeriod_in_fs))
            assert self.langevinPistonDecay >= 0.0
            langevinPistonDecay_in_ps = unit.Quantity(
                self.langevinPistonDecay, unit.picoseconds)
            langevinPistonDecay_in_fs = langevinPistonDecay_in_ps.value_in_unit(
                unit.femtoseconds)
            my_string += add_string_buffer("langevinPistonDecay", 
                                           str(langevinPistonDecay_in_fs))
            assert self.langevinPistonTemp >= 0.0
            my_string += add_string_buffer("langevinPistonTemp", 
                                           str(self.langevinPistonTemp))
        return my_string
            
class Periodic_boundary_conditions():
    """
    Settings involving the periodic boundaries for a NAMD calculation.
    
    Parameters:
    -----------
    PMEGridSpacing : float, Default 0.1
        The spacing of the PME grid in NAMD.
    wrapWater : str, Default "on"
        Upon output to a DCD, this function will wrap water molecules
        within the box.
    wrapAll : str, Default "off"
        Upon output to a DCD, this function will wrap all contiguous
        clusters of bonded atoms.
    wrapNearest : str, Default "on"
        Wrap coordinates to the nearest image to the origin.
    cellBasisVector1 : str or None, Default None
        If cellBasisVector1 is a str, then this setting is the first 
        box vector for the system.
    cellBasisVector2 : str or None, Default None
        Analogous to cellBasisVector1. All three of 
        cellBasisVector1/2/3 must be all defined or all None.
    cellBasisVector3 : str or None, Default None
        Analogous to cellBasisVector1. All three of 
        cellBasisVector1/2/3 must be all defined or all None.
    cellOrigin : str or None, Default None
        The origin of the unit cell.
    """
    def __init__(self):
        self.PMEGridSpacing = 0.1
        self.wrapWater = "on"
        self.wrapAll = "off"
        self.wrapNearest = "on"
        self.cellBasisVector1 = None
        self.cellBasisVector2 = None
        self.cellBasisVector3 = None
        self.cellOrigin = None
        return
    
    def fill_out_from_model(self, model):
        """
        Use the SEEKR Model() object settings to fill out this object
        with relevant quantities.
        """
        self.PMEGridSpacing = model.namd_settings.PMEGridSpacing
        return
    
    def assign_cell_basis_vectors(self, box_vectors, origin_vector=None):
        """
        Assign cell basis vectors from a box_vectors Quantity object
        such that NAMD can accept it as input within a NAMD input
        script.
        
        Parameters:
        -----------
        box_vectors : Quantity
            A 3x3 array Quantity that represents the box_vectors. This
            would be output from an parmed.Rst7.box_vectors object.
        origin_vector : Quantity, Default None
            A 1x3 array Quantity that represents the origin of the
            periodic box.
        """
        box_vectors_unitless = box_vectors.to_quantity().value_in_unit(
            unit.angstrom)
        ax = "{:10.6f}".format(box_vectors_unitless[0][0])
        ay = "{:10.6f}".format(box_vectors_unitless[0][1])
        az = "{:10.6f}".format(box_vectors_unitless[0][2])
        bx = "{:10.6f}".format(box_vectors_unitless[1][0])
        by = "{:10.6f}".format(box_vectors_unitless[1][1])
        bz = "{:10.6f}".format(box_vectors_unitless[1][2])
        cx = "{:10.6f}".format(box_vectors_unitless[2][0])
        cy = "{:10.6f}".format(box_vectors_unitless[2][1])
        cz = "{:10.6f}".format(box_vectors_unitless[2][2])
        self.cellBasisVector1 = "{0} {1} {2}".format(ax, ay, az)
        self.cellBasisVector2 = "{0} {1} {2}".format(bx, by, bz)
        self.cellBasisVector3 = "{0} {1} {2}".format(cx, cy, cz)
        if origin_vector is not None:
            ox = "{:10.6f}".format(origin_vector[0])
            oy = "{:10.6f}".format(origin_vector[1])
            oz = "{:10.6f}".format(origin_vector[2])
            self.cellOrigin = "{0}, {1}, {2}".format(ox, oy, oz)
        return
    
    def to_string(self):
        """
        Write all settings in this object to a format which can be
        provided to NAMD in an input file.
        """
        my_string = "\n# Periodic Boundary Conditions (PBC)\n"
        if self.PMEGridSpacing is not None:
            my_string += add_string_buffer("PME", "yes")
            assert self.PMEGridSpacing > 0.0
            PMEGridSpacing_in_nm = unit.Quantity(
                self.PMEGridSpacing, unit.nanometers)
            PMEGridSpacing_in_A = PMEGridSpacing_in_nm.value_in_unit(
                unit.angstroms)
            my_string += add_string_buffer("PMEGridSpacing", 
                                           str(PMEGridSpacing_in_A))
        my_string += add_string_buffer("wrapWater", self.wrapWater)
        my_string += add_string_buffer("wrapAll", self.wrapAll)
        my_string += add_string_buffer("wrapNearest", self.wrapNearest)
        if self.cellBasisVector1 is not None:
            assert self.cellBasisVector2 is not None
            assert self.cellBasisVector3 is not None
            my_string += add_string_buffer("cellBasisVector1", 
                                           self.cellBasisVector1)
            my_string += add_string_buffer("cellBasisVector2", 
                                           self.cellBasisVector2)
            my_string += add_string_buffer("cellBasisVector3", 
                                           self.cellBasisVector3)
            if self.cellOrigin is not None:
                cellOrig_str = ','.join(map(str, self.cellOrigin))
                my_string += add_string_buffer("cellOrigin", cellOrig_str)
        else:
            assert self.cellBasisVector2 is None
            assert self.cellBasisVector3 is None
            
        return my_string
            
class Namd_root():
    """
    An entire set of NAMD input settings.
    
    Attributes:
    -----------
    input_files : Input_files
        The settings for input files for a NAMD simulation.
    output_files : Output_files
        The settings for output files for a NAMD simulation.
    simulation_parameters : Simulation_parameters
        The settings for simulation setup.
    integrator_parameters : Integrator_parameters
        The settings involving the integrator.
    pressure_control : Pressure_control
        Settings used for controlling the system pressure.
    periodic_boundary_conditions : Periodic_boundary_conditions
        Settings used for the periodic boundaries.
    colvars_config_filename : str
        The name of the colvars (Collective variables) script for NAMD.
    """
    def __init__(self):
        self.input_files = None
        self.output_files = None
        self.simulation_parameters = None
        self.integrator_parameters = None
        self.pressure_control = None
        self.periodic_boundary_conditions = None
        self.colvars_config_filename = ""
        return
    
    def fill_out_from_model(self, model, anchor, output_filename):
        """
        Use the SEEKR Model() object settings to fill out this object
        with relevant quantities.
        """
        self.input_files = Input_files()
        self.input_files.fill_out_from_model(anchor)
        self.output_files = Output_files()
        self.output_files.fill_out_from_model(model, output_filename)
        self.simulation_parameters = Simulation_parameters()
        self.simulation_parameters.fill_out_from_model(model)
        self.integrator_parameters = Integrator_parameters()
        self.integrator_parameters.fill_out_from_model(model)
        if model.namd_settings.barostat is not None:
            self.pressure_control = Pressure_control()
            self.pressure_control.fill_out_from_model(model)
        self.periodic_boundary_conditions = Periodic_boundary_conditions()
        self.periodic_boundary_conditions.fill_out_from_model(model)
        self.colvars_config_filename = "colvars.script"
        return
    
    def to_string(self, header=""):
        """
        Write all settings in this object to a format which can be
        provided to NAMD in an input file.
        
        Parameters:
        -----------
        header : str, Default ""
            An optional header string to include in the NAMD input 
            script.
            
        Returns:
        --------
        my_string : str
            The string representing a NAMD input script containing all
            generic, non-colvars, non-milestoning portions.
        """
        my_string = header
        assert self.input_files is not None
        my_string += self.input_files.to_string()
        assert self.output_files is not None
        my_string += self.output_files.to_string()
        assert self.simulation_parameters is not None
        my_string += self.simulation_parameters.to_string()
        assert self.integrator_parameters is not None
        my_string += self.integrator_parameters.to_string()
        if self.pressure_control is not None:
            my_string += self.pressure_control.to_string()
        assert self.periodic_boundary_conditions is not None
        my_string += self.periodic_boundary_conditions.to_string()
        if self.colvars_config_filename:
            my_string += "\n# Collective Variables (Colvars)\n"
            my_string += add_string_buffer("colvars", "on")
            my_string += add_string_buffer(
                "colvarsConfig", self.colvars_config_filename)
        
        return my_string

def write_extended_system_file(filename, box_vectors, origin_vector=None):
    """
    Make a NAMD extended system file for this system
    """
    step = "0"
    box_vectors_unitless = box_vectors.value_in_unit(unit.angstrom)
    ax = "{:.3f}".format(box_vectors_unitless[0][0])
    ay = "{:.3f}".format(box_vectors_unitless[0][1])
    az = "{:.3f}".format(box_vectors_unitless[0][2])
    bx = "{:.3f}".format(box_vectors_unitless[1][0])
    by = "{:.3f}".format(box_vectors_unitless[1][1])
    bz = "{:.3f}".format(box_vectors_unitless[1][2])
    cx = "{:.3f}".format(box_vectors_unitless[2][0])
    cy = "{:.3f}".format(box_vectors_unitless[2][1])
    cz = "{:.3f}".format(box_vectors_unitless[2][2])
    if origin_vector is None:
        origin_vector = unit.Quantity([0.0,0.0,0.0], unit.angstrom)
    origin_unitless = origin_vector.value_in_unit(unit.angstrom)
    ox = "{:.3f}".format(origin_unitless[0])
    oy = "{:.3f}".format(origin_unitless[1])
    oz = "{:.3f}".format(origin_unitless[2])
    box_list = [step, ax, ay, az, bx, by, bz, cx, cy, cz, ox, oy, oz]
    with open(filename, "w") as f:
        f.write("# NAMD extended system configuration output file\n")
        f.write("#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y "\
                "o_z\n")
        f.write(" ".join(box_list))
        
    return

class Sim_namd():
    """
    Contain all the information necessary to run an MMVT SEEKR 
    calculation in OpenMM.
    
    Attributes:
    -----------
    namd_root : Namd_root()
        All non-SEEKR settings for NAMD input script.
    colvars_config : Colvars_config()
        All settings for the external Colvars script in NAMD.
    seekr_namd_settings : Seekr_namd_settings()
        All SEEKR-relevant settings for NAMD input script.
    header : str
        An informative string to put at the top of the input file.
    """
    def __init__(self):
        self.namd_root = Namd_root()
        self.colvars_config = None #Colvars_config()
        self.seekr_namd_settings = None #Seekr_namd_settings()
        self.header = "## MISSING HEADER ##"
        return
    
    def write_namd_script(self, model, anchor, restart_index, filename=None, 
                          base_filename=None):
        """
        Write out a NAMD script to a file.
        
        Parameters:
        -----------
        model : Model()
            The Model object for this system.
        anchor : Anchor()
            The Anchor object in whose directory to write the NAMD script.
        filename : str or None, Default None
            The absolute file path to write this NAMD script to. If left
            at None, a generic filename will be used.
        base_filename : str or None, Default None
            The name of the file (without the path) which will be the
            NAMD script. If base_filename is not None and filename is 
            None, the path will be added based on the anchor's 
            directory.
        """
        assert self.seekr_namd_settings.eval_stride\
            % self.namd_root.integrator_parameters.stepspercycle == 0, \
            "Eval stride must be a multiple of stepspercycle. eval_stride = "\
            "%s, stepspercycle = %s" % (
                self.seekr_namd_settings.eval_stride, 
                self.namd_root.integrator_parameters.stepspercycle)
              
        if filename is None:
            if base_filename is None:
                base_filename = NAMD_INPUT_FILENAME.format(restart_index)
            filename = os.path.join(model.anchor_rootdir , anchor.directory,
                                    anchor.production_directory, base_filename)
            
        my_string =  self.namd_root.to_string(header=self.header)
        my_string += self.seekr_namd_settings.to_string(
            model, anchor, self.namd_root.simulation_parameters.firsttimestep)
        with open(filename, "w") as f:
            f.write(my_string)
        if base_filename is None:
            base_filename = filename
        return base_filename
    
    def write_colvar_script(self, model, anchor, filename=None):
        """
        Write out a NAMD colvar script to a file.
        
        Parameters:
        -----------
        model : Model()
            The Model object for this system.
        anchor : Anchor()
            The Anchor object in whose directory to write the NAMD script.
        filename : str or None, Default None
            The absolute file path to write this NAMD colvar script to.
            If left at None, a generic filename will be used.
        
        """
        if filename is None:
            base_filename = self.namd_root.colvars_config_filename
            filename = os.path.join(model.anchor_rootdir , anchor.directory,
                                    anchor.production_directory, base_filename)
        
        colvar_script_string = self.colvars_config.to_string()
        with open(filename, "w") as f:
            f.write(colvar_script_string)
        
        return
        