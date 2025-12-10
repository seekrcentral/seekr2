"""
common_sim_sda

Base objects and routines for preparing and running SDA 
simulations.
"""

import xml.etree.ElementTree as ET
import os

import parmed
import numpy as np

from abserdes import Serializer

SDA_TRAJ_NAME = "traj"
SDA_COMPLEXES_NAME = "complexes"
SDA_INPUT_FILENAME = "sda.in"
APBS_INPUT_FILENAME = "apbs_input.xml"
REACTION_FILENAME = "p2.rxna"
SDA_RECEPTOR = "p1"
SDA_LIGAND = "p2"


class Ion():
    """
    Represent an Ion object in SDA for generating APBS grids and setting 
    ionic strength.
    
    Attributes:
    -----------
    radius : float
        The radius of the ion in units of Angstroms.
    charge : float
        The charge of the ion in units of e.
    conc : float
        The concentration of the ion in solution in units of moles per
        liter.
    """
    
    def __init__(self):
        self.radius = -1.0
        self.charge = -1.0
        self.conc = -1.0
        return
    
    def serialize(self, xmlIon):
        """
        Serialize this object to XML
        
        Parameters:
        -----------
        xmlIon : ElementTree.SubElement
            All sub elements get added to this root.
        """
        
        assert self.radius >= 0.0, "Ion radius must be set"
        xmlIonRadius = ET.SubElement(xmlIon, 'radius')
        xmlIonRadius.text = str(self.radius)
        xmlIonCharge = ET.SubElement(xmlIon, 'charge')
        xmlIonCharge.text = str(self.charge)
        assert self.conc >= 0.0, "Ion concentration must be set"
        xmlIonConc = ET.SubElement(xmlIon, 'conc')
        xmlIonConc.text = str(self.conc)
        return

class Solvent():
    """
    Parameters to represent the solvent within the BD simulation.
    
    Attributes:
    -----------
    debye_length : float
        The Debye length is a distance inversely related to the 
        strength and concentration of ions in solution.
    dielectric : float, Default 78.0
        The dielectric of solvent, relative to vacuum permittivity.
    relative_viscosity : float, Default 1.0
        Relative to water viscosity.
    kT : float
        Thermal energy relative to Boltzmann's constant times 298 K.
    desolvation_parameter : float, Default 1.0
        Factor that multiplies desolvation energy.
    ions : list
        A list of Ion() objects for APBS input
    """
    
    def __init__(self):
        self.dielectric = 78.0
        self.temperature = 298.15
        self.relative_viscosity = 1.0
        self.kT = -1.0
        self.desolvation_parameter = 1.0
        self.ions = []
        return
    
    def serialize(self, xmlSolvent, make_apbs_mode=True):
        """
        Serialize this object to XML
        
        Parameters:
        -----------
        xmlSolvent : ElementTree.SubElement
            All sub elements get added to this root.
        make_apbs_mode : bool
            Whether this object should be serialized for the
            make_apbs_inputs program.
        """
        
        if not make_apbs_mode:
            assert self.debye_length > 0.0, "Solvent Debye length must be assigned."
            xmlSolventDebye = ET.SubElement(xmlSolvent, 'debye_length')
            xmlSolventDebye.text = str(self.debye_length)
        assert self.dielectric > 0.0
        xmlSolventDielectric = ET.SubElement(xmlSolvent, 'dielectric')
        xmlSolventDielectric.text = str(self.dielectric)
        assert self.relative_viscosity > 0.0
        xmlSolventRelVisc = ET.SubElement(xmlSolvent, 'relative_viscosity')
        xmlSolventRelVisc.text = str(self.relative_viscosity)
        assert self.kT > 0.0
        xmlSolventKT = ET.SubElement(xmlSolvent, 'kT')
        xmlSolventKT.text = str(self.kT)
        assert self.desolvation_parameter >= 0.0
        xmlSolventDesolv = ET.SubElement(xmlSolvent, 'desolvation_parameter')
        xmlSolventDesolv.text = str(self.desolvation_parameter)
        if make_apbs_mode:
            xmlSolventIons = ET.SubElement(xmlSolvent, 'ions')
            for ion in self.ions:
                xmlIon = ET.SubElement(xmlSolventIons, 'ion')
                xmlIon.text = ion.serialize(xmlIon)
        
        return

class MainSDA():
    """
    Main SDA parameters (they don't belong to any group)
    
    Attributes:
    -----------
    dseed: integer, Default 0
        Seed number to initialize random forces. If 0, time clock is used.
    nrun: integer, Default 1
        Number of trajectories.
    timemax: float, Default 0
        Maximum time of trajectories in ps.
    probeb: float, Default 1.77
        Radius of the probe used to compute the exclusion volume grid
    probew: float, Default 1.4
        Radius of solvent probe used to calculate the solvent accessibilities 
        of solute atoms
    hexcl: float, default 0.1
        Exclusion grid spacions
    threshold: float, default 0.0
        Solvent accessible surface area (SASA) threshold to be used.
    epfct: float, default 0.5
        Factor to multiply to electrostatic potentials.
    edfct: float, default None
        Factor to multiply to electrostatid desolvation potentials.
    hdfct: float, default None
        Factor to multiply to hydrophobic desolvation potentials.
    restart: string
        If this string is not empty, BD simulations will be restarted.
        It is necessary to provide a complexes file.
    save_access: integer, Default 1
        Saves the SASA values calculated. It saves computational time if solutes are big.
    rboost: float, Default 1.0
        Boost distance (in Å) to move the solutes apart when more than novers moves 
        were not accepted due to overlaps.
    novers: integer, Default 150
        Number of overlaps allowed before making a boost. 
        Recommended value: 100-200, to avoid boosts.
    """

    def __init__(self):
        self.dseed = 0
        self.nrun = -1
        self.timemax = 0
        self.probep = 1.4
        self.probew = 1.2
        self.hexcl = 0.1
        self.threshold = 0.0
        self.epfct = 0.5
        self.edfct = 0.36
        self.hdfct = -0.013
        self.restart = None
        self.save_access = 1
        self.rboost = 1.0
        self.novers = 150

    def make_group(self):
        self.group = {}
        for attribute_name, value in self.__dict__.items():
            if attribute_name != "group" and value is not None:
                string = attribute_name + " = " + str(value) + "\n"
                self.group[attribute_name] = string


class Type_Calculation():
    """
    Calculation type to perform in BD simulations
    
    Attributes:
    -----------
    type : string (sda_2proteins, sdamm, sda_koff, sda_energy), Default None
        Indicates the type of simulation to run
    total solutes : integer, Default 2
        Total solutes to use.
    total grids : integer, Default 2
        Total grid types to use. Conformational ensamble counts as 1 type of grid.
    """

    def __init__(self):
        self.type = "sda_2proteins"
        self.total_solutes = 2
        self.total_grids = 2
        return
    
    def make_group(self):
        self.group = {}
        self.group["GROUP"] = "GROUP = Type_Calculation\n"
        for attribute_name, value in self.__dict__.items():
            if attribute_name != "group" and value is not None:
                string = "    " + attribute_name + " = " \
                    + str(value) + "\n"
                self.group[attribute_name] = string
        self.group["END_GROUP"] = "END GROUP\n"
        return

class ReactionCriteria():
    """
    Parameters for reaction criteria.
    
    Attributes:
    -----------
    computation: string
        Defines the type of calculations to perform during
    BD simulations. Options are:
            - association: calculates kon from trajectories that satisfy reaction criteria
            - docking: saves the encounter complexes that satisfy reaction criteria.
            - electron_transfer: computes electron transfer
            - all: save all the encounter complexes, reaction criteria is omited.
            - off_rc: disactivates reaction criteria.
    rxna12f: string
        Pathway and name of the reaction criteria file.
    et_sol1, et_sol2: string or None
        Pathways and names of the reaction criteria file 
        for the electron transfer calculations.
    dind: float, Default 0.
        Distance between contacts to be considered independent.
    nnnons: integer, Default 1
        Maximum number of independent contacts required to record encounter complexes.
    nwrec: integer, Default 1
        Windows distance number where start recording encounter complexes.
        The distance is defined as (win0+dwin*(nwrec-1))
    sdamd: integer, Default 1
        Saves the encounter complexes that are found in the nwrec windows distance.
    """    

    def __init__(self):
        self.computation = "association"
        self.rxna12f = REACTION_FILENAME
        self.et_sol1 = None
        self.et_sol2 = None
        self.dind = 0
        self.nnnons = 1
        self.nwrec = 1
        self.sdamd = 1
        return
    
    def make_group(self):
        self.group = {}
        self.group["GROUP"] = "GROUP = ReactionCriteria\n"
        for attribute_name, value in self.__dict__.items():
            if attribute_name != "group" and value is not None:
                string = "    " + attribute_name + " = " \
                    + str(value) + "\n"
                self.group[attribute_name] = string
        self.group["END_GROUP"] = "END GROUP\n"
        return
    

class RateCalculation():
    """
    Parameters for rate calculations.
    
    Attributes:
    -----------
    win0 : float, Default None
        Closest windows distance reaction
    nwin: integer, Default None
        Number of windows distances to monitore association rates
    dwin: float, Default 0.5
        Separation between windows distances
    nb_contact: integer, Default 1
        Number of independent contacts to be considered for rate calculations
    bootstrap: integer, Default 1
        If 1, association rates and first passage time are printed for every
        trajectory to perform later bootstrap calculations.
    fpt: integer, Default 0
        If 1, first passage time is printed.
    stop_traj: integer, Default 0
        If 1, the trajectory will stop when the most strict reaction criteria
        definition is met.
    analytical_correction: integer, Default 0
        If 1, the analytical correction for centrosymmetric forces is computed.
    """

    def __init__(self):
        self.win0 = 4.0
        self.nwin = 35
        self.dwin = 0.5
        self.nb_contact = 1
        self.bootstrap = 1
        self.fpt = 0
        self.stop_traj = 1
        self.analytical_correction = 0
        return

    def make_group(self):
        self.group = {}
        self.group["GROUP"] = "GROUP = RateCalculation\n"
        for attribute_name, value in self.__dict__.items():
            if attribute_name != "group" and value is not None:
                string = "    " + attribute_name + " = " \
                    + str(value) + "\n"
                self.group[attribute_name] = string
        self.group["END_GROUP"] = "END GROUP\n"
        return

class Analytical():
    """
    Parameters for rate calculations.
    
    Attributes:
    -----------
    h_analytic : float, Default 0.01
        defines the bin size (in Å) of the arrays of the precomputed 
        debye-hueckel values.
    debyeh: integer, Default 0
        if 1, activates Debye-Hueckel interactions for all solutes.
    dh_imgchg: integer, Default 0
        if 1 and a surface is present, the image-charge model used 
        for modelling electrostatic interactions with conducting surfaces 
        is extended to use debye-huckel interactions for all solutes whose 
        images lie outside of their electrostatic grids.
    ionic_strength: float or None, Default None
        Ionic strength (in M) for Debye-Hueckel interactions.
    hydrodynamic: integer, Default 0
        If 1, activates the mean-field hydrodynamics interactions.
    lvol_rcut: float or None, Default None
        Defines the local volume radius for hydrodynamics. lvol_rcut must
        be greater than vol_radius of any solute.
    surface_heigth: float or None, Default None
        Defines the heigth (in A) above a surface at which the image method 
        is applied for surface hydrodynamic interactions.
    hom_charged_surf: integer, Default 0
        If 1, activates the long-range Debye-Hueckel electrostatic treatment
        for a surface
    surface_fact: float or None, Default None
        Prefactor for the surface charge potential. If the layer is thin 
        (i.e.: one atom layer), surface_fact = 1. If the layer is thin,
        surface_fact = 0.5.
    gouy_chapman: integer, Default 0
        If 1, uses Gouy-Chapman treatment instead of Debye-Huckel for the long range 
        electrostatics of the surface.
    """

    def __init__(self):
        self.h_analytic = 0.01
        self.debyeh = 0
        self.dh_imgchg = 0
        self.ionic_strength = None
        self.hydrodynamic = 0
        self.lvol_rcut = None
        self.surface_heigth = None
        self.hom_charged_surf = 0
        self.surface_fact = None
        self.gouy_chapman = 0
        return
    
    def make_group(self):
        self.group = {}
        self.group["GROUP"] = "GROUP = Analytical\n"
        for attribute_name, value in self.__dict__.items():
            if attribute_name != "group" and value is not None:
                string = "    " + attribute_name + " = " \
                    + str(value) + "\n"
                self.group[attribute_name] = string
        self.group["END_GROUP"] = "END GROUP\n"
        return

class Solute(Serializer):

    """
    A solute for input to SDA
    
    Attributes:
    ----------

    type : string
        Indicates the type of molecule to be simulated (protein or \
        small organic compound).
    apbs_grid_dime : integer
        Size of the APBS grid dimensions.
    apbs_grid_spacing : float
        Space between the points conforming the grid.
    dielectric : float
        Dielectric constant of the solute.
    solute_grid : Solute_Grid()
        Solute grid object whose attributes will be used to build \
        an SDA input file.
    """

    def __init__(self):
        self.type = ""
        self.apbs_grid_dime = -1
        self.apbs_grid_spacing = -1
        self.dielectric = -1
        self.solute_grid = Solute_Grid()

class Solute_Grid(Serializer):

    """
    Parameters for rate calculations.
    
    Attributes:
    -----------
    nb_solute : integer, Default 1
        Number of solutes of this type.
    pdb_filename : string
        Pathway where the pdb file is located.
    diffusion_trans : float or None. Default None.
        Translational diffusion coefficient value in A^2/ps
    diffusion_rotat : float or None. Default None.
        Rotational diffusion coefficient value in rad/ps
    real_net_charge: float or None. Default None.
        Total net formal charge of the solute. 
        Used if the Debye field is activated.
    surface_charge_dens: float or None. Default None.
        If the solute is a surface, this option can be used to define
        an average charge density for surface. The interaction with other solutes 
        is then modelled using a Debye-Hueckel sphere approximation where 
        spherical charge solutes interact with a homogeneously charged surface.
    rotate : integer, default 1.
        Defines whether the solute can rotate or not. If 0,
        rotational diffusion coefficient is not used.
    surface: str or None, Default None.
        Defines if the solute is a surface.
    flex: integer or None, Default None
        Defines if a conformational ensamble is read for the solute.
    image_charge: integer or None, Default None
        Defines if the electrostatic image charge of this solute(s) must be computed. 
        It is relevant in the case that one solute is a metal surface.
    dh_radius: float or None, Default None
        Indicates the radius used to represent the solute cavity in the Debye-Hueckel 
        sphere charge model.
    vol_radius: float or None, Default None
        Indicates the radius used for calculating the local occupied volume fractions 
        of the mean-field hydrodynamic interaction model.
    list_conformation: string or None, Default None
        Plain text file that indicates the pathways and names of the different grids
        and population of each conformer of the solute.
    total_conf: integer or None, Default None
        total number of conformations of the solute
    method: string (random, minimum, metropolis) or None, Default None.
        method to be used for changing the conformations on the fly. Options:
            - random: choose a random conformation if nearest > 1. Energies are
            not evaluated.
            - minimum: the minimum energy between actual and the closest conformations
            is selected if nearest = 1
            - metropolis: Metropolis-Hastings algorithm is applied to choose a new
            conformation based on the population and the intermolecular energy with other
            solutes.
    initial_conf: integer, Default -1.
        selects the conformation at the beginning of each trajectory. If value is -1, the
        conformation will be selected randomly.
    frequency: float, Default 100.0
        Delay time between changes in conformation (in ps).
        If SDAMM, the delay is computed as a Gaussian random with average frequency
        and standard deviation std_frequency to avoid syncronisation between solutes.
    std_frequency: float, Default 0.0
        Standard deviation frequency (in ps).
    epf : string or None, Default None
        Pathway and name of the electrostatic potential grid file
    qef : string or None, Default None
        Pathway and name of the effective charges file
    edf : string or None, Default None
        Pathway and name of the electrostatic desolvation grid file
    hdf : string or None, Default None
        Pathway and name of the hydrophobic desolvation grid file
    lj_repf : string or None, Default None
        Pathway and name of the lennard-jones repulsion grid file
    """

    def __init__(self):
        self.nb_solute = 1
        self.pdb_filename = ""
        self.trans_diffusion = None
        self.rot_diffusion = None
        self.real_net_charge = None
        self.surface_charge_dens = None
        self.rotate = 1
        self.surface = 0
        self.flex = None
        self.image_charge = None
        self.dh_radius = None
        self.vol_radius = None
        self.list_conformation = None
        self.total_conf = None
        self.method = None
        self.initial_conf = None
        self.frequency = None
        self.std_frequency = None
        self.epf = ""
        self.qef = ""
        self.edf = ""
        self.hdf = ""
        self.lj_repf = None
        return
    
    def make_group(self):
        self.group = {}
        self.group["GROUP"] = "GROUP = Solute_Grid\n"
        for attribute_name, value in self.__dict__.items():
            if attribute_name != "group" and value is not None:
                string = "    " + attribute_name + " = " \
                    + str(value) + "\n"
                self.group[attribute_name] = string
        self.group["END_GROUP"] = "END GROUP\n"
        return


class Geometry():
    """
    Parameters to represent the geometry properties in where simulations are going to run.
    
    Attributes:
    -----------
    type : string (sphere, box, nambox), Default sphere
        Indicates whether a simulation runs in a box or in a sphere. If nambox, association rates
        are calculated using the NAM algorithm in a box.
    pbc : integer, Default 0
        If 1, periodic boundary conditions are applied when using a box.
    escape : integer, Default 1
        Indicates if a solute can escape. If the solute escapes, the trajectory will stop.
    surface : integer, Default 0
        Indicates if a surface is present.
    start_pos : float, default None
        Center-to-center distance (in amstrongs) of solutes at the beginning of a trajectory. 
        If the simulation runs in a box, this distance will be considered on the z-ax.
        If the simulation runs in a sphere or nambox, this distance is the b-surface.
    c : float, default None
        c-surface (in amstrongs) in case that simulations run in a sphere or nambox.
    <AX>min, <AX>max: float or None, default None
        <AX> = x, y, z. Min and max indicate the size of the ax.
        If escape is activated, the trajectory ends when the solute reaches zmax.
    """

    def __init__(self):
        self.type = "sphere"
        self.pbc = None
        self.surface = 0
        self.escape = 1
        self.start_pos = -1
        self.c = -1
        self.xmin = None
        self.xmax = None
        self.ymin = None
        self.ymax = None
        self.zmin = None
        self.zmax = None
        return

    def make_group(self):
        self.group = {}
        self.group["GROUP"] = "GROUP = Geometry\n"
        for attribute_name, value in self.__dict__.items():
            if attribute_name != "group" and value is not None:
                string = "    " + attribute_name + " = " \
                    + str(value) + "\n"
                self.group[attribute_name] = string
        self.group["END_GROUP"] = "END GROUP\n"
        return


class Timestep():
    """
    Parameters to represent the timestep configuration.
    
    Attributes:
    -----------
    variable : integer, Default 1
        If 1, timestep is variable, describing a linear change between 
        dt1 and dt2 when solutes are within swd1 and swd2 distance.
        If 0, just dt1 is considered.
    dt1 : float or None. Default None.
        Smallest timestep. Employed when solutes are closer than swd1.
    swd1 : float or None. Default None.
        Center-to-center distance at which timestep = dt1.
    dt2 : float or None. Default None.
        Longest timestep. Employed when solutes are beyond swd2
    swd2 : float or None. Default None.
        Center-to-center distance at which timestep = dt2.
    """

    def __init__(self):
        self.variable = 1
        self.dt1 = 1
        self.swd1 = 60
        self.dt2 = 20
        self.swd2 = 100
        return
    
    def make_group(self):
        self.group = {}
        self.group["GROUP"] = "GROUP = Timestep\n"
        for attribute_name, value in self.__dict__.items():
            if attribute_name != "group" and value is not None:
                string = "    " + attribute_name + " = " \
                    + str(value) + "\n"
                self.group[attribute_name] = string
        self.group["END_GROUP"] = "END GROUP\n"
        return
    

class Complexes():
    """
    Parameters to represent the timestep configuration.
    
    Attributes:
    -----------
    fcomplexes : string
        Pathway and name of the file where encounter complexes will be written.
    restart_complex : integer, Default 0.
        If 1, simulations will be restarted from the complexes file indicated
        in "restart" flag.
    binary_complex : integer, Default 0.
        If 1, complexes file will be in binary format.
    nb_complexes : integer, Default 500
        Number of encounter complexes to record. If this number is reached,
        the new encounter complexes with lower energy will replace the old
        encounter complexes with higher energy.
    rmsd_min : float, Default 1.0
        RMSD threshold (in A) between encounter complexes to consider them as
        the same encounter complex.
    ftrajectories : string
        Pathway and name of the trajectories file
    binary_trajectory : integer, Default 0.
        If 1, trajectories file will be in binary format.
    ntraj_rec : integer, Default -1
        Number of the trajectory to record. If -1, all trajectories
        are saved.
    freq_print : integer, Default 100
        Frequency (in steps) with which snapshots of the trajectory are saved.
    """

    def __init__(self):
        self.fcomplexes = ""
        self.restart_complex = 0
        self.binary_complex = 0
        self.nb_complexes = 500
        self.rmsd_min = 1.0
        self.ftrajectories = ""
        self.binary_trajectory = 0
        self.ntraj_rec = 1
        self.freq_print = 1000
        return
    
    def make_group(self):
        self.group = {}
        self.group["GROUP"] = "GROUP = Complexes\n"
        for attribute_name, value in self.__dict__.items():
            if attribute_name != "group" and value is not None:
                string = "    " + attribute_name + " = " \
                    + str(value) + "\n"
                self.group[attribute_name] = string
        self.group["END_GROUP"] = "END GROUP\n"
        return

class Input():
    """
    An SDA Input file object.
    
    Attributes:
    -----------
    nrun : int or None
        Number of trajectories.
    dseed : int
        Used to seed the random number generator.
    output : str
        Name of file used for outputting SDA results.
    timemax : float or None
        Maximum time per trajectory. If 0, trajectory ends when solute escapes.
    prot_data_grid_path : string
        Path of the protein's grid files.
    lig_data_grid_path : string
        Path of the ligand's grid files.
    trajectory_file : str or None
        Name of file used for outputting of trajectories.
    complexes_file: str or None
        Name of file used for outputting of encounter complexes.
    nb_complexes: int or None
        Number of encounter complexes to record
    ntraj: int or None
        Number of the trajectory to record. If -1, all are recorded.
    geometry: str or None
        Geometry used in the simulation (sphere, box, nambox)
    start_pos: float or None
        Position of the ligand at the beginning of a trajectory
    c: float or None
        c-surface
    <AX>min, <AX>max: float or None
        Initial and final position of box axes, being AX = x, y, z
    surface: int or None
        1 if a surface is present
    win0: float
        Closest windows distance
    nwin: integer
        Number of windows distances to monitor
    hydrodynamics : int or None, Default None
        If 1, mean-field hydrodynamics approximation will be used.
    """
    
    def __init__(self):
        self.nrun = -1
        self.timemax = 0
        self.dseed = 0
        self.type = "sda_2proteins"
        self.solutes = []
        self.solutes_path = ""
        self.total_solutes = 2
        self.total_grids = 2
        self.rxna12f = REACTION_FILENAME
        self.ftrajectories = SDA_TRAJ_NAME
        self.fcomplexes = SDA_COMPLEXES_NAME
        self.nb_complexes = 10000
        self.ntraj = 1
        self.geom_type = "sphere"
        self.start_pos = 100
        self.c = 200
        self.pbc = 0
        self.xmin = None
        self.xmax = None
        self.ymin = None
        self.ymax = None
        self.zmin = None
        self.zmax = None
        self.surface = 0
        self.rotate = 1
        self.win0 = 3.0
        self.nwin = 35
        self.dt_variable = 1
        self.dt1 = 1
        self.hydrodynamics = None
        return
    
    def make_input(self, make_analytical=False):
        """
        Creates the SDA input file.
        
        Parameters:
        -----------
        make_analytical : bool, Default False
            If True, the analytical group will be written into the 
            SDA input file.
        """
        
        self.main = MainSDA()
        self.main.nrun = self.nrun
        self.main.timemax = self.timemax
        self.main.make_group()

        self.type_calculation = Type_Calculation()
        self.type_calculation.type = self.type
        self.type_calculation.total_solutes = self.total_solutes
        self.type_calculation.total_grids = self.total_grids
        self.type_calculation.make_group()

        self.reaction_criteria = ReactionCriteria()
        self.reaction_criteria.make_group()

        self.rate_calculation = RateCalculation()
        self.rate_calculation.win0 = self.win0
        self.rate_calculation.nwin = self.nwin
        self.rate_calculation.make_group()

        self.solute_grids = []
        for solute in self.solutes:
            pqr_filename = os.path.basename(solute.pqr_filename)
            solute_name = '.'.join(pqr_filename.split(".")[:-1])
            pdb_filename = solute_name + "_noh.pdb"
            self.solute_grid = Solute_Grid()
            self.solute_grid.nb_solute = solute.solute_grid.nb_solute
            self.solute_grid.trans_diffusion = solute.solute_grid.trans_diffusion
            self.solute_grid.rot_diffusion = solute.solute_grid.rot_diffusion
            self.solute_grid.pdb_filename = os.path.join(self.solutes_path, \
                                                pdb_filename)
            if hasattr(solute, 'rotate'):
                self.solute_grid.rotate = solute.solute_grid.rotate
            if hasattr(solute, 'surface') and solute.solute_grid.surface == "yes":
                self.solute_grid.surface = 1
            self.solute_grid.epf = os.path.join(self.solutes_path, \
                                                solute_name + "_ep.grd")
            self.solute_grid.qef = os.path.join(self.solutes_path, \
                                                solute_name + ".echa")
            self.solute_grid.edf = os.path.join(self.solutes_path, \
                                                solute_name + "_ed.grd")
            self.solute_grid.hdf = os.path.join(self.solutes_path, \
                                                solute_name + "_hd.grd")
            if hasattr(solute, 'lj_repf'):
                self.solute_grid.lj_repf = os.path.join(self.solutes_path, \
                                                solute_name + "_ljrep.grd")
            self.solute_grid.make_group()
            self.solute_grids.append(self.solute_grid)

        if make_analytical:
            ionic_strength = 0
            for ion in self.solvent.ions:
                ionic_strength += ion.conc*ion.charge**2 / 2
            self.analytical = Analytical()
            self.analytical.make_group()
        else:
            self.analytical = None

        self.geometry = Geometry()
        self.geometry.type = self.geom_type
        self.geometry.start_pos = self.start_pos
        self.geometry.pbc = self.pbc
        if self.geometry.type == "sphere" or self.geometry.type == "nambox":
            self.geometry.c = self.c
        if self.geometry.type == "box" or self.geometry.type == "nambox":
            self.geometry.xmin = self.xmin
            self.geometry.xmax = self.xmax
            self.geometry.ymin = self.ymin
            self.geometry.ymax = self.ymax
            self.geometry.zmin = self.zmin
            self.geometry.zmax = self.zmax
        self.geometry.surface = self.surface   
        self.geometry.make_group()
        
        self.timestep = Timestep()
        self.timestep.variable = self.dt_variable
        self.timestep.dt1 = self.dt1
        if self.dt_variable is None:
            self.timestep.swd1 = None
            self.timestep.dt2 = None
            self.timestep.swd2 = None
        self.timestep.make_group()

        self.complexes = Complexes()
        self.complexes.fcomplexes = self.fcomplexes
        self.complexes.ftrajectories = self.ftrajectories
        self.complexes.make_group()

        return self
    
    def write(self, filename, make_apbs_mode=True):
        """
        Write the SDA input to plain text.
        
        Parameters:
        -----------
        filename : str
            The name of the file to write XML to.
        
        make_apbs_mode : bool
            Whether this object should be serialized for the
            make_apbs_inputs program.
        """
        
        root = self.make_input()
        with open(filename, 'w') as f:

            for parameter in root.type_calculation.group:            
                line = root.type_calculation.group[parameter]
                f.write(line)
            f.write("\n")

            for parameter in root.main.group:
                if parameter != "novers" and parameter != "rboost":            
                    line = root.main.group[parameter]
                    f.write(line)
            f.write("\n")

            for parameter in root.reaction_criteria.group:
                line = root.reaction_criteria.group[parameter]
                f.write(line)
            f.write("\n")

            for parameter in root.rate_calculation.group:
                line = root.rate_calculation.group[parameter]
                f.write(line)
            f.write("\n")

            if root.analytical is not None:
                for parameter in root.analytical.group:
                    line = root.analytical.group[parameter]
                    f.write(line)
                f.write("\n")

            for solute_grid in root.solute_grids:
                for parameter in solute_grid.group:
                    line = solute_grid.group[parameter]
                    f.write(line)
                f.write("\n")
            
            for parameter in root.geometry.group:
                line = root.geometry.group[parameter]
                f.write(line)
            f.write("\n")

            for parameter in root.timestep.group:
                line = root.timestep.group[parameter]
                f.write(line)
            f.write("\n")

            for parameter in root.complexes.group:
                line = root.complexes.group[parameter]
                f.write(line)
            f.write("\n")

            line = root.main.group["rboost"]
            f.write(line+"\n")
            line = root.main.group["novers"]
            f.write(line+"\n")

        return
    
class Reaction():
    """
    Reaction file class.
    """
    
    def __init__(self):
        self.atoms_rec = []
        self.atoms_lig = []
        return
    
    def make_input(self):
        """
        Creates an SDA reaction file.
        """
        
        rxna_lines = []
        for ghost_atom_rec, ghost_atom_lig in zip(self.atoms_rec, self.atoms_lig):
            reaction = "CNONS " + ghost_atom_rec + " |" + "{:>9.2f}".format(8) + \
            "| " + ghost_atom_lig
            rxna_lines.append(reaction)

        return rxna_lines
    
    def write(self, filename):
        """
        Write the SDA reaction file to a plain text.
        
        Parameters:
        -----------
        filename : str
            The name of the rxna output file.
        """
        
        root = self.make_input()
        with open(filename, 'w') as f:
            for line in root:
                f.write(line)
        return

class APBS_Input():
    def __init__(self):
        self.temperature = None
        self.solute = Solute()
        self.solvent = Solvent()
        self.solute_pqr = ""
        self.solute_name = ""

    def make_input(self):
        """
        Prepares APBS input files.
        """

        
        pqr_filename = os.path.basename(\
            self.solute.pqr_filename)
        solute_name = ".".join(pqr_filename.split(".")[:-1])
        grid = self.solute.apbs_grid_spacing
        dime = self.solute.apbs_grid_dime
        
        apbs_lines = []
        apbs_lines.append("read\n" \
                          +"  mol pqr "+pqr_filename+" \n" \
                          +"end \n" \
                          +"\n" \
                          +"elec name viz \n" \
                          +"  mg-manual \n" \
                          +"  dime " + str(dime) + " " + str(dime) + " " + str(dime) + "\n" \
                          +"  grid " + str(grid) + " " + str(grid) + " " + str(grid) + "\n" \
                          +"  gcent  mol 1 \n" \
                          +"  mol 1 \n" \
                          +"  lpbe \n" \
                          +"  bcfl sdh \n")
        for ion in self.solvent.ions:
            apbs_lines.append("  ion charge " + str(ion.charge) + \
                              " conc " + str(ion.conc) + \
                              " radius " + str(ion.radius) +"\n")
        apbs_lines.append("  pdie " + str(self.solute.dielectric) + "\n" \
                          +"  sdie " + str(self.solvent.dielectric) + "\n" \
                          +"  chgm spl0 \n" \
                          +"  srfm smol \n" \
                          +"  srad 0.0 \n" \
                          +"  swin 0.3 \n" \
                          +"  temp " + str(self.solvent.temperature) + "\n" \
                          +"  write pot uhbd " + solute_name + "_apbs\n" \
                          +"end \n" \
                          +"quit \n")
        return apbs_lines

    def write(self):

        """
        Write the APBS input file to a plain text.
        """

        apbs_lines = self.make_input()

        with open(self.solute_name+"_apbs.in", "w") as apbs_file:
            for line in apbs_lines:
                apbs_file.write(line)
            
        return

class SDA_grid_inputs():

    """
    Prepares SDA grid inputs.
    """

    def __init__(self):
        self.solute = Solute()
        self.temperature = None
        self.ionic_strength = None
        self.solute_name = ""
        self.solute_pdb = ""
        self.solvent = Solvent()
        self.ecm_input = ""
        self.ecm_lines = []
        self.tcha = ""
        self.echa = ""
        self.epf = ""
        self.edf_input = ""
        self.edf_lines = []
        self.hdf_input = ""
        self.hdf_lines = []
        self.lj_repf_input = ""
        self.lj_repf_lines = []
        self.edf = ""
        self.hdf = ""
        self.lj_repf = ""
        self.dime = 110
        self.grid_spacing = 0.5

    def make_inputs(self):
        
        """
        Create the grid and effective charge input files
        """

        self.ecm_lines.append("#------------------ pdb file name\n" \
                              +self.solute_pdb+"\n" \
                              +"#------------------ file with test charges for a molecule\n" \
                              +self.tcha+"\n" \
                              +"#------------------ grid file name\n" \
                              +self.epf+"\n" \
                              +"#------------------ probe, skin: expansion happens in [probe; probe+skin] interval\n" \
                              +"4.0, 3.0\n" \
                              +"#------------------ ionic strength, solvent dielectric\n" \
                              +str(self.ionic_strength)+", "+str(self.solvent.dielectric)+"\n" \
                              +"#------------------ file to write effective charges\n" \
                              +self.echa+"\n" \
                              +"#------------------ NEW: ADDED in  version flex, reg_charge, if 0 will not print echa, and then will need ecm_mkreglev for finalising effective charges\n" \
                              +"1.0")


        self.edf_lines.append("#------------------------------ h,ndimx,ndimy,ndimz\n" \
                              +str(self.grid_spacing)+", "+str(self.dime)+", "+str(self.dime)+", "+str(self.dime)+"\n" \
                              +"#------------------------------ iostr,epssol,rion\n" \
                              +"0.  "+str(self.solvent.dielectric)+"  1.5\n" \
                              +"#------------------------------ pfile\n" \
                              +self.solute_pdb+"\n" \
                              +"#------------------------------ efile, iform\n" \
                              +self.edf+"\n" \
                              +"0")
        
        self.hdf_lines.append("#------------------------------ h,ndimx,ndimy,ndimz\n" \
                              +str(self.grid_spacing)+", "+str(self.dime)+", "+str(self.dime)+", "+str(self.dime)+"\n" \
                              +"#------------------------------ a,b, grid-value\n" \
                              +"3.10, 4.35, 0.5\n" \
                              +"#------------------------------ pfile\n" \
                              +self.solute_pdb+"\n" \
                              +"#------------------------------ efile, iform\n" \
                              +self.hdf+"\n" \
                              +"0")

        if len(self.lj_repf_input) > 0:
            self.lj_repf_lines.append("#------------------------------ h,ndimx,ndimy,ndimz\n" \
                                      +str(self.grid_spacing)+", "+str(self.dime)+", "+str(self.dime)+", "+str(self.dime)+"\n" \
                                      +"#------------------------------ factor, nexp, fraction\n" \
                                      +"4096.d0, 6, 1.5d0\n" \
                                      +"#------------------------------ pfile\n" \
                                      +self.solute_pdb+"\n" \
                                      +"#------------------------------ efile, iform\n" \
                                      +self.lj_repf+"\n" \
                                      +"0")
        return


    def write(self):
        
        """
        Write the grid and effective charge input files
        """

        self.make_inputs()

        with open(self.ecm_input, "w") as ecm:
            for line in self.ecm_lines:
                ecm.write(line)

        with open(self.edf_input, "w") as ed:
            for line in self.edf_lines:
                ed.write(line)
        
        with open(self.hdf_input, "w") as hd:
            for line in self.hdf_lines:
                hd.write(line)

        if len(self.lj_repf_input) > 0:
            with open(self.lj_repf_input, "w") as lj_rep:
                for line in self.lj_repf_lines:
                    lj_rep.write(line)

        

class Atomic_parameters(Serializer):
    """
    An SDA atomic descriptions

    Attributes:
    --------------------
    vdw: float or None
        VdW radius of the atom.
    test_charge: float or None
        Test charge of the
    resname: string or None, Default None
        Name of the residue where the atom belongs
    type: string or none, Default None
        Atom name
    """

    def __init__(self):
        self.vdw = -1
        self.test_charge = -1
        self.resname = None
        self.type = None

class Hydropro():

    """
    Class to set the hydropro values

    Attributes:
    --------------------
    temperature: float
        Temperature of the simulation.
    pdbfile: str
        PDB file name containing the heavy atoms of the solute.
    """

    def __init__(self):
        self.temperature = -1
        self.pdbfile = ""

    def write_input(self, input_filename, small=False):

        """
        Write the HYDROpro input file to a plain text.
        
        Parameters:
        -----------
        input_filename : str
            The name of the HYDROpro input file to write.
        small : bool
            Boolean to indicate if the solute is a small
            organic compound. Necessary to set the AER value.
        """
        
        weight = self.get_MW(self.pdbfile)

        print("The Molecular Weight is %.2f \n" % weight)

        file = open(input_filename, 'w')
        file.write(self.pdbfile.split('.')[0].split('\\')[-1]+ "                        !Name of molecule\n")
        file.write(self.pdbfile.split('.')[0].split('\\')[-1]+ "                        !Name for output file\n")
        file.write(self.pdbfile.split('\\')[-1]+ "        !Strucutural (PBD) file\n")
        file.write("1               !Type of calculation\n")
        if small:
            file.write("1.2,            !AER, radius of primary elements\n")
        else:
            file.write("2.9,            !AER, radius of primary elements\n")
        file.write("-1,              !NSIG\n")

        self.temperature = self.temperature - 273.15
        file.write(str(self.temperature)+",            !T (temperature, centigrade)\n")

        if self.temperature >= 19.0 and self.temperature <= 21.0:
            file.write("0.01,           !ETA (Viscosity of the solvent in poises)\n")

        elif self.temperature >= 24.0 and self.temperature <= 26.0:
            file.write("0.0091,         !ETA (Viscosity of the solvent in poises)\n")

        else:
            file.write("0.0091,         !ETA (Viscosity of the solvent in poises)\n")

        file.write(str(weight)+",         !RM (Molecular weight)\n")
        psv = 0.730
        #float(raw_input("Enter the partial specific volume of molecule (cm3/g) : for example 0.702 for lysozyme protein\n"))
        file.write(str(psv)+ ",         !Partial specific volume, cm3/g\n")
        file.write("1.0,            !Solvent density, g/cm3\n")
        file.write("-1              !Number of values of Q\n")
        file.write("-1          !Number of intervals\n")
        file.write("0,              !Number of trials for MC calculation of covolume\n")
        file.write("1               !IDIF=1 (yes) for full diffusion tensors\n")
        file.write("*")
        file.close()
        print("The input file of HYDROpro is stored in the same directory with the name hydropro.dat")

    def read_output(self, output_filename):
        
        """
        Read the HYDROpro output file.
        
        output_filename:
        -----------
        filename : str
            Output file of the HYDROpro calculations, usually
            finishing in *-res.txt
        """

        filelines = open(output_filename,'r').readlines()
        trans = "Translational diffusion coefficient"
        rot = "Rotational diffusion coefficient"

        for line in filelines:
            if trans in line:
                #print("HYDROpro Result for Translational diffusion coefficient:\n")
                #print(line)
                trans_value= float(line[41:50])
                trans_value= trans_value* 10**4
                print("The value of Translational diffusion coefficient to be entered in SDA input file(unit: Ang^2/ps):   %f\n" %(trans_value))

            if rot in line:
                #print("HYDROpro Result for Rotational diffusion coefficient:\n")
                #print(line)
                rot_value= float(line[41:50])
                rot_value= rot_value/10**12
                print("The value of Rotational diffusion coefficient to be entered in SDA input file(unit: radian^2/ps):   %f\n" %(rot_value))
        
        return [trans_value, rot_value]

    def get_MW(self, pdbfilename):

        """
        Calculate the molecular weight of the solute.
        
        Parameters:
        -----------
        pdbfile: str
            PDB file name containing the heavy atoms of the solute.
        """

        AMU_dict = {
            'H': 1.00794,
            'C': 12.0107,
            'CL': 35.453,
            'N': 14.0067,
            'O': 15.9994,
            'P': 30.973762,
            'S': 32.065,
            'BR': 79.904,
            'I': 126.90447,
            'F': 18.9984032,
        }
        MW = 0
        filelines = open(pdbfilename, "r").readlines()
        for line in filelines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atomname = line[12:16].strip()
                element_one_letter = atomname[0]
                element_two_letters = atomname[:2]
                if element_two_letters in AMU_dict:
                    MW += AMU_dict[element_two_letters]
                elif element_one_letter in AMU_dict:
                    MW += AMU_dict[element_one_letter]
                else:
                    print("WARNING:", atomname, "not recognized." )
        return MW



def create_ghost_atom_from_atoms_center_of_mass(
        pqr_filename, atom_index_list, new_pqr_filename=None,
        center_molecule=True):
    """
    Add a ghost atom to a PQR file at the location of the center
    of mass of the atoms listed.
    
    Parameters:
    -----------
    pqr_filename : str
        The name of the PQR file where the ghost atom will be added
        into the file itself.
        
    atom_index_list : list
        A list of integers which are atom indices. The center of mass
        will be found for these atoms, and the ghost atom will be 
        located at that center of mass.
        
    new_pqr_filename : str or None
        If None, the pqr_filename will be overwritten. Otherwise, the
        new PQR with the ghost atom will be written to this file.
    """
    
    if new_pqr_filename is None:
        new_pqr_filename = pqr_filename
    pqr_struct = parmed.load_file(pqr_filename, skip_bonds=True)


    # Compute the center of mass for the selected group of atoms
    center_of_mass = np.array([[0., 0., 0.]])
    total_mass = 0.0
    for atom_index in atom_index_list:
        atom_pos = pqr_struct.coordinates[atom_index,:]
        atom_mass = pqr_struct.atoms[atom_index].mass
        if atom_mass == 0.0:
            atom_mass = 0.0001
        center_of_mass += atom_mass * atom_pos
        total_mass += atom_mass
    center_of_mass = center_of_mass / total_mass
    
    if center_molecule:
        print("calculating center of mass")
        # Compute the center of mass of the entire molecule to be transposed
        mol_center_of_mass = np.array([[0., 0., 0.]])
        mol_total_mass = 0.0
        step_range = max(len(pqr_struct.atoms) // 100, 1)
        
        pqr_coordinates = pqr_struct.coordinates
        for atom_index, atom in enumerate(pqr_struct.atoms):
            atom_pos = pqr_coordinates[atom_index,:]
            atom_mass = atom.mass
            if atom_mass == 0.0:
                atom_mass = 0.0001
            mol_center_of_mass += atom_mass * atom_pos
            mol_total_mass += atom_mass
            if atom_index % step_range == 0:
                print(f"{((atom_index + 1) / len(pqr_struct.atoms)) * 100:.0f}%")

        mol_center_of_mass = mol_center_of_mass / mol_total_mass
    ghost_atom = parmed.Atom(name="GHO", mass=0.0, charge=0.0, solvent_radius=0.0)
    ghost_structure = parmed.Structure()
    ghost_structure.add_atom(ghost_atom, "GHO", 1)
    ghost_structure.coordinates = np.array(center_of_mass)

    #pqr_complex = pqr_struct + ghost_structure
    pqr_complex = pqr_struct
    for residue in pqr_complex.residues:
        residue.chain = ""


    if center_molecule:
        print("centering molecule")
        step_range = max(len(pqr_complex.atoms) // 100, 1)
        new_coordinates = np.zeros(pqr_complex.coordinates.shape)
        pqr_coordinates = pqr_complex.coordinates
        for atom_index, atom in enumerate(pqr_struct.atoms):
            new_coordinates[atom_index,:] = pqr_coordinates[atom_index,:] \
                - mol_center_of_mass[0,:]
            
            atom.occupancy = atom.charge
            atom.bfactor = atom.solvent_radius

            if atom_index % step_range == 0:
                print(f"{((atom_index + 1) / len(pqr_complex.atoms)) * 100:.0f}%")
                
        pqr_complex.coordinates = new_coordinates
                

    
    ghost_index = len(pqr_complex.atoms)
    x, y, z = ghost_structure.coordinates[0][:]
    ghost_atom = "ATOM{:>7s} {:4s} {:3s}  {:>4s}{:>12.3f}{:>8.3f}{:>8.3f}".format(
    str(ghost_index), ghost_atom.name, ghost_atom.residue.name, str(ghost_atom.residue.idx), 
    x, y, z)


    pqr_complex.save(new_pqr_filename, overwrite=True, format='pdb')

    return ghost_atom