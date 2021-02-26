"""
sim_browndye2

Base objects and routines for preparing and running Browndye2 
simulations.
"""

import xml.etree.ElementTree as ET
from xml.dom import minidom
import os
import glob
import re

import parmed
import numpy as np

#import openmmvt.libraries.serializer.serializer as serializer

BROWNDYE_TRAJ_PREFIX = "traj"
BROWNDYE_INPUT_FILENAME = "input.xml"
BROWNDYE_APBS_INPUT_FILENAME = "apbs_input.xml"
BROWNDYE_RECEPTOR = "receptor"
BROWNDYE_LIGAND = "ligand"

class Ion():
    """
    Represent an Ion object in Browndye2 for generating APBS grids
    and Debye Length.
    
    Attributes:
    -----------
    radius : float
        The radius of the ion in units of Angstroms.
    charge : float
        The charge of the ion in units of 
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
        self.debye_length = -1.0
        self.dielectric = 78.0
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
        assert self.desolvation_parameter > 0.0
        xmlSolventDesolv = ET.SubElement(xmlSolvent, 'desolvation_parameter')
        xmlSolventDesolv.text = str(self.desolvation_parameter)
        if make_apbs_mode:
            xmlSolventIons = ET.SubElement(xmlSolvent, 'ions')
            for ion in self.ions:
                xmlIon = ET.SubElement(xmlSolventIons, 'ion')
                xmlIon.text = ion.serialize(xmlIon)
        
        return

class Time_step_tolerances():
    """
    Parameters to represent the limitations on the time step sizes
    within the BD simulation.
    
    Attributes:
    -----------
    force : float or None, Default None
        Dimensionless parameter that governs the splitting of the 
        adaptive time step in response to changes in force and torque.
    reaction : float or None, Default None
        Dimensionless parameter that governs the splitting of the 
        adaptive time step in response to changes in the reaction
        coordinate near a reaction.
    minimum_core_dt : float, default 0.0
        For a system with only rigid cores and frozen chains, this is
        the smallest possible time step size in picoseconds.
    minimum_chain_dt : float or None, default None
        For a system with unfrozen chains, this is the smallest 
        possible time step size in picoseconds.
    minimum_core_reaction_dt : float, default 0.0
        If a reaction coordinate is small, then the minimum step size
        given by minimum_core_dt can be overridden by the value given
        here.
    minimum_chain_reaction_dt : float or None, default None
        This has the same relation to minimum_chain_dt.
        
    """
    def __init__(self):
        self.force = None
        self.reaction = None
        self.minimum_core_dt = 0.0
        self.minimum_chain_dt = None
        self.minimum_core_reaction_dt = 0.0
        self.minimum_chain_reaction_dt = None
        return
    
    def serialize(self, xmlTime):
        """
        Serialize this object to XML
        
        Parameters:
        -----------
        xmlTime : ElementTree.SubElement
            All sub elements get added to this root.
        """
        if self.force is not None:
            assert self.force > 0.0
            xmlTimeForce = ET.SubElement(xmlTime, "force")
            xmlTimeForce.text = str(self.force)
        if self.reaction is not None:
            assert self.reaction > 0.0
            xmlTimeReaction = ET.SubElement(xmlTime, "reaction")
            xmlTimeReaction.text = str(self.reaction)
        assert self.minimum_core_dt >= 0.0
        xmlTimeMinCoreDt = ET.SubElement(xmlTime, "minimum_core_dt")
        xmlTimeMinCoreDt.text = str(self.minimum_core_dt)
        if self.minimum_chain_dt is not None:
            assert self.minimum_chain_dt >= 0.0
            xmlTimeMinChainDt = ET.SubElement(xmlTime, "minimum_chain_dt")
            xmlTimeMinChainDt.text = str(self.minimum_chain_dt)
        assert self.minimum_core_reaction_dt >= 0.0
        xmlTimeMinCoreRxnDt = ET.SubElement(xmlTime, "minimum_core_reaction_dt")
        xmlTimeMinCoreRxnDt.text = str(self.minimum_core_reaction_dt)
        if self.minimum_chain_reaction_dt is not None:
            assert self.minimum_chain_reaction_dt >= 0.0
            xmlTimeMinChainRxnDt = ET.SubElement(
                xmlTime, "minimum_chain_reaction_dt")
            xmlTimeMinChainRxnDt.text = str(self.minimum_chain_reaction_dt)
        return

class Electric_field():
    """
    In Browndye2, Electric_field tags define the various electric
    electric forces applied to a molecule
    
    Attributes:
    -----------
    grid_list : list, Default []
        A list of DX file names output from APBS.
    multipole_field : str or None, Default None
        A file of multipole extension outside grids.
    desolvation_field : str or None
        This block describes the desolvation field; it is just like the
        electric-field block but does not use the multipole field.
    eff_charges : str or None
        A file with effective charges and lumping information; output
        of lumped_charges.
    eff_charges_squared : str or None
        File like eff_charges above, but with the charges squared.
        
    """
    def __init__(self):
        self.grid_list = []
        self.multipole_field = None
        self.desolvation_field = None
        self.eff_charges = None
        self.eff_charges_squared = None
        self.copy = None
        return
    
    def serialize(self, xmlE):
        """
        Serialize this object to XML
        
        Parameters:
        -----------
        xmlE : ElementTree.SubElement
            All sub elements get added to this root.
        """
        for grid in self.grid_list:
            xmlGrid_i = ET.SubElement(xmlE, 'grid')
            xmlGrid_i.text = grid
        if self.multipole_field is not None:
            assert self.multipole_field
            xmlEMultipole = ET.SubElement(xmlE, "multipole_field")
            xmlEMultipole.text = self.multipole_field
        if self.desolvation_field is not None:
            assert self.desolvation_field
            xmlEDesolv = ET.SubElement(xmlE, "desolvation_field")
            xmlEDesolv.text = self.desolvation_field
        if self.eff_charges is not None:
            assert self.eff_charges
            xmlEEff = ET.SubElement(xmlE, "eff_charges")
            xmlEEff.text = self.eff_charges
        if self.eff_charges_squared is not None:
            assert self.eff_charges_squared
            xmlEEff2 = ET.SubElement(xmlE, "eff_charges_squared")
            xmlEEff2.text = self.eff_charges_squared
        if self.copy is not None:
            xmlECopy = ET.SubElement(xmlE, 'copy')
            xmlECopy.text = self.copy.serialize(xmlECopy)
        return

class Link():
    """
    In Browndye2, chains are connected to cores - defining the Link.
    
    Attributes:
    -----------
    core_name : str
        The name of the Core to link.
    core_residue : int
        The residue index within the Core to link to.
    chain_residue : int
        The residue index within the chain to link to.
    """
    def __init__(self):
        self.core_name = ""
        self.core_residue = -1
        self.chain_residue = -1
        return
        
    def serialize(self, xmlLink):
        """
        Serialize this object to XML
        
        Parameters:
        -----------
        xmlLink : ElementTree.SubElement
            All sub elements get added to this root.
        """
        assert self.core_name, "A name for the core must be provided"
        xmlLinkCoreName = ET.SubElement(xmlLink, "core_name")
        xmlLinkCoreName.text = self.core_name
        assert self.core_residue >= 0
        xmlLinkCoreResid = ET.SubElement(xmlLink, "core_residue")
        xmlLinkCoreResid.text = str(self.core_residue)
        assert self.chain_residue >= 0
        xmlLinkChainResid = ET.SubElement(xmlLink, "chain_residue")
        xmlLinkChainResid.text = str(self.chain_residue)
        return

class Chain():
    """
    In Browndye2, chains are flexible assemblies of spheres that can
    interact with each other and the cores with user-specified force
    fields.
    
    Attributes:
    -----------
    name : str
        The name of this Chain.
    atoms : str
        The filename of the atoms in this chain. When using the
        COFFDROP model, the atoms files for cores and chains
        must contain the coarse-grained "beads".
    link_list : list
        A list of Links, or connections between the chain and core. The
        name of the core, the number of the linking residue, and the
        number of the chain's linking residue are specified.
        
    """
    def __init__(self):
        self.name = ""
        self.atoms = ""
        self.link_list = []
        return
    
    def serialize(self, xmlChain):
        """
        Serialize this object to XML
        
        Parameters:
        -----------
        xmlChain : ElementTree.SubElement
            All sub elements get added to this root.
        """
        assert self.name, "A name for this chain must be provided"
        xmlChainName = ET.SubElement(xmlChain, "name")
        xmlChainName.text = self.name
        assert self.atoms, "An atoms file must be provided for this Chain"
        xmlChainAtoms = ET.SubElement(xmlChain, "atoms")
        xmlChainAtoms.text = self.atoms
        for link in self.link_list:
            xmlLink_i = ET.SubElement(xmlChain, 'link')
            xmlLink_i.text = link.serialize(xmlLink_i)
        return

class Core():
    """
    In Browndye2, cores are large rigid bodies composed of many smaller
    atoms.
    
    Attributes:
    -----------
    name : str
        The name to identify this core.
    atoms : str
        The PQRXML file of the atoms of this core, created by the 
        program pqr2xml.
    all_in_surface : str or None, Default "false"
        If this is true, then all of the atoms of the core are
        included, not just a shell made of the surface atoms.
    electric_field : Electric_field() or None
        Object containing APBS grids, desolvation grids, and other
        information pertaining to electrical forces.
    is_protein : str or None, Default "false"
        If this is true, then it affects how the effective charges are
        generated.
    dielectric : float, default 4.0
        The interior dielectric of the core.
    """
    def __init__(self):
        self.name = ""
        self.atoms = ""
        self.all_in_surface = "false"
        self.electric_field = Electric_field()
        self.is_protein = "false"
        self.dielectric = 4.0
        self.grid_spacing = 0.5
        return
    
    def serialize(self, xmlCore, make_apbs_mode=True):
        """
        Serialize this object to XML
        
        Parameters:
        -----------
        xmlCore : ElementTree.SubElement
            All sub elements get added to this root.
        make_apbs_mode : bool
            Whether this object should be serialized for the
            make_apbs_inputs program.
        """
        assert self.name, "A name must be assigned to a Core"
        xmlCoreName = ET.SubElement(xmlCore, "name")
        xmlCoreName.text = self.name
        assert self.atoms, "An atoms XML file must be provided for this Core"
        xmlCoreAtoms = ET.SubElement(xmlCore, "atoms")
        xmlCoreAtoms.text = self.atoms
        if self.all_in_surface is not None:
            assert self.all_in_surface.lower() in ["true","false"]
            xmlCoreAllInSurf = ET.SubElement(xmlCore, "all_in_surface")
            xmlCoreAllInSurf.text = self.all_in_surface.lower()
        if not make_apbs_mode:
            if self.electric_field is not None:
                xmlCoreElecField = ET.SubElement(xmlCore, 'electric_field')
                xmlCoreElecField.text = self.electric_field.serialize(
                    xmlCoreElecField)
        
        if self.is_protein is not None:
            assert self.is_protein.lower() in ["true","false"]
            xmlCoreIsProt = ET.SubElement(xmlCore, "is_protein")
            xmlCoreIsProt.text = self.is_protein.lower()
        assert self.dielectric > 0.0
        xmlCoreDielec = ET.SubElement(xmlCore, "dielectric")
        xmlCoreDielec.text = str(self.dielectric)
        if make_apbs_mode:
            assert self.grid_spacing > 0.0
            xmlCoreGridSpacing = ET.SubElement(xmlCore, "grid_spacing")
            xmlCoreGridSpacing.text = str(self.grid_spacing)
        return

class CoreCopy():
    """
    In Browndye2, cores may have all information copied to another
    core with a different translational/rotational position.
    
    Attributes:
    -----------
    parent : str
        The name of the core to copy from
    translation : list
        A list of 3 floats that describes the new translation.
    rotation : list
        A list of 9 floats that describes the new rotation.
    """
    def __init__(self):
        self.parent = ""
        self.translation = []
        self.rotation = []
        return
    
    def serialize(self, xmlCoreCopy):
        """
        Serialize this object to XML
        
        Parameters:
        -----------
        xmlCoreCopy : ElementTree.SubElement
            All sub elements get added to this root.
        """
        assert self.parent, "A name must be assigned to a CoreCopy"
        xmlCoreCpName = ET.SubElement(xmlCoreCopy, "parent")
        xmlCoreCpName.text = self.parent
        assert len(translation) == 3
        xmlCoreCpTrans = ET.SubElement(xmlCoreCopy, "translation")
        xmlCoreCpTrans.text = " ".join([str(i) for i in self.translation])
        assert len(translation) == 9
        xmlCoreCpRot = ET.SubElement(xmlCoreCopy, "rotation")
        xmlCoreCpRot.text = " ".join([str(i) for i in self.rotation])
        return

class Group():
    """
    In Browndye2, Cores and Chains can be assembled into groups. With
    two groups present, the second-order rate constant of the encounter
    of the groups may be computed.
    
    Cores may be connected by chains within a group. At this time, 
    however, SEEKR only supports a single Core in a group.
    
    Attributes:
    -----------
    name : str
        The name of this group
    core_list : list, Default []
        A list of Core() objects for this group
    chain_list : list, Default []
        A list of Chain() objects for this group
    """
    def __init__(self):
        self.name = ""
        self.core_list = []
        self.chain_list = []
        return
        
    def serialize(self, xmlGroup, make_apbs_mode=True):
        """
        Serialize this object to XML
        
        Parameters:
        -----------
        xmlCoreCopy : ElementTree.SubElement
            All sub elements get added to this root.
        """
        assert self.name, "A name must be assigned to a Group"
        xmlGroupName = ET.SubElement(xmlGroup, "name")
        xmlGroupName.text = self.name
        for core in self.core_list:
            xmlCore_i = ET.SubElement(xmlGroup, 'core')
            xmlCore_i.text = core.serialize(xmlCore_i, make_apbs_mode)
        for chain in self.chain_list:
            xmlChain_i = ET.SubElement(xmlGroup, 'chain')
            xmlChain_i.text = chain.serialize(xmlChain_i)
        return

class System():
    """
    A Browndye2 system description.
    
    Attributes:
    -----------
    force_field : str or None
        Accepts "molecular_mechanics" or "spline" as arguments to
        select the type of forcefield.
    parameters : str or None
        Parameters file for the force field. Optional for 
        "molecular_mechanics" with no chains.
    start_at_site : str or None
        if true, then all the molecular components are started 
        according to their atom input files; otherwise, the 
        second group is started on the b-sphere and the second-
        order rate constant is computed.
    b_radius : float or None
        This can be automatically computed, but if a value is given
        it overrides the automatic value.
    reaction_file : str
        Describes the reaction pathways and criteria.
    hydrodynamic_interactions : str or None
        Include hydrodynamic interactions.
    n_steps_between_hi_updates : int or None
        Number of steps between updates to HI.
    density_field : str or None
        DX file containing an external density from microscopy or
        tomography.
    atom_weights : str or None
        File of atomic weights for each atom name in the atoms file.
    solvent : Solvent() or None
        Solvent object with parameters.
    time_step_tolerances : Time_step_tolerances()
        Time_step_tolerances object with parameters.
    group_list : list
        A list of Group() objects for the system.
    """
    def __init__(self):
        self.force_field = None
        self.parameters = None
        self.start_at_site = "false"
        self.b_radius = None
        self.reaction_file = ""
        self.hydrodynamic_interactions = "true"
        self.n_steps_between_hi_updates = None
        self.density_field = None
        self.atom_weights = None
        self.solvent = Solvent()
        self.time_step_tolerances = Time_step_tolerances()
        self.group_list = []
        return
        
    def serialize(self, xmlSys, make_apbs_mode=True):
        """
        This function needs to be defined explicitly because the 
        serializer will not format this XML block the way that 
        Browndye2 will accept as input.
        
        Parameters:
        -----------
        xmlChain : ElementTree.SubElement
            All sub elements get added to this root.
        
        make_apbs_mode : bool
            Whether this object should be serialized for the
            make_apbs_inputs program.
        """
        if self.force_field is not None:
            xmlSysFF = ET.SubElement(xmlSys, "force_field")
            xmlSysFF.text = self.force_field
        if self.parameters is not None:
            xmlSysParam = ET.SubElement(xmlSys, "parameters")
            xmlSysParam.text = self.parameters
        if self.start_at_site is not None:
            assert self.start_at_site.lower() in ["true","false"]
            xmlSysStart = ET.SubElement(xmlSys, "start_at_site")
            xmlSysStart.text = self.start_at_site.lower()
        if self.b_radius is not None:
            xmlSysBRad = ET.SubElement(xmlSys, "b_radius")
            xmlSysBRad.text = str(self.b_radius)
        assert self.reaction_file
        xmlSysRxnFile = ET.SubElement(xmlSys, "reaction_file")
        xmlSysRxnFile.text = self.reaction_file
        if self.hydrodynamic_interactions is not None:
            assert self.hydrodynamic_interactions.lower() in ["true","false"]
            xmlSysHydro = ET.SubElement(xmlSys, "hydrodynamic_interactions")
            xmlSysHydro.text = self.hydrodynamic_interactions.lower()
        if self.n_steps_between_hi_updates is not None:
            assert self.n_steps_between_hi_updates > 0
            xmlSysNStepsHI = ET.SubElement(xmlSys, "n_steps_between_hi_updates")
            xmlSysNStepsHI.text = str(self.n_steps_between_hi_updates)
        if self.density_field is not None:
            xmlSysDens = ET.SubElement(xmlSys, "density_field")
            xmlSysDens.text = self.density_field
        if self.atom_weights is not None:
            xmlSysWeights = ET.SubElement(xmlSys, "atom_weights")
            xmlSysWeights.text = self.atom_weights
        xmlSysSolv = ET.SubElement(xmlSys, 'solvent')
        xmlSysSolv.text = self.solvent.serialize(xmlSysSolv, make_apbs_mode)
        xmlSysTime = ET.SubElement(xmlSys, 'time_step_tolerances')
        xmlSysTime.text = self.time_step_tolerances.serialize(xmlSysTime)
        for group in self.group_list:
            xmlGroup_i = ET.SubElement(xmlSys, 'group')
            xmlGroup_i.text = group.serialize(xmlGroup_i, make_apbs_mode)
        return
    
class Root():
    """
    An entire Browndye2 Input XML object.
    
    Attributes:
    -----------
    n_threads : int or None
        The number of threads used.
    seed : int
        Used to seed the random number generator
    output : str
        Name of file where output is sent (not trajectories).
    n_trajectories : int
        Total number of trajectories to be run.
    n_trajectories_per_output : int or None
        Number of time steps between updates to the output file.
    max_n_steps : int
        Number of time steps taken in a trajectory before giving up
        and declaring that it is "stuck".
    trajectory_file : str or None
        Prefix of files used for outputting of trajectories.
    n_steps_per_output : int or None
        Number of time steps between output to trajectory file.
    min_rxn_dist_file : str or None
        File to output minimum reaction coordinate for each trajectory.
    system : System()
        Information about the physical system and its motions are
        described inside this section.
    """
    def __init__(self):
        self.n_threads = 1
        self.seed = 11111113
        self.output = "results.xml"
        self.n_trajectories = -1
        self.n_trajectories_per_output = 1000
        self.max_n_steps = 1000000
        self.trajectory_file = BROWNDYE_TRAJ_PREFIX
        self.n_steps_per_output = 1000
        self.min_rxn_dist_file = None
        self.system = System()
        return
    
    def serialize(self, make_apbs_mode=True):
        """
        This function needs to be defined explicitly because the 
        serializer will not format this XML block the way that 
        Browndye2 will accept as input.
        
        Parameters:
        -----------
        xmlChain : ElementTree.SubElement
            All sub elements get added to this root.
        
        make_apbs_mode : bool
            Whether this object should be serialized for the
            make_apbs_inputs program.
        """
        xmlRoot = ET.Element('root')
        if self.n_threads is not None:
            assert self.n_threads > 0
            xmlNThreads = ET.SubElement(xmlRoot, "n_threads")
            xmlNThreads.text = str(self.n_threads)
        xmlSeed = ET.SubElement(xmlRoot, "seed")
        xmlSeed.text = str(self.seed)
        xmlOutput = ET.SubElement(xmlRoot, "output")
        xmlOutput.text = self.output
        assert self.n_trajectories > 0, "n_trajectories must be set."
        xmlNTraj = ET.SubElement(xmlRoot, "n_trajectories")
        xmlNTraj.text = str(self.n_trajectories)
        if self.n_trajectories_per_output is not None:
            assert self.n_trajectories_per_output > 0
            #if self.n_trajectories_per_output > self.n_trajectories:
            #    self.n_trajectories_per_output = self.n_trajectories
            xmlNTrajPerOut = ET.SubElement(xmlRoot, "n_trajectories_per_output")
            xmlNTrajPerOut.text = str(self.n_trajectories_per_output)
        assert self.max_n_steps > 0
        xmlMaxN = ET.SubElement(xmlRoot, "max_n_steps")
        xmlMaxN.text = str(self.max_n_steps)
        if self.trajectory_file is not None:
            assert self.trajectory_file
            xmlTrajFile = ET.SubElement(xmlRoot, "trajectory_file")
            xmlTrajFile.text = str(self.trajectory_file)
        if self.n_steps_per_output is not None:
            assert self.n_steps_per_output > 0
            xmlNStepsPerOut = ET.SubElement(xmlRoot, "n_steps_per_output")
            xmlNStepsPerOut.text = str(self.n_steps_per_output)
        if self.min_rxn_dist_file is not None:
            assert self.min_rxn_dist_file
            xmlMinRxnFile = ET.SubElement(xmlRoot, "min_rxn_dist_file")
            xmlMinRxnFile.text = str(self.min_rxn_dist_file)
        xmlSys = ET.SubElement(xmlRoot, "system")
        xmlSys.text = self.system.serialize(xmlSys, make_apbs_mode)
        return xmlRoot
    
    def write(self, filename, make_apbs_mode=True):
        """
        Write the Browndye2 calculation to XML.
        
        Parameters:
        -----------
        filename : str
            The name of the file to write XML to.
        
        make_apbs_mode : bool
            Whether this object should be serialized for the
            make_apbs_inputs program.
        """
        root = self.serialize(make_apbs_mode)
        xmlstr = minidom.parseString(ET.tostring(root)).toprettyxml(
            indent="   ")
        with open(filename, 'w') as f:
            f.write(xmlstr)
        return
    
class Reaction():
    """
    
    """
    def __init__(self):
        self.name = ""
        self.state_before = ""
        self.state_after = ""
        self.molecule0_group = ""
        self.molecule0_core = ""
        self.molecule1_group = ""
        self.molecule1_core = ""
        self.n_needed = 1
        self.pair_list = []
        
    def serialize(self, xmlRoot):
        """
        
        """
        assert self.name, "reaction name must be provided"
        xmlName = ET.SubElement(xmlRoot, "name")
        xmlName.text = self.name
        assert self.state_before, "reaction state_before must be provided"
        xmlStateBefore = ET.SubElement(xmlRoot, "state_before")
        xmlStateBefore.text = self.state_before
        assert self.state_after, "reaction state_after must be provided"
        xmlStateAfter = ET.SubElement(xmlRoot, "state_after")
        xmlStateAfter.text = self.state_after
        xmlCriterion = ET.SubElement(xmlRoot, "criterion")
        xmlMolecules = ET.SubElement(xmlCriterion, "molecules")
        xmlMolecule0 = ET.SubElement(xmlMolecules, "molecule0")
        assert self.molecule0_group != "" and self.molecule0_core != ""
        xmlMolecule0.text = self.molecule0_group + " " + self.molecule0_core
        xmlMolecule1 = ET.SubElement(xmlMolecules, "molecule1")
        assert self.molecule1_group != "" and self.molecule1_core != ""
        xmlMolecule1.text = self.molecule1_group + " " + self.molecule1_core
        assert self.n_needed >= 0
        xmlN_needed = ET.SubElement(xmlCriterion, "n_needed")
        xmlN_needed.text = str(self.n_needed)
        assert len(self.pair_list) > 0
        for pair in self.pair_list:
            xmlPair_i = ET.SubElement(xmlCriterion, "pair")
            xmlPair_i.text = pair.serialize(xmlPair_i)
        
        return
        
class Pair():
    """
    
    """
    def __init__(self):
        self.atom1_index = -1
        self.atom2_index = -1
        self.distance = -1.0
        return
    
    def serialize(self, xmlPair):
        """
        
        """
        assert self.atom1_index >= 0, "Pair atom1_index must be assigned."
        assert self.atom2_index >= 0, "Pair atom2_index must be assigned."
        xmlAtoms = ET.SubElement(xmlPair, "atoms")
        xmlAtoms.text = str(self.atom1_index) + " " + str(self.atom2_index)
        assert self.distance >= 0.0, "Pair distance must be assigned."
        xmlDistance = ET.SubElement(xmlPair, "distance")
        xmlDistance.text = str(self.distance)
        return
    
class Reaction_root():
    """
    
    """
    def __init__(self):
        self.first_state = ""
        self.reaction_list = []
        return
    
    def serialize(self):
        """
        
        """
        xmlRoot = ET.Element("roottag")
        assert self.first_state
        xmlFirstState = ET.SubElement(xmlRoot, "first_state")
        xmlFirstState.text = self.first_state
        assert len(self.reaction_list) > 0
        xmlReactions = ET.SubElement(xmlRoot, "reactions")
        for reaction in self.reaction_list:
            xmlReaction_i = ET.SubElement(xmlReactions, "reaction")
            xmlReaction_i.text = reaction.serialize(xmlReaction_i)
        return xmlRoot
    
    def write(self, filename):
        """
        Write the Browndye2 reaction file to XML.
        
        Parameters:
        -----------
        filename : str
            The name of the file to write XML to.
        """
        root = self.serialize()
        xmlstr = minidom.parseString(ET.tostring(root)).toprettyxml(
            indent="   ")
        with open(filename, 'w') as f:
            f.write(xmlstr)
        return

def add_ghost_atom_to_pqr_from_atoms_center_of_mass(
        pqr_filename, atom_index_list, new_pqr_filename=None):
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
    center_of_mass = np.array([[0., 0., 0.]])
    total_mass = 0.0
    for atom_index in atom_index_list:
        atom_pos = pqr_struct.coordinates[atom_index,:]
        atom_mass = pqr_struct.atoms[atom_index].mass
        center_of_mass += atom_mass * atom_pos
        total_mass += atom_mass
    center_of_mass = center_of_mass / total_mass
    ghost_atom = parmed.Atom(name="GHO", mass=0.0, charge=0.0, solvent_radius=0.0)
    ghost_structure = parmed.Structure()
    ghost_structure.add_atom(ghost_atom, "GHO", 1)
    ghost_structure.coordinates = np.array(center_of_mass)
    complex = pqr_struct + ghost_structure
    complex.save(new_pqr_filename, overwrite=True)
    ghost_index = len(complex.atoms)
    #ghost_index = complex.atoms[-1].number
    #assert complex.atoms[-1].name == "GHO"
    #print("ghost_index:", ghost_index)
    return ghost_index

def make_pqrxml(input_pqr_filename, browndye2_bin="", 
                output_xml_filename=None):
    """
    
    """
    if output_xml_filename is None:
        output_xml_filename = os.path.splitext(input_pqr_filename)[0] + ".xml"
    pqr2xml_binary = os.path.join(browndye2_bin, "pqr2xml")
    pqr2xml_command = pqr2xml_binary + " < " + input_pqr_filename + " > " \
        + output_xml_filename
    print("Running command:", pqr2xml_command)
    os.system(pqr2xml_command)
    assert os.path.exists(output_xml_filename), "Problem generating XML from "\
        "PQR file: file not created: %s" % output_xml_filename
    return output_xml_filename

def make_and_run_apbs(root, input_apbs_xml, browndye2_bin="", 
                      new_input_xml_base="input.xml"):
    """
    
    """
    assert root.system.solvent.ions is not None, \
        "Ions must be included for APBS calculations"
    assert len(root.system.solvent.ions) > 0, \
        "Ions must be included for APBS calculations"
    mol_name_list = []
    for group in root.system.group_list:
        mol_name = group.name
        mol_name_list.append(mol_name)
    
    bd_dir = os.path.dirname(input_apbs_xml)
    os.chdir(bd_dir)
    new_input_xml = os.path.join(bd_dir, new_input_xml_base)
    make_apbs_binary = os.path.join(browndye2_bin, "make_apbs_inputs")
    run_apbs_binary = os.path.join(browndye2_bin, "run_apbs_inputs")
    make_apbs_command = make_apbs_binary + " " + input_apbs_xml + " > " \
        + new_input_xml
    #std_out = subprocess.check_output(make_apbs_command, shell=True)
    print("Running command:", make_apbs_command)
    os.system(make_apbs_command)
    assert os.path.exists(new_input_xml), "Problem running make_apbs_input - "\
        "Output file not found: " + new_input_xml
    for mol_name in mol_name_list:
        apbs_input_glob = mol_name + "*.in"
        assert len(glob.glob(apbs_input_glob)), "Problem running " \
        "make_apbs_input - APBS input files not found: " + apbs_input_glob
    
    debye_length = -1.0
    with open(new_input_xml, "r") as f:
        for line in f.readlines():
            result = re.search("<debye_length> (.+?) </debye_length>", line)
            if result:
                debye_length = result.group(1)
                print("debye_length found:", debye_length)
                
    assert float(debye_length) > 0.0, "Problem: debye_length not generated."
    
    run_apbs_command = run_apbs_binary + " " + new_input_xml \
        + " > run_apbs_input.out"
    print("Running command:", run_apbs_command)
    os.system(run_apbs_command)
    for mol_name in mol_name_list:
        apbs_dx_glob = mol_name + "*.dx"
        assert len(glob.glob(apbs_dx_glob)) > 0, "Problem running " \
        "run_apbs_input - APBS DX file(s) not found: " + apbs_dx_glob
    
    return debye_length