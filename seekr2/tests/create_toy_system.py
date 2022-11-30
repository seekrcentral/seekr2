"""
create_toy_system.py

For creating systems of dummy particles that can be used to test
boundary and restraining forces, etc.

"""

import numpy as np

try:
    import openmm.unit as unit
except ImportError:
    import simtk.unit as unit

try:
    import openmm
except ImportError:
    import simtk.openmm as openmm
    
try:
    import openmm.app as openmm_app
except ImportError:
    import simtk.openmm.app as openmm_app
    
import seekr2.modules.mmvt_cvs.mmvt_closest_pair_cv as mmvt_closest_pair_cv

def make_toy_system_and_topology(num_particles):
    """
    Make openmm system and topology files for use with a toy system
    with n_particles.
    """
    system = openmm.System()
    topology = openmm_app.Topology()
    for i in range(num_particles):
        mass = 10.0 * openmm.unit.amu
        system.addParticle(mass)
        atom_name = "X{}".format(i)
        if not atom_name in openmm_app.element.Element._elements_by_symbol:
            element = openmm_app.element.Element(0, atom_name, atom_name, mass)
        else:
            element = openmm_app.element.Element._elements_by_symbol[atom_name]
        chain = topology.addChain()
        residue = topology.addResidue("UNK", chain)
        topology.addAtom(atom_name, element, residue)
        
    return system, topology

def make_toy_simulation_and_context(system, topology, positions):
    integrator = openmm.LangevinIntegrator(
        300.0*unit.kelvin, 1.0/unit.picoseconds, 0.002*unit.picoseconds)
    ref_platform = openmm.Platform.getPlatformByName('Reference')
    ref_properties = {}
    simulation = openmm_app.Simulation(
        topology, system, integrator, ref_platform, ref_properties)
    simulation.context.setPositions(positions)
    simulation.context.setVelocitiesToTemperature(300.0*unit.kelvin)
    return simulation

def assign_nonbonded_cv_info(cv, system, box_vectors, cutoff=None):
    """
    Update the cv with the correct nonbonded information.
    """
    cv.num_system_particles = system.getNumParticles()
    #forces = { force.__class__.__name__ : force for force in system.getForces() }
    #reference_force = forces['NonbondedForce']
    cv.exclusion_pairs = []
    #for index in range(reference_force.getNumExceptions()):
    #    [iatom, jatom, chargeprod, sigma, epsilon] \
    #        = reference_force.getExceptionParameters(index)
    #    cv.exclusion_pairs.append((iatom, jatom))
    
    if cutoff is None:
        cv.cutoff_distance = 0.4 * box_vectors.get_min_length()
    else:
        cv.cutoff_distance = cutoff
    print("cv.cutoff_distance:", cv.cutoff_distance)
        
    return