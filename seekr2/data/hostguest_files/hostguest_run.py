import openmm.app as app
import openmm as mm
import unit as unit
from sys import stdout
prmtop = app.AmberPrmtopFile('hostguest.parm7')
inpcrd = app.AmberInpcrdFile('hostguest.rst7')
system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1*unit.nanometer,
        constraints=app.HBonds)

input_files = ['hostguest_at1.9.pdb']
output_files = ['hostguest_1.9.pdb']

for input_file, output_file in zip(input_files, output_files):
    mypdb = app.PDBFile()

    integrator = mm.LangevinIntegrator(277.8*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)

    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}
    simulation = app.Simulation(prmtop.topology, system, integrator, platform, properties)
    simulation.context.setPositions(mypdb.positions)
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    print("inpcrd.boxVectors:", inpcrd.boxVectors)
    #simulation.minimizeEnergy()
    simulation.reporters.append(app.PDBReporter(output_file, 1000))
    simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True,
            potentialEnergy=True, temperature=True, volume=True))
    simulation.step(1000)
