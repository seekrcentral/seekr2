"""
test_common_sim_openmm.py
"""

import os

import pytest
import numpy as np
import mdtraj
try:
    import openmm.unit as unit
except ImportError:
    import simtk.openmm.unit as unit
from seekr2.modules import common_sim_openmm

def test_write_toy_pdb_file(toy_mmvt_model):
    anchor = toy_mmvt_model.anchors[1]
    out_file_name = os.path.join(
        toy_mmvt_model.anchor_rootdir, anchor.directory, 
        anchor.building_directory, "toy.pdb")
    system, topology = common_sim_openmm.make_toy_system_object(toy_mmvt_model)
    positions = np.array([[3.0, 4.0, 5.0]]) * unit.nanometers
    common_sim_openmm.write_toy_pdb_file(topology, positions, out_file_name)
    assert os.path.exists(out_file_name)
    traj = mdtraj.load(out_file_name)
    assert traj.xyz[0,0,0] == 3.0
    assert traj.xyz[0,0,1] == 4.0
    assert traj.xyz[0,0,2] == 5.0
    return

def test_make_toy_system_object(toy_mmvt_model):
    system, topology = common_sim_openmm.make_toy_system_object(toy_mmvt_model)
    force = system.getForce(0)
    assert force.getNumGroups() == 1
    assert force.getGlobalParameterName(0) == "k"
    return

def test_create_openmm_system(toy_mmvt_model, host_guest_mmvt_model):
    my_sim_openmm = common_sim_openmm.Common_sim_openmm()
    anchor = toy_mmvt_model.anchors[0]
    system, topology, positions, box_vectors, num_frames \
        = common_sim_openmm.create_openmm_system(
        my_sim_openmm, toy_mmvt_model, anchor, frame=0, load_state_file=None)
    assert num_frames == 2
    assert np.isclose(positions, np.array([[0.0, -0.7, 0.0]])\
                      *unit.nanometers).all()
    system2, topology2, positions2, box_vectors2, num_frames2 \
        = common_sim_openmm.create_openmm_system(
        my_sim_openmm, toy_mmvt_model, anchor, frame=1, load_state_file=None)
    assert np.isclose(positions2, np.array([[0.3, -0.7, 0.0]])\
                      *unit.nanometers).all()
    
    anchor2 = host_guest_mmvt_model.anchors[0]
    system3, topology3, positions3, box_vectors3, num_frames3 \
        = common_sim_openmm.create_openmm_system(
        my_sim_openmm, host_guest_mmvt_model, anchor2)
    

def test_add_barostat(host_guest_mmvt_model_npt):
    my_sim_openmm = common_sim_openmm.Common_sim_openmm()
    anchor = host_guest_mmvt_model_npt.anchors[0]
    system, topology, positions, box_vectors, num_frames \
        = common_sim_openmm.create_openmm_system(
        my_sim_openmm, host_guest_mmvt_model_npt, anchor, frame=0)
    forces = system.getForces()
    old_force_len = len(forces)
    common_sim_openmm.add_barostat(system, host_guest_mmvt_model_npt)
    forces2 = system.getForces()
    assert len(forces2) == old_force_len + 1
    barostat = forces2[-1]
    pressure = barostat.getDefaultPressure()
    assert pressure == 1.0 * unit.bar
    temperature = barostat.getDefaultTemperature()
    assert temperature == 298.15 * unit.kelvin
    
def test_add_platform_ref(toy_mmvt_model):
    my_sim_openmm = common_sim_openmm.Common_sim_openmm()
    anchor = toy_mmvt_model.anchors[0]
    system, topology, positions, box_vectors, num_frames \
        = common_sim_openmm.create_openmm_system(
        my_sim_openmm, toy_mmvt_model, anchor, frame=0)
    common_sim_openmm.add_platform(my_sim_openmm, toy_mmvt_model)
    assert my_sim_openmm.platform.getName() == "Reference"
    assert my_sim_openmm.properties == {}

@pytest.mark.needs_cuda
def test_add_platform_cuda(host_guest_mmvt_model):
    host_guest_mmvt_model.openmm_settings.cuda_platform_settings\
        .cuda_device_index = "2"
    my_sim_openmm = common_sim_openmm.Common_sim_openmm()
    anchor = host_guest_mmvt_model.anchors[0]
    system, topology, positions, box_vectors, num_frames \
        = common_sim_openmm.create_openmm_system(
        my_sim_openmm, host_guest_mmvt_model, anchor, frame=0)
    common_sim_openmm.add_platform(my_sim_openmm, host_guest_mmvt_model)
    assert my_sim_openmm.platform.getName() == "CUDA"
    assert my_sim_openmm.properties["CudaPrecision"] == "mixed"
    assert my_sim_openmm.properties["CudaDeviceIndex"] == "2"