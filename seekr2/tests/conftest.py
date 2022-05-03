"""
conftest.py

configurations for SEEKR2 tests
"""

import os
import pytest
import copy

import numpy as np

import seekr2.prepare as prepare
import seekr2.modules.common_base as base
import seekr2.tests.create_model_input as create_model_input

TEST_DIRECTORY = os.path.dirname(__file__)

@pytest.fixture(scope="session")
def host_guest_mmvt_model_input_persistent(tmpdir_factory):
    """
    Create a model object that is persistent across the tests in this file.
    """
    rootdir = tmpdir_factory.mktemp("hostguest_mmvt")
    host_guest_mmvt_model_input_persisent_obj \
        = create_model_input.create_host_guest_mmvt_model_input(
            rootdir, bd=True)
    return host_guest_mmvt_model_input_persisent_obj

@pytest.fixture()
def host_guest_mmvt_model_input(host_guest_mmvt_model_input_persistent):
    """
    Create a copy of the model input that is not persistent. But this 
    at least doesn't require us to generate an entirely new model 
    input.
    """
    host_guest_mmvt_model_input_obj = copy.deepcopy(
        host_guest_mmvt_model_input_persistent)
    return host_guest_mmvt_model_input_obj

@pytest.fixture(scope="session")
def host_guest_mmvt_model_persistent(host_guest_mmvt_model_input_persistent):
    """
    Create a model object that is persistent across the tests in this file.
    """
    os.chdir(TEST_DIRECTORY)
    host_guest_mmvt_model_obj, model_xml_path \
        = prepare.prepare(host_guest_mmvt_model_input_persistent, 
                          force_overwrite=False)
    model_dir = os.path.dirname(model_xml_path)
    host_guest_mmvt_model_obj.anchor_rootdir = os.path.abspath(model_dir)
    return host_guest_mmvt_model_obj

@pytest.fixture
def host_guest_mmvt_model(host_guest_mmvt_model_persistent):
    """
    Create a copy of the model that is not persistent. But this at least
    doesn't require us to generate an entirely new model
    """
    host_guest_mmvt_model = copy.deepcopy(host_guest_mmvt_model_persistent)
    return host_guest_mmvt_model

@pytest.fixture(scope="session")
def host_guest_mmvt_model_persistent_namd(tmpdir_factory, 
                                     host_guest_mmvt_model_input_persistent):
    """
    Create a model object that is persistent across the tests in this file.
    It uses the NAMD engine.
    """
    os.chdir(TEST_DIRECTORY)
    host_guest_mmvt_model_input_persistent.md_program = "namd"
    host_guest_mmvt_model_input_persistent.waterModel = "tip3p"
    rootdir = tmpdir_factory.mktemp("hostguest_mmvt_namd")
    host_guest_mmvt_model_input_persistent.root_directory = rootdir
    host_guest_mmvt_model_input_persistent.browndye_settings_input = None
    host_guest_mmvt_model_obj, model_xml_path \
        = prepare.prepare(host_guest_mmvt_model_input_persistent, 
                          force_overwrite=False)
    model_dir = os.path.dirname(model_xml_path)
    host_guest_mmvt_model_obj.anchor_rootdir = os.path.abspath(model_dir)
    host_guest_mmvt_model_obj.namd_settings.PMEGridSpacing = None
    return host_guest_mmvt_model_obj

@pytest.fixture
def host_guest_mmvt_model_namd(host_guest_mmvt_model_persistent_namd):
    """
    Create a copy of the model that is not persistent. But this at least
    doesn't require us to generate an entirely new model. It uses the NAMD
    engine.
    """
    host_guest_mmvt_model_namd = copy.deepcopy(
        host_guest_mmvt_model_persistent_namd)
    return host_guest_mmvt_model_namd

@pytest.fixture(scope="session")
def host_guest_mmvt_model_persistent_npt(tmpdir_factory, 
                                     host_guest_mmvt_model_input_persistent):
    """
    Create a model object that is persistent across the tests in this file.
    It uses the NAMD engine.
    """
    os.chdir(TEST_DIRECTORY)
    host_guest_mmvt_model_input_persistent.ensemble = "npt"
    rootdir = tmpdir_factory.mktemp("hostguest_mmvt_npt")
    host_guest_mmvt_model_input_persistent.root_directory = rootdir
    host_guest_mmvt_model_input_persistent.browndye_settings_input = None
    host_guest_mmvt_model_obj, model_xml_path \
        = prepare.prepare(host_guest_mmvt_model_input_persistent, 
                          force_overwrite=False)
    model_dir = os.path.dirname(model_xml_path)
    host_guest_mmvt_model_obj.anchor_rootdir = os.path.abspath(model_dir)
    return host_guest_mmvt_model_obj

@pytest.fixture
def host_guest_mmvt_model_npt(host_guest_mmvt_model_persistent_npt):
    """
    Create a copy of the model that is not persistent. But this at least
    doesn't require us to generate an entirely new model. It uses the NAMD
    engine.
    """
    host_guest_mmvt_model_npt = copy.deepcopy(
        host_guest_mmvt_model_persistent_npt)
    return host_guest_mmvt_model_npt

@pytest.fixture(scope="session")
def host_guest_mmvt_model_input_persistent_forcefield(tmpdir_factory):
    """
    Create a model object that is persistent across the tests in this file.
    Uses an OpenMM forcefield input.
    """
    rootdir = tmpdir_factory.mktemp("hostguest_mmvt")
    host_guest_mmvt_model_input_persisent_obj \
        = create_model_input.create_host_guest_mmvt_model_input(
            rootdir, bd=False, ff="forcefield")
    return host_guest_mmvt_model_input_persisent_obj

@pytest.fixture()
def host_guest_mmvt_model_input_forcefield(
        host_guest_mmvt_model_input_persistent_forcefield):
    """
    Create a copy of the model input that is not persistent. But this 
    at least doesn't require us to generate an entirely new model 
    input. Uses an OpenMM forcefield input.
    """
    host_guest_mmvt_model_input_obj = copy.deepcopy(
        host_guest_mmvt_model_input_persistent_forcefield)
    return host_guest_mmvt_model_input_obj

@pytest.fixture(scope="session")
def host_guest_mmvt_model_persistent_forcefield(
        host_guest_mmvt_model_input_persistent_forcefield):
    """
    Create a model object that is persistent across the tests in this file.
    Uses an OpenMM forcefield input.
    """
    os.chdir(TEST_DIRECTORY)
    host_guest_mmvt_model_obj, model_xml_path \
        = prepare.prepare(host_guest_mmvt_model_input_persistent_forcefield, 
                          force_overwrite=False)
    model_dir = os.path.dirname(model_xml_path)
    host_guest_mmvt_model_obj.anchor_rootdir = os.path.abspath(model_dir)
    return host_guest_mmvt_model_obj

@pytest.fixture
def host_guest_mmvt_model_forcefield(
        host_guest_mmvt_model_persistent_forcefield):
    """
    Create a copy of the model that is not persistent. But this at least
    doesn't require us to generate an entirely new model. 
    Uses an OpenMM forcefield input.
    """
    host_guest_mmvt_model = copy.deepcopy(host_guest_mmvt_model_persistent_forcefield)
    return host_guest_mmvt_model

@pytest.fixture(scope="session")
def host_guest_elber_model_input_persistent(tmpdir_factory):
    """
    Create a model object that is persistent across the tests in this file. 
    Uses an OpenMM forcefield input.
    """
    rootdir = tmpdir_factory.mktemp("hostguest_elber")
    host_guest_elber_model_input_persisent_obj \
        = create_model_input.create_host_guest_elber_model_input(
            rootdir, bd=True)
    return host_guest_elber_model_input_persisent_obj

@pytest.fixture()
def host_guest_elber_model_input(host_guest_elber_model_input_persistent):
    """
    Create a copy of the model input that is not persistent. But this 
    at least doesn't require us to generate an entirely new model 
    input.
    """
    host_guest_elber_model_input_obj = copy.deepcopy(
        host_guest_elber_model_input_persistent)
    return host_guest_elber_model_input_obj

@pytest.fixture(scope="session")
def host_guest_elber_model_persistent(host_guest_elber_model_input_persistent):
    """
    Create a model object that is persistent across the tests in this file.
    """
    os.chdir(TEST_DIRECTORY)
    host_guest_elber_model_obj, model_xml_path \
        = prepare.prepare(host_guest_elber_model_input_persistent, 
                          force_overwrite=False)
    model_dir = os.path.dirname(model_xml_path)
    host_guest_elber_model_obj.anchor_rootdir = os.path.abspath(model_dir)
    return host_guest_elber_model_obj

@pytest.fixture
def host_guest_elber_model(host_guest_elber_model_persistent):
    """
    Create a copy of the model that is not persistent. But this at least
    doesn't require us to generate an entirely new model
    """
    host_guest_elber_model = copy.deepcopy(host_guest_elber_model_persistent)
    return host_guest_elber_model

@pytest.fixture(scope="session")
def tiwary_mmvt_model_input_persistent(tmpdir_factory):
    """
    Create a model object that is persistent across the tests in this file.
    """
    rootdir = tmpdir_factory.mktemp("tiwary_mmvt")
    host_guest_mmvt_model_input_persisent_obj \
        = create_model_input.create_tiwary_mmvt_model_input(rootdir)
    return host_guest_mmvt_model_input_persisent_obj

@pytest.fixture()
def tiwary_mmvt_model_input(tiwary_mmvt_model_input_persistent):
    """
    Create a copy of the model input that is not persistent. But this 
    at least doesn't require us to generate an entirely new model 
    input.
    """
    tiwary_mmvt_model_input_obj = copy.deepcopy(
        tiwary_mmvt_model_input_persistent)
    return tiwary_mmvt_model_input_obj

@pytest.fixture(scope="session")
def tiwary_mmvt_model_persistent(tiwary_mmvt_model_input_persistent):
    """
    Create a model object that is persistent across the tests in this file.
    """
    os.chdir(TEST_DIRECTORY)
    tiwary_mmvt_model_obj, model_xml_path \
        = prepare.prepare(tiwary_mmvt_model_input_persistent, 
                          force_overwrite=False)
    model_dir = os.path.dirname(model_xml_path)
    tiwary_mmvt_model_obj.anchor_rootdir = os.path.abspath(model_dir)
    return tiwary_mmvt_model_obj

@pytest.fixture
def tiwary_mmvt_model(tiwary_mmvt_model_persistent):
    """
    Create a copy of the model that is not persistent. But this at least
    doesn't require us to generate an entirely new model
    """
    tiwary_mmvt_model = copy.deepcopy(tiwary_mmvt_model_persistent)
    return tiwary_mmvt_model

@pytest.fixture(scope="session")
def planar_mmvt_model_input_persistent(tmpdir_factory):
    """
    Create a model object that is persistent across the tests in this file.
    """
    rootdir = tmpdir_factory.mktemp("planar_mmvt")
    planar_mmvt_model_input_persisent_obj \
        = create_model_input.create_planar_mmvt_model_input(rootdir)
    return planar_mmvt_model_input_persisent_obj

@pytest.fixture()
def planar_mmvt_model_input(planar_mmvt_model_input_persistent):
    """
    Create a copy of the model input that is not persistent. But this 
    at least doesn't require us to generate an entirely new model 
    input.
    """
    planar_mmvt_model_input_obj = copy.deepcopy(
        planar_mmvt_model_input_persistent)
    return planar_mmvt_model_input_obj

@pytest.fixture(scope="session")
def planar_mmvt_model_persistent(planar_mmvt_model_input_persistent):
    """
    Create a model object that is persistent across the tests in this file.
    """
    os.chdir(TEST_DIRECTORY)
    planar_mmvt_model_obj, model_xml_path \
        = prepare.prepare(planar_mmvt_model_input_persistent, 
                          force_overwrite=False)
    model_dir = os.path.dirname(model_xml_path)
    planar_mmvt_model_obj.anchor_rootdir = os.path.abspath(model_dir)
    planar_mmvt_model_obj.openmm_settings.cuda_platform_settings = None
    planar_mmvt_model_obj.openmm_settings.reference_platform = True
    return planar_mmvt_model_obj

@pytest.fixture
def planar_mmvt_model(planar_mmvt_model_persistent):
    """
    Create a copy of the model that is not persistent. But this at least
    doesn't require us to generate an entirely new model
    """
    planar_mmvt_model = copy.deepcopy(planar_mmvt_model_persistent)
    return planar_mmvt_model

@pytest.fixture(scope="session")
def rmsd_mmvt_model_input_persistent(tmpdir_factory):
    """
    Create a model object that is persistent across the tests in this file.
    """
    rootdir = tmpdir_factory.mktemp("rmsd_mmvt")
    rmsd_mmvt_model_input_persisent_obj \
        = create_model_input.create_rmsd_mmvt_model_input(rootdir)
    return rmsd_mmvt_model_input_persisent_obj

@pytest.fixture()
def rmsd_mmvt_model_input(rmsd_mmvt_model_input_persistent):
    """
    Create a copy of the model input that is not persistent. But this 
    at least doesn't require us to generate an entirely new model 
    input.
    """
    rmsd_mmvt_model_input_obj = copy.deepcopy(
        rmsd_mmvt_model_input_persistent)
    return rmsd_mmvt_model_input_obj

@pytest.fixture(scope="session")
def rmsd_mmvt_model_persistent(rmsd_mmvt_model_input_persistent):
    """
    Create a model object that is persistent across the tests in this file.
    """
    os.chdir(TEST_DIRECTORY)
    rmsd_mmvt_model_obj, model_xml_path \
        = prepare.prepare(rmsd_mmvt_model_input_persistent, 
                          force_overwrite=False)
    model_dir = os.path.dirname(model_xml_path)
    rmsd_mmvt_model_obj.anchor_rootdir = os.path.abspath(model_dir)
    rmsd_mmvt_model_obj.openmm_settings.cuda_platform_settings = None
    rmsd_mmvt_model_obj.openmm_settings.reference_platform = True
    return rmsd_mmvt_model_obj

@pytest.fixture
def rmsd_mmvt_model(rmsd_mmvt_model_persistent):
    """
    Create a copy of the model that is not persistent. But this at least
    doesn't require us to generate an entirely new model
    """
    rmsd_mmvt_model = copy.deepcopy(rmsd_mmvt_model_persistent)
    return rmsd_mmvt_model

@pytest.fixture(scope="session")
def toy_mmvt_model_input_persistent(tmpdir_factory):
    """
    Create a model object that is persistent across the tests in this file.
    """
    rootdir = tmpdir_factory.mktemp("toy_mmvt")
    toy_mmvt_model_input_persisent_obj \
        = create_model_input.create_toy_mmvt_model_input(rootdir)
    return toy_mmvt_model_input_persisent_obj

@pytest.fixture()
def toy_mmvt_model_input(toy_mmvt_model_input_persistent):
    """
    Create a copy of the model input that is not persistent. But this 
    at least doesn't require us to generate an entirely new model 
    input.
    """
    toy_mmvt_model_input_obj = copy.deepcopy(
        toy_mmvt_model_input_persistent)
    return toy_mmvt_model_input_obj

@pytest.fixture(scope="session")
def toy_mmvt_model_persistent(toy_mmvt_model_input_persistent):
    """
    Create a model object that is persistent across the tests in this file.
    """
    os.chdir(TEST_DIRECTORY)
    toy_mmvt_model_obj, model_xml_path \
        = prepare.prepare(toy_mmvt_model_input_persistent, 
                          force_overwrite=False)
    model_dir = os.path.dirname(model_xml_path)
    toy_mmvt_model_obj.anchor_rootdir = os.path.abspath(model_dir)
    return toy_mmvt_model_obj

@pytest.fixture
def toy_mmvt_model(toy_mmvt_model_persistent):
    """
    Create a copy of the model that is not persistent. But this at least
    doesn't require us to generate an entirely new model
    """
    toy_mmvt_model = copy.deepcopy(toy_mmvt_model_persistent)
    return toy_mmvt_model

@pytest.fixture(scope="session")
def toy_elber_model_input_persistent(tmpdir_factory):
    """
    Create a model object that is persistent across the tests in this file.
    """
    rootdir = tmpdir_factory.mktemp("toy_mmvt")
    toy_elber_model_input_persisent_obj \
        = create_model_input.create_toy_elber_model_input(rootdir)
    return toy_elber_model_input_persisent_obj

@pytest.fixture()
def toy_elber_model_input(toy_elber_model_input_persistent):
    """
    Create a copy of the model input that is not persistent. But this 
    at least doesn't require us to generate an entirely new model 
    input.
    """
    toy_elber_model_input_obj = copy.deepcopy(
        toy_elber_model_input_persistent)
    return toy_elber_model_input_obj

@pytest.fixture(scope="session")
def toy_elber_model_persistent(toy_elber_model_input_persistent):
    """
    Create a model object that is persistent across the tests in this file.
    """
    os.chdir(TEST_DIRECTORY)
    toy_elber_model_obj, model_xml_path \
        = prepare.prepare(toy_elber_model_input_persistent, 
                          force_overwrite=False)
    model_dir = os.path.dirname(model_xml_path)
    toy_elber_model_obj.anchor_rootdir = os.path.abspath(model_dir)
    return toy_elber_model_obj

@pytest.fixture
def toy_elber_model(toy_elber_model_persistent):
    """
    Create a copy of the model that is not persistent. But this at least
    doesn't require us to generate an entirely new model
    """
    toy_elber_model = copy.deepcopy(toy_elber_model_persistent)
    return toy_elber_model

@pytest.fixture(scope="session")
def toy_multi_model_input_persistent(tmpdir_factory):
    """
    Create a model object that is persistent across the tests in this file.
    """
    rootdir = tmpdir_factory.mktemp("toy_multi")
    toy_multi_model_input_persisent_obj \
        = create_model_input.create_toy_multi_model_input(rootdir)
    return toy_multi_model_input_persisent_obj

@pytest.fixture()
def toy_multi_model_input(toy_multi_model_input_persistent):
    """
    Create a copy of the model input that is not persistent. But this 
    at least doesn't require us to generate an entirely new model 
    input.
    """
    toy_multi_model_input_obj = copy.deepcopy(
        toy_multi_model_input_persistent)
    return toy_multi_model_input_obj

@pytest.fixture(scope="session")
def toy_multi_model_persistent(toy_multi_model_input_persistent):
    """
    Create a model object that is persistent across the tests in this file.
    """
    os.chdir(TEST_DIRECTORY)
    toy_multi_model_obj, model_xml_path \
        = prepare.prepare(toy_multi_model_input_persistent, 
                          force_overwrite=False)
    model_dir = os.path.dirname(model_xml_path)
    toy_multi_model_obj.anchor_rootdir = os.path.abspath(model_dir)
    toy_multi_model_obj.anchors[0].starting_positions = np.array([[[-0.7, -0.7, 0.0]]])
    
    return toy_multi_model_obj

@pytest.fixture
def toy_multi_model(toy_multi_model_persistent):
    """
    Create a copy of the model that is not persistent. But this at least
    doesn't require us to generate an entirely new model
    """
    toy_multi_model = copy.deepcopy(toy_multi_model_persistent)
    return toy_multi_model


def compare_dicts(dict1, dict2):
    """
    Compare the values within two dictionaries and assert they are
    close.
    """
    for key1 in dict1:
        assert key1 in dict2
        assert np.isclose(dict1[key1], dict2[key1])
    return