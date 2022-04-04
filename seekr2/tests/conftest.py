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
def smoluchowski_mmvt_model_input_persistent(tmpdir_factory):
    """
    Create a model object that is persistent across the tests in this file.
    """
    rootdir = tmpdir_factory.mktemp("smoluchowski_mmvt")
    smoluchowski_mmvt_model_input_persisent_obj \
        = create_model_input.create_smoluchowski_mmvt_model_input(
            rootdir)
    return smoluchowski_mmvt_model_input_persisent_obj

@pytest.fixture()
def smoluchowski_mmvt_model_input(smoluchowski_mmvt_model_input_persistent):
    """
    Create a copy of the model input that is not persistent. But this 
    at least doesn't require us to generate an entirely new model 
    input.
    """
    smoluchowski_mmvt_model_input_obj = copy.deepcopy(
        smoluchowski_mmvt_model_input_persistent)
    return smoluchowski_mmvt_model_input_obj

@pytest.fixture(scope="session")
def smoluchowski_mmvt_model_persistent(
        smoluchowski_mmvt_model_input_persistent):
    """
    Create a model object that is persistent across the tests in this file.
    """
    os.chdir(TEST_DIRECTORY)
    smoluchowski_mmvt_model_obj, model_xml_path \
        = prepare.prepare(smoluchowski_mmvt_model_input_persistent, 
                          force_overwrite=False)
    model_dir = os.path.dirname(model_xml_path)
    smoluchowski_mmvt_model_obj.anchor_rootdir = os.path.abspath(model_dir)
    # make BD milestones
    smoluchowski_mmvt_model_obj.k_on_info = base.K_on_info()
    bd_milestone = base.BD_milestone()
    bd_milestone.index = 0
    bd_milestone.outer_milestone \
        = smoluchowski_mmvt_model_obj.anchors[-1].milestones[0]
    bd_milestone.inner_milestone \
        = smoluchowski_mmvt_model_obj.anchors[-2].milestones[0]
    smoluchowski_mmvt_model_obj.k_on_info.bd_milestones.append(bd_milestone)
    smoluchowski_mmvt_model_obj.browndye_settings = base.Browndye_settings()
    
    return smoluchowski_mmvt_model_obj

@pytest.fixture
def smoluchowski_mmvt_model(smoluchowski_mmvt_model_persistent):
    """
    Create a copy of the model that is not persistent. But this at least
    doesn't require us to generate an entirely new model
    """
    smoluchowski_mmvt_model = copy.deepcopy(smoluchowski_mmvt_model_persistent)
    return smoluchowski_mmvt_model

@pytest.fixture(scope="session")
def smoluchowski_elber_model_input_persistent(tmpdir_factory):
    """
    Create a model object that is persistent across the tests in this file.
    """
    rootdir = tmpdir_factory.mktemp("smoluchowski_elber")
    smoluchowski_elber_model_input_persistent_obj \
        = create_model_input.create_smoluchowski_elber_model_input(
            rootdir)
    return smoluchowski_elber_model_input_persistent_obj

@pytest.fixture(scope="session")
def smoluchowski_elber_model_persistent(
        smoluchowski_elber_model_input_persistent):
    """
    Create a model object that is persistent across the tests in this file.
    """
    os.chdir(TEST_DIRECTORY)
    smoluchowski_elber_model_obj, model_xml_path \
        = prepare.prepare(smoluchowski_elber_model_input_persistent, 
                          force_overwrite=False)
    model_dir = os.path.dirname(model_xml_path)
    smoluchowski_elber_model_obj.anchor_rootdir = os.path.abspath(model_dir)
    smoluchowski_elber_model_obj.k_on_info = base.K_on_info()
    bd_milestone = base.BD_milestone()
    bd_milestone.index = 0
    bd_milestone.outer_milestone \
        = smoluchowski_elber_model_obj.anchors[-1].milestones[1]
    bd_milestone.inner_milestone \
        = smoluchowski_elber_model_obj.anchors[-2].milestones[1]
    smoluchowski_elber_model_obj.k_on_info.bd_milestones.append(bd_milestone)
    smoluchowski_elber_model_obj.browndye_settings = base.Browndye_settings()
    return smoluchowski_elber_model_obj

@pytest.fixture
def smoluchowski_elber_model(smoluchowski_elber_model_persistent):
    """
    Create a copy of the model that is not persistent. But this at least
    doesn't require us to generate an entirely new model
    """
    smoluchowski_elber_model_obj = copy.deepcopy(
        smoluchowski_elber_model_persistent)
    return smoluchowski_elber_model_obj

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

def compare_dicts(dict1, dict2):
    """
    Compare the values within two dictionaries and assert they are
    close.
    """
    for key1 in dict1:
        assert key1 in dict2
        assert np.isclose(dict1[key1], dict2[key1])
    return