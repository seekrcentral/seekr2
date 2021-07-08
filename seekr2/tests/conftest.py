"""
conftest.py

configurations for all tests
"""

import os
import tempfile
import pytest
import copy

import seekr2.prepare as prepare
import seekr2.tests.create_model_input as create_model_input

TEST_DIRECTORY = os.path.dirname(__file__)

@pytest.fixture(scope="session")
def host_guest_mmvt_model_input_persisent(tmpdir_factory):
    """
    Create a model object that is persistent across the tests in this file.
    """
    rootdir = tmpdir_factory.mktemp("hostguest_mmvt")
    host_guest_mmvt_model_input_persisent_obj \
        = create_model_input.create_tryp_ben_mmvt_model_input(
            rootdir, bd=True)
    return host_guest_mmvt_model_input_persisent_obj

@pytest.fixture(scope="session")
def host_guest_mmvt_model_persistent(tmpdir_factory, host_guest_model_input):
    """
    Create a model object that is persistent across the tests in this file.
    """
    os.chdir(TEST_DIRECTORY)
    host_guest_mmvt_model_obj, model_xml_path \
        = prepare.generate_seekr2_model_and_filetree(
            host_guest_mmvt_model_input_persisent, force_overwrite=False)
    model_dir = os.path.dirname(model_xml_path)
    host_guest_mmvt_model_obj.anchor_rootdir = os.path.abspath(model_dir)
    return host_guest_mmvt_model_obj

@pytest.fixture
def host_guest_mmvt_model(tmpdir_factory, host_guest_mmvt_model_persistent):
    """
    Create a copy of the model that is not persistent. But this at least
    doesn't require us to generate an entirely new model
    """
    host_guest_mmvt_model = copy.deepcopy(host_guest_mmvt_model_persistent)
    return host_guest_mmvt_model

@pytest.fixture(scope="session")
def tryp_ben_mmvt_model_input_persistent(tmpdir_factory):
    """
    Create a model object that is persistent across the tests in this file.
    """
    rootdir = tmpdir_factory.mktemp("tryp_ben_mmvt")
    tryp_ben_mmvt_model_input_persistent_obj \
        = create_model_input.create_tryp_ben_mmvt_model_input(
            rootdir, bd=True)
    return tryp_ben_mmvt_model_input_persistent_obj

@pytest.fixture(scope="session")
def tryp_ben_mmvt_model_persistent(tmpdir_factory, 
                                   tryp_ben_mmvt_model_input_persistent):
    """
    Create a model object that is persistent across the tests in this file.
    """
    os.chdir(TEST_DIRECTORY)
    print("mmvt curdir:", os.getcwd())
    tryp_ben_mmvt_model_obj, model_xml_path \
        = prepare.generate_seekr2_model_and_filetree(
            tryp_ben_mmvt_model_input_persistent, force_overwrite=False)
    model_dir = os.path.dirname(model_xml_path)
    tryp_ben_mmvt_model_obj.anchor_rootdir = os.path.abspath(model_dir)
    return tryp_ben_mmvt_model_obj

@pytest.fixture
def tryp_ben_mmvt_model(tmpdir_factory, tryp_ben_mmvt_model_persistent):
    """
    Create a copy of the model that is not persistent. But this at least
    doesn't require us to generate an entirely new model
    """
    tryp_ben_mmvt_model_obj = copy.deepcopy(tryp_ben_mmvt_model_persistent)
    return tryp_ben_mmvt_model_obj

@pytest.fixture(scope="session")
def tryp_ben_elber_model_input_persistent(tmpdir_factory):
    """
    Create a model object that is persistent across the tests in this file.
    """
    rootdir = tmpdir_factory.mktemp("tryp_ben_elber")
    tryp_ben_elber_model_input_persistent_obj \
        = create_model_input.create_tryp_ben_elber_model_input(
            rootdir, bd=True)
    return tryp_ben_elber_model_input_persistent_obj

@pytest.fixture(scope="session")
def tryp_ben_elber_model_persistent(tmpdir_factory, 
                                   tryp_ben_elber_model_input_persistent):
    """
    Create a model object that is persistent across the tests in this file.
    """
    os.chdir(TEST_DIRECTORY)
    print("elber curdir:", os.getcwd())
    tryp_ben_elber_model_obj, model_xml_path \
        = prepare.generate_seekr2_model_and_filetree(
            tryp_ben_elber_model_input_persistent, force_overwrite=False)
    model_dir = os.path.dirname(model_xml_path)
    tryp_ben_elber_model_obj.anchor_rootdir = os.path.abspath(model_dir)
    return tryp_ben_elber_model_obj

@pytest.fixture
def tryp_ben_elber_model(tmpdir_factory, tryp_ben_elber_model_persistent):
    """
    Create a copy of the model that is not persistent. But this at least
    doesn't require us to generate an entirely new model
    """
    tryp_ben_elber_model_obj = copy.deepcopy(tryp_ben_elber_model_persistent)
    return tryp_ben_elber_model_obj