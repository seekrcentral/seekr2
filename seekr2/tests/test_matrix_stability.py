"""
test_matrix_stability.py

For systems with very slow kinetics, the conventional matrix solvers seem
to have trouble computing first passage times due to poor matrix conditioning.
Recreates the problem and test matrix solutions that will remain stable for
milestoning matrices with long timescales.
"""

import os
import pytest

import numpy as np

import seekr2.prepare as prepare
import seekr2.toy.smoluchowski as smol
import seekr2.tests.create_model_input as create_model_input
import seekr2.tests.test_analyze as test_analyze

def make_unstable_matrix_smoluchowski_mmvt_model(rootdir, num_input_anchors):
    """
    Create a milestoning model that will result in an unstable matrix
    """
    
    smol_input = create_model_input.create_smoluchowski_mmvt_model_input(
        rootdir, num_input_anchors+1)
    os.chdir(rootdir)
    smol_model, model_xml_path = prepare.prepare(
        smol_input, force_overwrite=False)
    model_dir = os.path.dirname(model_xml_path)
    smol_model.anchor_rootdir = os.path.abspath(model_dir)
    return smol_model

def make_mmvt_smol_model(rootdir, num_input_anchors):
    """
    Obtain two Smoluchowski models: One as an exact system to act as a 
    benchmark, and the second as a milestoning system with an unstable
    matrix.
    """
    model = make_unstable_matrix_smoluchowski_mmvt_model(rootdir, num_input_anchors)
    return model

def make_unstable_matrix_smoluchowski_elber_model(rootdir, num_input_anchors):
    """
    Create a milestoning model that will result in an unstable matrix
    """
    
    smol_input = create_model_input.create_smoluchowski_elber_model_input(
        rootdir, num_input_anchors+1)
    os.chdir(rootdir)
    smol_model, model_xml_path = prepare.prepare(
        smol_input, force_overwrite=False)
    model_dir = os.path.dirname(model_xml_path)
    smol_model.anchor_rootdir = os.path.abspath(model_dir)
    return smol_model

def make_elber_smol_model(rootdir, num_input_anchors):
    """
    Obtain two Smoluchowski models: One as an exact system to act as a 
    benchmark, and the second as a milestoning system with an unstable
    matrix.
    """
    model = make_unstable_matrix_smoluchowski_elber_model(rootdir, num_input_anchors)
    return model

def get_benchmark_time(potential_energy_function, num_input_anchors):
    """
    Find the exact solution for a Smoluchowski system for comparison with
    matrix results.
    """
    milestones = [1.0]
    absorbing_boundary = num_input_anchors
    calc = smol.SmoluchowskiCalculation1d(
        potential_energy_function, milestones=milestones, 
        absorbing_boundary=absorbing_boundary, n=401)
    outer_surface_flux = 4.0 * np.pi * absorbing_boundary**2 \
        * calc.regions[1].J_outward
    w_a = calc.potential_energy_function.evaluate_energy(calc.regions[1].a)
    inner_region_height = smol.expBetaW(w_a, calc.beta)
    outer_region_height = calc.regions[1].u_r_outward[0]
    inner_region_volume = calc.regions[0].partition_function \
        * (outer_region_height / inner_region_height)
    outer_region_volume = 1.0
    time = (inner_region_volume + outer_region_volume) / outer_surface_flux
    return time

def make_analytic_pi_alpha(num_regions, potential_energy_function):
    milestones = list(range(1, num_regions))
    absorbing_boundary = num_regions
    calc = smol.SmoluchowskiCalculation1d(
        potential_energy_function, milestones=milestones, 
        absorbing_boundary=absorbing_boundary)
    pi_alpha = np.zeros((1, len(calc.regions)))
    for i, region in enumerate(calc.regions):
        pi_alpha[0,i] = region.partition_function / calc.partition_function
    return pi_alpha

def compare_solvers_analytical_calculated(tmpdir, potential_energy_function):
    lowest = 5 # 55
    highest = 10 # 60
    interval = 5
    pre_equilibrium_approx = False
    print("i   benchmark_time  mmvt_time")
    for i in range(lowest, highest+interval, interval):
        mmvt_rootdir = os.path.join(
            tmpdir, "smol_mmvt_matrix_problem_{}".format(i))
        os.mkdir(mmvt_rootdir)
        elber_rootdir = os.path.join(
            tmpdir, "smol_elber_matrix_problem_{}".format(i))
        os.mkdir(elber_rootdir)
        num_input_anchors = i
        mmvt_model = make_mmvt_smol_model(mmvt_rootdir, num_input_anchors)
        analytic_pi_alpha = make_analytic_pi_alpha(i, potential_energy_function)
        elber_model = make_elber_smol_model(elber_rootdir, num_input_anchors)
        benchmark_time = get_benchmark_time(
            potential_energy_function, num_input_anchors)
        mmvt_model.k_on_info = None
        mmvt_time, dummy, mmvt_analysis = \
            test_analyze.make_smoluchowski_mmvt_analysis(
                mmvt_model, potential_energy_function, 
                pre_equilibrium_approx=pre_equilibrium_approx)
        elber_model.k_on_info = None
        #elber_time, dummy, elber_analysis = \
        #    test_analyze.make_smoluchowski_elber_analysis(
        #        elber_model, potential_energy_function)
        elber_time = 0.0
        
        #print(i, benchmark_time, elber_time, mmvt_time)
        assert np.isclose(benchmark_time, mmvt_time, rtol=1e-2)
    
    return


def test_make_method_comparison_data_for_plotting(tmpdir):
    potential_energy_function = smol.LinearPotentialEnergyFunction(a=1.0)
    compare_solvers_analytical_calculated(tmpdir, potential_energy_function)

@pytest.mark.skip()
def test_scipy_matrix_solver(tmpdir):
    #potential_energy_function = smol.FlatPotentialEnergyFunction()
    potential_energy_function = smol.LinearPotentialEnergyFunction(a=1.0)
    rootdir = os.path.join(tmpdir, "smol_matrix_problem")
    os.mkdir(rootdir)
    num_input_anchors = 5
    #model = make_mmvt_smol_model(rootdir, num_input_anchors)
    model = make_elber_smol_model(rootdir, num_input_anchors)
    benchmark_time = get_benchmark_time(
        potential_energy_function, num_input_anchors)
    #print("benchmark_time:", benchmark_time)
    model.k_on_info = None
    #mmvt_time, dummy1, my_analysis_mmvt = test_analyze.make_smoluchowski_mmvt_analysis(
    #    model, potential_energy_function)
    elber_time, dummy2, my_analysis_elber = test_analyze.make_smoluchowski_elber_analysis(
        model, potential_energy_function)
    #print(benchmark_time, elber_time, mmvt_time)
    assert np.isclose(benchmark_time, elber_time, rtol=1e-2)
    