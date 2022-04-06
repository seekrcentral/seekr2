# Add tests here for the entropic barrier and compare time scales

import numpy as np

import seekr2.run as run
import seekr2.analyze as analyze

def test_entropy_barrier_timescale(toy_mmvt_model):
    """
    Test entropy barrier system for if it recreates reasonable
    kinetics timescales.
    """
    num_steps = 1000000
    long_timescale_residence_time_in_ps = 1263.06
    toy_mmvt_model.openmm_settings.cuda_platform_settings = None
    toy_mmvt_model.openmm_settings.reference_platform = True
    toy_mmvt_model.openmm_settings.langevin_integrator.friction_coefficient \
        = 100.0
    toy_mmvt_model.calculation_settings.num_production_steps = num_steps
    toy_mmvt_model.calculation_settings.energy_reporter_interval = num_steps
    run.run(toy_mmvt_model, "any")
    analysis = analyze.analyze(toy_mmvt_model, num_error_samples=100, 
                               skip_checks=False)
    assert np.isclose(analysis.MFPTs[('anchor_0', 'bulk')], 
                      long_timescale_residence_time_in_ps,
                      rtol=0.5)