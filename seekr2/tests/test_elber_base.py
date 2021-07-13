"""
test_elber_base.py

Testing elber base objects
"""

def test_create_elber_model_input(tryp_ben_elber_model_input_persistent):
    """
    Just create the model input to see if it can be done without errors
    """
    assert tryp_ben_elber_model_input_persistent.calculation_type == "elber"
    return

def test_create_elber_model(tryp_ben_elber_model):
    """
    Create the model based on the model input to ensure that it works
    without errors.
    """
    assert tryp_ben_elber_model.calculation_type == "elber"
    return