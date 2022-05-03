"""
test_common_cv.py

Testing modules/common_cv.py
"""

import seekr2.modules.common_cv as common_cv

def test_assign_state_points_toy(toy_mmvt_model_input, toy_mmvt_model):
    stateA = common_cv.State_point()
    stateA.name = "stateA"
    stateA.location = [0.0, -0.7, 0.0]
    stateB = common_cv.State_point()
    stateB.name = "stateB"
    stateB.location = [0.0, 0.7, 0.0]
    toy_mmvt_model_input.cv_inputs[0].state_points = [stateA, stateB]
    common_cv.assign_state_points(toy_mmvt_model_input, toy_mmvt_model)
    assert toy_mmvt_model.anchors[0].name == "stateA"
    assert toy_mmvt_model.anchors[7].name == "stateB"
    return
    
def test_assign_state_points_hostguest(host_guest_mmvt_model_input, 
                                       host_guest_mmvt_model):
    bound_state = common_cv.State_point()
    bound_state.name = "bound"
    bound_state.location = 0.06
    bulk_state = common_cv.State_point()
    bulk_state.name = "bulk"
    bulk_state.location = 1.36
    host_guest_mmvt_model_input.cv_inputs[0].state_points \
        = [bound_state, bulk_state]
    common_cv.assign_state_points(host_guest_mmvt_model_input, 
                                  host_guest_mmvt_model)
    assert host_guest_mmvt_model.anchors[0].name == "bound"
    assert host_guest_mmvt_model.anchors[13].name == "bulk"
    return

