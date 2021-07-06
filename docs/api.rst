API Documentation
=================

Normally, one would run SEEKR2 directly on the Linux terminal. However, it is 
possible to perform an entire SEEKR2 calculation, along with all of its
auxiliary functions, from a custom python program by using SEEKR2's Application
Programming Interface (API).

Below is a sample Python script that runs an entire SEEKR2 calculation through 
the API::
  
  import seekr2
  my_directory = "~/test_api"
  my_model_input \
      = seekr2.tests.create_model_input.create_host_guest_mmvt_model_input(
          my_directory)
  my_model_input.temperature = 293.15
  my_model_input.cv_inputs[0].input_anchors[2].radius = 0.275
  # ... more model input modifications here
  model, xml_path = seekr2.prepare.generate_seekr2_model_and_filetree(
      my_model_input)
  model_directory = os.path.dirname(xml_path)
  model.anchor_rootdir = os.path.abspath(model_directory)
  seekr2.check.check_pre_simulation_all(model)
  seekr2.run.run(model, "any", min_b_surface_simulation_length=1000, 
      cuda_device_index="0")
  seekr2.analysis = analyze.analyze(model)
  seekr2.analysis.print_results()
  
For your own calculations, you will probably want to create your own 
function equivalent to the
**seekr2.tests.create_model_input.create_host_guest_mmvt_model_input**
function above. Examine this function directly to see how it is created (the
**create_host_guest_mmvt_model_input** function is located in 
**seekr2/seekr2/tests/create_model_input.py**).
  
This script is just intended to get developers started with the API. Please
consult the documentation within the code itself for details about the 
arguments to each of these functions.