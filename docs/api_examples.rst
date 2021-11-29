API Examples
============

An example of a script to run a SEEKR2 calculation can be found in Listing 1. 
This script defines a “model input” object, loads its parameters from an XML 
file, generates the full model, checks the model for potential problems, runs 
any and all of the MD and BD simulations, and analyzes the results.  

**Listing 1:** a sample Python script using the SEEKR2 API to perform a full 
SEEKR calculation. NOTE: this script should be run from the seekr2/seekr2/
directory to work as-is.::

  import seekr2
  model_input_filename = “data/sample_input_mmvt_openmm.xml"
  model_input = seekr2.prepare.common_prepare.Model_input()
  model_input.deserialize(model_input_filename, user_input=True)
  model, xml_path = seekr2.prepare.prepare(model_input)
  model.anchor_rootdir = os.path.dirname(xml_path)
  seekr2.modules.check.check_pre_simulation_all(model)
  seekr2.run.run(model, "any")
  analysis = seekr2.analyze.analyze(model)
  analysis.print_results()
  
The objects can be accessed/modified from the Python API at any step in the 
process. For instance, if one wishes to add or move milestones, one could add 
the lines shown in Listing 2 before the model is prepared. Or say that one 
wants to access the milestoning matrix directly, one could add the lines in 
listing 3 after analysis object is generated.  

**Listing 2:** a an additional snippet of sample Python code to use the SEEKR2 
API for first add an anchor to the calculation and generate the model. Then, 
to demonstrate that the radius can be changed after the model is generated, 
the anchor’s radius is changed and the model is regenerated.::

  from openmm import unit
  new_input_anchor = seekr2.modules.common_cv.Spherical_cv_anchor()
  new_input_anchor.radius = 0.1 * unit.nanometer
  model_input.cv_inputs[0].input_anchors.insert(1, new_input_anchor)
  model, xml_path = seekr2.prepare.prepare(model_input, force_overwrite=True)
  model_input.cv_inputs[0].input_anchors[1].radius = 0.11 * unit.nanometers
  model, xml_path = seekr2.prepare.prepare(model_input, force_overwrite=True)
  
**Listing 3:** The milestoning rate matrix “Q” or free energy profile are just 
two examples of many quantities that may be accessed directly from the Python 
analysis object.::

  print(analysis.main_data_sample.Q)
  print(analysis.free_energy_profile)

This script is just intended to give developers a glimpse of what is possible 
when using the API. Please consult the documentation within the code itself 
for details about the arguments to each of these functions.