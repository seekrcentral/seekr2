Prepare.py Inputs and Arguments
===============================

The prepare.py program takes a :doc:`Model Input File<model_input_files>` as 
input and generates the files necessary to run a SEEKR2 calculation.

usage::

  python prepare.py [-h] [-f] [-s] MODEL_INPUT_FILE
  
Required Arguments
------------------

**MODEL_INPUT_FILE**
  The name of a model input XML file for a SEEKR2 calculation. Samples model 
  files can be found in the following locations::

    seekr2/seekr2/data/sample_input_mmvt_openmm.xml
    seekr2/seekr2/data/sample_input_mmvt_namd.xml
    seekr2/seekr2/data/sample_input_elber_openmm.xml
    seekr2/seekr2/data/trypsin_benzamidine_files/input_tryp_ben_mmvt.xml
    seekr2/seekr2/data/trypsin_benzamidine_files/input_tryp_ben_elber.xml
  
  Also, an :doc:`entire page of documentation<model_input_files>` is dedicated 
  to explaining the inputs of the model input file.

Optional Arguments
------------------

**-h, --help**
  show help message and exit

**-f, --force_overwrite**
  Toggle whether to overwrite existing simulation 
  output files within any anchor that might have existed in an old model that 
  would be overwritten by generating this new model. If not toggled while
  attempting to overwrite an old model with a new model, this program will 
  throw an exception instead of performing any such overwrite.

**-s, --skip_checks**
  By default, pre-simulation checks will be run after 
  the preparation is complete, and if the checks fail, prepare.py will terminate
  and the SEEKR2 model will not be saved. This argument bypasses those checks 
  and allows the model to be generated regardless of check outcomes.
