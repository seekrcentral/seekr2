Analyze.py Inputs and Arguments
===============================

The analyze.py program takes the data generated from run.py, constructs a
milestoning model, and computes the kinetics and thermodynamics results.

usage::

  python analyze.py [-h] [-f] [-n NUM_ERROR_SAMPLES] [-p] [-d IMAGE_DIRECTORY] 
                  [-s] MDOEL_FILE
  
Required Arguments
------------------

**MODEL_FILE**
  Provide a path to a file named "model.xml" which can be found
  in the <root_directory> tag of the model input file provided to prepare.py.

Optional Arguments
------------------

**-h, --help**
  show help message and exit.

**-f, --force_warning**
  By default, missing statistics for any anchors will 
  generate fatal errors. Analyze.py will attempt to catch these situations early
  and abort gracefully. This option will instead raise a warning and attempt 
  the calculation anyway.

**-n NUM_ERROR_SAMPLES, --num_error_samples NUM_ERROR_SAMPLES**
  Specify the number of error samples to generate for estimating 
  error/uncertainty of computed values. Default: 1000.
  
**-S STRIDE_ERROR_SAMPLES, --stride_error_samples STRIDE_ERROR_SAMPLES**
  Specify the number of strides betweensaved error samples. An argument of 
  None automatically assigns the quantity at the number of milestones in the 
  model squared. Default: None
  
**-K SKIP_ERROR_SAMPLES, --skip_error_samples SKIP_ERROR_SAMPLES**
  Specify the number of error samples. An argument of None automatically 
  assigns the quantity at ten times the number of milestones in the model 
  squared. Default: None.

**-d IMAGE_DIRECTORY, --image_directory IMAGE_DIRECTORY**
  Define the directory where all plots and images will be saved. By default, 
  graphics will be saved to the 'images_and_plots/' directory in the model's 
  anchor root directory.

**-s, --skip_checks**
  By default, post-simulation checks will be run before 
  the analysis is started, and if the checks fail, the analysis will not proceed.
  This argument bypasses those checks and allows the analysis to proceed anyways.
  
**-t MINIMUM_TIME, --minimum_time MINIMUM_TIME**
  A user may wish to skip an amount of simulation time for each anchor before 
  counting the transitions for milestoning analysis. Enter the time (in ps) to 
  skip a portion of the production simulations when performing analysis.
  
**-T MAXIMUM_TIME, --maximum_time MAXIMUM_TIME**
  A user may wish to stop the analysis of simulation time for each anchor at 
  a particular time. Enter the time (in ps) at which to end the analysis at 
  a given anchor if the simulation time exceeds it.
