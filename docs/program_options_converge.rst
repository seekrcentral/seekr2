Converge.py Inputs and Arguments
================================

The converge.py program takes the data generated from run.py and provides the
user with information about simulation progress, including the lengths of 
simulations, the number of transitions (bounces) against milestones, and 
convergence information.

In addition, the run.py program can be configured to use the output of
converge.py's API to run simulations to convergence.

usage::

  python converge.py [-h] [-s K_ON_STATE] [-d IMAGE_DIRECTORY] [-c CUTOFF] 
                   [-m MINIMUM_ANCHOR_TRANSITIONS] [-p] MODEL_FILE

  
Required Arguments
------------------

**MODEL_FILE**
    Provide a path to a file named "model.xml" which can be found
    in the <root_directory> tag of the model input file provided to prepare.py.

Optional Arguments
------------------

-h, --help            show help message and exit.

-s K_ON_STATE, --k_on_state K_ON_STATE
                      Define the bound state (anchor index) that will be used 
                      to compute k-on. If left blank, and k-on statistics 
                      exist, then the first end state in the model will be 
                      chosen by default.

-d IMAGE_DIRECTORY, --image_directory IMAGE_DIRECTORY
                      Define the directory where all plots and images will 
                      be saved. By default, graphics will be saved to the 
                      'images_and_plots/' directory in the model's anchor root 
                      directory.

-c CUTOFF, --cutoff CUTOFF
                      The minimum convergence that must be achieved before 
                      concluding that the calculations have converged for a 
                      given anchor.

-m MINIMUM_ANCHOR_TRANSITIONS, --minimum_anchor_transitions MINIMUM_ANCHOR_TRANSITIONS
                       Enter a minimum number of transitions that must be 
                       observed per milestone in a given MD anchor as a 
                       criteria for the simulations to be considered converged.
