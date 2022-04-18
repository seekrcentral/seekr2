Running Calculations in SEEKR2
==============================

This page provides a cursory overview of a typical SEEKR2 calculation.

SEEKR2 calculations have three distinct stages:

1. Prepare
2. Run
3. Analyze

These are the main components, although there can also be more optional
components to a SEEKR calculation.

Appropriately, in SEEKR2, the programs that accomplish each of these stages 
are named:

1. prepare.py
2. run.py
3. analyze.py

Prepare
-------

In order to run the prepare stage, SEEKR2 typically needs a "model input" XML 
file. Some examples of input files can be found in the following locations::

  seekr2/seekr2/data/sample_input_mmvt_openmm.xml
  seekr2/seekr2/data/sample_input_mmvt_namd.xml
  seekr2/seekr2/data/sample_input_elber_openmm.xml
  
Inside these files are entries like temperature, type of calculations, system
parameters, and starting atomic positions. This file also defines the 
milestone locations and shapes. More information about the structure of
model input files can be found in the 
:doc:`Model Input Files<model_input_files>` documentation.

In addition, a number of input parameter, topology, and structure files must
be provided which provide the simulation engines with forcefield information,
starting atomic positions for each of the anchors (for MD), and waterbox 
dimensions. Examples of such files can also be found in the seekr2/seekr2/data/
directory::

  seekr2/seekr2/data/hostguest_files

The prepare.py script typically takes one argument - the model input XML file. 
For example::

  python prepare.py data/sample_input_mmvt_openmm.xml

The script will construct a directory and file tree located at the path 
specified by the <root_directory> tag in the model input XML file.

More arguments to prepare.py can be found by running prepare.py with the "-h" 
argument, or by consulting the :doc:`Program Options<program_options>` 
documentation.

Additionally, more systems and sample input files can be found in the 
seekr2_systems repository: https://github.com/seekrcentral/seekr2_systems.git.

Run
---

The prepare stage will generate a filetree of all necessary files and 
directories at the location specified by the <root_directory> tag of the
model input file. Within this directory is a file named *model.xml*. The
path to this file will be the argument of most programs used in SEEKR2
in the run stage and beyond.

*NOTE: the model.xml file is not intended to be modified directly without
re-running prepare.py on the model input file.*

The run stage is started with the run.py script::

  python run.py any /path/to/root_directory/model.xml
  
The word "any" is the INSTRUCTION argument to run.py, and instructs run.py to
run any unfinished MD or BD simulations. One may use "any_md" or "any_bd" as 
inputs to INSTRUCTION if one wishes to simulation only unfinished MD or 
unfinished BD simulations, respectively. There are many more inputs to 
INSTRUCTION and arguments to run.py, and they are all listed in 
:doc:`Program Options<program_options>`.

Following the run.py command, simulations will run until completion or 
interruption. SEEKR saves checkpoints for MD and BD simulations, so if an 
interruption occurs, one may re-enter the run.py command and the calculation
will pick up where it left off.

One may check on the progress and convergence of the calculation with the
converge.py program::

  python converge.py /path/to/root_directory/model.xml

The program converge.py will also generate convergence plots and images,
which may be found in the "plots_and_images" subfolder of <root_directory>.

More arguments to both run.py and converge.py can be found by running 
either program with the "-h" argument, or by consulting the 
:doc:`Program Options<program_options>` documentation.

Analyze
-------

The final stage of a typical SEEKR2 calculation involves analyzing the results
of the simulations to construct the kinetics and thermodynamics of the 
process in question.

The analyze stage is launches with the analyze.py script::

  python analyze.py /path/to/root_directory/model.xml

The analyze.py script constructs the milestoning model, populates it with the
probabilities and times from the simulations, and computes other quantities, 
like the error margins.

Upon completion, interesting quantities will be printed to the screen, and
more images will be saved in the "plots_and_images" subfolder of 
<root_directory>, including a free energy profile (potential of mean force) of
each milestone.