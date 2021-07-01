Running Calculations in SEEKR2
==============================

This page provides a cursory overview of a typical SEEKR2 calculation.

SEEKR2 calculations have three distinct stages:

1. Preparation
2. Running
3. Analysis

These are the main components, although there can also be more optional
components to a SEEKR calculation.

Appropriately, in SEEKR2, the programs that accomplish each of these stages 
are named:

1. prepare.py
2. run.py
3. analyze.py

Preparation
-----------

In order to run the preparation stage, SEEKR2 needs a "model input" XML file. 
Some examples of input files can be found in the following locations::

  seekr2/seekr2/data/sample_input_mmvt_openmm.xml
  seekr2/seekr2/data/sample_input_mmvt_namd.xml
  seekr2/seekr2/data/sample_input_elber_openmm.xml
  seekr2/seekr2/data/trypsin_benzamidine_files/input_tryp_ben_mmvt.xml
  seekr2/seekr2/data/trypsin_benzamidine_files/input_tryp_ben_elber.xml
  
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

The prepare.py script typically takes one argument - the input XML file. For
example:

``python prepare.py data/sample_input_mmvt_openmm.xml``

The script will construct a directory and file tree located at the path 
specified by the <root_directory> tag in the model input XML file.

More arguments to prepare.py can be found by running prepare.py with the "-h" 
argument, or by consulting the :doc:`Program Options<program_options>` 
documentation.