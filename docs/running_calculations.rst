Running Calculations in SEEKR2
==============================

SEEKR2 calculations have three distinct stages:

1. Preparation
2. Running
3. Analysis

These are the main components, although there can also be more optional
components to a SEEKR calculation.

Appropriately, the programs that accomplish each of these stages are named:

1. prepare.py
2. run.py
3. analyze.py

Preparation
-----------

In order to run the preparation stage, SEEKR2 needs an input XML file. Examples
of input files can be found in the following locations::

  seekr2/seekr2/data/sample_input_mmvt_openmm.xml
  seekr2/seekr2/data/sample_input_mmvt_namd.xml
  seekr2/seekr2/data/sample_input_elber_openmm.xml
  seekr2/seekr2/data/trypsin_benzamidine_files/input_tryp_ben_mmvt.xml
  seekr2/seekr2/data/trypsin_benzamidine_files/input_tryp_ben_elber.xml
  
Inside these files are entries like temperature, type of calculations, system
parameters, and starting atomic positions. This file also defines the 
milestone locations and shapes.

The prepare.py script typically takes one argument - the input XML file. For
example:

``python prepare.py data/sample_input_mmvt_openmm.xml``