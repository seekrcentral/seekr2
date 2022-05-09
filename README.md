seekr2
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/seekrcentral/seekr2/workflows/CI/badge.svg)](https://github.com/seekrcentral/seekr2/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/seekrcentral/seekr2/branch/master/graph/badge.svg)](https://codecov.io/gh/seekrcentral/seekr2/branch/master)


[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![python](https://img.shields.io/badge/python-3.8-blue.svg)](https://www.python.org/)


```
##########################################      ########
   ###   ######  ######  ### ###  #####      ####      ####
  #####   ##  #   ##  #   ##  #    ##  #    ##           ###
  ##      ##      ##      ## #     ##  #     ##           ##
   ##     ####    ####    ###      ##  #                 ###
    ##    ##      ##      ###     #####                ####
  #####   ##  #   ##  #   ## #     ##  #            ###
   ###   ######  ######  ###  ##  ###  #         ## 
##########################################    ##
         #######   #  #      ###           ##            #
       ##       #########       ##       ##  ####    #####
         ###      #  #   #######        ##       #######

Simulation-Enabled Estimation of Kinetic Rates - Version 2
```

## Overview
**Fast** and **versatile** multiscale milestoning to compute molecular 
thermodynamics and kinetics.

Prepare and run milestoning calculations in the OpenMM, NAMD, and/or Browndye2
simulation engines to obtain the kinetics and thermodynamics of molecular 
processes such as: ligand-receptor binding and unbinding, membrane 
permeability, internal molecular dynamics, and many other potential 
situations.

This README is only a quickstart guide to get SEEKR2 up and running as soon as
possible. To see more detailed **instructions** and **tutorials**, please see 
https://seekr2.readthedocs.io/en/latest or the docs/ subfolder.

## Quick Install

### Dependencies
Many of the dependencies for SEEKR2 will be installed alongside SEEKR2, but
some must be installed separately, and are installed first, before SEEKR2.

#### OpenMM (recommended)

OpenMM is recommended for the molecular dynamics (MD) stage of SEEKR2. SEEKR2 
also needs the SEEKR2 OpenMM Plugin in order to use OpenMM for MD simulations.

The easiest, quickest way to install the SEEKR2 OpenMM Plugin is to use
Conda. If you don't already have Conda installed, Download Conda with 
Python version 3.8 from 
https://conda.io/en/latest/miniconda.html and run the downloaded script and 
fill out the prompts. 

With Conda working, install the SEEKR2 OpenMM Plugin:

```
conda install -c conda-forge seekr2_openmm_plugin
```

Alternatively, NAMD2 may be used for MD if desired. See the NAMD2 section
below for installation instructions.

#### Browndye2 (optional; required for k-on calculations)

SEEKR2 needs Browndye2 if Brownian dynamics (BD) simulations will be run 
(necessary for k-on calculations). Please see (https://browndye.ucsd.edu/) 
for Browndye2 installation instructions.

#### NAMD2 (optional, required if not using OpenMM)

If OpenMM is not desirable or available for the MD simulations, SEEKR2 may 
use NAMD2 in order to run MD simulations. NAMD2 is frequently already 
available on shared computing resources like clusters and supercomputers.
Not all SEEKR2 options may be available using NAMD2.

If you wish to install NAMD2 yourself, please see 
http://www.ks.uiuc.edu/Research/namd/ for instructions to install NAMD.

### Install SEEKR2

Once the dependencies are installed, we may install SEEKR2. First, clone this 
repository and install the package:

```
git clone https://github.com/seekrcentral/seekr2.git
cd seekr2
python setup.py install
```

### Testing SEEKR2 (Optional)
To test SEEKR2, run the following command in the seekr2/ directory:

```
python setup.py test
```

Additional continuous integration tests may be run from the Python scripts in 
the seekr2/seekr2/continuous_integration/ directory.

## Run

A SEEKR2 calculation needs a "Model Input File" to run. Several examples may
be found in seekr2/seekr2/data. Execute the following commands within the 
seekr2/seekr2/ directory to run a sample calculation on the host-guest system:

(It is assumed that both Browndye2 and the the SEEKR2 OpenMM Plugin have
already been installed).

```
python prepare.py data/sample_input_mmvt_openmm.xml
python run.py any ~/test_mmvt_openmm/model.xml
python converge.py ~/test_mmvt_openmm/model.xml
python analyze.py ~/test_mmvt_openmm/model.xml
```

### Important Options and Hints

* In general, SEEKR2 programs can be run with the '-h' argument to see all
available options. Please see https://seekr2.readthedocs.io/en/latest for a
detailed description of programs and options.

* the **run.py** program accepts an INSTRUCTION as its first argument, and
the following options are available:
  * "any" - run any MD or BD calculation that still needs to be run.
  * "any_md" - run any MD calculation that still needs to be run.
  * "any_bd" - run any BD calculation that still needs to be run.
  * INTEGER - run the MD anchor whose index is INTEGER. Examples: "0", "1", etc.
  * "b_surface" - run BD simulations starting at the b-surface.

* When an option must be modified, it is usually necessary to make the change
in the **input file** and not the **model file**, and then re-run 
prepare.py to make a new model file. In some exceptions, the model file may be
directly (and carefully) modified.



## Authors and Contributors

The following people have contributed directly to the coding and validation
efforts of SEEKR2 (listed an alphabetical order of last name). 
Thanks also to everyone who has helped or will help improve this project by 
providing feedback, bug reports, or other comments.

* Rommie Amaro (principal investigator)
* Ilker Deveci (developer)
* Hilliary Frank (contributor)
* Ben Jagger (developer)
* Anand Ojha (developer)
* Andy Stokely (developer)
* Natalie Timms (contributor)
* Lane Votapka (lead developer)
* Jeff Wagner (contributor)

### Citing SEEKR2

If you use SEEKR2, please cite the following paper:

* Votapka, L. W.; Stokely, A. M.; Ojha, A. A.; Amaro, R. E. SEEKR2: Versatile Multiscale Milestoning Utilizing the OpenMM Molecular Dynamics Engine. J. Chem. Inf. Mod. In Review. https://doi.org/10.33774/chemrxiv-2021-pplfs

You may also optionally cite the following papers:

* Votapka, L. W.; Jagger, B. R.; Heyneman, A. L.; Amaro, R. E. SEEKR: Simulation Enabled Estimation of Kinetic Rates, A Computational Tool to Estimate Molecular Kinetics and Its Application to Trypsin–Benzamidine Binding. J. Phys. Chem. B 2017, 121 (15), 3597–3606. https://doi.org/10.1021/acs.jpcb.6b09388. 

* Jagger, B. R.; Ojha, A. A.; Amaro, R. E. Predicting Ligand Binding Kinetics Using a Markovian Milestoning with Voronoi Tessellations Multiscale Approach. J. Chem. Theory Comput. 2020. https://doi.org/10.1021/acs.jctc.0c00495. 

* Jagger, B. R.; Lee, C. T.; Amaro, R. E. Quantitative Ranking of Ligand Binding Kinetics with a Multiscale Milestoning Simulation Approach. J. Phys. Chem. Lett. 2018, 9 (17), 4941–4948. https://doi.org/10.1021/acs.jpclett.8b02047. 

* Votapka LW, Amaro RE (2015) Multiscale Estimation of Binding Kinetics Using Brownian Dynamics, Molecular Dynamics and Milestoning. PLOS Computational Biology 11(10): e1004381. https://doi.org/10.1371/journal.pcbi.1004381

### Copyright

Copyright (c) 2021, Lane Votapka


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.5.
