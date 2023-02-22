SEEKR2: Multiscale Milestoning Calculations of Kinetics and Thermodynamics
==========================================================================

:Release: |release|
:Date: |today|

**SEEKR2** Simulation Enabled Estimation of Kinetic Rates version 2
is a package to perform simulations at multiple scales for the calculation of
kinetics and thermodynamics quantities. 

Most typically, one would use this package to obtain the kinetics and 
thermodynamics of binding and unbinding between a small molecule and 
biomolecular receptor. Also typically, this package performs molecular 
dynamics (MD) and Brownian dynamics (BD) simulations of the binding process 
within multiple regions of phase space, utilizing expensive, but highly 
accurate MD for regions where close molecular interaction would require high 
accuracy, but utilizing faster, but more approximate BD for regions of greater 
intermolecular distance. This way the advantages of each simulation regime is 
augmented and their disadvantages are reduced. The regimes are then combined 
using the theory of milestoning such that the kinetics and thermodynamics of the 
process under interest is obtained.

Two versions of milestoning theory are available to use in SEEKR2:

1. **Markovian milestoning with Voronoi Tesselations (MMVT)**:  Vanden-Eijnden, E.; Venturoli, M. Markovian Milestoning with Voronoi Tessellations. J. Chem. Phys. 2009, 130 (19), 194101. https://doi.org/10.1063/1.3129843

2. **The original formulation of milestoning theory (Elber milestoning)**: Faradjian, A. K.; Elber, R. Computing Time Scales from Reaction Coordinates by Milestoning. J. Chem. Phys. 2004, 120 (23), 10880–10889. https://doi.org/10.1063/1.1738640

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation
   running_calculations
   tutorial
   model_input_files
   program_options
   api
   faq

Cite SEEKR2
===========

If you wish to cite SEEKR2, please cite the following paper:

* Votapka, L. W.; Stokely, A. M.; Ojha, A. A.; Amaro, R. E. SEEKR2: Versatile Multiscale Milestoning Utilizing the OpenMM Molecular Dynamics Engine. J. Chem. Inf. Mod. 2022 62 (13), 3253-3262. DOI: 10.1021/acs.jcim.2c00501

You may also optionally cite one or more of the following papers:

* Votapka, L. W.; Jagger, B. R.; Heyneman, A. L.; Amaro, R. E. SEEKR: Simulation Enabled Estimation of Kinetic Rates, A Computational Tool to Estimate Molecular Kinetics and Its Application to Trypsin–Benzamidine Binding. J. Phys. Chem. B 2017, 121 (15), 3597–3606. https://doi.org/10.1021/acs.jpcb.6b09388. 

* Jagger, B. R.; Ojha, A. A.; Amaro, R. E. Predicting Ligand Binding Kinetics Using a Markovian Milestoning with Voronoi Tessellations Multiscale Approach. J. Chem. Theory Comput. 2020. https://doi.org/10.1021/acs.jctc.0c00495. 

* Jagger, B. R.; Lee, C. T.; Amaro, R. E. Quantitative Ranking of Ligand Binding Kinetics with a Multiscale Milestoning Simulation Approach. J. Phys. Chem. Lett. 2018, 9 (17), 4941–4948. https://doi.org/10.1021/acs.jpclett.8b02047. 

* Votapka LW, Amaro RE (2015) Multiscale Estimation of Binding Kinetics Using Brownian Dynamics, Molecular Dynamics and Milestoning. PLOS Computational Biology 11(10): e1004381. https://doi.org/10.1371/journal.pcbi.1004381


Getting Involved
================

Please report **bugs** or **enhancement requests** through the `Issue 
Tracker`_. 

.. _Issue Tracker: https://github.com/seekrcentral/seekr2/issues

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`