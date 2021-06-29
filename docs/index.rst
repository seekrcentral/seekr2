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

Getting Involved
================

Please report **bugs** or **enhancement requests** through the `Issue 
Tracker`_. 

.. _Issue Tracker: https://github.com/seekrcentral/seekr2/issues

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   running_calculations
   tutorial
   api
   



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`