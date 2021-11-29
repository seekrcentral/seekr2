Frequently Asked Questions
==========================

.. contents:: Contents
   :depth: 3
  
What is SEEKR2?
---------------
Simulation Enabled Estimation of Kinetic Rates is, simultaneously, a collection
of scientific software programs, and also a technical approach to estimate the
kinetics and thermodynamics of certain molecular processes, especially ligand-
receptor binding, using Milestoning theory. SEEKR2 is the second iteration of
the software, with better speed, capabilities, and usability than the earlier
version.

What can SEEKR2 do?
-------------------
At present, the main utility of SEEKR2 is to estimate the kinetics of binding
and unbinding (represented by the quantities k-on and k-off), in addition to
the thermodynamics of binding (represented by the Gibbs free energy of 
binding). 

Testing and preliminary calculations indicate that SEEKR2 can 
consistently obtain binding/unbinding rate constants within an order of 
magnitude of the experimental quantity, and sometimes much better, even for 
systems with very long residence times. SEEKR2 can also usually obtain Gibbs
free energies of binding within ~1 kcal/mol of experiment. In addition, the
SEEKR method has been able to rank compounds by their residence time and
affinity to a receptor with good accuracy.

More details about SEEKR2's accuracy and performance can be obtained from the
publications on the :doc:`Index page<index>`.

Development is ongoing, and more capabilities for SEEKR2 are being actively
pursued. The kinetics and thermodynamics of intramolecular motion (such as 
hinge or pocket opening/closing), intersite transfer, peptide folding, and
processes can be estimated using SEEKR2. Publications of these applications
are forthcoming.

What if I doubt SEEKR2's usefulness/validity?
---------------------------------------------



How can I use SEEKR2 for my research?
-------------------------------------


How does one get a benchmark of a SEEKR2 calculation?
-----------------------------------------------------
If the SEEKR2 run.py program terminates successfully, the last thing it will
print is a benchmark (in ns/day) for an MD anchor. Therefore, if you wish to
obtain a benchmark, consider running a short SEEKR2 calculation by passing a
relatively small number (like 10000) to the "-t" argument of run.py.


How many CPUs are optimal for a SEEKR2 calculation?
---------------------------------------------------
In NAMD mode, using multiple CPUs is likely to make the calculation go faster.
However, in OpenMM mode, the CUDA platform is designed to use only one CPU,
so using multiple CPUs in OpenMM mode is not likely to make the calculation go
any faster.


How are the convergence plots calculated in SEEKR2?
---------------------------------------------------



How does SEEKR2 compute the convergence value for the individual anchors?
-------------------------------------------------------------------------
