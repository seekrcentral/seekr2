Analyze Results
===============

Once all MD and BD simulations are complete, then one must run the "analysis"
stage, which uses the statistics from the simulations, along with milestoning
theory, to obtain kinetic and thermodynamic quantities.

Analyze
-------

One may run a generic analysis with the following command:

``python analyze.py /path/to/directory/model.xml``

If insufficient statistics exist within any anchors, the program should throw
an error. If this happens, then more simulations must be run for those anchors
before any analysis may be performed. (Note: sufficient statistics to perform
an analysis does not necessarily indicate convergence. Refer to "Converge"
section below to converge your simulations and maximize accuracy).

Under some circumstances, one may wish to adjust the number of Monte Carlo
matrices sampled, which are used to estimate error intervals:

``python analyze.py /path/to/directory/model.xml -n 2000``

Or, if one wishes to leave out error calculations entirely:

``python analyze.py /path/to/directory/model.xml -n 0``

In some unusual circumstances, the default MMVT theory will be unable to
accurately solve an ill-conditioned matrix, which might occur with very long-
timescale processes. In this case, the pre-equilibrium approximation can be
used in place of MMVT to approximate kinetics:

``python analyze.py /path/to/directory/model.xml -p``

Converge
--------

In many cases, it can be very helpful to analyze the convergence of computed
kinetics quantities and/or anchor statistics. A generic convergence analysis 
can be computed using the following command:

``python converge.py /path/to/directory/model.xml``

All convergence plots will be saved to the "images_and_plots/" directory.

One may also want to see which, if any, anchors have converged so that no more
simulations need to be run. To do this, one may do an RMSD convergence analysis:

``python converge.py /path/to/directory/model.xml -r``

The program will print which anchors have been converged or not.

Many other options exist for this program, which can be found using the "-h"
argument.