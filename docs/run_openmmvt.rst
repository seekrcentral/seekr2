Running OpenMMVT
================

Once you have successfully created an input XML file, you can prepare and run
the OpenMMVT calculations.

All of these commands are run with the assumption that the user is in the
openmmvt/openmmvt/ directory.

Also, most openmmvt programs have been equipped with a '-h' argument for
easy reference to usage.

Prepare
-------

Let us assume that your input XML file is named "my_input.xml"

``python prepare_1d_spherical.py /path/to/my_input.xml``

This will create a file tree in the location defined by the <root_directory>
tag in the input XML. Looking inside this directory, you should see the 
anchor directories, and if you provided BrownDye inputs, you should also see
at least one bd_milestone directory and a b_surface directory. Most 
importantly, notice the model.xml file, which contains all information that 
future stages of the calculation will use.

In a pinch, this model.xml file can be modified manually, although if possible,
it is recommended to rerun prepare_1d_spherical.py.

This entire directory can be copied to different file systems as-is, and if 
that file system has OpenMMVT installed, then you should be able to just run 
any subsequent calculations by providing the path to model.xml when required.

Once the directory containing model.xml and the anchor directories has been
created or transferred to the computer where calculations will be performed,
the calculation is ready to run.

Run With OpenMM
---------------

Assuming that OpenMM and the OpenMMVT plugin have been installed, one may run
calculations for the anchors.

For an anchor where all parameter, topology, coordinates, and box_vectors are
defined, one may run the following command.

``python runner_openmm.py 1 /path/to/directory/model.xml``

Replace "1" with the index of the anchor you wish to run. This is the easiest
way to run an anchor.

However, if parameter, topology, coordinates, or box_vectors are missing for 
a given anchor, then it is possible to save a state that has crossed into an
adjacent anchor.

``python runner_openmm.py 1 /path/to/directory/model.xml -s``

This will save OpenMM state files in anchor_1/prod/states. States are named
as "openmmvt_BOUNDARY_COUNT" where BOUNDARY is an internal index of a Voronoi
cell's boundary (for spherical VT's, '1' is generally the downward
boundary, while '2' is the upward boundary), and COUNT is incremented with 
every bounce, therefore larger COUNT numbers are occurring later in the 
simulation.

Then load the saved state into the next anchor over.

``python runner_openmm.py 0 /path/to/directory/model.xml -l /path/to/directory/
anchor_1/prod/states/openmmvt_1_1``

All anchors must be run before proceeding to the analysis stage. To do this
one could create a loop to run all anchors (assuming that all parameters, 
topology, coordinates, and box_vectors have been defined for all anchors).

``for i in 0 1 2 3 4 5 6 7 8; do python runner_openmm.py $i 
/path/to/directory/model.xml; done``

This, of course, assumed there were 9 anchors (0 thru 8).

Again, use the "-h" argument to see additional options or arguments.

Run With NAMD
-------------

If NAMD has been installed, then the MD stages can be run with NAMD if the
"<md_program>" option has been set to "namd".

Simulations may be run for anchor 1 using a command like the following:

``python runner_namd.py 1 /path/to/directory/model.xml``

The same state save/load capabilities mentioned in the "Run With OpenMM" 
section above can be used in NAMD.

If a specific version or build of NAMD is desired besides the generic "namd2"
command, one may use the "-n" option. For instance:

``python runner_namd.py 1 /path/to/directory/model.xml -n "charmrun namd2"``

One may also add additional NAMD arguments using the (-a) argument:

``python runner_namd.py 1 /path/to/directory/model.xml -a "++n 4 ++ppn 3"``

Consult the `NAMD manual <http://www.ks.uiuc.edu/Research/namd/2.14/ug/>`_ for
detailed information about NAMD commands and arguments.

All anchors must be run before proceeding to the analysis stage. To do this
one could create a loop to run all anchors (assuming that all parameters, 
topology, coordinates, and box_vectors have been defined for all anchors).

``for i in 0 1 2 3 4 5 6 7 8; do python runner_namd.py $i 
/path/to/directory/model.xml; done``

This, of course, assumed there were 9 anchors (0 thru 8).

Again, use the "-h" argument to see additional options or arguments.

Run With BrownDye
-----------------

In order to compute k-ons, BrownDye2 must be installed, and the input XML file
must be configured to enable BrownDye runs.

First, the b_surface simulations must be run:

``python runner_browndye2 b_surface /path/to/directory/model.xml``

Once this is complete, one may run the BD Milestones:

``python runner_browndye2 0 /path/to/directory/model.xml``

Where 0 is the index of the BD Milestone.

Again, use the "-h" argument to see additional options or arguments.