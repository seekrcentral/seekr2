Preparing Input File
====================

OpenMMVT calculations require the preparation of an input XML file. Samples
of such files can be found in the openmmvt/openmmvt/data folder named 
sample_input.xml and sample2_input.xml.

It is recommended that you copy one of these to another location and adapt
them for your own purposes.

``cp sample_input.xml my_input.xml``

Opening the file "my_input.xml" in a text editor, one can see a number of
adjustible parameters. 

Warning: it is usually safe to adjust the *text contents* of an XML tag, but
one should not change any of the *tags* or *tag attributes*. For instance, 
the attribute labeled as ``class='Model_input'`` should not be changed, but 
the 277.8 in ``<temperature>277.8</temperature>`` can be safely changed.

Adjustible parameters
---------------------
Note: Most parameters (excluding file names or paths) are case-insensitive.

* **temperature** - The temperature (in Kelvin) of all starting velocity 
  distributions and thermostats.
* **pressure** - The pressure (in bar) of the barostat if the "npt" option is 
  chosen for **ensemble**, ignored otherwise.
* **ensemble** - Options include "nvt" for constant volume/temperature 
  simulations, and "npt" for constant pressure/temperature simulations.
* **root_directory** - A path to directory which will be created, if it does
  not yet exist, where all calculation files will be written.
* **md_program** - Which program to use for MD. Options include "openmm" or 
  "namd".
* **md_output_frequency** - The number of MD timesteps between outputs to MD-
  stage DCD trajectory writes, restart backups, and state data reports.
* **md_steps_per_anchor** - The default number of MD timesteps that should be
  run per Voronoi tesselation (VT) cell (anchor). This value can be adjusted
  later on a per-anchor basis if needed.
* **constraints** - Impose bond or angle constraints. For OpenMM,
  the following options are available: "None", "HBonds", "AllBonds", or 
  "AllAngles". For NAMD, the following options are available: "None", "Water",
  or "HBonds". This option will allow one to increase simulation timestep with
  minimal loss of simulation accuracy (typically).
* **rigidWater** - For OpenMM, this argument will be passed when creating a 
  System, and if "True", will constrain the water bond lengths and angles. When
  using NAMD, this setting has no effect.
* **hydrogenMass** - In OpenMM, a number (in AMU) may be placed inside this 
  tag to allow hydrogen mass repartitioning (HMR). If a different mass for 
  hydrogens are chosen, then OpenMM will automatically repartition mass to the 
  hydrogen nuclei from the bonded heavy atom. If left empty, this tag will use
  the default hydrogen mass. If using NAMD, this setting currently has no 
  effect.
* **timestep** - The time step size (in picoseconds) for the MD simulation.
* **cv_input** - Settings for collective variables (CV) are used to construct
  the Voronoi tesselation and the resulting milestoning model.
  
  * **group1** - The indices of the atoms in the system whose center of mass
    will be one end of the distance collective variable. In a receptor-ligand 
    binding calculation, this would probably be a list of atom indices whose 
    center of mass would define the receptor's binding site.
  * **group2** - The indices of the atoms in the system whose center of mass
    will be the *other* end of the distance CV. In a receptor-ligand binding
    calculation, this would probably be a list of atom indices whose center of
    mass would define the ligand's center of mass - and probably would include
    all atoms within the ligand.
    
Preparing anchors
-----------------
All anchors must be defined within the "<input_anchors>" section. You will need
as many <input_anchor> sections as there are anchors - one <input_anchor> tag
per anchor.

* **radius** - The radius (in nanometers) of the *center* of the Voronoi cell 
  (that center is also officially known as the Anchor).
  This is an important distinction - these will not be the locations of the 
  milestones (edges of the Voronoi cells). The milestones will be placed,
  automatically, halfway between each pair of anchor radii. Example, if your
  lowest anchor has a radius of 0.1, and the next one has a radius of 0.2, then
  the milestone will be placed at a radius of 0.15 - halfway in between.
* **starting_amber_params** - Parameters for AMBER input files.

  * **prmtop_filename** - A path to a .prmtop or .parm7 file for this anchor.
  * **inpcrd_filename** - A path to a .inpcrd or .rst7 file for this anchor.
  * **box_vectors** - Box vectors may be entered manually here. If left empty,
    then the box vectors in the inpcrd file will be used.
  * **pdb_coordinates_filename** - A path to a .pdb file which can provide the
    starting positions for each atom. If left empty, then the coordinates of
    the inpcrd file will be used.
    
* **bound_state** - If set to "True", then OpenMMVT will assume that if the
  system exists within the VT cell defined by this anchor, then it is in the
  bound state. At least one anchor should have this quantity set to "True".
* **bulk_anchor** - If set to "True", then OpenMMVT will assume that if the 
  system exists within the VT cell defined by this anchor, then it is in the
  "bulk" state, or it has escaped from the binding site. Exactly one anchor
  should have this quantity set to "True" - probably the anchor with the 
  largest radius.

Browndye settings
-----------------
If a k-on calculation is desired, then you must supply setting for BrownDye. If
the <browndye_settings> section is left empty, then no BrownDye can be run, and
k-on cannot be obtained.

* **binary_directory** - A path to the location of all the BrownDye executables
* **receptor_pqr_filename** - A path to the PQR file that represents the 
  receptor
* **ligand_pqr_filename** - A path to the PQR file that represents the ligand.
* **apbs_grid_spacing** - The resolution of the APBS grid.
* **receptor_indices** - The atom indices whose center of mass will be the
  center of the binding site. This is likely to be an identical list to the
  <group1> list in the <cv_input> section mentioned above.
* **ligand_indices** - The atom indices whose center of mass will be the
  center of the ligand molecule. This is likely to be an identical list to the
  <group2> list in the <cv_input> section mentioned above.
* **ions** - Define which ions will exist in the APBS calculation.

  * **ion** - Create an ion block for each type of ion present.
  
    * **radius** - The radius of the ion (in Angstroms).
    * **charge** - The charge of the ion (in proton charge).
    * **conc** - The concentration of the ion in the system (in moles per 
      liter).
    
* **num_bd_milestone_steps** - The number of simulations to run per point on
  the first hitting point distribution.
* **num_b_surface_steps** - The number of simulations to start from the 
  b-surface.
* **n_threads** - The number of threads (cores) to use in the BrownDye 
  calculation.
