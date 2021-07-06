Model Input Files
=================

Typical SEEKR2 calculations must begin with a model input XML file. (Although
this is not necessary if the :doc:`API<api>` is used.) This page describes the
parameters of the model input file so that users may adapt the sample inputs
to their own systems.

<calculation_type>
  The type of milestoning model to employ. At present,
  the only available options are 1) "mmvt" for Markovian milestoning with
  Voronoi Tesselations, and 2) "elber" for the original version of milestoning
  proposed by Elber et. al.
  
<calculation_settings>
  A block of settings specific to the setting
  chosen for the **<calculation_type>** tag. The **<calculation_settings>** 
  variable will be different for different settings of the 
  **<calculation_type>**.
  
  If **<calculation_type>** is set to "mmvt", these are the tags within
  **<calculation_settings>**.
  
  <md_output_interval>
    The interval (in MD timesteps) between outputs
    of simulation state information, trajectory frames, and restart checkpoints.
  
  **<md_steps_per_anchor>**
    The total number of MD timesteps that will be
    run per anchor (assuming optional convergence criteria are met). This 
    parameter can also be later adjusted with the "-t" argument in run.py for
    individual anchors, if a different number of steps is desired than what is
    entered here in the model input file.
    
  If **<calculation_type>** is set to "elber", these are the tags within
  **<calculation_settings>**.
  
  **<temperature_equil_progression>**
    If a temperature equilibration step is
    desired, one may enter a list of floats representing temperatures (in
    K) and the MD simulation engine will run this series of temperatures
    for **<num_temperature_equil_steps>** number of timesteps per temperature
    step.
    
  **<num_temperature_equil_steps>**
    If a temperature equilibration step is
    desired, one must enter a number of timesteps here, along with a 
    progression of temperatures in **<temperature_equil_progression>**
    representing the number of steps to run at each temperature.
    
  **<umbrella_force_constant>**
    Since Elber milestoning requires a harmonic
    force to restrain the system close to a given milestone, this parameter
    represents the force constant (in units of kJ per mole per nm^2).
  
  **<fwd_rev_interval>**
    In Elber milestoning, the umbrella stage generates
    an equilibrium distribution of conformations, which are used to spawn
    reversal-stage, and eventually forward-stage, simulations. This parameter
    represents the number of umbrella stage timesteps between spawning reversal-
    stage simulations.
    
  **<rev_output_interval>**
    The interval (in MD timesteps) between outputs
    of simulation state information, trajectory frames, and restart checkpoints
    in the reversal stage.
    
  **<fwd_output_interval>**
    The interval (in MD timesteps) between outputs
    of simulation state information, trajectory frames, and restart checkpoints
    in the forward stage.
    
**<temperature>**
  The temperature (in K) to use for *all* stages of the 
  SEEKR2 calculation - including starting temperatures, thermostated 
  temperatures, analysis, etc.
  
**<pressure>**
  If the **<ensemble>** is set to "npt" (constant pressure and
  temperature), this parameter defines the target pressure (in bar) for the 
  simulation barostat. If **<ensemble>** is set to "nvt" (constant volume and
  constant temperature), this parameter is ignored.
  
**<ensemble>**
  Defines the ensemble of the MD simulations. Options at this 
  time include "nvt" (constant volume and temperature), "npt" (constant pressure
  and temperature).
  
**<root_directory>**
  A filesystem path to a directory name where the SEEKR2
  calculation files will be written. The model.xml file will be written just
  inside this directory.
  
**<md_program>**
  Options for this parameter include "openmm" to use the 
  OpenMM engine to run the MD simulations, and "namd" to run the MD simulations
  using the NAMD2 engine.
  
**<constraints>**
  Whether any bonds or angles should be constrained in the MD
  simulation. Options include:
  
  #. "none" - No bonds or angles are constrained
  #. "hbonds" - Only bonds involving hydrogen atoms are constrained.
  #. "allbonds" - All bonds are constrained and rigid.
  #. "hangles" - All bonds are constrained, and so are angles involving
     hydrogen atoms.
     
**<rigid_water>** 
  Whether water molecules will have a rigid angle. Options
  include "True" or "False".
  
**<hydrogen_mass>**
  In order to use hydrogen mass repartitioning (HMR), set
  this parameter to a new mass (in AMU) to move from a bonded heavy atom to the
  hydrogen atom.
  
**<timestep>**
  The MD timestep (in ps). (Will be automatically converted to fs for NAMD 
  runs).

**<nonbonded_cutoff>**
  The nonbonded cutoff (in nm) to use in MD.

**<cv_inputs>**
  Settings defining the collective variables (CV) and their
  associated anchors and milestones will be defined within this block.
  
  **<cv_input>**
    One of the input structures for a CV. The type of CV is 
    defined by the "class" attribute in the XML tag.
    
    **<group#>**
      CVs are functions of atomic positions, and in general,
      involve functions of centers of masses of groups of molecules. Depending
      on how many inputs there are to the CV function, a given CV will have
      a certain number of atomic groups. These parameters are lists of atom
      indices in the input structures (starting from zero).
      
    **<input_anchors>**
      The "input anchors" are not actual anchor objects, 
      but are used as inputs to construct the anchors in the model. This is a
      block of such objects
      
      **<input_anchor>**
        One of the anchor inputs used to create model 
        anchors. Different CVs can have different anchor input types that 
        contain different attributes. The "class" attribute of the input_anchor
        matches to a particular CV.
        
Spherical Collective Variable
-----------------------------

A "spherical" collective variable is defined as, simply, the distance between
two center of masses of groups of atoms. This CV is called "spherical" because
the surfaces located at constant values of the CV create spherical milestones.

The following are attributes of spherical anchor inputs:

**<radius>**
  The radius (in nm) of the *anchor* (between the milestones). 
  The significance of this variables comes from Voronoi tesselation (VT) 
  definitions where the milestones are drawn exactly halfway between anchor 
  points. The strict definition of surfaces in Voronoi tesselations don't need
  to be enforced. In the absence of *<lower_milestone_radius>* and 
  *<upper_milestone_radius>*, the *<radius>* parameter will be used to define 
  the milestones - the milestones will lie exactly halfway between adjacent 
  anchor radii. Even if *<lower_milestone_radius>* and 
  *<upper_milestone_radius>* are defined, this parameter is required and 
  should be entered since SEEKR2 uses it to keep track of changes to the model.
  
**<lower_milestone_radius>**
  If the lower milestone of this anchor should
  not lie exactly halfway between this anchor's radius and the lower anchor's
  radius, then enter a radius (in nm) where the lower milestone should be.
  
**<upper_milestone_radius>**
  If the upper milestone of this anchor should
  not lie exactly halfway between this anchor's radius and the upper anchor's
  radius, then enter a radius (in nm) where the upper milestone should be.
  
**<starting_amber_params>**
  This block contains parameters and starting
  atomic positions and box vectors for a system that uses the Amber forcefield.
  
  **<prmtop_filename>**
    Enter a file path for the Amber parameter/topology
    file (format: .prmtop or .parm7)

  **<pdb_coordinates_filename>**
    Enter a path to a PDB which contains the
    starting atomic positions for this anchor.
  
    The PDB file may also contain a CRYST line that defines box vectors
    
    ``CRYST1   40.142   40.329   32.472  90.00  90.00  90.00 P 1``        

  **<box_vectors>**
    Optionally enter a box vector object for this anchor to
    define the simulation box vectors. If left empty, the box vectors will be
    taken from the "CRYST" line in the PDB file (if it exists).
    
    XML for a box vector object looks like (units are in nm)

.. code-block::

  <box_vectors class"Box_vectors">
    <ax type="float">4.0251</ax>
    <ay type="float">0.0</ay>
    <az type="float">0.0</az>
    <bx type="float">0.0</bx>
    <by type="float">4.044</by>
    <bz type="float">0.0</bz>
    <cx type="float">0.0</cx>
    <cy type="float">0.0</cy>
    <cz type="float">3.2561</cz>
  </box_vectors>
  

**<bound_state>**
  Mark this anchor as a "bound state" This affects k-on
  calculations primarily. Options include "True" or "False".
  
**<bulk_anchor>**
  Mark this anchor as a "bulk anchor" which represents, in
  effect, the completely dissociated state. This state should have no 
  structure or parameter information assigned, and is primarily used to
  calculate the k-off. The bulk anchor may coincide with the locations of
  BD milestones.
  
Browndye Settings
-----------------

The **<browndye_settings_input>** section of the model input file may be 
optionally left blank, but it must be filled out if BD simulations and k-on 
calculations are desired.

**<binary_directory>**
  The directory containing Browndye2 programs (binaries). This
  parameter may be an empty string if the Browndye2 programs have been placed
  in the PATH environmental variable.
  
**<receptor_pqr_filename>**
  A file path to a PQR file representing the 
  receptor.
  
**<ligand_pqr_filename>**
  A file path to a PQR file representing the ligand.

**<apbs_grid_spacing>**
  The space (in Angstroms) between grid points in an 
  APBS calculation.
  
**<receptor_indices>**
  The atom indices (numbering starting from zero) 
  defining the binding site in the PQR file specified by the
  **<receptor_pqr_filename>**. The center of mass of these atoms will be taken
  to define the center of the binding site.
  
**<ligand_indices>**
  The atom indices (numbering starting from zero) 
  defining the center of the molecule in the PQR file specified by the
  **<ligand_pqr_filename>**. The center of mass of these atoms will be taken
  to define the center of the ligand molecule.
  
**<ions>**
  A block of Ion objects that will be used in the APBS calculations
  and, by extension, the BD calculations.
  
  **<ion>**
    An object representing an ion. You should format these inputs
    the same as would be input into APBS.
  
    **<radius>**
      The radius of the ion (in Angstroms).
    
    **<charge>**
      The charge of the ion (in proton charge "e")
    
    **<conc>**
      The concentration of the electrolyte (in moles per liter)
    
**<num_bd_milestone_trajectories>**
  The number of BD trajectories to run
  per encounter complex extracted from the b-surface simulations
  
**<num_b_surface_trajectories>**
  The total number of trajectories to run for
  the b-surface simulations.
  
**<max_b_surface_trajs_to_extract>**
  The maximum number of encounter 
  complexes extracted from the b-surface simulations to run the BD milestones
  simulations.
  
**<n_threads>**
  The number of CPUs to use in the Browndye2 calculations.