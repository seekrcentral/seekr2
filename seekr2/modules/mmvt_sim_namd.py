"""
mmvt_sim_namd.py

Base objects and routines for preparing NAMD simulations
to run for MMVT calculations.
"""

import os
import datetime

import parmed

import seekr2.modules.common_sim_namd as common_sim_namd
from seekr2.modules.common_sim_namd import add_string_buffer

dt_date = datetime.datetime.now()
SEEKR_MMVT_PROD_HEADER = "#" + "*"*64 + "\n" \
    + "#  Name SEEKR MMVT - Production Stage\n" \
    + "#  Date " + dt_date.strftime("%Y-%b-%d") + "\n" \
    + "#  Generated with the SEEKR2 automatic NAMD input file generator\n" \
    + "#" + "*"*64 + "\n"

class Colvars_config():
    """
    The settings used to fill out the NAMD colvars script.
    
    Attributes:
    -----------
    colvarstrajfrequency : int, Default 0
        The frequency to write the values of each colvar. The default
        value of 0 will not write the trajectory.
    colvarsrestartfrequency : int, Default 0
        The frequency to write a colvars restart file. The default
        value of 0 will not write the restarts.
    colvars_string_list : list
        A list of strings which are multiline colvar commands which
        define the collective variables to monitor. These strings are
        written to the colvars script.
    """
    
    def __init__(self):
        self.colvarstrajfrequency = 0
        self.colvarsrestartfrequency = 0
        self.colvars_string_list = []
        return
    
    def fill_out_from_model(self, model):
        """
        Use the SEEKR Model() object settings to fill out this object
        with relevant quantities.
        """
        
        self.colvars_string_list = []
        for cv in model.collective_variables:
            cv_string = cv.make_namd_colvar_string()
            self.colvars_string_list.append(cv_string)
        return
    
    def to_string(self):
        """
        Convert all settings in this object to a format which can be
        provided to NAMD colvar input file.
        """
        
        my_string = ""
        assert self.colvarstrajfrequency >= 0
        my_string += add_string_buffer("colvarstrajfrequency", 
                                       str(self.colvarstrajfrequency))
        assert self.colvarsrestartfrequency >= 0
        my_string += add_string_buffer("colvarsrestartfrequency", 
                                       str(self.colvarsrestartfrequency))
        
        for cv_string in self.colvars_string_list:
            my_string += cv_string
        
        return my_string
    
    def write_file(self, filename):
        """
        Convert all settings in this object to a string and then write
        a NAMD colvars script file.
        
        Parameters:
        -----------
        filename : str
            The name of the colvars script file to write to.
        """
        
        my_string = self.to_string()
        with open(filename, "w") as f:
            f.write(my_string)
        return
    
class Seekr_namd_settings():
    """
    The settings within the NAMD script that define variables used
    by SEEKR.
    
    Attributes:
    -----------
    max_steps : int
        The total number of steps to run within the NAMD simulation.
    eval_stride : int
        How many steps to take between evaluations of the colvars and
        possible bounce events.
    save_state : bool
        Whether to save the system state upon collision with a boundary.
    save_one_state_for_all_boundaries : bool
        Whether to save at least one state for each of the boundaries
        in the simulation.
    check_state_interval : float
        If one state is being saved for all boundaries, define an 
        interval at which to check each of the boundaries to see
        if a state is saved for it. If so, then stop saving states.
    """
    
    def __init__(self):
        self.max_steps = -1
        self.eval_stride = -1
        self.save_state = False
        self.save_one_state_for_all_boundaries = False
        self.check_state_interval = 1000
        return
        
    def fill_out_from_model(self, model):
        """
        Use the SEEKR Model() object settings to fill out this object
        with relevant quantities.
        """
        
        self.max_steps = model.calculation_settings.num_production_steps
        self.eval_stride = model.namd_settings.eval_stride
        return
        
    def to_string(self, model, anchor, starting_incubation_steps=0):
        """
        Convert all settings in this object to a format which can be
        provided to the SEEKR portion of the NAMD script.
        """
        
        my_string = "\n# SEEKR variables\n"
        my_string += add_string_buffer("set crossing", "0")
        my_string += add_string_buffer("set whoami", "none")
        my_string += add_string_buffer("set incubation_steps", "0")
        assert self.max_steps >= 0
        my_string += add_string_buffer("set max_steps", str(self.max_steps))
        assert self.eval_stride >= 0
        my_string += add_string_buffer("set EVAL_STRIDE", str(self.eval_stride))
        my_string += add_string_buffer("set CURR_ANCHOR", str(anchor.index))
        
        for milestone in anchor.milestones:
            assert milestone.alias_index >= 0
            assert milestone.neighbor_anchor_index >= 0
            anchor_index_set = "set ANCHOR_{0}".format(milestone.alias_index)
            my_string += add_string_buffer(anchor_index_set, 
                                       str(milestone.neighbor_anchor_index))
            assert milestone.index >= 0
            milestone_index_set = "set MILESTONE_{0}".format(
                milestone.alias_index)
            my_string += add_string_buffer(milestone_index_set, 
                                       str(milestone.index))
        my_string += add_string_buffer("set rev_last", "0")
        my_string += add_string_buffer("set total_bounces", "0")
        my_string += add_string_buffer(
            "set save_one_state_for_all_boundaries", 
            self.save_one_state_for_all_boundaries)
        my_string += add_string_buffer("set check_state_interval", 
                                       self.check_state_interval)
        if self.save_one_state_for_all_boundaries:
            my_string += add_string_buffer("set save_state_prefix", 
                                           "states/namdmmvt")
        else:
            my_string += add_string_buffer("set save_state_prefix", 
                                           "\"\"")
        my_string += add_string_buffer("set save_all_states", self.save_state)
        alias_index_list = []
        for milestone in anchor.milestones:
            alias_index_list.append(str(milestone.alias_index))
        alias_index_list_str = " ".join(alias_index_list)
        routine_head_string = """\n# SEEKR routines\n
proc assert condition {
    if {![uplevel 1 expr $condition]} {
        return -code error "assertion failed: $condition"
    }
}
assert {[expr "$check_state_interval %% $EVAL_STRIDE"] == 0}
proc all_boundaries_have_state {myglob alias_indices} {
    set orig_alias_indices $alias_indices
    set all_state_files [glob -nocomplain $myglob]
    foreach state_file $all_state_files {
        set state_alias_index [lindex [split [file tail $state_file] "_"] 1]
        assert {[lsearch -exact $orig_alias_indices $state_alias_index] >= 0}
        set alias_index [lsearch -exact $alias_indices $state_alias_index]
        if {$alias_index >= 0} {
            set alias_indices [lreplace $alias_indices $alias_index $alias_index]
        }
    }
    if {[llength $alias_indices] == 0} {
        return True
    } else {
        return False
    }
}
checkpoint
for {set stepnum %d} {$stepnum < $max_steps} {incr stepnum $EVAL_STRIDE} {
    if {([expr "$stepnum %% $check_state_interval"] == 0) && """\
    """($save_one_state_for_all_boundaries == True)} {
        set myglob "${save_state_prefix}*"
        set result [all_boundaries_have_state $myglob "%s"] 
        if {$result == True} {
            set save_state_prefix ""
            set save_one_state_for_all_boundaries False
        }
    }
    run $EVAL_STRIDE""" % (starting_incubation_steps, alias_index_list_str)
        my_string += routine_head_string
        
        done_cvs = []
        
        for milestone in anchor.milestones:
            cv = model.collective_variables[milestone.cv_index]
            if cv.index not in done_cvs:
                my_string += "\n    set cv_val_{0} [cv colvar "\
                        "collective_variable_{0} value]".format(cv.index)
                done_cvs.append(cv.index)
        
        for i, milestone in enumerate(anchor.milestones):
            cv = model.collective_variables[milestone.cv_index]
            eval_str = cv.get_namd_evaluation_string(milestone)
            anchor_index = "ANCHOR_{0}".format(milestone.alias_index)
            milestone_index = "MILESTONE_{0}".format(
                milestone.alias_index)
            
            if i == 0:
                my_string += "\n    if "
            else:
                my_string += "\n    } elseif "
            
            bounce_string1 = """{ %s } {
        if {$rev_last > 100} {
            print "trajectory stuck"
            abort
        }
        if {$stepnum == 0} {
            puts "MMVT simulation bouncing on first step: the system will """\
            """be trapped behind a boundary. Check and revise MMVT """\
            """boundary definitions and/or atomic positions."
            abort
        }
        puts "SEEKR: Cell Collision: current: $CURR_ANCHOR, new: $%s, """\
        """stepnum: $stepnum"
        if {$whoami != $%s } {
            puts "SEEKR: Milestone Transition: anchor: $CURR_ANCHOR, """\
            """source: $whoami, destination: $%s, stepnum: $stepnum, """\
            """incubation steps: $incubation_steps"
            set incubation_steps 0
        }""" % (eval_str, anchor_index, milestone_index, milestone_index)
            save_state_string = """
        if {(${save_state_prefix} != "") || ($save_all_states == True)} {
            output ${save_state_prefix}_%d_${total_bounces}
        }""" \
             % milestone.alias_index
            bounce_string2 = """
        incr total_bounces
        revert
        rescalevels -1
        checkpoint
        incr rev_last
        set whoami $%s
""" % milestone_index
            my_string += bounce_string1
            my_string += save_state_string
            my_string += bounce_string2
        
        routine_tail_string = """    } else {
        checkpoint
        set rev_last 0
    }
    incr incubation_steps $EVAL_STRIDE
} """
        my_string += routine_tail_string
        return my_string
    
def create_sim_namd(model, anchor, output_filename):
    """
    Take all relevant model and anchor information and generate
    the necessary NAMD scripts to run the simulation.
    
    Parameters
    ----------
    model : Model()
        The Model which contains the anchor and all settings for
        the simulation that is about to be run.
        
    anchor : Anchor()
        The anchor object that this NAMD simulation applies to.
        
    output_filename : str
        The base name of various output files produced by NAMD.
        
    Returns
    -------
    sim_namd : Sim_namd()
        An object to store all the information needed to run a NAMD
        SEEKR calculation.
    """
    
    sim_namd = common_sim_namd.Sim_namd()
    box_vectors = None
    if anchor.amber_params is not None:
        box_vectors = anchor.amber_params.box_vectors
        
    elif anchor.charmm_params is not None:
        raise Exception("Charmm systems not yet implemented")
    
    else:
        print("anchor.index", anchor.index)
        raise Exception("No Amber or Charmm input settings detected.")
    
    sim_namd.namd_root = common_sim_namd.Namd_root()
    sim_namd.namd_root.fill_out_from_model(model, anchor, output_filename)
    if box_vectors is not None:
        sim_namd.namd_root.periodic_boundary_conditions\
            .assign_cell_basis_vectors(box_vectors)
    sim_namd.colvars_config = Colvars_config()
    sim_namd.colvars_config.fill_out_from_model(model)
    sim_namd.seekr_namd_settings = Seekr_namd_settings()
    sim_namd.seekr_namd_settings.fill_out_from_model(model)
    return sim_namd
