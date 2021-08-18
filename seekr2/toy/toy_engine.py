"""
toy_engine.py


"""

bd_results_string = \
"""<rates>
  <solvent>
    <kT> 1 </kT>
    <debye_length> 1.0 </debye_length>
    <dielectric> 78 </dielectric>
    <vacuum_permittivity> 0.000142 </vacuum_permittivity>
    <water_viscosity> 0.243 </water_viscosity>
    <relative_viscosity> 1 </relative_viscosity>
  </solvent>
    <hydrodynamic_interactions> 0 </hydrodynamic_interactions>
  <time_step_tolerances>
    <minimum_core_dt> 0.2 </minimum_core_dt>
    <minimum_core_rxn_dt> 0.05 </minimum_core_rxn_dt>
    <minimum_chain_dt> 0 </minimum_chain_dt>
    <minimum_chain_rxn_dt> 0 </minimum_chain_rxn_dt>
  </time_step_tolerances>
  <molecule_info>
    <b_radius> {b_radius} </b_radius>
    <b_reaction_rate> {b_reaction_rate} </b_reaction_rate>
    <n_reactions> 1 </n_reactions>
  </molecule_info>
  <reactions>
    <n_trajectories> {n_trajectories} </n_trajectories>
    <stuck> 0 </stuck>
    <escaped> {escaped} </escaped>
    <completed>
      <name> {name1} </name>
      <n> {n1} </n>
    </completed>
    <completed>
      <name> {name2} </name>
      <n> {n2} </n>
    </completed>
  </reactions>
  <n_bd_steps> 0 </n_bd_steps>
</rates>
"""

import numpy as np

def write_b_surface_output_file(output_file_name, k_on_src, b_transition_probs,
                                milestone_transition_probs):
    CONV_FACTOR = 602000000.0
    b_reaction_rate = k_on_src / CONV_FACTOR
    N_TRAJS = 1000000000
    for key in b_transition_probs:
        if key == "escaped":
            escaped = int(b_transition_probs["escaped"] * N_TRAJS)
        else:
            key1 = key
            name1 = "b_%s" % key
            n1 = int(b_transition_probs[key] * N_TRAJS)
            
    for key in milestone_transition_probs:
        if key == "escaped":
            escaped2 = int(milestone_transition_probs["escaped"] * n1)
            pass
        else:
            name2 = "%s_%s" % (key1, key)
            n2 = int(milestone_transition_probs[key] * n1)
    
    new_result = bd_results_string.format(
        b_radius=0.0, b_reaction_rate=b_reaction_rate, n_trajectories=N_TRAJS,
        escaped=escaped+escaped2, name1=name1, n1=n1, name2=name2, n2=n2)
    #print("new_result:", new_result)
    with open(output_file_name, "w") as f:
        f.write(new_result)
    
    return

class BrownianParticle():
    def __init__(self, mass, dimensions, position, velocity, diffusion,
                 potential_energy_function):
        self.mass = mass
        self.dimensions = dimensions
        assert len(position.shape) == dimensions
        self.position = position
        self.velocity = velocity
        self.diffusion = diffusion
        self.potential_energy_function = potential_energy_function
        return



class ToyIntegrator():
    """
    Base class for toy integrators.
    """
    def __init__(self):
        pass

class SmoluchowskiSphericalMMVTIntegrator(ToyIntegrator):
    """
    An integrator that generates MMVT output files based on a Smoluchowski
    region.
    """
    def __init__(self, calc, index, output_file_name, style="openmm"):
        self.calc = calc
        self.index = index
        self.style = style.lower()
        self.output_file_name = output_file_name
        self._generate_output_file_header()
        self.in_state = 1
        self.time = 0.0
        self.bounce_counter = 0
        self.timestep = 0.002
        self.N_i_j_alpha_dict, self.R_i_alpha_dict, self.N_alpha_beta_dict, self.T_alpha \
            = self.calc.regions[self.index].produce_mmvt_statistics(self.index)
        return
        
    def _generate_output_file_header(self):
        if self.style=="openmm":
            header = "#\"Bounced boundary ID\",\"bounce index\",\"total time (ps)\n"
        elif self.style=="namd":
            header = "# NAMD TEST OUTPUT\n"
        else:
            header = "UNKNOWN OUTPUT STYLE"
        with open(self.output_file_name, "w") as f:
            f.write(header)
        return
    
    def _write_bounce_to_output_file(self, line):
        with open(self.output_file_name, "a") as f:
            f.write(line)
        return
    """ # TODO: remove if the other version is good enough
    def _write_random_transition(self, starting_step, num_steps):
        # TODO: perhaps this algorithm can be improved by not sampling times
        #  for both transitions and bounces, but rather using the number
        #  of transitions per bounce to decide when a bounce becomes a 
        #  transition?
        MAX_ITER = 1e9
        step_time = 0.0
        N_i_j_alpha_dict, R_i_alpha_dict, N_alpha_beta_dict, T_alpha \
            = self.calc.regions[self.index].produce_mmvt_statistics(self.index)
        available_transition_dict = {}
        available_transitions = []
        transition_probabilities = []
        total_transitions_out_of_state = 0
        for key in N_i_j_alpha_dict:
            if key[0] == self.in_state:
                available_transition_dict[key[1]] = N_i_j_alpha_dict[key]
                total_transitions_out_of_state += N_i_j_alpha_dict[key]
        
        probability_sum = sum(available_transition_dict.values())
        for key in available_transition_dict:
            available_transitions.append(key)
            transition_probabilities.append(available_transition_dict[key]/probability_sum)
        
        if len(available_transitions) == 0:
            time_in_this_transition = 1e9
            avg_time_between_self_bounces = T_alpha / N_alpha_beta_dict[self.in_state]
        else:
            next_state = np.random.choice(a=np.array(available_transitions), 
                                          p=np.array(transition_probabilities))
            avg_transition_time = R_i_alpha_dict[self.in_state] / total_transitions_out_of_state
            time_in_this_transition = np.random.exponential(avg_transition_time)
            avg_time_between_self_bounces = R_i_alpha_dict[self.in_state] / N_alpha_beta_dict[self.in_state]
            
        # sample times in while loop until time_between_transitions is exceeded
        counter = 0
        while True:
            time_in_this_self_bounce = np.random.exponential(avg_time_between_self_bounces)
            step_time += time_in_this_self_bounce
            if step_time > time_in_this_transition:
                break
            if self.style=="openmm":
                bounce_str = "{},{},{}\n".format(self.in_state, self.bounce_counter, self.time+step_time)
                
            #print(bounce_str)
            self._write_bounce_to_output_file(bounce_str)
            self.bounce_counter += 1
            if self.bounce_counter >= starting_step + num_steps:
                return
            if counter > MAX_ITER: 
                raise Exception("Max iterations reached.")
            counter += 1
        
        self.time += time_in_this_transition   
        self.in_state = next_state 
        if self.style=="openmm":
            bounce_str = "{},{},{}\n".format(self.in_state, self.bounce_counter, self.time)
        #print(bounce_str)
        self._write_bounce_to_output_file(bounce_str)
        self.bounce_counter += 1
        return
    """
    
    def _write_random_transition(self, starting_step, num_steps):
        # TODO: perhaps this algorithm can be improved by not sampling times
        #  for both transitions and bounces, but rather using the number
        #  of transitions per bounce to decide when a bounce becomes a 
        #  transition?
        MAX_ITER = 1e9
        step_time = 0.0
        
        available_transition_dict = {}
        available_transitions = []
        transition_probabilities = []
        total_transitions_out_of_state = 0
        for key in self.N_i_j_alpha_dict:
            if key[0] == self.in_state:
                available_transition_dict[key[1]] = self.N_i_j_alpha_dict[key]
                total_transitions_out_of_state += self.N_i_j_alpha_dict[key]
        
        probability_sum = sum(available_transition_dict.values())
        for key in available_transition_dict:
            available_transitions.append(key)
            transition_probabilities.append(available_transition_dict[key]/probability_sum)
        
        if len(available_transitions) == 0:
            time_in_this_transition = 1e9
            avg_time_between_self_bounces = self.T_alpha / self.N_alpha_beta_dict[self.in_state]
            bounces_per_transition = 1e9
        else:
            next_state = np.random.choice(a=np.array(available_transitions), 
                                          p=np.array(transition_probabilities))
            avg_transition_time = self.R_i_alpha_dict[self.in_state] / total_transitions_out_of_state
            time_in_this_transition = np.random.exponential(avg_transition_time)
            avg_time_between_self_bounces = self.R_i_alpha_dict[self.in_state] / self.N_alpha_beta_dict[self.in_state]
            avg_bounces_per_transition = self.N_alpha_beta_dict[self.in_state] / total_transitions_out_of_state
            bounces_per_transition = np.random.exponential(avg_bounces_per_transition)
            
        # sample times in while loop until time_between_transitions is exceeded
        counter = 1
        while True:
            if counter > bounces_per_transition:
                break
            time_in_this_self_bounce = np.random.exponential(avg_time_between_self_bounces)
            step_time += time_in_this_self_bounce
            
            if self.style=="openmm":
                bounce_str = "{},{},{}\n".format(self.in_state, self.bounce_counter, self.time+step_time)
            elif self.style=="namd":
                if self.in_state == 1:
                    new_anchor = self.index - 1
                    new_milestone = self.index - 1
                elif self.in_state == 2:
                    new_anchor = self.index + 1
                    new_milestone = self.index
                    
                if new_anchor == -1:
                    new_anchor = 1
                    new_milestone = self.index
                    
                template="SEEKR: Cell Collision: current: {}, new: {}, stepnum: {}\n"
                bounce_str = template.format(self.index, new_anchor, int((self.time+step_time)/self.timestep))
                if self.bounce_counter == 0:
                    template="SEEKR: Milestone Transition: anchor: {}, source: {}, destination: {}, stepnum: {}, incubation steps: {}\n"
                    bounce_str += template.format(
                        self.index, "none", new_milestone, 
                        int((self.time+step_time)/self.timestep), 
                        int(step_time/self.timestep))
            
            self._write_bounce_to_output_file(bounce_str)
            self.bounce_counter += 1
            
            if self.bounce_counter >= starting_step + num_steps:
                return
            
            if counter > MAX_ITER:
                raise Exception("Max iterations exceeded.")
            counter += 1
        
        self.time += step_time   
        prev_state = self.in_state
        self.in_state = next_state 
        if self.style=="openmm":
            bounce_str = "{},{},{}\n".format(self.in_state, self.bounce_counter, self.time)
        elif self.style=="namd":
            if self.in_state == 1:
                new_anchor = self.index - 1
                old_milestone = self.index
                new_milestone = self.index - 1
            elif self.in_state == 2:
                new_anchor = self.index + 1
                old_milestone = self.index - 1
                new_milestone = self.index
            if self.bounce_counter == 0:
                old_milestone = "none"
            template="SEEKR: Cell Collision: current: {}, new: {}, stepnum: {}\n"\
                     +"SEEKR: Milestone Transition: anchor: {}, source: {}, destination: {}, stepnum: {}, incubation steps: {}\n"
            bounce_str = template.format(self.index, new_anchor, int((self.time)/self.timestep),
                                         self.index, old_milestone, new_milestone, int((self.time)/self.timestep), int(step_time/self.timestep))
        self._write_bounce_to_output_file(bounce_str)
        #print("bounce_str:", bounce_str)
        self.bounce_counter += 1
        return
        
    def step(self, number):
        starting_step = self.bounce_counter
        while self.bounce_counter < number:
            self._write_random_transition(starting_step, number)
            
        return
    
class SmoluchowskiSphericalElberIntegrator(ToyIntegrator):
    """
    An integrator that generates MMVT output files based on a Smoluchowski
    region.
    """
    def __init__(self, calc, index, output_file_name, style="openmm"):
        self.calc = calc
        self.index = index
        self.style = style.lower()
        self.output_file_name = output_file_name
        self._generate_output_file_header()
        self.in_state = index
        self.time = 0.0
        self.bounce_counter = 0
        self.timestep = 1.0
        self.elberN_ij, self.elberR_i = self.calc.produce_elber_statistics()
        return
        
    def _generate_output_file_header(self):
        if self.style=="openmm":
            header = "#\"Bounced boundary ID\",\"bounce index\",\"total time (ps)\n"
        else:
            header = "UNKNOWN OUTPUT STYLE"
        with open(self.output_file_name, "w") as f:
            f.write(header)
        return
    
    def _write_bounce_to_output_file(self, line):
        with open(self.output_file_name, "a") as f:
            f.write(line)
        return
    
    def _write_random_transition(self, starting_step, num_steps):
        # TODO: perhaps this algorithm can be improved by not sampling times
        #  for both transitions and bounces, but rather using the number
        #  of transitions per bounce to decide when a bounce becomes a 
        #  transition?
        
        available_transition_dict = {}
        available_transitions = []
        transition_probabilities = []
        total_transitions_out_of_state = 0
        for key in self.elberN_ij:
            if key[0] == self.in_state:
                available_transition_dict[key[1]] = self.elberN_ij[key]
                total_transitions_out_of_state += self.elberN_ij[key]
        
        probability_sum = sum(available_transition_dict.values())
        for key in available_transition_dict:
            available_transitions.append(key)
            transition_probabilities.append(available_transition_dict[key]/probability_sum)
        
        next_state = np.random.choice(a=np.array(available_transitions), 
                                        p=np.array(transition_probabilities))
        avg_transition_time = self.elberR_i[self.in_state]
        step_time = np.random.exponential(avg_transition_time)
        
        transition_alias = next_state - self.in_state + 2
        self.time += step_time   
        
        if self.style=="openmm":
            bounce_str = "{},{},{}\n".format(transition_alias, self.bounce_counter, step_time)
        #print(bounce_str)
        self._write_bounce_to_output_file(bounce_str)
        self.bounce_counter += 1
        return
        
    def step(self, number):
        starting_step = self.bounce_counter
        while self.bounce_counter < number:
            self._write_random_transition(starting_step, number)
            
        return
    
