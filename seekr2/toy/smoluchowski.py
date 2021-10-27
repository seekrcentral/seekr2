"""
smoluchowski.py

Define the object used by the Toy engine's Smoluchowski integrator, as well
as some of SEEKR2's tests.
"""

from collections import defaultdict

import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as plt

def expBetaW(w, beta):
    return np.exp(-beta * w)

class FlatPotentialEnergyFunction():
    """
    A completely flat potential energy surface.
    """
    def __init__(self):
        pass
    
    def evaluate_force(self, variables):
        return 0.0
    
    def evaluate_energy(self, variables):
        return 0.0
    
class LinearPotentialEnergyFunction():
    """
    A linear potential energy surface.
    """
    def __init__(self, a=1.0):
        self.a = a
    
    def evaluate_force(self, variable):
        return -self.a
    
    def evaluate_energy(self, variable):
        return self.a * variable
    
class QuadraticPotentialEnergyFunction():
    """
    A quadratic potential energy surface.
    """
    def __init__(self, a=1.0):
        self.a = a
    
    def evaluate_force(self, variable):
        return -2.0 * self.a * variable
    
    def evaluate_energy(self, variable):
        return self.a * variable**2
    
class CoulombicPotentialEnergyFunction():
    """
    A quadratic potential energy surface.
    """
    def __init__(self, q1q2=1.0, a=1.0):
        self.q1q2 = q1q2
        self.a = a
    
    def evaluate_force(self, variable):
        return 0.0 # TODO: fill out if necessary
    
    def evaluate_energy(self, variable):
        r = variable
        if r < self.a: # make a flat-bottomed potential
            r = self.a
        return self.q1q2 / (4.0 * np.pi * r)

class System():
    def __init__(self, particles, potential_energy_function):
        self.particles = particles
        self.potential_energy_function = potential_energy_function
        return
    
class Integrator():
    def __init__(self):
        pass

#TODO: make linear version?

class SmoluchowskiSphericalRegion():
    def __init__(self, a, b, potential_energy_function, diffusion, beta, n=101):
        self.a = a
        self.b = b
        self.potential_energy_function = potential_energy_function
        self.diffusion = diffusion
        self.beta = beta
        self.n = n
        self.h = (self.b-self.a)/(self.n-1)
        self.r_s = np.zeros(self.n) #np.arange(self.a, self.b+self.h, self.h)
        for i in range(n):
            self.r_s[i] = self.a + i*self.h
        assert self.r_s.shape[0] == self.n
        self.u_r_outward = None
        self.J_outward = None
        self.u_r_inward = None
        self.J_inward = None
        self.partition_function = self.make_partition_function()
        self.outward_A1 = None
        self.inward_A1 = None
        self.compute_outward_flux()
        self.compute_inward_flux()
        self.weight = 0.0
        self.N_alpha_down = None
        self.N_alpha_up = None
        self.N_alpha_down_up = 0.0
        self.N_alpha_up_down = 0.0
        self.R_alpha = None
        self.make_mmvt_statistics()
        
        return
    
    def compute_outward_flux(self):
        r_s_forward = self.r_s
        self.u_r_outward = np.zeros(self.n)
        
        # Compute denominator
        denominator_vals = np.zeros(self.n)
        for i, r in enumerate(r_s_forward):
            jacobian = r**2
            w = self.potential_energy_function.evaluate_energy(r)
            denominator_vals_prime = np.zeros(self.n - i)
            for j, r_prime in enumerate(r_s_forward[i:]):
                if np.isclose(r_prime, 0.0):
                    denominator_vals_prime[j] = 0.0
                    continue
                
                w2 = self.potential_energy_function.evaluate_energy(r_prime)
                denominator_vals_prime[j] = 1.0 / (r_prime**2 * expBetaW(w2, self.beta))
            
            denominator_val_prime = simps(denominator_vals_prime, dx=self.h)
            denominator_vals[i] = jacobian*expBetaW(w, self.beta)*denominator_val_prime
        
        outer_integral = simps(denominator_vals, dx=self.h)
        denominator = 4.0*np.pi*outer_integral
        self.outward_A1 = 1.0/denominator
        if np.isclose(denominator, 0.0):
            print("smoluchowski.py: denominator is nearly zero:", denominator)
            
        # Compute numerator
        for i, r in enumerate(r_s_forward):
            w = self.potential_energy_function.evaluate_energy(r)
            numerator_vals_forward = np.zeros(self.n - i)
            for j, r_prime in enumerate(r_s_forward[i:]):
                if np.isclose(r_prime, 0.0):
                    numerator_vals_forward[j] = 0.0
                    continue
                w2 = self.potential_energy_function.evaluate_energy(r_prime)
                numerator_vals_forward[j] = 1.0 / (r_prime**2 * expBetaW(w2, self.beta))
                
            numerator_val_forward = simps(numerator_vals_forward, dx=self.h)
            self.u_r_outward[i] = expBetaW(w, self.beta)*numerator_val_forward / denominator
        
        self.J_outward = self.diffusion / (self.b**2 * denominator)
        return
    
    def compute_inward_flux(self):
        r_s_forward = self.r_s #np.arange(self.a, self.b+self.h, self.h)
        self.u_r_inward = np.zeros(self.n)
        if np.isclose(self.a, 0.0):
            self.J_inward = 0.0
            for i, r in enumerate(r_s_forward):
                w = self.potential_energy_function.evaluate_energy(r)
                self.u_r_inward[i] = expBetaW(w, self.beta) / self.partition_function
            return
        
        # Compute denominator
        denominator_vals = np.zeros(self.n)
        for i, r in enumerate(r_s_forward):
            jacobian = r**2
            w = self.potential_energy_function.evaluate_energy(r)
            denominator_vals_prime = np.zeros(i+1)
            for j, r_prime in enumerate(r_s_forward[:i+1]):
                if np.isclose(r_prime, 0.0):
                    denominator_vals_prime[j] = 0.0
                    continue
                w2 = self.potential_energy_function.evaluate_energy(r_prime)
                denominator_vals_prime[j] = 1.0 / (r_prime**2 * expBetaW(w2, self.beta))
            
            denominator_val_prime = simps(denominator_vals_prime, dx=self.h)
            denominator_vals[i] = jacobian*expBetaW(w, self.beta)*denominator_val_prime
        
        denominator = 4.0*np.pi*simps(denominator_vals, dx=self.h)
        self.inward_A1 = 1.0/denominator
        if np.isclose(denominator, 0.0):
            print("smoluchowski.py denominator is nearly zero:", denominator)
        
        for i, r in enumerate(r_s_forward):
            w = self.potential_energy_function.evaluate_energy(r)
            numerator_vals_backward = np.zeros(i+1)
            for j, r_prime in enumerate(r_s_forward[:i+1]):
                if np.isclose(r_prime, 0.0):
                    numerator_vals_backward[j] = 0.0
                    continue
                w2 = self.potential_energy_function.evaluate_energy(r_prime)
                numerator_vals_backward[j] = 1.0 / (r_prime**2 * expBetaW(w2, self.beta))
                
            numerator_val_backward = simps(numerator_vals_backward, dx=self.h)
            self.u_r_inward[i] = expBetaW(w, self.beta)*numerator_val_backward / denominator
            
        self.J_inward = self.diffusion / (self.a**2 * denominator)
        return
    
    def check_solution(self, u_q, a_boundary_value=None, b_boundary_value=None):
        # make sure u_q values are solutions to the Smol equation
        r_s_forward = self.r_s #np.arange(self.a, self.b+self.h, self.h)
        for i, r in enumerate(r_s_forward):
            if i == 0 or i == 1 or i == r_s_forward.shape[0]-1 or i == r_s_forward.shape[0]-2:
                continue
            r_plus = r + self.h
            r_plus2 = r + 2.0*self.h
            r_minus = r - self.h
            r_minus2 = r - 2.0*self.h
            w_center = self.potential_energy_function.evaluate_energy(r)
            w_plus2 = self.potential_energy_function.evaluate_energy(r_plus2)
            w_minus2 = self.potential_energy_function.evaluate_energy(r_minus2)
            u_r_inward_center = u_q[i]
            u_r_inward_plus = u_q[i+1]
            u_r_inward_plus2 = u_q[i+2]
            u_r_inward_minus = u_q[i-1]
            u_r_inward_minus2 = u_q[i-2]
            grad_u_plus = (u_r_inward_plus2 - u_r_inward_center) / (2.0*self.h)
            grad_u_minus = (u_r_inward_center - u_r_inward_minus2) / (2.0*self.h)
            grad_w_plus = (w_plus2 - w_center) / (2.0*self.h)
            grad_w_minus = (w_center - w_minus2) / (2.0*self.h)
            
            term1 = r_plus**2*(grad_u_plus + self.beta*u_r_inward_plus*grad_w_plus)
            term2 = r_minus**2*(grad_u_minus + self.beta*u_r_inward_minus*grad_w_minus)
            
            smol_eq = (term1 - term2) / (2.0*self.h)
            if not np.isclose(smol_eq, 0.0, atol=0.001):
                print("Smol equation far from zero: {}".format(smol_eq))
        
        # Check boundary values
        if a_boundary_value is not None:
            if not np.isclose(u_q[0], a_boundary_value, atol=0.001):
                print("Boundary a far from condition: was {}, expected {}".format(u_q[0], a_boundary_value))
                
        if b_boundary_value is not None:
            if not np.isclose(u_q[-1], b_boundary_value, atol=0.001):
                print("Boundary b far from condition: was {}, expected {}".format(u_q[-1], b_boundary_value))
                
        # Check area under curve
        area_vals = np.zeros(self.n)
        for i, r in enumerate(r_s_forward):
            jacobian = r**2
            area_vals[i] = jacobian * u_q[i]
        
        area = 4.0*np.pi*simps(area_vals, dx=self.h)
        if not np.isclose(area, 1.0, atol=0.001):
            print("Area far from 1.0: {}".format(area))
            
        return
    
    def make_partition_function(self):
        r_s_forward = self.r_s #np.arange(self.a, self.b+self.h, self.h)
        area_vals = np.zeros(self.n)
        for i, r in enumerate(r_s_forward):
            jacobian = r**2
            w = self.potential_energy_function.evaluate_energy(r)
            area_vals[i] = jacobian*expBetaW(w, self.beta)
            
        area = 4.0*np.pi*simps(area_vals, dx=self.h)
        return area
    
    """
   def make_mmvt_statistics(self):
        w_a = self.potential_energy_function.evaluate_energy(self.a)
        w_b = self.potential_energy_function.evaluate_energy(self.b)
        #self.N_alpha_down = self.a**2 * expBetaW(w_a, self.beta)
        #self.N_alpha_up = self.b**2 * expBetaW(w_b, self.beta)
        self.T_alpha = self.partition_function
        #self.k_alpha_down = self.N_alpha_down / self.T_alpha
        #self.k_alpha_up = self.N_alpha_up / self.T_alpha
        # I denotes current - flux J integrated over area
        I_outward_weighted = 4.0 * np.pi * self.b**2 * self.J_outward
        I_inward_weighted = 4.0 * np.pi * self.a**2 * self.J_inward
        
        denominator = (I_outward_weighted + I_inward_weighted)
        if np.isclose(I_inward_weighted, 0.0):
            self.R_alpha_up_down = 0.0
            self.N_alpha_up_down = 0.0
        else:
            self.R_alpha_up_down = 1.0 / I_inward_weighted
            self.N_alpha_up_down = 1.0 #I_inward_weighted / denominator
            
        self.N_alpha_down_up = 1.0 #I_outward_weighted / denominator
        self.R_alpha = 1.0 / denominator
        self.R_alpha_down_up = 1.0 / I_outward_weighted
        return
    """
    
    # fix this better tomorrow...
    
    def make_mmvt_statistics(self):
        w_a = self.potential_energy_function.evaluate_energy(self.a)
        w_b = self.potential_energy_function.evaluate_energy(self.b)
        #self.N_alpha_down = self.a**2 * expBetaW(w_a, self.beta)
        #self.N_alpha_up = self.b**2 * expBetaW(w_b, self.beta)
        self.T_alpha = 1000.0 #self.partition_function
        #self.k_alpha_down = self.N_alpha_down / self.T_alpha
        #self.k_alpha_up = self.N_alpha_up / self.T_alpha
        # I denotes current - flux J integrated over area
        I_outward_weighted = 4.0 * np.pi * self.b**2 * self.J_outward
        I_inward_weighted = 4.0 * np.pi * self.a**2 * self.J_inward
        
        total_flux = (I_outward_weighted + I_inward_weighted)
        if np.isclose(I_inward_weighted, 0.0):
            self.R_alpha_up_down = 0.0
            self.N_alpha_up_down = 0.0
            
            #time_per_transition_down_up = 1.0 / I_outward_weighted
            #self.N_alpha_down_up = self.T_alpha / time_per_transition_down_up #1.0 
            #self.R_alpha = #1.0 / denominator
            self.R_alpha_down_up = self.T_alpha
        else:
            time_per_transition_up_down = 1.0 / I_inward_weighted
            time_per_transition_down_up = 1.0 / I_outward_weighted
            time_per_transition_sum = time_per_transition_up_down + time_per_transition_down_up
            time_fraction_up_down = time_per_transition_up_down / time_per_transition_sum
            time_fraction_down_up = time_per_transition_down_up / time_per_transition_sum
            self.R_alpha_up_down = self.T_alpha * time_fraction_up_down
            self.N_alpha_up_down =  self.T_alpha / time_per_transition_sum
            
            self.R_alpha_down_up = self.T_alpha * time_fraction_down_up
            self.N_alpha_down_up = self.T_alpha / time_per_transition_sum
            self.R_alpha = None
            
        return
    
    def plot_functions(self):
        r_s_forward = self.r_s # np.arange(self.a, self.b+self.h, self.h)
        
        fig, ax = plt.subplots()
        ax.plot(r_s_forward, self.u_r_outward)
        plt.ylabel("$u_r$")
        plt.xlabel("r")
        plt.title("Outward flux concentration")
        plt.show()
        
        fig, ax = plt.subplots()
        ax.plot(r_s_forward, self.u_r_inward)
        plt.ylabel("$u_r$")
        plt.xlabel("r")
        plt.title("Inward flux concentration")
        plt.show()
        return
    
    def produce_mmvt_statistics(self, i):
        T_alpha = self.T_alpha
        N_backwards = self.N_alpha_down
        N_forwards = self.N_alpha_up
        R_i_backwards = self.R_alpha_up_down
        R_i_forwards = self.R_alpha_down_up
        N_ij_backwards = self.N_alpha_up_down
        N_ij_forwards = self.N_alpha_down_up
        
        N_i_j_alpha_dict = defaultdict(int)
        R_i_alpha_dict = defaultdict(float)
        N_alpha_beta_dict = defaultdict(int)
        new_time_factor = (R_i_forwards + R_i_backwards) / T_alpha
        new_T_alpha = new_time_factor * T_alpha
        if i == 0:
            N_alpha_beta_dict[1] = N_forwards #new_time_factor
            R_i_alpha_dict[1] = new_T_alpha
        else:
            N_i_j_alpha_dict[(1, 2)] = N_ij_forwards
            N_i_j_alpha_dict[(2, 1)] = N_ij_backwards
            R_i_alpha_dict[1] = R_i_forwards
            R_i_alpha_dict[2] = R_i_backwards
            N_alpha_beta_dict[1] = N_backwards # * new_time_factor
            if N_forwards is not None:
                N_alpha_beta_dict[2] = N_forwards # * new_time_factor
            else:
                N_alpha_beta_dict[2] = N_backwards
        
        return N_i_j_alpha_dict, R_i_alpha_dict, N_alpha_beta_dict, new_T_alpha
        
"""
class SmoluchowskiSphericalMMVTIntegrator(Integrator):
    def __init__(self, temperature):
        self.boltzmanns_constant
        self.temperature = temperature
        self.beta = 1.0 / (self.boltzmanns_constant * self.temperature)
        return
    
    def prepare(self, system, inner_boundary_position, outer_boundary_position):
        assert len(system.particles) == 1
        particle = system.particles[0]
        assert particle.dimensions == 1
        
class Simulation():
    def __init__(self):
        pass
"""

class SmoluchowskiCalculation1d():
    def __init__(self, potential_energy_function, milestones, 
                 absorbing_boundary, reflecting_boundary=0.0, diffusion=1.0,
                 beta=1.0, n=101):
        self.potential_energy_function = potential_energy_function
        self.milestones = milestones
        self.absorbing_boundary = absorbing_boundary
        self.reflecting_boundary = reflecting_boundary
        self.regions = []
        self.diffusion = diffusion
        self.beta = beta
        self.n = n
        self.partition_function = 0.0
        self.make_regions()
        self.make_elber_statistics()
        self.make_state_transition_rates()
        
        return
    
    def make_regions(self):
        first_region = SmoluchowskiSphericalRegion(
            self.reflecting_boundary, self.milestones[0], 
            self.potential_energy_function, self.diffusion, self.beta, self.n)
        
        self.regions.append(first_region)
        self.partition_function += first_region.partition_function
        
        for i, milestone in enumerate(self.milestones[:-1]):
            region = SmoluchowskiSphericalRegion(
                self.milestones[i], self.milestones[i+1], 
                self.potential_energy_function, self.diffusion, self.beta, 
                self.n)
            
            self.regions.append(region)
            self.partition_function += region.partition_function
            
        last_region = SmoluchowskiSphericalRegion(
            self.milestones[-1], self.absorbing_boundary,
            self.potential_energy_function, self.diffusion, self.beta, self.n)
        self.regions.append(last_region)
        self.partition_function += last_region.partition_function
        
        for region in self.regions:
            region.weight = region.partition_function / self.partition_function
        return
    
    def make_state_transition_rates(self):
        #self.downward_rates = []
        #self.upward_rates = []
        N_MAGNITUDE = 100.0 # 100.0
        prev_rate = 0.0
        for i in range(len(self.regions)-1):
            a = self.regions[i].partition_function
            b = self.regions[i+1].partition_function
            if self.regions[i].N_alpha_down_up > 1.0:
                down_rate = N_MAGNITUDE*self.regions[i].N_alpha_down_up
            else:
                down_rate = 10000.0
            up_rate_for_prev = (b)*down_rate - prev_rate
            up_rate = (b/a)*down_rate - prev_rate/a
            #self.downward_rates.append(down_rate)
            #self.upward_rates.append(up_rate)
            self.regions[i].N_alpha_up = up_rate
            self.regions[i+1].N_alpha_down = down_rate
            prev_rate = -up_rate_for_prev + b*down_rate
            
        #print("self.upward_rates:", self.upward_rates)
        #print("self.downward_rates:", self.downward_rates)
        return
            
    def make_elber_statistics(self):
        self.n_list = []
        self.t_list = []
        for i, region in enumerate(self.regions[:-1]):
            lower_region = region
            upper_region = self.regions[i+1]
            if i == 0:
                n_up = 1.0 #upper_region.u_r_inward[0]
                self.n_list.append((n_up,))
            else:
                n_down = 1.0 #lower_region.u_r_inward[-1]
                n_up = 1.0 #upper_region.u_r_inward[0]
                self.n_list.append((n_down, n_up))
            time = 1.0 / (lower_region.J_inward + upper_region.J_outward)
            self.t_list.append(time)
            
    def produce_elber_statistics(self):
        elberN_ij = defaultdict(float)
        elberR_i = defaultdict(float)
        num_milestones = len(self.regions)
        for i, region1 in enumerate(self.regions):
            if i == 0:
                region2 = self.regions[i+1]
                elberN_ij[(0,1)] = 1.0
                # need to make sure that u and exp(-beta*W) match up
                #  on the edge.
                w_b1 = self.potential_energy_function.evaluate_energy(region1.b)
                region1_edge_value = expBetaW(w_b1, region1.beta)
                region2_edge_value = region2.u_r_outward[0]
                region2_current_I = 4.0 * np.pi * region2.b**2 * region2.J_outward
                elberR_i[0] = (region1.partition_function * region2_edge_value/region1_edge_value + 1.0) / (region2_current_I)
                
            elif i == num_milestones-1:
                elberN_ij[(num_milestones-1, num_milestones-2)] = 1.0
                region1_current_I = 4.0 * np.pi * region1.a**2 * region1.J_inward / region1_edge_value
                elberR_i[i] = 1.0 / region1_current_I
                    
            else:
                region2 = self.regions[i+1]
                region1_edge_value = region1.u_r_inward[-1]
                region2_edge_value = region2.u_r_outward[0]
                region1_current_I = 4.0 * np.pi * region1.a**2 * region1.J_inward / region1_edge_value
                region2_current_I = 4.0 * np.pi * region2.b**2 * region2.J_outward / region2_edge_value
                new_total_volume = 1/region1_edge_value + 1/region2_edge_value
                
                total_flux = region2_current_I + region1_current_I
                elberN_ij[(i,i+1)] = region2_current_I / total_flux
                elberN_ij[(i,i-1)] = region1_current_I / total_flux 
                elberR_i[i] = new_total_volume / total_flux
                
        return elberN_ij, elberR_i
            
def make_smoluchowski_calculation_from_model(model, potential_energy_function,
                                             beta=1.0, diffusion=1.0):
    milestones = []
    milestone_dict = {}
    for i, anchor in enumerate(model.anchors):
        for milestone in anchor.milestones:
            milestone_dict[milestone.index] = milestone.variables["radius"]
            
    for i in range(len(milestone_dict)):
        milestones.append(milestone_dict[i])
    
    absorbing_boundary = milestones.pop()
    calc = SmoluchowskiCalculation1d(
        potential_energy_function, milestones=milestones, 
        absorbing_boundary=absorbing_boundary, beta=beta, diffusion=diffusion)
    return calc

if __name__ == "__main__":
    potential_energy_function = FlatPotentialEnergyFunction()
    #potential_energy_function = LinearPotentialEnergyFunction()
    #potential_energy_function = QuadraticPotentialEnergyFunction(a=0.1)
    smol = SmoluchowskiSphericalRegion(1.0, 2.0, potential_energy_function, 1.0, 1.0)
    
    #print("smol.k_alpha_down: ", smol.k_alpha_down)
    #print("smol.k_alpha_up: ", smol.k_alpha_up)
    print("smol.N_alpha_down_up: ", smol.N_alpha_down_up)
    print("smol.N_alpha_up_down: ", smol.N_alpha_up_down)
    print("smol.R_alpha: ", smol.R_alpha)
    
    
    print("checking inward")
    smol.check_solution(smol.u_r_inward)
    print("checking_outward")
    smol.check_solution(smol.u_r_outward)
        
    #smol.plot_functions()
    
    calc = SmoluchowskiCalculation1d(
        potential_energy_function, milestones=[1.0, 2.0, 3.0], 
        absorbing_boundary=4.0)
    
    for i, region in enumerate(calc.regions):
        print("region: {}".format(i))
        print("  region.k_alpha_down: ", region.k_alpha_down)
        print("  region.k_alpha_up: ", region.k_alpha_up)
        print("  region.N_alpha_down_up: ", region.N_alpha_down_up)
        print("  region.N_alpha_up_down: ", region.N_alpha_up_down)
        print("  region.R_alpha: ", region.R_alpha)
        
    print("calc.n_list:", calc.n_list)
    print("calc.t_list:", calc.t_list)