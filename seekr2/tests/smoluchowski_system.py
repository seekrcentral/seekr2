"""
Generate a system whose values are computed using the Smoluchowski
equation for exact verification of the analysis methods.

The time-dependent Smoluchowski Equation is:

du/dt = -dJ/dx = d(D(z) * du/dx -v(x,t) * u)/dx

Where u = u(x,t) is the time-dependent concentration, 
J = J(x,t) is a time-dependent flux, D(z) is a local
diffusivity, and v(x,t) is flow velocity.

A number of assumptions are made in this program, including
that D(x) is a constant D, that the flow velocity depends only
on the potential of mean force W(x) on the system: 

v(x) = -D * dW(x)/dx * (1/kT)

Where k is Boltzmann's constant, and T is temperature. The value
1/kT is often abbreviated "beta".

In addition, if concentration and flux is constant in time, 
then the time-independent Smoluchowski equation may be used.

0 = d(D * (du/dx + beta * dW/dx * u))/dx

In 1D, the time-independent equation above can be fairly 
easily solved using simple ordinary differential equation (ODE) 
methods.

In contrast, the time-dependent system may also be solved, but
with more involved, trickier partial differential equation (PDE)
solution methods.

For the time-independent situation, the general solution becomes:

u(x) = A1 * exp(-beta*W(x))*(integrate_0^x (exp(beta*W(x') / D) dx'
    + A2 * exp(-beta*W(x))

Where A1 and A2 are constants of integration, and may be found
by applying boundary conditions. 

This script performs both ODE and PDE to obtain useful information
about the system.

"""

from math import exp, sin, cos, log
from collections import defaultdict

import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as plt

import seekr2.analyze as analyze
import seekr2.modules.common_base as base
import seekr2.modules.mmvt_base as mmvt_base

beta = 0.1
max_time = 10.0
time_points = 101
k = (max_time)/(time_points-1)
max_coeffs = 100

def expW_constant(x, q=1.0):
    """ A Boltzman distribution of a constant PMF."""
    return 1.0 / q

def expW_linear(x, q=1.0):
    """A Boltzmann distribution of a linearly-increasing PMF."""
    return exp(-beta * abs(x)) / q

def expW_quadratic(x, q=1.0):
    """A Boltzmann distribution of a quadratic PMF function."""
    return exp(-beta * x**2) / q

def get_partition_function(func, a, b, n):
    """
    Given a Boltzmann distribution function, and a domain, compute the
    partition function of the distribution.
    """
    h = (b-a)/(n-1)
    x_s = np.arange(a, b+h, h)
    func_vals = np.zeros(n)
    for i, x in enumerate(x_s):
        func_vals[i] = func(x, 1.0)
    return simps(func_vals, dx=h)

class Smoluchowski():
    """
    A 1D system whose dynamics are described by the Smoluchowski
    equation, depending on a potential of mean force (PMF) function:
    W(x), and a constant diffusivity D.
    """
    def __init__(self, a, b, expWfunc, n=101, D=1.0):
        self.a = a
        self.b = b
        self.span = b-a
        self.n = n
        self.h = (b-a)/(n-1)
        self.expWfunc = expWfunc
        self.expWq = get_partition_function(self.expWfunc, a, b, n)
        self.D = D
        self.coeffs = []
        self.x_s = np.arange(self.a, self.b+self.h, self.h)
        self.u_q_forward = 1.0
        self.u_q_backward = 1.0
        self.u_x_forward = None
        self.u_x_backward = None
        self.J_backward = None
        self.J_forward = None
        for m in range(1,max_coeffs+1):
            A_m = self.find_fourier_coeff(self.expWfunc, m)
            self.coeffs.append(A_m)
        self.reflect_lower = False
        self.fill_out_flux_conc()
        
        return
    
    def fill_out_flux_conc(self):
        """
        For this Smoluchowski domain, calculate the concentration given
        that flux will be flowing
        """
        h = (self.b-self.a)/(self.n-1)
        x_s_forward = np.arange(self.a, self.b+h, h)
        x_s_backward = np.arange(self.b+h, self.a, -h)
        #x_s = np.arange(self.a, self.b+h, h)
        self.u_x_forward = np.zeros(self.n)
        self.u_x_backward = np.zeros(self.n)
        
        denominator_vals = np.zeros(self.n)
        for i, x_double_prime in enumerate(x_s_forward):
            denominator_vals[i] = 1.0/self.expWfunc(x_double_prime)
            
        denominator = simps(denominator_vals, dx=h)
        
        for i, x in enumerate(x_s_forward):
            numerator_vals_forward = np.zeros(self.n - i)
            for j, x_prime in enumerate(x_s_forward[i:]):
                numerator_vals_forward[j] = 1.0/self.expWfunc(x_prime)
            integrated = simps(numerator_vals_forward, dx=h)
            self.u_x_forward[i] = (self.expWfunc(x) / self.expWfunc(self.a)) * integrated / denominator
        
        self.u_q_forward = simps(self.u_x_forward, dx=h)
        self.J_forward = self.D * (1.0 / self.expWfunc(self.a)) / denominator
        
        for i, x in enumerate(x_s_forward):
            numerator_vals_backward = np.zeros(i+1)
            for j, x_prime in enumerate(x_s_forward[:i+1]):
                numerator_vals_backward[j] = 1.0/self.expWfunc(x_prime)
            integrated = simps(numerator_vals_backward, dx=h)
            self.u_x_backward[i] = (self.expWfunc(x) / self.expWfunc(self.b)) * integrated / denominator
            
        self.u_q_backward = simps(self.u_x_backward, dx=h)
        self.J_backward = self.D * (1.0 / self.expWfunc(self.b)) / denominator
        return
    
    def find_fourier_coeff(self, func, m, a=None, b=None, offset=0.0,
                           partition_function=None):
        """
        Find the Fourier coefficients needed to recreate the Fourier
        series of the equilibrium distribution.
        """
        if a is None:
            a = self.a
        if b is None:
            b = self.b
        h = (b-a)/(self.n-1)
        span = b-a
        x_s = np.arange(a, b+h, h)
        func_vals = np.zeros(self.n)
        if partition_function is None:
            partition_function = self.expWq
        for i, x in enumerate(x_s):
            func_vals[i] = func(x+offset, partition_function) \
                * sin(np.pi*m*x/span)
            
        coeff = 2.0*simps(func_vals, dx=h)/span
        return coeff
    
    def get_total_time(self):
        total_time = 0.0
        if self.reflect_lower:
            coeffs = []
            a = 0
            b = 2.0*self.b
            partition_function = get_partition_function(
                self.expWfunc, a, b, self.n)
            for m in range(1,max_coeffs+1):
                A_m = self.find_fourier_coeff(
                    self.expWfunc, m, a=a, b=b, offset=-self.b, 
                    partition_function=partition_function)
                coeffs.append(A_m)
        else:
            a = self.a
            b = self.b
            coeffs = self.coeffs
        span = b-a
        for j, A_m in enumerate(coeffs):
            m = j+1
            total_time -= A_m * (cos(np.pi*m*b/span)-cos(np.pi*m*a/span)) / m**3
        
        total_time *= span**3 / (np.pi**3*self.D)
        
        return total_time
    
    def compute_MMVT_kinetics_quantities(self):
        """
        Compute quantities that may be used for MMVT calculations.
        """
        J_0_flux_fraction = self.expWfunc(self.a, self.expWq)
        J_span_flux_fraction = self.expWfunc(self.b, self.expWq)
        total_flux = J_0_flux_fraction + J_span_flux_fraction
        J_0_flux_fraction /= total_flux
        J_span_flux_fraction /= total_flux
        
        T_alpha = self.get_total_time()
        N_backwards = J_0_flux_fraction
        N_forwards = J_span_flux_fraction
        k_backwards = N_backwards / T_alpha
        k_forwards = N_forwards / T_alpha
        
        R_i_forwards = self.u_q_forward / self.J_forward
        R_i_backwards = self.u_q_backward / self.J_backward
        if self.reflect_lower:
            N_ij_forwards = 0.0
            N_ij_backwards = 0.0
        else:
            N_ij_forwards = 1.0
            N_ij_backwards = 1.0
        
        return k_backwards, k_forwards, T_alpha, N_backwards, N_forwards, \
            R_i_backwards, R_i_forwards, N_ij_backwards, N_ij_forwards

def make_smol_model(tmp_path, num_anchors, intervals):
    basename = "testmmvt.dat"
    mymodel = base.Model()
    mymodel.temperature = 300.0
    mymodel.calculation_type = "mmvt"
    mymodel.anchor_rootdir = tmp_path
    mymodel.num_anchors = num_anchors+1
    mymodel.num_milestones = mymodel.num_anchors - 1
    
    # TEMPORARY: toy system will eventually have its own settings
    mymodel.openmm_settings = base.Openmm_settings()
    
    intervals.append(1.0)
    
    for index in range(mymodel.num_anchors):
        boundary1 = float(index)
        boundary2 = float(index)+intervals[index]
        x0 = boundary1 + 0.5
        milestones = []
        anchor = mmvt_base.MMVT_anchor()
        if index == 0:
            milestone_count = 1
            milestone1 = base.Milestone()
            milestone1.index = 0
            milestone1.neighbor_anchor_index = 1
            milestone1.alias_index = 1
            milestone1.cv_index = 0
            milestones.append(milestone1)
            anchor.bulkstate = False
            end_state = True
        elif index > 0 and index < mymodel.num_anchors-1:
            milestone_count = 2
            milestone1 = base.Milestone()
            milestone1.index = index-1
            milestone1.neighbor_anchor_index = index-1
            milestone1.alias_index = 1
            milestone1.cv_index = 0
            milestone2 = base.Milestone()
            milestone2.index = index
            milestone2.neighbor_anchor_index = index+1
            milestone2.alias_index = 2
            milestone2.cv_index = 0
            milestones.append(milestone1)
            milestones.append(milestone2)
            anchor.bulkstate = False
            end_state = False
        elif index == mymodel.num_anchors-1:
            milestone_count = 1
            milestone1 = base.Milestone()
            milestone1.index = index-1
            milestone1.neighbor_anchor_index = index-1
            milestone1.alias_index = 1
            milestone1.cv_index = 0
            milestones.append(milestone1)
            anchor.bulkstate = True
            end_state = False
        
        anchor.name = "anchor_%d" % index
        anchor.index = index
        anchor.directory = anchor.name
        anchor.md_mmvt_output_glob = basename
        anchor.md = True
        anchor.endstate = end_state
        anchor.milestones = milestones
        mymodel.anchors.append(anchor)
    return mymodel

if __name__ == "__main__":
    print("The following code computes the mean first passage time")
    print("for a quadratic PMF from an innermost milestone at 1.0, ")
    print("a reflecting boundary at 0.0, and an absorbing boundary")
    print("at 10.0.")
    #func = expW_constant
    #func = expW_linear
    func = expW_quadratic
    
    n = 101
    D = 0.02
    
    a1 = 0.0
    b1 = 1.0
    a2 = 1.0
    b2 = 10.0
    smol1 = Smoluchowski(a1, b1, func, n=n, D=D)
    smol2 = Smoluchowski(a2, b2, func, n=n, D=D)
    q1 = smol1.expWq
    q2 = smol2.expWq
    k_backwards, k_forwards, T_alpha, N_backwards, N_forwards, R_i_backwards, \
        R_i_forwards, N_ij_backwards, N_ij_forwards \
        = smol2.compute_MMVT_kinetics_quantities()
    
    J2 = q2 / (R_i_forwards + R_i_backwards)
    time = R_i_forwards + q1/J2
    print("time:", time)