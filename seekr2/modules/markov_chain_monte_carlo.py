"""
markov_chain_monte_carlo.py

Algorithms used to sample distributions of rate matrices for Elber and MMVT
milestoning calculations.

Citation: Noe, F. "Probability Distributions of molecular observables
computed from Markov models." J. Chem. Phys. 2008, 128, No. 244103.
Distribution is:  p(Q|N) = p(Q)p(N|Q)/p(N) = 
    p(Q) PI(q_ij**N_ij * exp(-q_ij * Ri))
"""
import random
from copy import deepcopy

import numpy as np
from scipy.stats import expon
    
def irreversible_stochastic_matrix_algorithm_sample(Q, N, R):
    """
    Sample an additional matrix Q from the distribution using a Monte-Carlo
    method. Use a modified version of Algorithm 1 from Noe 2008.
    """
    # step 1 - Initialize the matrix (already done)
    Qnew = deepcopy(Q)
    # step 2.1 - Generate uniform random variables: i,j members of {1,...,m},
    #  delta, r between 0 and 1
    m = Q.shape[0]
    i = random.choice(range(m))
    j = random.choice(range(m))
    if Q[i,j] == 0.0: return Qnew
    if Q[i,i] == 0.0: return Qnew
    if N[i,j] == 0: return Qnew
    if R[i,0] == 0: return Qnew
    # delta must be sampled such that when it's subtracted from q_ij, the 
    #  result doesn't go below zero
    delta = Q[i,j] * (1.0 - expon.rvs())
    random_uniform = random.random()  
    
    # step 2.2 - copy Q to a new matrix - already done
    
    # step 2.3 - Nonreversible element shift
    # prior probability: P(Q)
    log_prior_probability_old = -Qnew[i,j]
    log_prior_probability_new = -Qnew[i,j] + delta
    log_p_Q_old = N[i,j] * np.log(Qnew[i,j])  - Qnew[i,j] * R[i,0]
    log_p_Q_new = N[i,j] * np.log(Qnew[i,j] - delta) \
        - (Qnew[i,j] - delta) * R[i,0]
    p_acc =  log_p_Q_new - log_p_Q_old + log_prior_probability_new \
        - log_prior_probability_old
        
    if np.log(random_uniform) <= p_acc: 
        #log(r) can be directly compared to 
        # log-likelihood acceptance, p_acc
        Qnew[i,i] = Qnew[i,i] + delta
        Qnew[i,j] = Qnew[i,j] - delta
        
    return Qnew