"""
test_markov_chain_monte_carlo.py

Unit tests for the MCMC sampling algorithms to estimate error bars in 
milestoning calculations.
"""

import numpy as np
import matplotlib.pyplot as plt

import seekr2.modules.common_base as base
import seekr2.modules.mmvt_cvs.mmvt_cv_base as mmvt_cv_base
import seekr2.modules.common_analyze as common_analyze
import seekr2.modules.mmvt_analyze as mmvt_analyze
import seekr2.modules.markov_chain_monte_carlo as markov_chain_monte_carlo

def make_2d_test_plots(v1_data_points, v2_data_points, X, Y, Z, max_x, max_y):
    """
    Make the plots showing similarities between MCMC procedure and
    analytical distribution.
    """
    fig, axs = plt.subplots(nrows=1, ncols=2)
    #im = axs[0].imshow(Z, vmin=abs(Z).min(), vmax=abs(Z).max(), extent=[0, 1, 0, 1],
    #               cmap=plt.cm.jet)
    #im.set_interpolation('bilinear')
    p = axs[0].pcolor(X, Y, Z, cmap=plt.cm.jet, vmin=abs(Z).min(), 
                      vmax=abs(Z).max())
    axs[0].set_title("Analytical distribution")
    axs[0].set_xlabel("$q_{0,1}$")
    axs[0].set_ylabel("$q_{1,0}$")
        
    axs[1].hist2d(v1_data_points, v2_data_points, bins=50, 
               range=[[0.0, max_x],[0.0, max_y]], 
               cmap=plt.cm.jet)
    axs[1].set_title("Sampled using MCMC procedure")
    axs[1].set_xlabel("$q_{0,1}$")
    axs[1].set_ylabel("$q_{1,0}$")
    plt.show()
    
    return

def compute_2d_distribution_avg_std(bounds_x, bounds_y, Z):
    """
    Given a 2d probability distribution, compute its average and standard
    deviation.
    """
    n = 0.0
    sum_x = 0.0
    sum_y = 0.0
    for i, x in enumerate(bounds_x):
        for j, y in enumerate(bounds_y):
            value = Z[j,i]
            sum_x += x * value
            sum_y += y * value
            n += value
            
    avg_x = sum_x / n
    avg_y = sum_y / n
    
    sum_x2 = 0.0
    sum_y2 = 0.0
    for i, x in enumerate(bounds_x):
        for j, y in enumerate(bounds_y):
            value = Z[j,i]
            sum_x2 += (x-avg_x)**2 * value
            sum_y2 += (y-avg_y)**2 * value
    
    std_x = np.sqrt(sum_x2 / n)
    std_y = np.sqrt(sum_y2 / n)
    
    return avg_x, avg_y, std_x, std_y

def mcmc_algorithm_any_2x2_common(algorithm, Q, N, R, max_x, max_y, 
                                       make_plot=False):
    """
    
    """
    num = 100000
    stride = 4
    skip = 40
    q1_distribution = []
    q2_distribution = []
    for counter in range(num * (stride) + skip):
        Qnew = algorithm(Q, N, R)
        Q = Qnew
        q1_distribution.append(Q[0,1])
        q2_distribution.append(Q[1,0])
    
    n = 100
    bounds_x = np.linspace(0, max_x, n)
    bounds_y = np.linspace(0, max_y, n)
    X,Y = np.meshgrid(bounds_x, bounds_y)
    distribution_2d_histogram = np.zeros((n, n))
    for i, x in enumerate(bounds_x):
        for j, y in enumerate(bounds_y):
            # fill out distribution here
            value = x**N[0,1] * np.exp(-x*R[0,0] - Q[0,1]) \
                * y**N[1,0] * np.exp(-y*R[1,0] - Q[1,0])
            distribution_2d_histogram[j,i] = value
                
    # TODO: check std. dev's
    avg_q1, avg_q2, std_q1, std_q2 = compute_2d_distribution_avg_std(
        bounds_x, bounds_y, distribution_2d_histogram)
    expected_avg_q1 = np.mean(q1_distribution)
    expected_avg_q2 = np.mean(q2_distribution)
    expected_std_q1 = np.std(q1_distribution)
    expected_std_q2 = np.std(q2_distribution)
    assert np.isclose(avg_q1, expected_avg_q1, rtol=0.5, atol=0.01)
    assert np.isclose(avg_q2, expected_avg_q2, rtol=0.5, atol=0.01)
    assert np.isclose(std_q1, expected_std_q1, rtol=0.5, atol=0.01)
    assert np.isclose(std_q2, expected_std_q2, rtol=0.5, atol=0.01)
    
    # make plots
    if make_plot:
        make_2d_test_plots(q1_distribution, q2_distribution, X, Y, 
                           distribution_2d_histogram, max_x, max_y)
    return

def test_mcmc_algorithm_1_2x2_elber():
    """
    Test a simple 2x2 Elber system to ensure that the correct matrices
    are generated and sampled. Optionally generate a plot showing this.
    """
    max_x = 3.0 * 12.0 / 500.0
    max_y = 3.0 * 30.0 / 150.0
    N = np.array([[0, 12], [30, 0]])
    t = np.array([[500],[150]])
    R = np.array([[500], [150]])
    model = base.Model()
    model.num_milestones = 2
    data_sample = common_analyze.Data_sample(model)
    data_sample.N_ij = N
    data_sample.R_i = t
    data_sample.compute_rate_matrix()
    Q = data_sample.Q
    assert Q[0,0] == -12.0 / 500.0
    assert Q[0,1] == 12.0 / 500.0
    assert Q[1,0] == 30.0 / 150.0
    assert Q[1,1] == -30.0 / 150.0
    algorithm = markov_chain_monte_carlo\
            .irreversible_stochastic_matrix_algorithm_sample
    mcmc_algorithm_any_2x2_common(algorithm, Q, N, R, max_x, max_y, 
                                       make_plot=False)
    
    return

def test_mcmc_3x3_mmvt(tmpdir_factory):
    """
    Use the same statistics to generate both MMVT and Elber rate matrices.
    Then use the data to sample matrices and compare the distributions to
    ensure that the MMVT and Elber matrix sampling methods are consistent.
    """
    num = 10000 #100000
    stride = 4
    skip = 90
    n_anchors = 4
    n_milestones = 3
    
    # generate data to feed directly into MMVT_data_sample()
    model = base.Model()
    model.num_milestones = n_milestones
    model.num_anchors = n_anchors
    anchor0 = mmvt_cv_base.MMVT_toy_anchor()
    anchor0.index = 0
    milestone_0_0 = base.Milestone()
    milestone_0_0.index = 0
    milestone_0_0.neighbor_anchor_index = 1
    milestone_0_0.alias_index = 1
    anchor0.milestones = [milestone_0_0]
    
    anchor1 = mmvt_cv_base.MMVT_toy_anchor()
    anchor1.index = 1
    milestone_1_0 = base.Milestone()
    milestone_1_0.index = 0
    milestone_1_0.neighbor_anchor_index = 0
    milestone_1_0.alias_index = 1
    milestone_1_1 = base.Milestone()
    milestone_1_1.index = 1
    milestone_1_1.neighbor_anchor_index = 2
    milestone_1_1.alias_index = 2
    anchor1.milestones = [milestone_1_0, milestone_1_1]
    
    anchor2 = mmvt_cv_base.MMVT_toy_anchor()
    anchor2.index = 2
    milestone_2_0 = base.Milestone()
    milestone_2_0.index = 1
    milestone_2_0.neighbor_anchor_index = 1
    milestone_2_0.alias_index = 1
    milestone_2_1 = base.Milestone()
    milestone_2_1.index = 2
    milestone_2_1.neighbor_anchor_index = 3
    milestone_2_1.alias_index = 2
    anchor2.milestones = [milestone_2_0, milestone_2_1]
    
    anchor3 = mmvt_cv_base.MMVT_toy_anchor()
    anchor3.index = 3
    milestone_3_0 = base.Milestone()
    milestone_3_0.index = 2
    milestone_3_0.neighbor_anchor_index = 2
    milestone_3_0.alias_index = 1
    anchor3.milestones = [milestone_3_0]
    
    model.anchors = [anchor0, anchor1, anchor2, anchor3]
    
    # MMVT stats
    N_alpha_beta = {(0,1):12, (1,0):12,
                    (1,2):12, (2,1):12,
                    (2,3):6,  (3,2):6}
    k_alpha_beta = {(0,1):20.0, (1,0):10.0,
                    (1,2):10.0,  (2,1):(40.0/3.0),
                    (2,3):(20.0/3.0), (3,2):20.0}
    N_i_j_alpha = [{},
                   {(0,1):4, (1,0):4}, 
                   {(1,2):2, (2,1):2},
                   {}]
    R_i_alpha_total = [{0: 1.2},
                       {0: 1.2, 1:1.2},
                       {1: 1.2, 2:0.6},
                       {2: 0.6}]
    T_alpha_total = [1.2,
                     2.4,
                     1.8,
                     0.6]
    
    main_data_sample = mmvt_analyze.MMVT_data_sample(
            model, N_alpha_beta, k_alpha_beta, N_i_j_alpha, 
            R_i_alpha_total, T_alpha_total)
    main_data_sample.calculate_pi_alpha()
    main_data_sample.fill_out_data_quantities()
    main_data_sample.compute_rate_matrix()
    main_data_sample.calculate_thermodynamics()
    main_data_sample.calculate_extra_thermodynamics()
    main_data_sample.calculate_kinetics()
    mmvt_Q = main_data_sample.Q
    
    # Elber stats
    N_i_j = {(0,1): 4, (1,0): 4, (1,2): 2, (2,1): 2}
    R_i = {0: 2.4, 1: 2.4, 2: 1.2}
    elber_N = np.array([[0, 4, 0],
                        [4, 0, 2],
                        [0, 2, 0]])
    elber_R = np.array([[2.4],
                        [2.4],
                        [1.2]])
    
    elber_Q = np.zeros((n_milestones, n_milestones))
    for i in range(n_milestones):
        for j in range(n_milestones):
            key = (i,j)
            if key in N_i_j:
                elber_Q[i,j] = N_i_j[key] / R_i[i]
    
    for i in range(n_milestones):
        elber_Q[i,i] = -np.sum(elber_Q[i,:])
    
    # Make sure the two models make identical matrices
    assert np.isclose(mmvt_Q, elber_Q).all()
    
    # Now compare the distributions of both of them
    
    # MMVT matrix sampler
    data_sample_list, p_i_error, free_energy_profile_err, \
        free_energy_anchors_err, MFPTs_error, k_off_error, k_ons_error \
        = mmvt_analyze.monte_carlo_milestoning_error(
        main_data_sample, num=num, stride=stride, skip=skip, verbose=True)
    
    mmvt_q1_distribution = []
    mmvt_q2_distribution = []
    for data_sample in data_sample_list:
        mmvt_q1_distribution.append(data_sample.Q[0,1])
        mmvt_q2_distribution.append(data_sample.Q[1,0])
    
    # Elber matrix sampler
    elber_q1_distribution = []
    elber_q2_distribution = []
    for counter in range(num * (stride) + skip):
        #if verbose: print("MCMC stepnum: ", counter)
        elber_Qnew = markov_chain_monte_carlo\
            .irreversible_stochastic_matrix_algorithm_sample(
                elber_Q, elber_N, elber_R)
        if counter > skip and counter % stride == 0:
            elber_q1_distribution.append(elber_Q[0,1])
            elber_q2_distribution.append(elber_Q[1,0])
        
        elber_Q = elber_Qnew
    
    assert np.isclose(np.mean(elber_q1_distribution), 
                      np.mean(mmvt_q1_distribution), rtol=0.5, atol=0.01)
    assert np.isclose(np.std(elber_q1_distribution), 
                      np.std(mmvt_q1_distribution), rtol=0.5, atol=0.01)
    assert np.isclose(np.mean(elber_q2_distribution), 
                      np.mean(mmvt_q2_distribution), rtol=0.5, atol=0.01)
    assert np.isclose(np.std(elber_q2_distribution), 
                      np.std(mmvt_q2_distribution), rtol=0.5, atol=0.01)
    
    return
    
    