"""
test_analyze.py

Testing analyze.py
"""

import os
from collections import defaultdict

import numpy as np

import seekr2.modules.common_analyze as common_analyze
import seekr2.modules.mmvt_analyze as mmvt_analyze
import seekr2.analyze as analyze
import seekr2.modules.common_base as base
import seekr2.modules.mmvt_base as mmvt_base
import seekr2.tests.smoluchowski_system as smoluchowski

this_dir = os.path.dirname(os.path.realpath(__file__))

test_output_filename = os.path.join(this_dir, 
                                    "data/test_analyze_outputfile.txt")

def test_read_output_file():
    N_i_j_alpha, R_i_alpha_list, R_i_alpha_average, \
        R_i_alpha_std_dev, R_i_alpha_total, N_alpha_beta, \
        T_alpha_list, T_alpha_average, T_alpha_std_dev, \
        T_alpha_total, existing_lines \
            = mmvt_analyze.openmm_read_output_file_list(
                [test_output_filename])
    
    N_i_j_alpha_dict1 = N_i_j_alpha
    R_i_alpha_dict1 = R_i_alpha_total
    N_alpha_beta_dict1 = N_alpha_beta
    T_alpha1 = T_alpha_total
    #N_i_j_alpha_dict1, R_i_alpha_dict1, N_alpha_beta_dict1, T_alpha1 = \
    #    analyze.openmm_read_output_file_list([test_output_filename])
    
    N_i_j_alpha_dict2 = {(1, 2): 52, (2, 1): 52}
    R_i_alpha_dict2 = {1: 1658.696, 2: 198.912}
    N_alpha_beta_dict2 = {1: 2423, 2: 98}
    T_alpha2 = 1954.760
    
    for key in N_i_j_alpha_dict1:
        assert key in N_i_j_alpha_dict2
        assert np.isclose(N_i_j_alpha_dict1[key], N_i_j_alpha_dict2[key])
        
    for key in R_i_alpha_dict1:
        assert key in R_i_alpha_dict2
        assert np.isclose(R_i_alpha_dict1[key], R_i_alpha_dict2[key])
        
    for key in N_alpha_beta_dict1:
        assert key in N_alpha_beta_dict2
        assert np.isclose(N_alpha_beta_dict1[key], N_alpha_beta_dict2[key])
    
    assert np.isclose(T_alpha1, T_alpha2)
    
    N_i_j_alpha, R_i_alpha_list, R_i_alpha_average, \
        R_i_alpha_std_dev, R_i_alpha_total, N_alpha_beta, \
        T_alpha_list, T_alpha_average, T_alpha_std_dev, \
        T_alpha_total, existing_lines \
            = mmvt_analyze.openmm_read_output_file_list([test_output_filename, 
                                                    test_output_filename], 
                                                   skip_restart_check=True)
    
    N_i_j_alpha_dict1 = N_i_j_alpha
    R_i_alpha_dict1 = R_i_alpha_total
    N_alpha_beta_dict1 = N_alpha_beta
    T_alpha1 = T_alpha_total
    #N_i_j_alpha_dict1, R_i_alpha_dict1, N_alpha_beta_dict1, T_alpha = \
    #    analyze.openmm_read_output_file_list([test_output_filename, 
    #                                   test_output_filename])
        
    for key in N_i_j_alpha_dict1:
        assert key in N_i_j_alpha_dict2
        assert np.isclose(N_i_j_alpha_dict1[key], 2*N_i_j_alpha_dict2[key], 
                          rtol=0.01)
        
    for key in N_alpha_beta_dict1:
        assert key in N_alpha_beta_dict2
        assert np.isclose(N_alpha_beta_dict1[key], 2*N_alpha_beta_dict2[key], 
                          rtol=0.01)
    
    return

def test_minor2d():
    A = np.array([[1,2,3],[4,5,6],[7,8,9]])
    B = np.array([[1,3],[7,9]])
    C = np.array([[1,2],[4,5]])
    D = np.array([[2,8],[3,9]])
    assert common_analyze.minor2d(A, 1, 1).all() == B.all()
    assert common_analyze.minor2d(A, 2, 2).all() == C.all()
    assert common_analyze.minor2d(A, 1, 0).all() == D.all()
    return
    
def test_minor1d():
    A = np.array([1,2,3])
    B = np.array([1,3])
    C = np.array([2,3])
    D = np.array([1,2])
    assert common_analyze.minor1d(A, 1).all() == B.all()
    assert common_analyze.minor1d(A, 0).all() == C.all()
    assert common_analyze.minor1d(A, 2).all() == D.all()
    return

def make_fake_output_file_osc(anchor, tmp_path, timestep=1.0):
    num_steps = 50
    
    mmvt_output_filename = os.path.join(
                    tmp_path, anchor.name, "prod", 
                    "%s%d.%s" % (mmvt_base.OPENMMVT_BASENAME, 1, 
                                         mmvt_base.OPENMMVT_EXTENSION))
    with open(mmvt_output_filename, "w") as f:
        if anchor.index == 0:
            for i in range(num_steps+1):
                line = "%d,%d,%f\n" % (1, i, i*timestep)
                f.write(line)
            
        else:
            for i in range(num_steps+1):
                if (i % 2) == 0:
                    line = "%d,%d,%f\n" % (2, i, i*timestep)
                    f.write(line)
                else:
                    line = "%d,%d,%f\n" % (1, i, i*timestep)
                    f.write(line)
    return

def make_fake_output_file2(anchor, tmp_path, ups=1, downs=9, timestep=1.0):
    num_steps = 50
    total = ups + downs
    
    mmvt_output_filename = os.path.join(
                    tmp_path, anchor.name, "prod", 
                    "%s%d.%s" % (mmvt_base.OPENMMVT_BASENAME, 1, 
                                         mmvt_base.OPENMMVT_EXTENSION))
    with open(mmvt_output_filename, "w") as f:
        if anchor.index == 0:
            for i in range(num_steps+1):
                line = "%d,%d,%f\n" % (1, i, i*timestep)
                f.write(line)
                
        else:
            for i in range(num_steps+1):
                if (i % total) < ups:
                    line = "%d,%d,%f\n" % (2, i, i*timestep)
                    f.write(line)
                else:
                    line = "%d,%d,%f\n" % (1, i, i*timestep)
                    f.write(line)
    return



"""
def make_smol_calculation(tmp_path, func=None):
    num_anchors = 10
    D = 0.01
    interval = 1.0
    n = 101
    
    intervals = []
    for i in range(num_anchors):
        intervals.append(interval)
    
    if func is None:
        func = smoluchowski.expW_constant
    
    q_s = np.zeros(num_anchors)
    mymodel = smoluchowski.make_smol_model(tmp_path, num_anchors, intervals)
    my_analysis = analyze.Analysis(mymodel)
    elberN_ij = defaultdict(float)
    elberR_i = defaultdict(float)
    smols = []
    for i, anchor in enumerate(mymodel.anchors[:-1]):
        a = interval*i
        b = interval*(i+1)
        smol = smoluchowski.Smoluchowski(a, b, func, n=n, D=D)
        q_s[i] = smol.expWq
        if i == 0:
            smol.reflect_lower = True
        k_backwards, k_forwards, T_alpha, N_backwards, N_forwards, \
            R_i_backwards, R_i_forwards, N_ij_backwards, N_ij_forwards \
            = smol.compute_MMVT_kinetics_quantities()
        
        N_i_j_alpha_dict = defaultdict(int)
        R_i_alpha_dict = defaultdict(float)
        N_alpha_beta_dict = defaultdict(int)
        new_time_factor = (R_i_forwards + R_i_backwards) / T_alpha
        new_T_alpha = new_time_factor * T_alpha
        if i == 0:
            N_alpha_beta_dict[1] = new_time_factor
            R_i_alpha_dict[1] = new_T_alpha
        else:
            N_i_j_alpha_dict[(1, 2)] = N_ij_forwards
            N_i_j_alpha_dict[(2, 1)] = N_ij_backwards
            R_i_alpha_dict[1] = R_i_forwards
            R_i_alpha_dict[2] = R_i_backwards
            N_alpha_beta_dict[1] = N_backwards * new_time_factor
            N_alpha_beta_dict[2] = N_forwards * new_time_factor
        
        anchor_stats = mmvt_analyze.MMVT_anchor_statistics(alpha=i)
        anchor_stats.N_i_j_alpha = N_i_j_alpha_dict
        anchor_stats.R_i_alpha_total = R_i_alpha_dict
        anchor_stats.R_i_alpha_std_dev = R_i_alpha_dict
        anchor_stats.R_i_alpha_list = {}
        for key in anchor_stats.R_i_alpha_total:
            anchor_stats.R_i_alpha_list[key] = []
        anchor_stats.N_alpha_beta = N_alpha_beta_dict
        anchor_stats.T_alpha_total = new_T_alpha
        anchor_stats.T_alpha_std_dev = new_T_alpha
        for key in N_alpha_beta_dict:
            anchor_stats.k_alpha_beta[key] = N_alpha_beta_dict[key] \
                / new_T_alpha

        #    N_i_j_alpha_dict, R_i_alpha_dict, N_alpha_beta_dict, new_T_alpha, 
        #    alpha=i)
        # FIll out values here...
        my_analysis.anchor_stats_list.append(anchor_stats)
        smols.append(smol)
        
    for i, anchor in enumerate(mymodel.anchors[:-1]):
        smol1 = smols[i]
        if i == 0:
            smol2 = smols[i+1]
            elberN_ij[(0,1)] = 1.0
            # need to make sure that u and exp(-beta*W) match up
            #  on the edge.
            smol1_edge_value = smol1.expWfunc(smol1.b, q=smol1.expWq)
            elberR_i[0] = (smol2.u_q_forward + (1.0/smol1_edge_value)) / (smol2.J_forward)
        elif i == mymodel.num_milestones-1:
            elberN_ij[(mymodel.num_milestones-1,mymodel.num_milestones-2)] = 1.0
            elberR_i[mymodel.num_milestones-1] = (smol1.u_q_backward) / (smol1.J_backward)
        else:
            smol2 = smols[i+1]
            elberN_ij[(i,i+1)] = smol2.J_forward / (smol2.J_forward + smol1.J_backward)
            elberN_ij[(i,i-1)] = smol1.J_backward / (smol2.J_forward + smol1.J_backward)
            elberR_i[i] = (smol2.u_q_forward + smol1.u_q_backward) / (smol2.J_forward + smol1.J_backward)
    
    my_analysis.mmvt_check_anchor_stats()
    
    #my_analyze._calculate_equilibrium_probability()
    #my_analyze._calculate_overall_statistics()
    #my_analysis.extract_data()
    my_analysis.fill_out_data_samples()
    my_analysis.main_data_sample.pi_alpha = np.zeros(mymodel.num_anchors)
    for i, anchor in enumerate(mymodel.anchors[:-1]):
        my_analysis.main_data_sample.pi_alpha[i] = q_s[i] / np.sum(q_s)
    my_analysis.fill_out_data_samples()
    my_analysis.process_data_samples()
    my_analysis.main_data_sample.Q = np.zeros((mymodel.num_milestones, 
                       mymodel.num_milestones), dtype=np.longdouble)
    elberQ = np.zeros((mymodel.num_milestones, 
                       mymodel.num_milestones), dtype=np.longdouble)
    for i in range(mymodel.num_milestones):
        for j in range(mymodel.num_milestones):
            if my_analysis.main_data_sample.R_i[i] == 0.0:
                my_analysis.main_data_sample.Q[i,j] = 0.0
            else:
                my_analysis.main_data_sample.Q[i,j] \
                    = my_analysis.main_data_sample.N_ij[i,j] \
                    / my_analysis.main_data_sample.R_i[i]
                if elberR_i[i] > 0.0:
                    elberQ[i,j] = elberN_ij[i,j] / elberR_i[i]
                
    for i in range(mymodel.num_milestones):
        my_analysis.main_data_sample.Q[i][i] = \
            -np.sum(my_analysis.main_data_sample.Q[i])
        elberQ[i][i] = -np.sum(elberQ[i])
    
    #my_analyze._rate_mat_to_prob_mat()
    #print("my_analyze.Q:", my_analyze.Q)
    #print("elberQ:", elberQ)
    #print("my_analyze.K:", my_analyze.K)
    #my_analyze.calculate_kinetics()
    my_analysis.main_data_sample.calculate_kinetics()
    mmvt_time = my_analysis.main_data_sample.MFPTs[(0,"bulk")]
    #print("mmvt_time:", mmvt_time)
    my_analysis.main_data_sample.Q = elberQ
    my_analysis.main_data_sample.calculate_kinetics()
    elber_time = my_analysis.main_data_sample.MFPTs[(0,"bulk")]
    #print("elber_time:", elber_time)
    
    a1 = 0.0
    b1 = interval
    a2 = interval
    b2 = interval*num_anchors
    smol1 = smoluchowski.Smoluchowski(a1, b1, func, n=n, D=D)
    smol2 = smoluchowski.Smoluchowski(a2, b2, func, n=n, D=D)
    q1 = smol1.expWq
    q2 = smol2.expWq
    k_backwards, k_forwards, T_alpha, N_backwards, N_forwards, R_i_backwards, \
        R_i_forwards, N_ij_backwards, N_ij_forwards \
        = smol2.compute_MMVT_kinetics_quantities()
    
    J2 = q2 / (R_i_forwards + R_i_backwards)
    correct_time = R_i_forwards + q1/J2
    #print("correct_time:", correct_time)
    print("Time predicted by Elber:", elber_time, "Time predicted by MMVT:", 
          mmvt_time, "Exact time:", correct_time)
    
    ""
    x_s = np.arange(0.0, num_anchors, interval)
    func_vals1 = np.zeros(num_anchors)
    func_vals2 = np.zeros(num_anchors)
    print("q_s:", q_s)
    for i, x in enumerate(x_s):
        print("i:", i, "my_analyze.pi_alpha[i]:", my_analyze.pi_alpha[i], "q_s[i]:", q_s[i] / np.sum(q_s))
        func_vals1[i] = my_analyze.pi_alpha[i]
        func_vals2[i] = q_s[i] / np.sum(q_s)
    
    plt.plot(x_s, func_vals1, "g", x_s, func_vals2, "r")
    plt.show()
    ""
    return mmvt_time, elber_time, correct_time
    

def test_smoluchowski_solution_flat_1(tmp_path):
    print("Constant PMF:")
    mmvt_time, elber_time, true_time = make_smol_calculation(tmp_path)
    assert np.isclose(mmvt_time, true_time, rtol=0.001)
    assert np.isclose(elber_time, true_time, rtol=0.001)
    
    print("linear PMF:")
    func = smoluchowski.expW_linear
    mmvt_time, elber_time, true_time = make_smol_calculation(tmp_path, func)
    assert np.isclose(mmvt_time, true_time, rtol=0.001)
    assert np.isclose(elber_time, true_time, rtol=0.001)
    
    print("quadratic PMF:")
    func = smoluchowski.expW_quadratic
    mmvt_time, elber_time, true_time = make_smol_calculation(tmp_path, func)
    assert np.isclose(mmvt_time, true_time, rtol=0.001)
    assert np.isclose(elber_time, true_time, rtol=0.001)
"""