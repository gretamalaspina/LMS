from sklearn import cluster
import numpy as np
from scipy.sparse import linalg as sla
import metis
import pypardiso
import random

from scipy import sparse as sp
from scipy.sparse.linalg import spsolve


def clustering(adj_list, K=2, splitting = 'metis'):
    print('BEGIN SPLITTING')
    if splitting == 'spectral':
        cl_labels = cluster.spectral_clustering(Adj, n_clusters=K, n_init = 5, eigen_solver='arpack', eigen_tol=0.005)

    elif splitting == 'metis':
        graph = metis.adjlist_to_metis(adj_list)
        weight, cl_labels = metis.part_graph(graph, nparts=K, ubvec=[1.1], recursive=False, objtype = 'cut', ncuts = 10)

    else:
        raise Exception('Unrecognized splitting method.')

    #cl_labels is a vector of length N where N is the number of points
    #if j = cl_labels[i] then point i belogns to the j-th subset of points

    cl_variables = [[] for k in range(K)]
    #cl_variables[j] will hold the list of variabels in the j-th subset

    cl_var_dic = [{} for k in range(K)]
    #cl_var_dic[j] handles the variables numbering restricted to the j-th subset
    #if s = cl_var_dic[j][i] then the i-th global variable is the s-th variable in subset j


    sizes = np.zeros(K)
    for i, l in enumerate(cl_labels):
        cl_variables[l].extend([2 * i, 2 * i + 1])
        cl_var_dic[l][2*i] = int(sizes[l])
        cl_var_dic[l][2*i+1] = int(sizes[l]+1)
        sizes[l] += 2


    print('clusters sizes:', [len(cl_variables[k]) for k in range(K)])
    print('SPLITTING FINISHED')
    return cl_variables, cl_labels, cl_var_dic


def split_LM_step(N_pt, K, R, cl_J, compl_J, mu, cl_var, cl_eq, compl_eq, RHS_mat='B'):
    #
    # approximatly solves the linear system (J.TJ + muI)d = -J.TR exploiting the splitting into suproblem
    # :param N_pt: total number of points
    # :param K: number of suproblems we split into
    # :param R: residual vector
    # :param cl_J[k]: Jacobian of the observations in subproblem k
    # :param compl_J[k]: Jacobian wrt the variables in cl_var[k] of the observations in compl_eq
    # :param mu: current value of the damping parameter
    # :param cl_var[k]: variables in subproblem k
    # :param cl_eq[k]: observations in subproblem k
    # :param compl_eq: observations between different subproblems
    # :param RHS_mat: matrix used for the RHS correction
    #
    # :return d direction, approximate solution of (J.TJ + muI)d = -J.TR
    #         cl_G[k] gradient wrt the variables in subproblem k
    #
    #
    # if RHS_mat = 'none' we compute d by solving the linear systems corresponding to the subproblems independently
    #     that is, we solve Hd = -g
    # if RHS_mat = 'B' we apply a RHS correction before solving the independent systems
    #     that is, we solve Hd = -g + beta*B
    #
    # cl_G[k] is the subvector of g = J.TR that only involves derivatives wrt variables in subproblem k
    # cl_H[k] is the submatrix of (J.TJ + muI) that only involves derivatives wrt variables in subproblem k
    # B[k][s] is the submatrix of (J.TJ + muI) that involves mixed derivatives in subproblems k and s
    #
    cl_H = [[] for k in range(K)]
    cl_H_fac = [[] for k in range(K)]
    compl_eq_ind = [eq[0] for eq in compl_eq]
    cl_G = [[] for k in range(K)]

    d = np.zeros(2 * N_pt)
    for k in range(K):
        cl_H[k] = cl_J[k].T @ cl_J[k] + compl_J[k].T @ compl_J[k] + mu * sp.eye(np.shape(cl_J[k])[1])
        cl_G[k] = cl_J[k].T @ R[cl_eq[k]] + compl_J[k].T @ R[compl_eq_ind]

    if RHS_mat == 'none':
        print('   no RHS correction')

        # solve the K linear systems independently
        for k in range(K):
            d[cl_var[k]] = pypardiso.spsolve(cl_H[k],-cl_G[k])
        return d, cl_G

    elif RHS_mat == 'B':
        print('   using B as RHS correction, with 1 coefficient')


        # we need to solve several linear systems with the same coefficient matrix cl_H[k] so we store the factorization
        for k in range(K):
            cl_H_fac[k] = pypardiso.factorized(cl_H[k])

        B = [[[] for k in range(K)] for k in range(K)]
        for i in range(K):
            for j in range(i + 1, K):
                B[i][j] = compl_J[i].T @ compl_J[j]
                B[j][i] = B[i][j].T


        # compute the correction coefficient beta
        u = [np.zeros(len(cl_var[k])) for k in range(K)]
        v = [np.zeros(len(cl_var[k])) for k in range(K)]
        w = [np.zeros(len(cl_var[k])) for k in range(K)]
        v_aux = [np.zeros(len(cl_var[k])) for k in range(K)]

        w_aux = [np.zeros(len(cl_var[k])) for k in range(K)]
        for k in range(K):
            for j in range(K):
                if j != k:
                    u[k] += B[k][j] @ cl_G[j]
            w_aux[k] = cl_H_fac[k](u[k])
            v_aux[k] = cl_H_fac[k](cl_G[k])

        for k in range(K):
            for j in range(K):
                if j != k:
                    v[k] += B[k][j] @ v_aux[j]
                    w[k] += B[k][j] @ w_aux[j]
        uw = [[] for k in range(K)]
        uw_p = 0
        uwv_p = 0
        for k in range(K):
            uw[k] = u[k] + w[k]
            uw_p += uw[k] @ uw[k]
            uwv_p += uw[k] @ v[k]

        beta = uwv_p / uw_p #correction coefficient


        # solve the K linear systems, with corrected RHS
        for k in range(K):
            d[cl_var[k]] = cl_H_fac[k](beta * u[k] - cl_G[k])

    else:
        raise Exception('Unrecognized RHS_mat choice. Available options: none, B ')


    return d, cl_G


def split_indices(observations, K, cl_var_dic, cl_eq, eq_compl, ind_anchors = []):
    """
    Does the same thing as evaluations.find_indices but separately for
         cl_eq (that is, separately for each subproblem )
         eq_compl (that is, for the observations between different subproblems)
    """
    print('   splitting indices')
    cl_indices = [[] for k in range(K)]
    cl_ptr = [[] for k in range(K)]
    obs_variables = observations[0]
    compl_indices = [[] for k in range(K)]
    compl_ptr = [[0] for k in range(K)]

    #cl_eq
    for k in range(K):
        cl_indices[k] = [cl_var_dic[k][v] for j in cl_eq[k] for v in obs_variables[j]]
        cl_ptr[k] = [0]
        for j in cl_eq[k]:
            eq = obs_variables[j]
            cl_ptr[k].append(cl_ptr[k][-1] + len(eq))

    #eq_compl
    for obs in eq_compl:
        if len(obs) == 3:
            j,k0,k1 = obs
            count_k0 = 0
            count_k1 = 0
            if observations[-2][j] == 'TD':
                if obs_variables[j][0] not in ind_anchors:
                    compl_indices[k0].append(cl_var_dic[k0][obs_variables[j][0]])
                    count_k0+=1
                if obs_variables[j][1] not in ind_anchors:
                    compl_indices[k0].append(cl_var_dic[k0][obs_variables[j][1]])
                    count_k0 += 1
                if obs_variables[j][2] not in ind_anchors:
                    compl_indices[k1].append(cl_var_dic[k1][obs_variables[j][2]])
                    count_k1 += 1
                if obs_variables[j][3] not in ind_anchors:
                    compl_indices[k1].append(cl_var_dic[k1][obs_variables[j][3]])
                    count_k1 += 1
                compl_ptr[k0].append(compl_ptr[k0][-1] + count_k0)
                compl_ptr[k1].append(compl_ptr[k1][-1] + count_k1)

            elif observations[-2][j] == 'BE':
                if obs_variables[j][0] not in ind_anchors:
                    compl_indices[k0].append(cl_var_dic[k0][obs_variables[j][0]])
                    count_k0 += 1
                if obs_variables[j][1] not in ind_anchors:
                    compl_indices[k0].append(cl_var_dic[k0][obs_variables[j][1]])
                    count_k0 += 1

                if obs_variables[j][2] not in ind_anchors:
                    compl_indices[k1].append(cl_var_dic[k1][obs_variables[j][2]])
                    count_k1 += 1
                if obs_variables[j][3] not in ind_anchors:
                    compl_indices[k1].append(cl_var_dic[k1][obs_variables[j][3]])
                    count_k1 += 1
                compl_ptr[k0].append(compl_ptr[k0][-1] + count_k0)
                compl_ptr[k1].append(compl_ptr[k1][-1] + count_k1)

            elif observations[-2][j] == 'EQ':
                if obs_variables[j][0] not in ind_anchors:
                    compl_indices[k0].append(cl_var_dic[k0][obs_variables[j][0]])
                    compl_ptr[k0].append(compl_ptr[k0][-1] + 1)
                if obs_variables[j][1] not in ind_anchors:
                    compl_indices[k1].append(cl_var_dic[k1][obs_variables[j][1]])
                    compl_ptr[k1].append(compl_ptr[k1][-1] + 1)

            else:
                print(observations[-2][j])
                raise Exception('splitting line 204: len(obs)==3')
            for k in range(K):
                if k!=k0 and k!=k1:
                    compl_ptr[k].append(compl_ptr[k][-1])
        else:
            j,k0,k1,k2 = obs
            count_k0 = 0
            count_k1 = 0
            count_k2 = 0
            if obs_variables[j][0] not in ind_anchors:
                compl_indices[k0].append(cl_var_dic[k0][obs_variables[j][0]])
                count_k0+=1
            if obs_variables[j][1] not in ind_anchors:
                compl_indices[k0].append(cl_var_dic[k0][obs_variables[j][1]])
                count_k0+=1
            if obs_variables[j][2] not in ind_anchors:
                compl_indices[k1].append(cl_var_dic[k1][obs_variables[j][2]])
                count_k1+=1
            if obs_variables[j][3] not in ind_anchors:
                compl_indices[k1].append(cl_var_dic[k1][obs_variables[j][3]])
                count_k1+=1
            if obs_variables[j][4] not in ind_anchors:
                compl_indices[k2].append(cl_var_dic[k2][obs_variables[j][4]])
                count_k2+=1
            if obs_variables[j][5] not in ind_anchors:
                compl_indices[k2].append(cl_var_dic[k2][obs_variables[j][5]])
                count_k2+=1
            if k0 == k1:
                compl_ptr[k0].append(compl_ptr[k0][-1] + count_k0+count_k1)
                compl_ptr[k2].append(compl_ptr[k2][-1] + count_k2)
            if k0 == k2:
                compl_ptr[k0].append(compl_ptr[k0][-1] + count_k0+count_k2)
                compl_ptr[k1].append(compl_ptr[k1][-1] + count_k1)
            if k1 == k2:
                compl_ptr[k0].append(compl_ptr[k0][-1] + count_k0)
                compl_ptr[k2].append(compl_ptr[k2][-1] + count_k2+count_k1)
            if k0!= k1 and k0!= k2 and k1!=k2:
                compl_ptr[k0].append(compl_ptr[k0][-1] + count_k0)
                compl_ptr[k2].append(compl_ptr[k2][-1] + count_k2)
                compl_ptr[k1].append(compl_ptr[k1][-1] + count_k1)
            for k in range(K):
                if k!=k0 and k!=k1 and k!=k2:
                    compl_ptr[k].append(compl_ptr[k][-1])


    print('   finished splitting indices')



    return cl_indices, cl_ptr, compl_indices, compl_ptr



