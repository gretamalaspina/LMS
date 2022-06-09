import numpy as np
import re
import itertools as it
import numpy as np
import re
import itertools as it
import os.path

def setup_problem(folderpath, PO_file, weights):
    """
    :param folderpath: folder containing the problem files
           PO_file: file containing the coordinate observations
           weights: set of weights to be used for the creation of the underlying network
    :return observations: list containing the observations
                for the i-th observation we have
                observations[0][i] = variables involved
                observations[1][i] = observed value
                observations[2][i] = observation kind (e.g. 'TD', 'AN', 'PO', ...)
                observations[3][i] = standard deviation
            N_pt: number of points
            adj_list: adjacency list of the underlying network
            X0: initial guess, defined from the coordinate observations
    """

    print('STARTING SETUP')

    # read the PO file and set the initial guess
    f = open(folderpath + PO_file)
    init_strings = f.read()
    X0, list_pt, N_pt, st_devs = read_variables(init_strings)


    if os.path.isfile(folderpath + 'project.obs'):
        f = open(folderpath + 'project.obs')
    else:
        f = open(folderpath + 'project_sk1.obs')
    obs_string = f.read()
    observations, N_TD, N_AN, N_BE, N_PL, N_EQ, N_pt, list_pt, adj_list, new_points_list, new_points_dic = parse_observations_file(
        obs_string, list_pt, N_pt, weights,  observations = [[],[],[],[]], count_obs = np.zeros(5))



    # if there are any equality observations, set the initial guess for the involved points to be the same
    list_EQ = [observations[0][i] for i, k in enumerate(observations[2]) if k == 'EQ']
    set_equal_start = 1
    if set_equal_start == 1:
        for eq in list_EQ:
            X0[eq[1]] = X0[eq[0]]

    # add coordinates observations using the values that we got from the PO file (same as initial guess)
    observations, N_PO = coordinates_observations(X0, observations, st_devs)


    print('problem with ', len(X0), 'variables and ', len(observations[0]), 'observations')

    print('SETUP FINISHED')

    return observations, N_pt, adj_list, X0

def read_variables(var_string):
    X = []
    st_devs = []
    list_pt = {}
    count_pt = -1
    for line in var_string.splitlines()[7:]:
        elements = re.split('\s+', line)
        if len(elements)>=2:
            p = elements[1]
            count_pt += 1
            list_pt[p] = count_pt
            X.append(float(elements[2][:-1]))
            X.append(float(elements[3][:-1]))
            if len(elements[4])>1:
                # if standard dev for the coordinate obs is available in the data, we use the availables value
                st_devs.append(float(elements[4]))
                st_devs.append(float(elements[4]))
            else:
                # otherwise, use default value 20
                st_devs.append(20.0)
                st_devs.append(20.0)


    return X, list_pt, count_pt+1, st_devs


def parse_observations_file(observations_string, list_pt, count_pt, weights=[1,3], observations = [[],[],[],[]], count_obs = np.zeros(5)):

    obs_variables,obs_value,obs_kind,obs_variance = observations
    count_TD, count_AN, count_BE, count_PL, count_EQ = count_obs
    adj_list = [[] for i in range(count_pt)]
    new_points_dic = {}
    new_points_list = []

    for line_n, line in enumerate(observations_string.splitlines()):
        elements = re.split('\s+', line)
        if line.find('TD') > 0:
            # point-point distance observation
            pt_eq = [elements[1], elements[2]] #involved points
            val = float(elements[4]) #observed value
            for pt in pt_eq:
                #check if point already in the list of variables and add it if it's not
                if pt not in list_pt:
                    new_points_list.append(pt)
                    new_points_dic[pt] = 'unknown'
            equation, count_pt, list_pt, adj_list = new_equation(2, val, pt_eq, count_pt, list_pt, adj_list)
            obs_value.append(val)
            obs_variance.append(float(elements[5]))
            obs_variables.append(equation[:4])
            obs_kind.append('TD')
            count_TD += 1


        if line.find("AN") > 0:
            # angle observation
            pt_eq = [elements[1], elements[2], elements[3]]
            val = float(elements[5])*np.pi/200 #convert the observed value to radiants
            for pt in pt_eq:
                if pt not in list_pt:
                    new_points_list.append(pt)
                    new_points_dic[pt] = 'unknown'

            equation, count_pt, list_pt, adj_list= new_equation(3, val, pt_eq, count_pt, list_pt, adj_list, weights)
            # equation.append(float(elements[6]))
            obs_variables.append(equation[:6])
            obs_kind.append('AN')

            obs_value.append(val)
            obs_variance.append(float(elements[6])*200)
            # the *100 is not a mistake: the mistake was in generating the data, here we are just taking that into account
            #needs to be changed if we correct the generator
        if line.find("BE") > 0:
            #bearing observation
            pt_eq = [elements[1], elements[2]]
            val = float(elements[4])
            for pt in pt_eq:
                if pt not in list_pt:
                    new_points_list.append(pt)
                    new_points_dic[pt] = 'unknown'
            equation, count_pt, list_pt, adj_list = new_equation(2, val, pt_eq, count_pt, list_pt, adj_list,
                                                                 weights)
            obs_value.append(val)
            obs_variance.append(np.cos(float(elements[5])))
            obs_variables.append(equation[:4])
            obs_kind.append('BE')
            count_BE += 1

        if line.find("CL") > 0:
            pt_eq = [elements[2], elements[1], elements[3]]
            for pt in pt_eq:
                if pt not in list_pt:
                    new_points_list.append(pt)
                    new_points_dic[pt] = 'unknown'

            val = 0
            equation, count_pt, list_pt, adj_list = new_equation(3, val, pt_eq, count_pt, list_pt, adj_list, weights)
            count_PL += 1
            obs_value.append(val)
            obs_variance.append(float(0.015))
            obs_variables.append(equation[:6])
            obs_kind.append('PL')

        if line.find("CH") > 0:
            pt_eq = [elements[2], elements[1], elements[3]]
            for pt in pt_eq:
                if pt not in list_pt:
                    new_points_list.append(pt)
                    new_points_dic[pt] = 'unknown'

            val = 0
            equation, count_pt, list_pt, adj_list = new_equation(3, val, pt_eq, count_pt, list_pt, adj_list, weights)
            count_PL += 1
            obs_value.append(val)
            obs_variance.append(float(0.15))
            obs_variables.append(equation[:6])
            obs_kind.append('PL')

            pt_eq = [elements[1], elements[2]]
            val = float(elements[5])
            equation, count_pt, list_pt,adj_list = new_equation(2, val, pt_eq, count_pt, list_pt, adj_list, weights)
            equation.append(float(elements[6]))

            obs_value.append(val)
            obs_variance.append(0.015)
            obs_variables.append(equation[:4])
            obs_kind.append('TD')
            count_TD += 1

        if line.find('EQ') > 0:
            pt_eq = [elements[0], elements[1]]
            equation, count_pt, list_pt, adj_list = new_equation(2, 0, pt_eq, count_pt, list_pt, adj_list)
            equation.append(0.1)
            count_EQ += 1
            obs_value.append(0)
            obs_variance.append(float(1e-04))
            obs_variables.append([equation[0],equation[2]])
            obs_kind.append('EQ')
            for pt in pt_eq:
                if pt not in list_pt:
                    new_points_list.append(pt)
                    new_points_dic[pt] = 'unknown'

            equation.append(0.1)
            count_EQ += 1
            obs_value.append(0)
            obs_variance.append(float(1e-06))
            obs_variables.append([equation[1], equation[3]])
            obs_kind.append('EQ')
            if pt_eq[0] in new_points_dic:
                new_points_dic[pt_eq[0]] = pt_eq[1]
            if pt_eq[1] in new_points_dic:
                new_points_dic[pt_eq[1]] = pt_eq[0]





    observations = [obs_variables,obs_value,obs_kind,obs_variance]

    return (observations, count_TD, count_AN, count_BE, count_PL, count_EQ, count_pt, list_pt, adj_list, new_points_list, new_points_dic)


def coordinates_observations(X0, observations, st_devs):
    obs_variables,obs_value,obs_kind,obs_variance = observations
    N_PO = np.shape(X0)[0]
    for i in range(N_PO):
        st_dev = st_devs[i]
        obs_variables.append([i])
        obs_value.append(X0[i])
        obs_kind.append('PO')
        obs_variance.append(st_dev)
    observations = [obs_variables,obs_value,obs_kind,obs_variance]

    return observations, N_PO






def new_equation(n, val, pt_eq, count_pt, list_pt, adj_list, weights = [1,1,3]):
    equation = []
    for k in range(n):
        p = pt_eq[k]
        if p in list_pt:
            ind_pt = list_pt[p]
        else:
            print('new point')

            list_pt[p] = count_pt
            ind_pt = count_pt
            adj_list.append([])
            count_pt += 1
        equation.extend([2*ind_pt, 2*ind_pt+1])
    equation.append(val)
    w1, w2, w3 = weights
    if n == 2 and val!=0:
        i = list_pt[pt_eq[0]]
        j = list_pt[pt_eq[1]]
        adj_list[i].append([j,w1])
        adj_list[j].append([i, w1])
    elif n == 2 and val==0:
        i = list_pt[pt_eq[0]]
        j = list_pt[pt_eq[1]]
        adj_list[i].append([j, w3])
        adj_list[j].append([i, w3])

    else:
        i = list_pt[pt_eq[0]]
        j = list_pt[pt_eq[1]]
        k = list_pt[pt_eq[2]]

        adj_list[i].append([j, w2])
        adj_list[j].append([i,w2])

        adj_list[i].append([k,w2])
        adj_list[k].append([i, w2])

        adj_list[k].append([j, w2])
        adj_list[j].append([k, w2])

    return equation, count_pt, list_pt, adj_list


def read_initial_points(string, list_pt):
    X = np.zeros(2 * np.shape(list_pt)[0])
    for line in string.splitlines():
        elements = re.split(';', line)
        point = int(re.split('p', elements[0])[1])
        if point in list_pt:
            j = list_pt[point]
            X[2 * j] = elements[1]
            X[2 * j + 1] = elements[2]
    return X
