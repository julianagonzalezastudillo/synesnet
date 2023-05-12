import os.path
import scipy.io as sio
import numpy as np
from scipy import io
import networkx as nx
from net_metrics import local_laterality
import bct
from net.core import coreness
from copy import deepcopy


class bcolors:
    WARNING = '\033[93m'
    ENDC = '\033[0m'


path = '/Users/juliana.gonzalez/ownCloud/graph_analysis/'
net_path = '/Users/juliana.gonzalez/ownCloud/graph_analysis/net_metrics/'
if not os.path.exists(net_path):
    os.makedirs(net_path)

n_name = []
with open(path + "BN_Atlas_246_LUT_reoriented.txt", "r") as filestream:
    for line in filestream:
        n_name.append(line.split(",")[0])

for sub in np.arange(1, 35):
    sub_file = 'CorrMatrix_Subject{0}.mat'.format(str(sub).zfill(3))
    print(sub_file)

    # load correlation
    fc_file = path + 'symetrical_corr_mat/' + sub_file
    fc = io.loadmat(fc_file)

    # connectivity matrix
    # work with 3 possibilities of connectivity matrices
    # 1. corr = [-1, 1]
    # Xfc = deepcopy(fc['CorrMatrix'])
    # net_file_Xfc = net_path + '_'.join(('net_metrics', 'Subject{0}.mat'.format(str(sub).zfill(3))))

    # 2. abs(corr) = [0, 1]
    # Xfc_abs = abs(deepcopy(fc['CorrMatrix']))
    # net_file_abs = net_path + '_'.join(('net_metrics', 'Subject{0}'.format(str(sub).zfill(3)), '_abs.mat'))

    # 3. corr[corr>=0] = [0, 1]  # wining matrices !!!
    Xfc_thr = deepcopy(fc['CorrMatrix'])
    Xfc_thr[Xfc_thr <= 0] = 0
    net_file_thr = net_path + '_'.join(('net_metrics', 'Subject{0}'.format(str(sub).zfill(3)), 'thr.mat'))

    # for X, net_file in zip([Xfc, Xfc_abs, Xfc_thr], [net_file_Xfc, net_file_abs, net_file_thr]):
    for X, net_file in zip([Xfc_thr], [net_file_thr]):
        np.fill_diagonal(X, 0)

        # network metrics
        G = nx.from_numpy_matrix(X)  # to nx format
        strength = np.array([v for k, v in G.degree(weight = 'weight')])  # strength
        lat = local_laterality(X, n_name)  # laterality
        geff = bct.efficiency_wei(X)
        leff = bct.efficiency_wei(X, local=True)
        C, isCore = coreness(X)  # coreness
        # sw = nx.sigma(X)

        # save
        Xnet = {'strength': strength,
                'laterality': lat,
                'global_efficiency': geff,
                'local_efficiency': leff,
                'coreness': C,
                'isCore': isCore
                }

        sio.savemat(net_file, Xnet)



