"""
=================================
            SYNESNET
=================================
This module is design to load data.
"""

from copy import deepcopy
import numpy as np
import os.path
from scipy import io

path = '/Users/juliana.gonzalez/ownCloud/graph_analysis/'
node_file = os.path.join(path, 'BN_Atlas_246_LUT_reoriented.txt')
results_file = os.path.join(path, 'resultsROI_Subject001_Condition001.mat')


def load_net_metrics(net_path, metric, corr_type, num_sub, idx_select):
    xnet = np.array([
        io.loadmat(os.path.join(net_path, f'net_metrics_Subject{str(sub).zfill(3)}{corr_type}'))[metric][0]
        for sub in range(1, num_sub + 1)])
    x_ = np.zeros(np.shape(xnet))
    x_[:, idx_select] = xnet[:, idx_select]
    return x_


def load_fc():
    # Load all connectivity matrices
    num_sub = len(os.listdir(os.path.join(path, 'symetrical_corr_mat')))
    fc_all = []
    for sub in np.arange(num_sub) + 1:
        # load correlation
        sub_file = 'CorrMatrix_Subject{0}.mat'.format(str(sub).zfill(3))
        fc_file = path + 'symetrical_corr_mat/' + sub_file
        fc = io.loadmat(fc_file)

        # 3. corr[corr>=0] = [0, 1]  # winning matrices !!!
        fc_thr = deepcopy(fc['CorrMatrix'])
        fc_thr[fc_thr <= 0] = 0
        np.fill_diagonal(fc_thr, 0)

        # Append matrices
        fc_all.append(fc_thr)

    fc_syn = np.array(fc_all[:17])
    fc_ctr = np.array(fc_all[17:])

    return fc_syn, fc_ctr


def load_node_names():
    # nodes names
    n_name, n_name_full = [], []
    with open(node_file, "r") as filestream:
        for line in filestream:
            name, name_full = line.split(",")[:2]
            n_name.append(name)
            n_name_full.append(name_full.strip())

    return n_name, n_name_full


def load_xyz():
    xyz = io.loadmat(results_file)['xyz'][0][:-1]
    xyz = np.array([x[0] for x in xyz])
    return xyz
