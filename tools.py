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
import scipy.io as sio

path = '/Users/juliana.gonzalez/ownCloud/graph_analysis/'
node_file = os.path.join(path, 'BN_Atlas_246_LUT_reoriented.txt')
results_file = os.path.join(path, 'resultsROI_Subject001_Condition001.mat')

# number of subjects
num_sub = len(os.listdir(os.path.join(path, 'symetrical_corr_mat')))


def load_net_metrics(net_path, metric, corr_type='_thr', num_sub=num_sub, idx_select=slice(None)):
    xnet = np.array([
        io.loadmat(os.path.join(net_path, f'net_metrics_Subject{str(sub).zfill(3)}{corr_type}'))[metric][0]
        for sub in range(1, num_sub + 1)])

    # Create an empty matrix of subjects and nodes and replace the significant nodes with its values
    x = np.zeros(np.shape(xnet))
    x[:, idx_select] = xnet[:, idx_select]

    # Put inf and nan values to zero
    x[np.isnan(x) | np.isinf(x)] = 0

    # Split synesthetic and control subjects
    x_syn = x[:17, :]
    x_ctr = x[17:, :]

    return x_syn, x_ctr


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

    # Split synesthetic and control subjects
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


def save_mat_file(X, xyz, rgb_values, n_name, X_name, plot_path):
    # select index for left and right hemisphere
    lh_ind = [index for index, element in enumerate(np.array(n_name)) if element.endswith('_L')]
    rh_ind = [index for index, element in enumerate(np.array(n_name)) if element.endswith('_R')]

    for ind, side in zip([lh_ind, rh_ind, range(len(X))], ('_lh', '_rh', '')):
        Xvalues = {
            'Xnet': X[ind],
            'xyz': xyz[ind],
            'color': rgb_values[ind],
            'names': np.array(n_name)[np.nonzero(X[ind])[0] + ind[0]],  # to mention only the significant nodes
            'names_idx': np.nonzero(X[ind])[0]
        }

        nodes_file = os.path.join(plot_path, f'{X_name}{side}.mat')
        sio.savemat(nodes_file, Xvalues)
        print(nodes_file)

