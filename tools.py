"""
=================================
            SYNESNET
=================================
This module is design to load data.
"""

from copy import deepcopy
import numpy as np
import scipy.io as sio
from scipy import io
from config import DATA_DIR, NET_DIR, NODE_FILE, INFO_FILE, N_SUB, CORR_TYPE


def load_net_metrics(
    metric, corr_type=CORR_TYPE, net_path=NET_DIR, n_sub=N_SUB, idx_select=slice(None)
):
    """
    Load metrics computed and saved previously.
    :param metric: "coreness" or "strength"
    :param corr_type: "_thr" for metrics computed on positive connectivity matrices.
    :param net_path: path where network metric file is saved
    :param n_sub: total number of subjects (34)
    :param idx_select: selection of nodes to load
    :return:
        x_syn: matrix with vectors for selected metric for each synesthets subject
        x_ctr: matrix with vectors for selected metric for each control subject
    """
    xnet = np.array(
        [
            io.loadmat(net_path / f"net_metrics_Subject{str(sub).zfill(3)}{corr_type}")[
                metric
            ][0]
            for sub in range(1, n_sub + 1)
        ]
    )

    # Create an empty matrix of subjects and nodes and replace the significant nodes with its values
    x = np.zeros(np.shape(xnet))
    x[:, idx_select] = xnet[:, idx_select]

    # Put inf and nan values to zero
    x[np.isnan(x) | np.isinf(x)] = 0

    # Split synesthetic and control subjects
    x_syn = x[:17, :]
    x_ctr = x[17:, :]

    return x_syn, x_ctr


def load_fc(n_sub=N_SUB):
    """
    Load all connectivity matrices.
    :param n_sub: total number of subjects (34)
    :return:
        fc_syn: FC matrices for 17 synethets subjects
        fc_ctr: FC matrices for 17 control subjects
    """
    fc_all = []
    for sub in np.arange(n_sub) + 1:
        # load correlation
        sub_file = "CorrMatrix_Subject{0}.mat".format(str(sub).zfill(3))
        fc_file = DATA_DIR / "symetrical_corr_mat/" / sub_file
        fc = io.loadmat(fc_file)

        # connectivity matrix
        # work with 3 possibilities of connectivity matrices
        # 1. corr = [-1, 1]
        # Xfc = deepcopy(fc['CorrMatrix'])
        # net_file_Xfc = net_path + '_'.join(('net_metrics', 'Subject{0}.mat'.format(str(sub).zfill(3))))

        # 2. abs(corr) = [0, 1]
        # Xfc_abs = abs(deepcopy(fc['CorrMatrix']))
        # net_file_abs = net_path + '_'.join(('net_metrics', 'Subject{0}'.format(str(sub).zfill(3)), '_abs.mat'))

        # 3. corr[corr>=0] = [0, 1]  # winning matrices !!!
        fc_thr = deepcopy(fc["CorrMatrix"])
        fc_thr[fc_thr <= 0] = 0
        np.fill_diagonal(fc_thr, 0)

        # Append matrices
        fc_all.append(fc_thr)

    # Split synesthetic and control subjects
    fc_syn = np.array(fc_all[:17])
    fc_ctr = np.array(fc_all[17:])

    return fc_syn, fc_ctr


def load_node_names():
    """
    Load nodes names from file.
    :return:
        n_name: node names abbreviation
        n_name_full: node names complete
    """
    n_name, n_name_full = [], []
    with open(NODE_FILE, "r") as filestream:
        for line in filestream:
            name, name_full = line.split(",")[:2]
            n_name.append(name)
            n_name_full.append(name_full.strip())

    return n_name, n_name_full


def load_xyz():
    """
    Load 3D positions from file.
    :return:
        xyz: matrix with 3D positions [246 x 3]
    """
    xyz = io.loadmat(INFO_FILE)["xyz"][0][:-1]
    xyz = np.array([x[0] for x in xyz])
    return xyz


def save_mat_file(X, xyz, rgb_values, n_name, X_name, plot_path):
    """
    Save .mat to do 3D brain plots.
    :param X: vector with nodes values
    :param xyz: 3D nodes positions
    :param rgb_values: 4D matrix with nodes colors [246 x 4]
    :param n_name: nodes names
    :param X_name: file name
    :param plot_path: path to save
    :return: save .mat
    """
    # select index for left and right hemisphere
    lh_ind = [
        index
        for index, element in enumerate(np.array(load_node_names()[0]))
        if element.endswith("_L")
    ]
    rh_ind = [
        index
        for index, element in enumerate(np.array(load_node_names()[0]))
        if element.endswith("_R")
    ]

    for ind, side in zip([lh_ind, rh_ind, range(len(X))], ("_lh", "_rh", "")):
        Xvalues = {
            "Xnet": X[ind],
            "xyz": xyz[ind],
            "color": rgb_values[ind],
            "names": np.array(n_name)[
                np.nonzero(X[ind])[0] + ind[0]
            ],  # to mention only the significant nodes
            "names_idx": np.nonzero(X[ind])[0],
        }

        sio.savemat(plot_path / f"{X_name}{side}.mat", Xvalues)
        print(plot_path / f"{X_name}{side}.mat")
