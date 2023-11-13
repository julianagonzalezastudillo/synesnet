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
    # Load all connectivity matrices
    fc_all = []
    for sub in np.arange(n_sub) + 1:
        # load correlation
        sub_file = "CorrMatrix_Subject{0}.mat".format(str(sub).zfill(3))
        fc_file = DATA_DIR / "symetrical_corr_mat/" / sub_file
        fc = io.loadmat(fc_file)

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
    # nodes names
    n_name, n_name_full = [], []
    with open(NODE_FILE, "r") as filestream:
        for line in filestream:
            name, name_full = line.split(",")[:2]
            n_name.append(name)
            n_name_full.append(name_full.strip())

    return n_name, n_name_full


def load_xyz():
    xyz = io.loadmat(INFO_FILE)["xyz"][0][:-1]
    xyz = np.array([x[0] for x in xyz])
    return xyz


def save_mat_file(X, xyz, rgb_values, n_name, X_name, plot_path):
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
