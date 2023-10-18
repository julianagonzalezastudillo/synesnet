from copy import deepcopy
import numpy as np
import os.path
from scipy import io

path = '/Users/juliana.gonzalez/ownCloud/graph_analysis/'
node_file = os.path.join(path, 'BN_Atlas_246_LUT_reoriented.txt')
results_file = os.path.join(path, 'resultsROI_Subject001_Condition001.mat')


def load_net_metrics(net_path, metric, corr_type, num_sub, idx_select):
    Xnet = np.array([
        io.loadmat(os.path.join(net_path, f'net_metrics_Subject{str(sub).zfill(3)}{corr_type}'))[metric][0]
        for sub in range(1, num_sub + 1)])
    X_ = np.zeros(np.shape(Xnet))
    X_[:, idx_select] = Xnet[:, idx_select]
    return X_


def load_fc():
    # Load all connectivity matrices
    num_sub = len(os.listdir(os.path.join(path, 'symetrical_corr_mat')))
    Xfc_all = []
    for sub in np.arange(num_sub) + 1:
        # load correlation
        sub_file = 'CorrMatrix_Subject{0}.mat'.format(str(sub).zfill(3))
        fc_file = path + 'symetrical_corr_mat/' + sub_file
        fc = io.loadmat(fc_file)

        # 3. corr[corr>=0] = [0, 1]  # winning matrices !!!
        Xfc_thr = deepcopy(fc['CorrMatrix'])
        Xfc_thr[Xfc_thr <= 0] = 0
        np.fill_diagonal(Xfc_thr, 0)

        # Append matrices
        Xfc_all.append(Xfc_thr)

    Xfc_syn = np.array(Xfc_all[:17])
    Xfc_ctr = np.array(Xfc_all[17:])

    return Xfc_syn, Xfc_ctr


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
