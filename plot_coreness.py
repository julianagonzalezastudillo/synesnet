#%% PLOT
import os.path
import numpy as np
from scipy import io
import pandas as pd
import matplotlib.pyplot as plt
from os import path
import sys
sys.path.append(path.abspath('../netviz'))
from viz.netlocal import plot_3d_local_metric
from tools import load_net_metrics, load_xyz, load_node_names


path = '/Users/juliana.gonzalez/ownCloud/graph_analysis/'
net_path = os.path.join(path, 'net_metrics', 'new')
plot_path = os.path.join(os.getcwd(), 'plots', 'glb', 'new')
strength_stat_file = os.path.join(os.getcwd(), 'plots', 'glb', 'strength_thr_t-val.mat')
coreness_stat_file = os.path.join(plot_path, 'coreness_norm_by_rand_conserve_strenght_distribution_thr_t-val_selection.mat')

# CONSTANTS
SELECTION = True  # True: select base on strength significance
corr_type = '_thr'
metric = 'coreness'

# number of subjects
num_sub = len(os.listdir(os.path.join(path, 'symetrical_corr_mat')))

# nodes nodes positions and names
xyz = load_xyz()
n_name, n_name_full = load_node_names()

# select index for left and right hemisphere
lh_ind = [index for index, element in enumerate(np.array(n_name)) if element.endswith('_L')]
rh_ind = [index for index, element in enumerate(np.array(n_name)) if element.endswith('_R')]

# Open strength t-val file
if SELECTION:
    strength_stat = io.loadmat(strength_stat_file)
    idx_select = strength_stat['names_idx'][0]
    coreness_stat = io.loadmat(coreness_stat_file)
    idx_select_coreness = coreness_stat['names_idx'][0]
    idx_select = list(set(idx_select) & set(idx_select_coreness))
else:
    idx_select = slice(None)

# Load net metrics for all subjects
Xnet = load_net_metrics(net_path, metric, corr_type, num_sub, idx_select)

# Create an empty matrix of subjects and nodes and replace the significant nodes with its values
X_ = np.zeros([np.shape(Xnet)[0], len(n_name)])  # 34x246
X_[:, idx_select] = Xnet[:, idx_select]
Xnet = X_

# Split synesthetic and control subjects
Xnet[np.isnan(Xnet) | np.isinf(Xnet)] = 0
Xnet_syn = Xnet[:17, :]
Xnet_ctr = Xnet[17:, :]

# Mean across subjects
Xnet_syn_mean = Xnet_syn.mean(axis=0)
Xnet_ctr_mean = Xnet_ctr.mean(axis=0)

# Get max and min across all subjects without node selection
Xnet_ = load_net_metrics(net_path, metric, corr_type, num_sub, slice(None))
Xnet_syn_mean_ = np.mean(Xnet_[:17, :], axis=0)
Xnet_ctr_mean_ = np.mean(Xnet_[17:, :], axis=0)
X_max = np.max([Xnet_syn_mean_, Xnet_ctr_mean_])
X_min = np.min([Xnet_syn_mean_, Xnet_ctr_mean_])
norm = plt.Normalize(vmin=X_min, vmax=X_max)
kwargs = {"norm": norm}

# Create DataFrame for results
pd.set_option("display.precision", 2)
df = pd.DataFrame({'node': np.array(n_name)[idx_select],
                   'node_complete': np.array(n_name_full)[idx_select],
                   'node_idx': idx_select,
                   f'{metric}_synes': Xnet_syn.mean(axis=0)[idx_select],
                   f'{metric}_ctr': Xnet_ctr.mean(axis=0)[idx_select],
                   })

for X, X_name in zip([Xnet_syn_mean, Xnet_ctr_mean],
                     [f'{metric}{corr_type}_mean_syn',
                      f'{metric}{corr_type}_mean_ctr']):

    # Plot and generate scatter info to plot in matlab
    X_size = abs(pow(X, 2) / pow(X_max, 2)) * 80
    fig, ax, scatter, cbar = plot_3d_local_metric(X_size, X, xyz, n_name, return_scatter=True, **kwargs)
    X_name += '_selection' if SELECTION else ''
    # plt.savefig(os.path.join(plot_path, f'{X_name}.png'), transparent=True)
    plt.show()

    cmap = scatter.get_cmap()
    rgb_values = cmap(norm(X))

    # save .mat to plot 3D brain with matlab
    # for each hemisphere individually
    for ind, side in zip([lh_ind, rh_ind, range(len(X))], ('_lh', '_rh', '')):
        Xvalues = {'Xnet': X[ind],
                   'xyz': xyz[ind],
                   'color': rgb_values[ind],
                   'names': np.array(n_name)[np.nonzero(X[ind])[0]+ind[0]],  # to mention only the significant nodes
                   'names_idx': np.nonzero(X[ind])[0]
                   }

        nodes_file = os.path.join(plot_path, f'{X_name}{side}.mat')
        # sio.savemat(nodes_file, Xvalues)
        print(nodes_file)
