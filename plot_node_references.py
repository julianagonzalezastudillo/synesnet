#%% PLOT
import os.path
import numpy as np
from scipy import io
from scipy import stats
import pandas as pd
import statsmodels.stats.multitest as smt
import scipy.io as sio
from statsmodels.stats.multitest import multipletests
import statsmodels
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap
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
P_VAL = 0.05
SELECTION = True  # True: select base on strength significance
corr_type = '_thr'
metric_list = ['coreness_norm_by_rand_conserve_strenght_distribution']

cmap = plt.colormaps['Spectral'].reversed()

p_val_type = 'p-val_corrected' if SELECTION else 'p-val'

# number of subjects
num_sub = len(os.listdir(os.path.join(path, 'symetrical_corr_mat')))

# nodes nodes positions and names
xyz = load_xyz()
n_name, n_name_full = load_node_names()

# Open strength t-val file
strength_stat = io.loadmat(strength_stat_file)
idx_select = strength_stat['names_idx'][0]
idx_select = idx_select if SELECTION else slice(None)

for net_key in metric_list:
    # Load net metrics for all subjects
    Xnet = load_net_metrics(net_path, net_key, corr_type, num_sub, idx_select)

    # Create an empty matrix of subjects and nodes and replace the significant nodes with its values
    X_ = np.zeros([np.shape(Xnet)[0], len(n_name)])  # 34x246
    X_[:, idx_select] = Xnet[:, idx_select]
    Xnet = X_
    lh_ind = [index for index, element in enumerate(np.array(n_name)) if element.endswith('_L')]
    rh_ind = [index for index, element in enumerate(np.array(n_name)) if element.endswith('_R')]

    # Split synesthetic and control subjects
    Xnet[np.isnan(Xnet) | np.isinf(Xnet)] = 0
    Xnet_syn = Xnet[:17, :]
    Xnet_ctr = Xnet[17:, :]

    # Perform t-test
    t_val, p_val = stats.ttest_ind(Xnet_syn, Xnet_ctr)
    t_val[np.isnan(t_val)] = 0  # Replace NaNs with 0

    # Perform FDR correction
    p_val = p_val[idx_select]
    rejected, corrected_p_values, _, _ = multipletests(p_val, method='fdr_bh')
    statsmodels.stats.multitest.fdrcorrection(p_val, alpha=P_VAL, method='indep', is_sorted=False)

    # Corrected by hand FDR (Benjamini-Hochberg)
    rank = np.arange(len(p_val))+1
    rank_idx = np.argsort(p_val)
    N = len(p_val)
    p_val_corrected_ = np.sort(p_val)*N/rank
    p_val_corrected = np.zeros(len(p_val))
    p_val_corrected[rank_idx] = p_val_corrected_

    # Create DataFrame for results
    pd.set_option("display.precision", 2)
    df = pd.DataFrame({'node': np.array(n_name)[idx_select],
                       'node_complete': np.array(n_name_full)[idx_select],
                       'node_idx': idx_select,
                       f'{net_key}_synes': Xnet_syn.mean(axis=0)[idx_select],
                       f'{net_key}_ctr': Xnet_ctr.mean(axis=0)[idx_select],
                       't-val': np.array(t_val)[idx_select],
                       'p-val': p_val,
                       'p-val_corrected': p_val_corrected
                       })

    # for t-values, keep significant t-values
    X_ = np.zeros(len(n_name))  # 246
    X_[idx_select] = np.where((df[p_val_type] > P_VAL), 0, df['t-val'])
    X_t_val = X_

    X = X_t_val

    # Plot and generate scatter info to plot in matlab
    fig, ax, scatter, cbar = plot_3d_local_metric(X, X, xyz, n_name, cmap=cmap, return_scatter=True)
    plt.show()
    cmap = scatter.get_cmap()
    rgb_values = cmap(X)

    # Create numeric references for selected nodes (strength + coreness)
    df_sorted = df.sort_values(by='p-val_corrected').reset_index(drop=True)
    df_sorted['Order'] = df_sorted.index + 1  # +1 to start indexing from 1

    # for each hemisphere individually
    for ind, side in zip([lh_ind, rh_ind], ('lh', 'rh')):
        mask = (np.array(df_sorted["node_idx"]) >= ind[0]) & (np.array(df_sorted["node_idx"]) < ind[-1])
        Xvalues = {'Xnet': np.zeros(np.shape(X))[ind],
                   'xyz': xyz[ind],
                   'color': rgb_values[ind],
                   'names': list(df_sorted["Order"][mask]),
                   'names_idx': list(df_sorted["node_idx"][mask] - ind[0])
                   }
        nodes_file = os.path.join(plot_path, f'reference_numbers_{side}.mat')
        # sio.savemat(nodes_file, Xvalues)
