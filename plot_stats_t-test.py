"""
=================================
            SYNESNET
=================================
Plot statistical analysis performed in net_stats_t-test.py
results in .png and also save for complementary 3D oplot in matlab
"""

import os.path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from os import path
import pandas as pd
import sys
from viz.netlocal import plot_3d_local_metric
from tools import load_xyz, load_node_names, save_mat_file


sys.path.append(path.abspath('../netviz'))
path = '/Users/juliana.gonzalez/ownCloud/graph_analysis/'
plot_path = os.path.join(os.getcwd(), 'plots', 'glb', 'new')

# CONSTANTS
SELECTION = True  # True: select base on strength significance
P_VAL = 0.05
corr_type = '_thr'
metric_list = ['strength', 'coreness']
cmap = plt.colormaps['Spectral'].reversed()

# nodes nodes positions and names
xyz = load_xyz()
n_name, n_name_full = load_node_names()

# Open strength t-val file
stats_file = os.path.join(path, 'results', 'stats_results.csv')
df_stats = pd.read_csv(stats_file)

for metric in metric_list:
    X_name = f'{metric}{corr_type}_t-val'
    p_val_type = 'p-val' if metric == 'strength' else 'p-val_corrected'
    idx_select = np.array(df_stats['node_idx'][(df_stats['metric'] == metric) & (df_stats[p_val_type] < P_VAL)])
    t_vals = np.array(df_stats['t-val'][(df_stats['metric'] == metric) & (df_stats[p_val_type] < P_VAL)])

    # Create an empty matrix of subjects and nodes and replace the significant nodes with its values
    X = np.zeros(np.shape(n_name))
    X[idx_select] = t_vals

    # Set colormap parameters
    norm = colors.TwoSlopeNorm(vmin=-max(abs(X)), vmax=max(abs(X)), vcenter=0)
    X_max = None
    kwargs = {"norm": norm}

    # Plot and generate scatter info to plot in matlab
    X_size = abs(pow(X, 2) / max(abs(pow(X, 2)))) * 80
    fig, ax, scatter, cbar = plot_3d_local_metric(X_size, X, xyz, n_name, cmap=cmap, return_scatter=True, **kwargs)
    # plt.savefig(os.path.join(plot_path, f'{X_name}.png'), transparent=True)
    plt.show()

    cmap = scatter.get_cmap()
    rgb_values = cmap(norm(X))

    # save .mat to plot 3D brain with matlab
    # for each hemisphere individually
    save_mat_file(X, xyz, rgb_values, n_name, X_name, plot_path)

