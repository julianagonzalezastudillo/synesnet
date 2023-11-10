"""
=================================
            SYNESNET
=================================
Plot coreness results in .png and also save for complementary 3D oplot in matlab
"""

import os.path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from viz.netlocal import plot_3d_local_metric
from tools import load_net_metrics, load_xyz, load_node_names, save_mat_file
import seaborn as sns
from matplotlib.ticker import MaxNLocator
from config import DATA_DIR, PLOT_PATH, P_VAL, CORR_TYPE
from pathlib import Path


PLOT_DIR = Path(PLOT_PATH)
PLOT_DIR.mkdir(parents=True, exist_ok=True)

# CONSTANTS
SELECTION = False  # True: select base on strength significance
metric = 'coreness'

# nodes nodes positions and names
xyz = load_xyz()
n_name, n_name_full = load_node_names()

if SELECTION:
    # Open strength t-val file
    stats_file = os.path.join(DATA_DIR, 'results', 'stats_results.csv')
    df_stats = pd.read_csv(stats_file)
    idx_select = np.array(df_stats['node_idx'][(df_stats['metric'] == metric) & (df_stats['p-val_corrected'] < P_VAL)])
else:
    idx_select = slice(None)

# Load net metrics for all subjects
Xnet_syn, Xnet_ctr = load_net_metrics(metric, idx_select=idx_select)

# Mean across subjects
Xnet_syn_mean = Xnet_syn.mean(axis=0)
Xnet_ctr_mean = Xnet_ctr.mean(axis=0)

# Get max and min across all subjects without node selection
Xnet_syn_, Xnet_ctr_ = load_net_metrics(metric, idx_select=slice(None))
Xnet_syn_mean_ = np.mean(Xnet_syn_, axis=0)
Xnet_ctr_mean_ = np.mean(Xnet_ctr_, axis=0)
X_max = np.max([Xnet_syn_mean_, Xnet_ctr_mean_])
X_min = np.min([Xnet_syn_mean_, Xnet_ctr_mean_])
norm = plt.Normalize(vmin=X_min, vmax=X_max)
kwargs = {"norm": norm}

for X, X_name in zip([Xnet_syn_mean, Xnet_ctr_mean],
                     [f'{metric}{CORR_TYPE}_mean_syn',
                      f'{metric}{CORR_TYPE}_mean_ctr']):

    # Plot and generate scatter info to plot in matlab
    X_size = abs(pow(X, 2) / pow(X_max, 2)) * 80
    fig, ax, scatter, cbar = plot_3d_local_metric(X_size, X, xyz, n_name, return_scatter=True, **kwargs)
    X_name += '_selection' if SELECTION else ''
    # plt.savefig(PLOT_DIR / f"{X_name}.png", transparent=True)
    plt.show()

    cmap = scatter.get_cmap()
    rgb_values = cmap(norm(X))

    # save .mat to plot 3D brain with matlab
    # for each hemisphere individually
    save_mat_file(X, xyz, rgb_values, n_name, X_name, PLOT_DIR)

#%% Plot coreness distribution
fig, ax = plt.subplots(figsize=(7, 6), dpi=300)

# BLUE
sns.histplot(Xnet_ctr_mean_, kde=True, ax=ax, label='KDE', color='#203FB6', element="step", linewidth=3)
kde_line = ax.lines[-1]  # The KDE line is the last line in the axis
kde_line.set_linewidth(3)

# RED
sns.histplot(Xnet_syn_mean_, kde=True, ax=ax, label='KDE', color='#F62336', element="step", linewidth=3)
kde_line = ax.lines[-1]
kde_line.set_linewidth(3)

plt.xlim(X_min, X_max)
# plt.xlabel("coreness")
plt.ylabel(" ")
ax.yaxis.set_major_locator(MaxNLocator(5))
ax.xaxis.set_major_locator(MaxNLocator(5))

# Remove the top and right axes
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

# save
# plt.savefig(PLOT_DIR / "coreness_distribution.png", bbox_inches='tight', pad_inches=0, dpi=300, transparent=True)
plt.show()
