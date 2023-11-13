"""
=================================
            SYNESNET
=================================
Plot FC differences for the 5 nodes with significant difference between synesthets and controls.
"""

import os.path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as colors
import pandas as pd
import scipy.io as sio
from viz.netlocal import plot_3d_local_metric
import seaborn as sns
from tools import load_fc, load_xyz, load_node_names, save_mat_file
from config import DATA_DIR, PLOT_DIR, P_VAL


# Open strength t-val file
stats_file = os.path.join(DATA_DIR, "results", "stats_results.csv")
df_stats = pd.read_csv(stats_file)
idx_select = np.array(
    df_stats["node_idx"][
        (df_stats["metric"] == "coreness") & (df_stats["p-val_corrected"] < P_VAL)
    ]
)

# nodes positions and names
xyz = load_xyz()
n_name, full_n_name = load_node_names()

# Load all connectivity matrices
Xfc_syn, Xfc_ctr = load_fc()

# Reshape all matrices
Xfc_syn_mean = Xfc_syn.mean(axis=0)
Xfc_ctr_mean = Xfc_ctr.mean(axis=0)

# save matrices for plot
Xfc = Xfc_syn_mean[idx_select] - Xfc_ctr_mean[idx_select]
Xfc_diff_mean = np.mean(Xfc_syn_mean[idx_select] - Xfc_ctr_mean[idx_select], axis=0)

# plot each node in
for X, node_idx in zip(
    np.vstack((Xfc, Xfc_diff_mean)),
    np.append(idx_select, 247),
):
    # Plot and generate scatter info to plot in matlab
    X_size = abs(pow(X, 2) / max(abs(pow(X, 2)))) * 80
    norm = colors.TwoSlopeNorm(vmin=-max(abs(X)), vmax=max(abs(X)), vcenter=0)
    kwargs = {"norm": norm}
    fig, ax, scatter, cbar = plot_3d_local_metric(
        X_size, X, xyz, n_name, return_scatter=True, **kwargs
    )
    # plt.savefig(PLOT_DIR / f"Xfc_mean.png"), transparent=True)
    plt.show()

# save .mat to plot 3D brain with matlab
# for each hemisphere individually
cmap = scatter.get_cmap()
rgb_values = cmap(norm(X))
save_mat_file(X, xyz, rgb_values, n_name, "Xfc_mean", PLOT_DIR)

print(
    f"mean connectivity difference: {np.mean(Xfc_syn_mean[idx_select] - Xfc_ctr_mean[idx_select])}"
)

# save matrices of FC difference
file_name = PLOT_DIR / "Xfc_syn-Xfc_ctr(5nodes).mat"
X_diff = Xfc_syn_mean[idx_select] - Xfc_ctr_mean[idx_select]
X_ = {"diff_correlation_syn_ctr_five_nodes": X_diff}
# sio.savemat(file_name, X_)

file_name = PLOT_DIR / "Xfc_syn-Xfc_ctr_mean(5nodes).mat"
X_diff_mean = np.mean(Xfc_syn_mean[idx_select] - Xfc_ctr_mean[idx_select], axis=0)
X_ = {"mean_diff_correlation_syn_ctr_five_nodes": X_diff_mean}
# sio.savemat(file_name, X_)

file_name = PLOT_DIR / "Xfc_syn_sub(5nodes).mat"
X_sub = Xfc_syn[:, idx_select, :]
X_ = {"correlation_syn_five_nodes": X_sub}
# sio.savemat(file_name, X_)

file_name = PLOT_DIR / "Xfc_ctr_sub(5nodes).mat"
X_sub = Xfc_ctr[:, idx_select, :]
X_ = {"correlation_ctr_five_nodes": X_sub}
# sio.savemat(file_name, X_)


# %% Plot connectivity matrices for nodes with significant difference in terms of coreness
sns.set_theme(style="white")
cmap_colors = ["white", "#FFEC4A", "#F62336", "#80121B"]
positions = np.linspace(0, 1, len(cmap_colors))
cmap = LinearSegmentedColormap.from_list(
    "custom_colormap", list(zip(positions, cmap_colors))
)

xticklabels = np.array(full_n_name)[idx_select]
yticklabels = np.array(full_n_name)[idx_select]

X_max = np.max(
    [Xfc_syn_mean[idx_select][:, idx_select], Xfc_ctr_mean[idx_select][:, idx_select]]
)
X_min = np.min(
    [Xfc_syn_mean[idx_select][:, idx_select], Xfc_ctr_mean[idx_select][:, idx_select]]
)

# Draw the heatmap with the mask and correct aspect ratio
for X, sub_type in zip(
    [Xfc_syn_mean, Xfc_ctr_mean],
    ["SYN", "CTR"],
):
    f, ax = plt.subplots(figsize=(12, 10))
    heatmap = sns.heatmap(
        X[idx_select][:, idx_select],
        cmap=cmap,
        vmax=X_max,
        vmin=X_min,
        linewidths=0.0,
        cbar_kws={"shrink": 0.9},
        xticklabels=xticklabels,
        yticklabels=yticklabels,
    )

    # Rotate the yticklabels
    plt.xticks(rotation=30, fontsize=16)
    plt.yticks(rotation=30, fontsize=16)
    ax.set_aspect("equal")
    plt.tight_layout()
    cbar = heatmap.collections[0].colorbar
    cbar.ax.tick_params(labelsize=16)  # Adjust the font size as needed
    plt.title(sub_type, fontsize=24)
    plt.show()
