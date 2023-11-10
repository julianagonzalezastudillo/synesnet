"""
=================================
            SYNESNET
=================================
Plot statistical analysis performed in net_stats_t-test.py for strength and coreness.
Save results in .png and also save for complementary 3D plot in matlab.
Also generate .mat file to plot reference nodes in matlab.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
from viz.netlocal import plot_3d_local_metric
from tools import load_xyz, load_node_names, save_mat_file
from config import DATA_DIR, PLOT_DIR, P_VAL, CORR_TYPE


# CONSTANTS
SELECTION = True  # True: select base on strength significance
metric_list = ["strength", "coreness"]
cmap = plt.colormaps["Spectral"].reversed()

# nodes nodes positions and names
xyz = load_xyz()
n_name, n_name_full = load_node_names()

# Open strength t-val file
stats_file = DATA_DIR / "results" / "stats_results.csv"
df_stats = pd.read_csv(stats_file)

for metric in metric_list:
    X_name = f"{metric}{CORR_TYPE}_t-val"
    p_val_type = "p-val" if metric == "strength" else "p-val_corrected"
    idx_select = np.array(
        df_stats["node_idx"][
            (df_stats["metric"] == metric) & (df_stats[p_val_type] < P_VAL)
        ]
    )
    t_vals = np.array(
        df_stats["t-val"][
            (df_stats["metric"] == metric) & (df_stats[p_val_type] < P_VAL)
        ]
    )

    # Create an empty matrix of subjects and nodes and replace the significant nodes with its values
    X = np.zeros(np.shape(n_name))
    X[idx_select] = t_vals

    # Set colormap parameters
    norm = colors.TwoSlopeNorm(vmin=-max(abs(X)), vmax=max(abs(X)), vcenter=0)
    X_max = None
    kwargs = {"norm": norm}

    # Plot and generate scatter info to plot in matlab
    X_size = abs(pow(X, 2) / max(abs(pow(X, 2)))) * 80
    fig, ax, scatter, cbar = plot_3d_local_metric(
        X_size, X, xyz, n_name, cmap=cmap, return_scatter=True, **kwargs
    )
    # plt.savefig(PLOT_DIR / f"{X_name}.png", transparent=True)
    plt.show()

    cmap = scatter.get_cmap()
    rgb_values = cmap(norm(X))

    # save .mat to plot 3D brain with matlab
    # for each hemisphere individually
    # save_mat_file(X, xyz, rgb_values, n_name, X_name, PLOT_DIR)

# Plot node references
# the 16 significant nodes order as in the results table
ref_name = "reference_numbers"

# Filter significant nodes directly using boolean indexing
significant_coreness_nodes = (
    df_stats[df_stats["metric"] == "coreness"]
    .sort_values(by="p-val_corrected")["node_idx"]
    .values
)

# Create an array X_all with 1s for significant nodes
X_all = np.zeros(np.shape(X))
X_all[significant_coreness_nodes] = 0.1

# Replace n_names by numbers in the table
for i, idx in enumerate(significant_coreness_nodes):
    n_name[idx] = str(i + 1)

# White nodes
rgb_values = np.ones(np.shape(rgb_values))

# save_mat_file(X_all, xyz, rgb_values, n_name, ref_name, PLOT_DIR)
