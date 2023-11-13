"""
=================================
            SYNESNET
=================================
Plot coreness difference between synesthets and controls in .png and also save for complementary 3D to plot in matlab
"""

import numpy as np
import matplotlib.pyplot as plt
from viz.netlocal import plot_3d_local_metric
from tools import load_net_metrics, load_xyz, load_node_names, save_mat_file
from config import PLOT_DIR, CORR_TYPE


# CONSTANTS
X_name = f"coreness{CORR_TYPE}_diff_syn-ctr"

# nodes nodes positions and names
xyz = load_xyz()
n_name, n_name_full = load_node_names()

# Load net metrics for all subjects
Xnet_syn, Xnet_ctr = load_net_metrics("coreness")

# Mean across subjects
Xnet_syn_mean = Xnet_syn.mean(axis=0)
Xnet_ctr_mean = Xnet_ctr.mean(axis=0)

# Coreness difference between synesthets and controls
X = Xnet_syn_mean - Xnet_ctr_mean

# Set plot parameters
X_max = np.max(abs(X))
norm = plt.Normalize(vmin=-X_max, vmax=X_max)
kwargs = {"norm": norm}
X_size = abs(pow(X, 2) / pow(X_max, 2)) * 80
fig, ax, scatter, cbar = plot_3d_local_metric(
    X_size, X, xyz, n_name, return_scatter=True, **kwargs
)
# plt.savefig(PLOT_DIR / f"{X_name}.png"), transparent=True)
plt.show()

# Get parameters for matlab plot
cmap = scatter.get_cmap()
rgb_values = cmap(norm(X))

# save .mat to plot 3D brain with matlab
# for each hemisphere individually
# save_mat_file(X, xyz, rgb_values, n_name, X_name, PLOT_DIR)
