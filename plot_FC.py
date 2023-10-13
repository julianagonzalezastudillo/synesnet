import os.path
from scipy import io
import numpy as np
from copy import deepcopy
from tools import load_net_metrics
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as colors
import scipy.io as sio
from viz.netlocal import plot_3d_local_metric


path = '/Users/juliana.gonzalez/ownCloud/graph_analysis/'
net_path = os.path.join(path, 'net_metrics')
strength_stat_file = os.path.join(os.getcwd(), 'plots', 'glb', 'strength_thr_t-val.mat')
node_file = os.path.join(path, 'BN_Atlas_246_LUT_reoriented.txt')
results_file = os.path.join(path, 'resultsROI_Subject001_Condition001.mat')

corr_type = "_thr"
num_sub = len(os.listdir(os.path.join(path, 'symetrical_corr_mat')))
idx_select = [5, 70, 128, 133, 9, 64, 68, 25, 10]  # obtained from previous results

# nodes names
n_name, n_name_full = [], []
with open(node_file, "r") as filestream:
    for line in filestream:
        name, name_full = line.split(",")[:2]
        n_name.append(name)
        n_name_full.append(name_full.strip())

# nodes positions
xyz = io.loadmat(results_file)['xyz'][0][:-1]
xyz = np.array([x[0] for x in xyz])

# Load net metrics for all subjects
Xnet = load_net_metrics(net_path, "coreness", corr_type, num_sub, slice(None))  # [sub x nodes]

# Split synesthetic and control subjects
# Xnet[np.isnan(Xnet) | np.isinf(Xnet)] = 0
Xnet_syn = Xnet[:17, :]
Xnet_ctr = Xnet[17:, :]
Xnet_syn_mean = Xnet[:17, :].mean(axis=0)
Xnet_ctr_mean = Xnet[17:, :].mean(axis=0)

# Load all connectivity matrices
Xfc_all = []
for sub in np.arange(num_sub)+1:
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

# Reshape all matrices
Xfc_all = np.concatenate([mat.reshape(1, *np.shape(Xfc_thr)) for mat in Xfc_all], axis=0)
Xfc_syn = Xfc_all[:17]
Xfc_ctr = Xfc_all[17:]
Xfc_syn_mean = Xfc_all[:17].mean(axis=0)
Xfc_ctr_mean = Xfc_all[17:].mean(axis=0)

# save matrices for plot
Xfc = Xfc_syn_mean[idx_select] - Xfc_ctr_mean[idx_select]
Xfc_max = np.max(Xfc)
Xfc_min = np.min(Xfc)

cmap_colors = ['#11205E', '#203FB6', '#86CAFF', 'white', '#FFEC4A', '#F62336', '#80121B']
positions = np.linspace(0, 1, len(cmap_colors))
cmap = LinearSegmentedColormap.from_list('custom_colormap', list(zip(positions, cmap_colors)))

# plot each node in
for X, node_idx in zip(np.vstack((Xfc, np.mean(Xfc_syn_mean[idx_select] - Xfc_ctr_mean[idx_select], axis=0))),
                       np.append(idx_select, 247)):
    X_size = abs(pow(X, 2) / max(abs(pow(X, 2)))) * 80

    norm = colors.TwoSlopeNorm(vmin=-max(abs(X)), vmax=max(abs(X)), vcenter=0)
    kwargs = {"norm": norm}
    fig, ax, scatter, cbar = plot_3d_local_metric(X_size, X, xyz, n_name, cmap=cmap, return_scatter=True, **kwargs)
    plt.show()

    cmap = scatter.get_cmap()
    rgb_values = cmap(norm(X))

    file_name = f"{node_idx}_Xfc_mean.mat"
    # file_name = f"Xfc_mean.mat"
    X_ = {'Xnet': X,
           'xyz': xyz,
           'color': rgb_values,
           'names': np.array(n_name),  # to mention only the significant nodes
           'names_idx': np.nonzero(X)
           }
    # sio.savemat(file_name, X_)
    # print(f"{node_idx}: {n_name_full[node_idx]}")


# save matrices of FC difference
file_name = "Xfc_syn-Xfc_ctr(9nodes).mat"
idx_select = [5, 70, 128, 133, 9, 64, 68, 25, 10]  # obtained from previous results
X_diff = Xfc_syn_mean[idx_select] - Xfc_ctr_mean[idx_select]
X_ = {'diff_correlation_syn_ctr_nine_nodes': X_diff}
# sio.savemat(file_name, X_)

file_name = "Xfc_syn-Xfc_ctr_mean(9nodes).mat"
X_diff_mean = np.mean(Xfc_syn_mean[idx_select] - Xfc_ctr_mean[idx_select], axis=0)
X_ = {'mean_diff_correlation_syn_ctr_nine_nodes': X_diff_mean}
# sio.savemat(file_name, X_)

file_name = "Xfc_syn_sub(9nodes).mat"
X_sub = Xfc_syn[:, idx_select, :]
X_ = {'correlation_syn_nine_nodes': X_sub}
# sio.savemat(file_name, X_)

file_name = "Xfc_ctr_sub(9nodes).mat"
X_sub = Xfc_ctr[:, idx_select, :]
X_ = {'correlation_ctr_nine_nodes': X_sub}
# sio.savemat(file_name, X_)