import os.path
from scipy import io
import numpy as np
from copy import deepcopy
from tools import load_net_metrics
import matplotlib.pyplot as plt

path = '/Users/juliana.gonzalez/ownCloud/graph_analysis/'
net_path = os.path.join(path, 'net_metrics')
strength_stat_file = os.path.join(os.getcwd(), 'plots', 'glb', 'strength_thr_t-val.mat')
node_file = os.path.join(path, 'BN_Atlas_246_LUT_reoriented.txt')

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

# mean excluding zero
Xfc_syn_sum_along_axis0 = np.sum(Xfc_syn, axis=0)
Xfc_syn_non_zero_count = np.count_nonzero(Xfc_syn, axis=0)
Xfc_syn_mean_non_zero = np.divide(Xfc_syn_sum_along_axis0, Xfc_syn_non_zero_count, where=Xfc_syn_non_zero_count != 0)

Xfc_ctr_sum_along_axis0 = np.sum(Xfc_ctr, axis=0)
Xfc_ctr_non_zero_count = np.count_nonzero(Xfc_ctr, axis=0)
Xfc_ctr_mean_non_zero = np.divide(Xfc_ctr_sum_along_axis0, Xfc_ctr_non_zero_count, where=Xfc_ctr_non_zero_count != 0)


for Xnet_mean, Xfc_mean, plot_name in zip([Xnet_syn_mean, Xnet_ctr_mean], [Xfc_syn_mean, Xfc_ctr_mean], ["SYN", "CTR"]):
    fig, axs = plt.subplots(3, 3, figsize=(12, 12))
    # Loop through node_idx and populate the subplots
    for row in range(3):
        for col in range(3):
            node_idx = idx_select[row * 3 + col]

            # Scatter plot
            idx_non_zero = np.where(Xfc_mean[node_idx] != 0)
            # axs[row, col].scatter(Xnet_mean, Xfc_[:, node_idx, :], alpha=0.5)
            x = Xnet_mean[idx_non_zero]
            y = Xfc_mean[node_idx][idx_non_zero]
            axs[row, col].scatter(x, y, alpha=0.5)

            # Calculate the correlation coefficient
            correlation = np.corrcoef(x, y)[0, 1]

            # Create a correlation line
            fit = np.polyfit(x, y, 2)
            coefficients = np.polyfit(x, y, 2)
            correlation_line = np.poly1d(fit)

            # Plot the correlation line
            # axs[row, col].plot(Xnet_mean, correlation_line(Xnet_mean), color='red',
            #                    label=f'Correlation Line (r = {correlation:.2f})')
            axs[row, col].plot(np.sort(x), correlation_line(np.sort(x)),
                               label=f'(r = {correlation:.2f})', color='red')

            # Set labels and title
            axs[row, col].set_xlabel('coreness')
            axs[row, col].set_ylabel('correlation')
            axs[row, col].set_title(f'{n_name_full[node_idx]}, C={Xnet_mean[node_idx]:.2f}')

            # Add a legend
            axs[row, col].legend()

    fig.suptitle(f"{plot_name}", fontsize=16)
    plt.tight_layout()

    # Show the subplots
    plt.show()
