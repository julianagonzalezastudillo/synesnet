import os.path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tools import load_fc, load_net_metrics, load_node_names
from config import DATA_DIR, PLOT_PATH, P_VAL


metric = "coreness"

# Open stats file
stats_file = os.path.join(DATA_DIR, 'results', 'stats_results.csv')
df_stats = pd.read_csv(stats_file)
idx_select = np.array(df_stats['node_idx'][(df_stats['metric'] == metric) & (df_stats['p-val_corrected'] < P_VAL)])

# nodes positions and names
n_name, full_n_name = load_node_names()

# Load net metrics for all subjects
Xnet_syn, Xnet_ctr = load_net_metrics(metric)  # [sub x nodes]
Xnet_syn_mean = Xnet_syn.mean(axis=0)
Xnet_ctr_mean = Xnet_ctr.mean(axis=0)

# Load connectivity matrices for all subjects
Xfc_syn, Xfc_ctr = load_fc()
Xfc_syn_mean = Xfc_syn.mean(axis=0)
Xfc_ctr_mean = Xfc_ctr.mean(axis=0)

# mean excluding zero
Xfc_syn_sum_along_axis0 = np.sum(Xfc_syn, axis=0)
Xfc_syn_non_zero_count = np.count_nonzero(Xfc_syn, axis=0)
Xfc_syn_mean_non_zero = np.divide(Xfc_syn_sum_along_axis0, Xfc_syn_non_zero_count, where=Xfc_syn_non_zero_count != 0)

Xfc_ctr_sum_along_axis0 = np.sum(Xfc_ctr, axis=0)
Xfc_ctr_non_zero_count = np.count_nonzero(Xfc_ctr, axis=0)
Xfc_ctr_mean_non_zero = np.divide(Xfc_ctr_sum_along_axis0, Xfc_ctr_non_zero_count, where=Xfc_ctr_non_zero_count != 0)

ncols = 2
nrows = 3
for Xnet_mean, Xfc_mean, plot_name in zip([Xnet_syn_mean, Xnet_ctr_mean], [Xfc_syn_mean, Xfc_ctr_mean], ["SYN", "CTR"]):
    fig, axs = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
    # Loop through node_idx and populate the subplots
    for node_idx, row, col in zip(idx_select, np.sort(np.tile(np.arange(nrows), ncols)),
                                  np.tile(np.arange(ncols), nrows)):
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
        axs[row, col].set_title(f'{full_n_name[node_idx]}, C={Xnet_mean[node_idx]:.2f}')

        # Add a legend
        axs[row, col].legend()

    fig.suptitle(f"{plot_name}", fontsize=16)
    plt.tight_layout()

    # Show the subplots
    plt.show()
