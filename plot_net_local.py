#%% PLOT
import os.path
import numpy as np
from scipy import io
from scipy import stats
import pandas as pd
import statsmodels.stats.multitest as smt
import scipy.io as sio
import math
from statsmodels.stats.multitest import multipletests
import statsmodels
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from os import path
import sys
sys.path.append(path.abspath('../netviz'))
# from net3D import nodes_3D
from viz.netlocal import plot_3d_local_metric

path = '/Users/juliana.gonzalez/ownCloud/graph_analysis/'
net_path = path + 'net_metrics/'
num_sub = len(os.listdir(path + 'symetrical_corr_mat'))

# nodes positions
results_file = path + 'resultsROI_Subject001_Condition001.mat'
results = io.loadmat(results_file)
xyz_all = results['xyz'][0][:-1]
xyz_all = np.array([xyz_all[i][0][:] for i in np.arange(np.size(xyz_all))])

# nodes names
n_name_all = []
n_name_full_all = []
with open(path + "BN_Atlas_246_LUT_reoriented.txt", "r") as filestream:
    for line in filestream:
        n_name_all.append(line.split(",")[0])
        n_name_full_all.append(line.split(",")[1].strip())

corr_type = '_thr'
metric_list = ['coreness_norm_strength_selection']
# metric_list = ['strength']

for net_key in metric_list:
    # Load node names
    if net_key.split('_')[-1] == 'selection':
        # Open strength t-val file
        strength_stat_file = os.path.join(os.getcwd(), 'plots', 'glb', 'strength_thr_t-val.mat')
        strength_stat = io.loadmat(strength_stat_file)
        idx_select = strength_stat['names_idx'][0]
        n_name = list(np.array(n_name_all)[idx_select])
        n_name_full = list(np.array(n_name_full_all)[idx_select])
        xyz = xyz_all[idx_select]
    else:
        n_name = n_name_all
        n_name_full = n_name_full_all
        xyz = xyz_all

    lh_ind = [index for index, element in enumerate(n_name) if element.endswith('_L')]
    rh_ind = [index for index, element in enumerate(n_name) if element.endswith('_R')]

    # Load net metrics for all subjects
    Xnet = np.array([io.loadmat(net_path + f'net_metrics_Subject{str(sub).zfill(3)}{corr_type}')[net_key][0]
                     for sub in range(1, num_sub + 1)])

    # Split synesthetic and control subjects
    Xnet_syn = Xnet[:17, :]
    Xnet_ctr = Xnet[17:, :]

    # Perform t-test
    t_val, p_val = stats.ttest_ind(Xnet_syn, Xnet_ctr)
    t_val = [0 if math.isnan(x) else x for x in t_val]  # put to 0 the t-val=nan

    # Perform FDR correction
    rejected, corrected_p_values, _, _ = multipletests(p_val, method='fdr_bh')
    statsmodels.stats.multitest.fdrcorrection(p_val, alpha=0.05, method='indep', is_sorted=False)

    # Corrected by hand FDR (Benjamini-Hochberg)
    rank = np.arange(len(p_val))+1
    rank_idx = np.argsort(p_val)
    N = len(p_val)
    p_val_corrected_ = p_val_corrected = np.sort(p_val)*N/rank
    p_val_corrected[rank_idx] = p_val_corrected_

    # Mean across subjects
    Xnet_syn_mean = Xnet_syn.mean(axis=0)
    Xnet_ctr_mean = Xnet_ctr.mean(axis=0)

    # Create DataFrame for results
    pd.set_option("display.precision", 2)
    df = pd.DataFrame({'node': n_name,
                       'node_complete': n_name_full,
                       f'{net_key}_synes': Xnet_syn.mean(axis=0),
                       f'{net_key}_ctr': Xnet_ctr.mean(axis=0),
                       't-val': t_val,
                       'p-val': p_val,
                       'p-val_corrected': p_val_corrected})

    # for t-values, keep significant t-values
    X_t_val = np.where((df['p-val'] > 0.05), 0, df['t-val'])  # t-val higher than 0.05

    # print Hub list
    print('-' * 100)
    print(df[['node', f'{net_key}_synes']].sort_values(by=f'{net_key}_synes', key=abs).tail(30).to_string())
    print(f'Average synes {net_key}: {Xnet_syn.mean(axis=0).mean()}')

    print('-' * 100)
    print(df[['node', f'{net_key}_ctr']].sort_values(by=f'{net_key}_ctr', key=abs).tail(30).to_string())
    print(f'Average ctr {net_key}: {Xnet_ctr.mean(axis=0).mean()}')

    print('-' * 100)
    print(df.loc[df['p-val'] < 0.05, ["node", "node_complete", "p-val", "p-val_corrected",  "t-val"]].sort_values(by='p-val'))

    for X, X_name in zip([Xnet_syn_mean, Xnet_ctr_mean, X_t_val],
                         [f'{net_key}{corr_type}_mean_syn',
                          f'{net_key}{corr_type}_mean_ctr',
                          f'{net_key}{corr_type}_t-val']):
        norm_ = norm_ = colors.TwoSlopeNorm(vmin=-max(abs(X)), vmax=max(abs(X)), vcenter=0)
        kwargs = {"norm": norm_} if X_name == f'{net_key}{corr_type}_t-val' else {}

        fig, ax, scatter, cbar = plot_3d_local_metric(X, xyz, n_name, return_scatter=True, **kwargs)
        plot_name = f'{X_name}.png'
        plt.savefig(os.getcwd() + '/plots/' + plot_name, transparent=True)
        plt.show()

        # save .mat to plot 3D brain with matlab
        cmap = scatter.get_cmap()
        if X_name == f'{net_key}{corr_type}_t-val':
            norm = plt.Normalize(vmin=-abs(X).max(), vmax=abs(X).max())
        else:
            norm = plt.Normalize(vmin=X.min(), vmax=X.max())
        rgb_values = cmap(norm(X))
        Xvalues = {'Xnet': X,
                   'xyz': xyz,
                   'color': rgb_values,
                   'names': np.array(n_name)[np.nonzero(X_t_val)],  # to mention only the significant nodes
                   'names_idx': np.nonzero(X_t_val)
                   }
        nodes_file = os.path.join(os.getcwd(), 'plots', 'glb', f'{X_name}.mat')
        sio.savemat(nodes_file, Xvalues)

        # fo each hemisphere individually
        for ind, side in zip([lh_ind, rh_ind], ('lh', 'rh')):
            Xvalues = {'Xnet': X[ind],
                       'xyz': xyz[ind],
                       'color': rgb_values[ind],
                       'names': np.array(n_name)[np.nonzero(X_t_val[ind])[0]+ind[0]],  # to mention only the significant nodes
                       'names_idx': np.nonzero(X_t_val[ind])[0]
                       }
            nodes_file = os.path.join(os.getcwd(), 'plots', 'glb', f'{X_name}_{side}.mat')
            sio.savemat(nodes_file, Xvalues)

        #TODO save colorbar
        # cmap = cbar.cmap
        # norm = cbar.norm
        # num_colors = 256  # Number of colors in the colorbar
        # colorbar_rgb = [cmap(norm(i))[:3] for i in range(num_colors)]
        # np.savetxt(os.getcwd() + '/plots/glb/colorbar.txt', colorbar_rgb, fmt='%f', delimiter='\t')

        # Plot the colorbar in Python
        # fig, ax = plt.subplots()
        # cbar = plt.colorbar(sm, ax=ax)
        # plt.show()

#%% Plot ROIs from the literature
functional_activation = np.array([32, 54, 73, 78])-1
literature = np.array([20, 32, 36, 54, 61, 62, 67, 69, 71, 73])-1
selection = literature
X_ = np.zeros(len(n_name_all))
X_[:] = 0
X_[selection] = 1
norm = plt.Normalize(vmin=-abs(X_).max(), vmax=abs(X_).max())
rgb_values = cmap(norm(X_))
Xvalues = {'Xnet': X_,
           'xyz': xyz_all,
           'color': rgb_values,
           'names': np.array(n_name_all)[selection],  # to mention only the significant nodes
           'names_idx': np.nonzero(X_)
           }
nodes_file = os.path.join(os.getcwd(), 'plots', 'glb', 'literature.mat')
# sio.savemat(nodes_file, Xvalues)