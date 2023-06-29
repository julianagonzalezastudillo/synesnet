#%% PLOT
import os.path
import numpy as np
from scipy import io
from scipy import stats
import pandas as pd
import statsmodels.stats.multitest as smt
import scipy.io as sio
import math

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
xyz = results['xyz'][0]
xyz = np.array([xyz[i][0][:] for i in np.arange(np.size(xyz))])

# nodes names
n_name = []
n_name_full = []
with open(path + "BN_Atlas_246_LUT_reoriented.txt", "r") as filestream:
    for line in filestream:
        n_name.append(line.split(",")[0])
        n_name_full.append(line.split(",")[1].strip())

lh_ind = [index for index, element in enumerate(n_name) if element.endswith('_L')]
rh_ind = [index for index, element in enumerate(n_name) if element.endswith('_R')]

corr_type = '_thr'
metric_list = ['strength', 'coreness_norm']
for net_key in metric_list:
    # Load net metrics for all subjects
    Xnet = np.array([io.loadmat(net_path + f'net_metrics_Subject{str(sub).zfill(3)}{corr_type}')[net_key][0]
                     for sub in range(1, num_sub + 1)])

    # Split synesthetic and control subjects
    Xnet_syn = Xnet[:17, :]
    Xnet_ctr = Xnet[17:, :]

    # Perform t-test
    t_val, p_val = stats.ttest_ind(Xnet_syn, Xnet_ctr)
    t_val = [0 if math.isnan(x) else x for x in t_val]  # put to 0 the t-val=nan

    # Mean across subjects
    Xnet_syn_mean = Xnet_syn.mean(axis=0)
    Xnet_ctr_mean = Xnet_ctr.mean(axis=0)

    pd.set_option("display.precision", 2)
    # Create DataFrame for results
    df = pd.DataFrame({'node': n_name,
                       'node_complete': n_name_full,
                       f'{net_key}_synes': Xnet_syn.mean(axis=0),
                       f'{net_key}_ctr': Xnet_ctr.mean(axis=0),
                       't-val': t_val,
                       'p-val': p_val})

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
    print(df.loc[df['p-val'] < 0.05, ["node", "node_complete", "p-val", "t-val"]].sort_values(by='p-val'))

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

