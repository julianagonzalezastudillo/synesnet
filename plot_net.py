#%% PLOT
import os.path
import numpy as np
from scipy import io
from scipy import stats
import pandas as pd
import statsmodels.stats.multitest as smt
import scipy.io as sio

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
metric_list = ['strength']
for net_key in metric_list:
    # load net metrics all subjects
    Xnet = []
    for sub in np.arange(num_sub)+1:
        net_file = net_path + '_'.join(('net_metrics', 'Subject' + str(sub).zfill(3) + '{0}'.format(corr_type)))
        net_values = io.loadmat(net_file)
        Xnet.append(net_values[net_key][0])

    Xnet = np.array(Xnet)

    # first 17 subjects correspond synesthetic subjects
    Xnet_syn = Xnet[:17, :]
    Xnet_ctr = Xnet[17:, :]

    # t-test
    n_t_val = []
    n_p_val = []
    for n_idx in np.arange(len(n_name)):
        t_val, p_val = stats.ttest_ind(Xnet_syn[:, n_idx], Xnet_ctr[:, n_idx])
        n_t_val = np.append(n_t_val, t_val)
        n_p_val = np.append(n_p_val, p_val)

    # Bonferroni correction
    alpha = 0.05
    n_tests = len(n_name)
    # adjusted_alpha = alpha / n_tests
    # reject, p_corrected, _, _ = smt.multipletests(n_p_val, alpha=adjusted_alpha, method='bonferroni')
    # for i in range(len(n_name)):
    #     if reject[i]:
    #         print(f"Electrode {n_name[i]}: p-value = {p_corrected[i]:.4f} (rejected)")
    #     else:
    #         print(f"Electrode {n_name[i]}: p-value = {p_corrected[i]:.4f}")
    #
    rejected, corrected_p_values, _, _ = smt.multipletests(n_p_val, method='fdr_bh')
    # Print the rejected hypotheses and their corresponding corrected p-values
    for i, p_value, corrected_p_value in zip(range(len(n_name)), n_p_val, corrected_p_values):
        if rejected[i]:
            print(
                f"Electrode {n_name[i]}: Reject null hypothesis (p-value: {p_value}, corrected p-value: {corrected_p_value})")


    # print Hub list
    # Mean across subjects
    Xnet_syn_mean = Xnet_syn.mean(axis=0)
    Xnet_ctr_mean = Xnet_ctr.mean(axis=0)

    pd.set_option("display.precision", 2)
    df = pd.DataFrame()
    df['node'] = n_name
    df['node_complete'] = n_name_full
    df['{0}_synes'.format(net_key)] = Xnet_syn_mean
    df['{0}_ctr'.format(net_key)] = Xnet_ctr_mean
    df['t-val'] = n_t_val
    df['p-val'] = n_p_val

    print('-'*100)
    print(df[['node', '{0}_synes'.format(net_key)]].sort_values(by='{0}_synes'.format(net_key), key=abs).tail(30).to_string())
    print('Average synes {0}: {1}'.format(net_key, Xnet_syn.mean(axis=0).mean()))

    print('-'*100)
    print(df[['node', '{0}_ctr'.format(net_key)]].sort_values(by='{0}_ctr'.format(net_key), key=abs).tail(30).to_string())
    print('Average ctr {0}: {1}'.format(net_key, Xnet_ctr.mean(axis=0).mean()))

    print('-'*100)
    print(df.loc[df['p-val'] < 0.05, ["node", "node_complete", "p-val", "t-val"]].sort_values(by='p-val'))

    # Plot 3D nodes
    # normalize by highest values
    # X_max = np.max([Xnet_syn_mean, Xnet_ctr_mean])
    # Xnet_syn_mean_norm = Xnet_syn_mean/X_max
    # Xnet_ctr_mean_norm = Xnet_ctr_mean/X_max
    # X_max_norm = np.max([Xnet_syn_mean_norm, Xnet_ctr_mean_norm])
    # cmap = colors.LinearSegmentedColormap.from_list('YlBu',
    # ['#FCEEA5', '#EEF8DF', '#BDE2EE', '#80B7D6', '#4A79B7', '#313595'])

    for X, X_name in zip([Xnet_syn_mean, Xnet_ctr_mean],
                         ['{0}{1}_mean_syn'.format(net_key, corr_type),
                          '{0}{1}_mean_ctr'.format(net_key, corr_type)]):
        fig, ax, scatter, cbar = plot_3d_local_metric(X, xyz, n_name, return_scatter=True)
        plot_name = '{0}.png'.format(X_name)
        plt.savefig(os.getcwd() + '/plots/' + plot_name, transparent=True)
        plt.show()

        # save .mat to plot 3D brain with matlab
        cmap = scatter.get_cmap()
        norm = plt.Normalize(vmin=X.min(), vmax=X.max())
        rgb_values = cmap(norm(X))
        Xvalues = {'Xnet': X,
                   'xyz': xyz,
                   'color': rgb_values,
                   }
        nodes_file = os.getcwd() + '/plots/glb/{0}.mat'.format(X_name)
        sio.savemat(nodes_file, Xvalues)

        for ind, side in zip([lh_ind, rh_ind], ('lh', 'rh')):
            Xvalues = {'Xnet': X[ind],
                       'xyz': xyz[ind],
                       'color': rgb_values[ind],
                       }
            nodes_file = os.getcwd() + '/plots/glb/{0}_{1}.mat'.format(X_name, side)
            sio.savemat(nodes_file, Xvalues)

    # for t-values, keep significant t-values
    X_t_val = np.where((df['p-val'] > 0.05), 0, df['t-val'])  # t-val higher than 0.05
    X_name = '{0}_t-val'.format(net_key)
    fig, ax, scatter, cbar = plot_3d_local_metric(X_t_val, xyz, n_name, return_scatter=True)
    plot_name = '{0}{1}_t-val.png'.format(net_key, corr_type)
    plt.savefig(os.getcwd() + '/plots/' + plot_name, transparent=True )
    plt.show()

    # save .mat to plot 3D brain with matlab
    import matplotlib.cm as cm
    cmap = scatter.get_cmap()
    norm = plt.Normalize(vmin=X_t_val.min(), vmax=X_t_val.max())
    rgb_values = cmap(norm(X_t_val))

    # save names
    Xnames = {'names': np.array(n_name)[np.nonzero(X_t_val)],
              'names_idx': np.nonzero(X_t_val)}
    nodes_names = os.getcwd() + '/plots/glb/{0}_names.mat'.format(net_key)
    sio.savemat(nodes_names, Xnames)

    Xvalues = {'Xnet': X_t_val,
               'xyz': xyz,
               'color': rgb_values,
               }
    nodes_file = os.getcwd() + '/plots/glb/{0}_t-val.mat'.format(net_key)
    sio.savemat(nodes_file, Xvalues)

    for ind, side in zip([lh_ind, rh_ind], ('lh', 'rh')):
        Xvalues = {'Xnet': X_t_val[ind],
                   'xyz': xyz[ind],
                   'color': rgb_values[ind],
                   }
        nodes_file = os.getcwd() + '/plots/glb/{0}_t-val_{1}.mat'.format(net_key, side)
        sio.savemat(nodes_file, Xvalues)

        # save names
        n_name_h = np.array(n_name)[ind]
        Xnames = {'names': n_name_h[np.nonzero(X_t_val[ind])],
                  'names_idx': np.nonzero(X_t_val[ind])}
        nodes_names = os.getcwd() + '/plots/glb/{0}_names_{1}.mat'.format(net_key, side)
        sio.savemat(nodes_names, Xnames)

        cmap = cbar.cmap
        norm = cbar.norm
        num_colors = 256  # Number of colors in the colorbar
        colorbar_rgb = [cmap(norm(i))[:3] for i in range(num_colors)]
        np.savetxt(os.getcwd() + '/plots/glb/colorbar.txt', colorbar_rgb, fmt='%f', delimiter='\t')

        # Plot the colorbar in Python
        fig, ax = plt.subplots()
        cbar = plt.colorbar(sm, ax=ax)
        plt.show()

