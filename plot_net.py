#%% PLOT
import os.path
import numpy as np
from scipy import io
from scipy import stats
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os import path
import sys
sys.path.append(path.abspath('../netviz'))
from net3D import nodes_3D


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
# n_name_full = []
with open(path + "BN_Atlas_246_LUT_reoriented.txt", "r") as filestream:
    for line in filestream:
        n_name.append(line.split(",")[0])
        # n_name_full.append(line.split(",")[1].strip())

for net_key in ['strength', 'laterality']:
    # for net_key in ['strength']:
    # load net metrics all subjects
    Xnet = []
    for sub in np.arange(num_sub)+1:
        net_file = net_path + '_'.join((net_key, 'Subject' + str(sub).zfill(3)))
        net_values = io.loadmat(net_file)
        Xnet.append(net_values[net_key][0])

    # Xnet = abs(np.array(Xnet))  # take absolute value of connectivity matrices (undirected networks)
    Xnet = np.array(Xnet)

    # first 17 subjects correspond synesthetic subjects
    Xnet_syn = Xnet[:17, :]
    Xnet_ctr = Xnet[17:, :]

    #%% t-test
    n_t_val = []
    n_p_val = []
    for n_idx in np.arange(len(n_name)):
        t_val, p_val = stats.ttest_ind(Xnet_syn[:, n_idx], Xnet_ctr[:, n_idx])
        n_t_val = np.append(n_t_val, t_val)
        n_p_val = np.append(n_p_val, p_val)

    #%% print Hub list
    # Mean across subjects
    Xnet_syn_mean = Xnet_syn.mean(axis=0)
    Xnet_ctr_mean = Xnet_ctr.mean(axis=0)

    pd.set_option("display.precision", 2)
    df = pd.DataFrame()
    df['node'] = n_name
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
    print(df.loc[df['p-val'] < 0.05, ["node", "p-val"]].sort_values(by='p-val'))

    #%% Plot 3D nodes
    # normalize by highest values
    X_max = np.max([Xnet_syn_mean, Xnet_ctr_mean])
    # Xnet_syn_mean_norm = Xnet_syn_mean/X_max
    # Xnet_ctr_mean_norm = Xnet_ctr_mean/X_max
    # X_max_norm = np.max([Xnet_syn_mean_norm, Xnet_ctr_mean_norm])

    for X, X_name in zip([Xnet_syn_mean, Xnet_ctr_mean],
                         ['{0}_mean_syn'.format(net_key),
                          '{0}_mean_ctr'.format(net_key)]):
        fig = plt.figure(figsize=(12, 10), dpi=400)
        gs = gridspec.GridSpec(1, 20, wspace=0.5)
        ax1 = fig.add_subplot(gs[:, 0:19], projection='3d', title=X_name)
        ax2 = fig.add_subplot(gs[:, 19])

        ax_1, ax_2 = nodes_3D(X, xyz, n_name, ax1, ax2, nodes_val_max=X_max)

        plot_name = '{0}_abs.png'.format(X_name)
        plt.savefig(os.getcwd() + '/plots/' + plot_name, transparent=True)
        plt.show()

    #%% for t-values
    # significant t-values
    X_t_val = np.where((df['p-val'] > 0.05), 0, df['t-val'])  # t-val higher than 0.05

    X_name = '{0}_t-val'.format(net_key)
    fig = plt.figure(figsize=(12, 10), dpi=400)
    gs = gridspec.GridSpec(1, 20, wspace=0.5)
    ax1 = fig.add_subplot(gs[:, 0:19], projection='3d', title=X_name)
    ax2 = fig.add_subplot(gs[:, 19])

    ax_1, ax_2 = nodes_3D(X_t_val, xyz, n_name, ax1, ax2, nodes_val_max=None)

    plot_name = '{0}_abs.png'.format(X_name)
    plt.savefig(os.getcwd() + '/plots/' + plot_name, transparent=True)
    plt.show()