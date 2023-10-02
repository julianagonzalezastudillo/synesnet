#%% PLOT
import os.path
import numpy as np
from scipy import io
from scipy import stats
import pandas as pd
import statsmodels.stats.multitest as smt
import scipy.io as sio
from statsmodels.stats.multitest import multipletests
import statsmodels
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap
from os import path
import sys
sys.path.append(path.abspath('../netviz'))
from viz.netlocal import plot_3d_local_metric


def load_net_metrics(net_path, net_key, corr_type, num_sub, idx_select):
    Xnet = np.array([
        io.loadmat(os.path.join(net_path, f'net_metrics_Subject{str(sub).zfill(3)}{corr_type}'))[net_key][0]
        for sub in range(1, num_sub + 1)])
    X_ = np.zeros([np.shape(Xnet)[0], len(n_name)])
    X_[:, idx_select] = Xnet[:, idx_select]
    return X_


path = '/Users/juliana.gonzalez/ownCloud/graph_analysis/'
net_path = os.path.join(path, 'net_metrics')
node_file = os.path.join(path, 'BN_Atlas_246_LUT_reoriented.txt')
results_file = os.path.join(path, 'resultsROI_Subject001_Condition001.mat')
strength_stat_file = os.path.join(os.getcwd(), 'plots', 'glb', 'strength_thr_t-val.mat')
coreness_stat_file = os.path.join(os.getcwd(), 'plots', 'glb', 'coreness_norm_by_rand_conserve_strenght_distribution_thr_t-val_selection.mat')

# CONSTANTS
P_VAL = 0.05
SELECTION = True  # True: select base on strength significance
corr_type = '_thr'
metric_list = ['coreness']
# ['coreness']
# ['coreness_norm']  # ['strength', 'coreness_norm_strength_selection']
# ['coreness_norm_by_rand_conserve_strenght_distribution']
binarize_coreness_size = True
cmap = plt.colormaps['Spectral'].reversed()

if SELECTION:
    p_val_type = 'p-val_corrected'
else:
    p_val_type = 'p-val'

# number of subjects
num_sub = len(os.listdir(os.path.join(path, 'symetrical_corr_mat')))

# nodes positions
xyz = io.loadmat(results_file)['xyz'][0][:-1] 
xyz = np.array([x[0] for x in xyz])

# nodes names
n_name, n_name_full = [], []
with open(node_file, "r") as filestream:
    for line in filestream:
        name, name_full = line.split(",")[:2]
        n_name.append(name)
        n_name_full.append(name_full.strip())

# Open strength t-val file
strength_stat = io.loadmat(strength_stat_file)
idx_select = strength_stat['names_idx'][0]
idx_select = idx_select if SELECTION else slice(None)

for net_key in metric_list:

    if SELECTION and net_key == 'coreness':
        coreness_stat = io.loadmat(coreness_stat_file)
        idx_select_coreness = coreness_stat['names_idx'][0]
        idx_select = list(set(idx_select) & set(idx_select_coreness))

    # Load net metrics for all subjects
    Xnet = load_net_metrics(net_path, net_key, corr_type, num_sub, idx_select)

    # Create an empty matrix of subjects and nodes and replace the significant nodes with its values
    X_ = np.zeros([np.shape(Xnet)[0], len(n_name)])  # 34x246
    X_[:, idx_select] = Xnet[:, idx_select]
    Xnet = X_
    lh_ind = [index for index, element in enumerate(np.array(n_name)) if element.endswith('_L')]
    rh_ind = [index for index, element in enumerate(np.array(n_name)) if element.endswith('_R')]

    # Split synesthetic and control subjects
    Xnet[np.isnan(Xnet) | np.isinf(Xnet)] = 0
    Xnet_syn = Xnet[:17, :]
    Xnet_ctr = Xnet[17:, :]

    # Perform t-test
    t_val, p_val = stats.ttest_ind(Xnet_syn, Xnet_ctr)
    t_val[np.isnan(t_val)] = 0  # Replace NaNs with 0

    # Perform FDR correction
    p_val = p_val[idx_select]
    rejected, corrected_p_values, _, _ = multipletests(p_val, method='fdr_bh')
    statsmodels.stats.multitest.fdrcorrection(p_val, alpha=P_VAL, method='indep', is_sorted=False)

    # Corrected by hand FDR (Benjamini-Hochberg)
    rank = np.arange(len(p_val))+1
    rank_idx = np.argsort(p_val)
    N = len(p_val)
    p_val_corrected_ = np.sort(p_val)*N/rank
    p_val_corrected = np.zeros(len(p_val))
    p_val_corrected[rank_idx] = p_val_corrected_

    # Mean across subjects
    Xnet_syn_mean = Xnet_syn.mean(axis=0)
    Xnet_ctr_mean = Xnet_ctr.mean(axis=0)

    # Create DataFrame for results
    pd.set_option("display.precision", 2)
    df = pd.DataFrame({'node': np.array(n_name)[idx_select],
                       'node_complete': np.array(n_name_full)[idx_select],
                       f'{net_key}_synes': Xnet_syn.mean(axis=0)[idx_select],
                       f'{net_key}_ctr': Xnet_ctr.mean(axis=0)[idx_select],
                       't-val': np.array(t_val)[idx_select],
                       'p-val': p_val,
                       'p-val_corrected': p_val_corrected
                       })

    # for t-values, keep significant t-values
    X_ = np.zeros(len(n_name))  # 246
    X_[idx_select] = np.where((df[p_val_type] > P_VAL), 0, df['t-val'])
    X_t_val = X_

    # print Hub list
    print('-' * 100)
    print(df[['node', f'{net_key}_synes']].sort_values(by=f'{net_key}_synes', key=abs).tail(30).to_string())
    print(f'Average synes {net_key}: {Xnet_syn.mean(axis=0).mean()}')

    print('-' * 100)
    print(df[['node', f'{net_key}_ctr']].sort_values(by=f'{net_key}_ctr', key=abs).tail(30).to_string())
    print(f'Average ctr {net_key}: {Xnet_ctr.mean(axis=0).mean()}')

    print('-' * 100)
    if SELECTION:
        print(df.sort_values(by='p-val_corrected'))
    else:
        print(df.loc[df['p-val'] < P_VAL, ["node", "node_complete", f'{net_key}_synes', f'{net_key}_ctr', "p-val", "p-val_corrected", "t-val"]]
              .sort_values(by='p-val'))

    for X, X_name in zip([Xnet_syn_mean, Xnet_ctr_mean, X_t_val],
                         [f'{net_key}{corr_type}_mean_syn',
                          f'{net_key}{corr_type}_mean_ctr',
                          f'{net_key}{corr_type}_t-val']):

        # Colormap for t-val
        if X_name == f'{net_key}{corr_type}_t-val':
            cmap = plt.colormaps['Spectral'].reversed()
            norm = colors.TwoSlopeNorm(vmin=-max(abs(X)), vmax=max(abs(X)), vcenter=0)
            X_max = None
        elif X_name == 'coreness_thr_mean_syn' or X_name == 'coreness_thr_mean_ctr':
            # cmap_colors = ['#203FB6', '#008AFB', '#B3D4FF', 'white', '#FFEC4A', '#FF6611', '#F62336']
            cmap_colors = ['#11205E', '#203FB6', '#86CAFF', 'white', '#FFEC4A', '#F62336', '#80121B']
            positions = np.linspace(0, 1, len(cmap_colors))
            cmap = LinearSegmentedColormap.from_list('custom_colormap', list(zip(positions, cmap_colors)))
            norm = plt.Normalize(vmin=0, vmax=1)
            X_max = np.max([Xnet_syn_mean, Xnet_ctr_mean])
        else:
            norm = plt.Normalize(vmin=X.min(), vmax=X.max())
            X_max = None
        kwargs = {"norm": norm}

        # Plot and generate scatter info to plot in matlab
        if X_max is None:
            X_size = abs(pow(X, 2) / max(abs(pow(X, 2)))) * 80
        else:
            X_size = abs(pow(X, 2) / pow(X_max, 2)) * 80
        fig, ax, scatter, cbar = plot_3d_local_metric(X_size, X, xyz, n_name, cmap=cmap,
                                                      return_scatter=True, **kwargs)
        # plt.savefig(os.path.join(os.getcwd(), 'plots', f'{X_name}.png'), transparent=True)
        plt.show()

        cmap = scatter.get_cmap()
        rgb_values = cmap(norm(X))

        # save .mat to plot 3D brain with matlab
        Xvalues = {'Xnet': X,
                   'xyz': xyz,
                   'color': rgb_values,
                   'names': np.array(n_name)[np.nonzero(X_t_val)],  # to mention only the significant nodes
                   'names_idx': np.nonzero(X_t_val)
                   }
        X_name += '_selection' if SELECTION else ''
        nodes_file = os.path.join(os.getcwd(), 'plots', 'glb', f'{X_name}.mat')
        sio.savemat(nodes_file, Xvalues)

        # for each hemisphere individually
        for ind, side in zip([lh_ind, rh_ind], ('lh', 'rh')):
            Xvalues = {'Xnet': X[ind],
                       'xyz': xyz[ind],
                       'color': rgb_values[ind],
                       'names': np.array(n_name)[np.nonzero(X_t_val[ind])[0]+ind[0]],  # to mention only the significant nodes
                       'names_idx': np.nonzero(X_t_val[ind])[0]
                       }

            nodes_file = os.path.join(os.getcwd(), 'plots', 'glb', f'{X_name}_{side}.mat')
            sio.savemat(nodes_file, Xvalues)
            print(nodes_file)

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
X_ = np.zeros(len(n_name))
X_[:] = 0
X_[selection] = 1
norm = plt.Normalize(vmin=-abs(X_).max(), vmax=abs(X_).max())
rgb_values = cmap(norm(X_))
Xvalues = {'Xnet': X_,
           'xyz': xyz,
           'color': rgb_values,
           'names': np.array(n_name)[selection],  # to mention only the significant nodes
           'names_idx': np.nonzero(X_)
           }
# nodes_file = os.path.join(os.getcwd(), 'plots', 'glb', 'literature.mat')
# sio.savemat(nodes_file, Xvalues)
