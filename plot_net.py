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

# load net metrics all subjects
Xnet = []
for sub in np.arange(num_sub)+1:
    net_file = net_path + '_'.join(('strength', 'Subject' + str(sub).zfill(3)))
    strength = io.loadmat(net_file)
    Xnet.append(strength['strength'][0])

# Xnet = abs(np.array(Xnet))  # take absolute value of connectivity matrices (undirected networks)
Xnet = np.array(Xnet)

# first 17 subjects correspond synesthetic subjects
Xnet_syn = Xnet[:17, :]
Xnet_ctr = Xnet[17:, :]
Xnet_mean = Xnet_ctr.mean(axis=0)  # mean across subjects

#%% t-test
n_t_val = []
n_p_val = []
for n_idx in np.arange(len(n_name)):
    t_val, p_val = stats.ttest_ind(Xnet_syn[:, n_idx], Xnet_ctr[:, n_idx])
    n_t_val = np.append(n_t_val, t_val)
    n_p_val = np.append(n_p_val, p_val)

#%% print Hub list
pd.set_option("display.precision", 2)
df = pd.DataFrame()
df['node'] = n_name
df['strength_synes'] = Xnet_syn.mean(axis=0)
df['strength_ctr'] = Xnet_ctr.mean(axis=0)
df['t-val'] = n_t_val
df['p-val'] = n_p_val

print('-'*100)
print(df[['node', 'strength_synes']].sort_values(by='strength_synes', key=abs).tail(15).to_string())
print('Average synes strength: {0}'.format(Xnet_syn.mean(axis=0).mean()))

print('-'*100)
print(df[['node', 'strength_ctr']].sort_values(by='strength_ctr', key=abs).tail(15).to_string())
print('Average ctr strength: {0}'.format(Xnet_ctr.mean(axis=0).mean()))

print('-'*100)
print(df.loc[df['p-val'] < 0.05, ["node", "p-val"]])


#%% Plot 3D nodes
fig = plt.figure(figsize=(12, 10), dpi=400)
gs = gridspec.GridSpec(1, 20, wspace=0.5)
ax1 = fig.add_subplot(gs[:, 0:19], projection='3d')
ax2 = fig.add_subplot(gs[:, 19])

# X = Xnet_syn.mean(axis=0)
X = np.where((df['p-val'] > 0.05), 0, df['t-val'])  # t-val higher than 0.05
ax_1, ax_2 = nodes_3D(X, xyz, n_name, ax1, ax2)

plt.tight_layout(pad=0.5)
plot_name = 'strength_t-val.png'
plt.savefig(os.getcwd() + '/plots/' + plot_name, transparent=True)
plt.show()