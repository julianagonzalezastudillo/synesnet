#%% PLOT
import os.path
import numpy as np
from scipy import io
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
from os import path
import sys
sys.path.append(path.abspath('../netviz'))


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

# transform metric to dataframe
df = pd.DataFrame()
sub_list = np.arange(num_sub)+1
df['sub'] = sub_list
df['sub_type'] = ['synes']*17 + ['ctr']*17
for net_key in ['global_efficiency']:
    # load net metrics all subjects
    metric = []
    for sub in sub_list:
        net_file = net_path + '_'.join(('net_metrics', 'Subject' + str(sub).zfill(3)))
        Xnet = io.loadmat(net_file)
        metric.append(Xnet[net_key][0])

    metric = np.array(metric)
    df[net_key] = metric

#%% plot
fig, ax = plt.subplots()
ax = sns.catplot(data=df, x="sub_type", y="global_efficiency", kind="box", plot_kws=dict(alpha=0.3))
ax = sns.swarmplot(data=df, x="sub_type", y="global_efficiency", color=".25")
plt.show()
