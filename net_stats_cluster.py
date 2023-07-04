import numpy as np
from scipy import io
import os.path
import glob
from os import listdir
import matplotlib.pyplot as plt
from scipy import signal, stats, fft, sparse
# from Topomap_separate import *
from statsmodels.stats.anova import AnovaRM
# import matlab.engine
from sklearn.cluster import SpectralClustering
from sklearn.metrics import pairwise_distances
from collections import Counter
# import community
from scipy.stats import wilcoxon
import networkx as nx
from scipy.spatial.distance import cdist
import mne


# Calculate t-statistic using Wilcoxon signed-rank test
def perm_t_stat(X):
    mean = X.mean(0)
    t_values, pv = stats.wilcoxon(X)
    return t_values


# nodes positions
path = '/Users/juliana.gonzalez/ownCloud/graph_analysis/'
path = '/Users/juliana.gonzalez/ownCloud/graph_analysis/'
net_path = path + 'net_metrics/'
num_sub = len(os.listdir(path + 'symetrical_corr_mat'))

# nodes names
n_name = []
n_name_full = []
with open(path + "BN_Atlas_246_LUT_reoriented.txt", "r") as filestream:
    for line in filestream:
        n_name.append(line.split(",")[0])
        n_name_full.append(line.split(",")[1].strip())

results_file = path + 'resultsROI_Subject001_Condition001.mat'
results = io.loadmat(results_file)
xyz = results['xyz'][0]
pos = np.array([xyz[i][0][:] for i in np.arange(np.size(xyz))])
pos = pos[:-1]

# Calculate pairwise distances between node positions
distances = cdist(pos, pos)
threshold = 20
# threshold distance in meters
adjacency = sparse.csr_matrix((distances <= threshold).astype(int))
# min = 2.3  # np.min(distances[np.nonzero(distances)])
# max = 166.6  # np.max(distances)

n_observations = 17  # number of subjects
pval = 0.05  # significance level
df = n_observations - 1  # degrees of freedom for the test
thresh = stats.t.ppf(1 - pval / 2, df)  # threshold for t-values

# Load net metrics for all subjects
net_key = 'coreness_norm'
Xnet = np.array([io.loadmat(net_path + f'net_metrics_Subject{str(sub).zfill(3)}{"_thr"}')[net_key][0]
                 for sub in range(1, num_sub + 1)])

# Split synesthetic and control subjects
Xnet_syn = Xnet[:17, :]
Xnet_ctr = Xnet[17:, :]
X = [Xnet_syn, Xnet_ctr]

# Perform permutation-based cluster test
T_obs_2, clusters_2, p_values_2, _ = mne.stats.permutation_cluster_1samp_test(Xnet_syn, threshold=thresh, n_permutations=1024,
                                                                        adjacency=adjacency, tail=1, n_jobs=-1)

print(clusters_2)  # Print the clusters found
print(p_values_2)  # Print the p-values associated with the clusters

T = np.zeros(np.shape(Xnet))
List_Clusters_relevant = []

# print clusters
# clusters = []
# for i in range(len(clusters_2)):
#     cluster = [n_name_full[list(n[0])[0]] for n in i]
#     print(cluster)
#     clusters = clusters.append([n_name_full[list(n[0])[0]] for n in clusters_2])
# print(clusters)

# Identify clusters with p-values below the threshold (0.001)
for i in range(len(p_values_2)):
    if p_values_2[i] < 0.05:
        List_Clusters_relevant.append(clusters_2[i])

print(List_Clusters_relevant)  # Print the relevant clusters

# Mark the electrodes belonging to relevant clusters in the T array
for k in range(len(List_Clusters_relevant)):
    for cluster in range(len(List_Clusters_relevant[k])):
        print(np.array(n_name_full)[List_Clusters_relevant[k][cluster]])
        T[List_Clusters_relevant[k][cluster]] = 1

T_obs_2 = T * T_obs_2  # Apply the T array to the T_obs_2 values
