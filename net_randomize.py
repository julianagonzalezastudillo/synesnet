import os.path
import numpy as np
from scipy import io
import bct
from copy import deepcopy
import time

t = time.time()
# path = '/network/lustre/iss02/aramis/users/juliana.gonzalez/synesnet/'
path = '/Users/juliana.gonzalez/ownCloud/graph_analysis/'
num_sub = len(os.listdir(path + 'symetrical_corr_mat'))

#%% Create a directory to save the adjacency matrices
rand_path = path + 'rand_mat'
os.makedirs(rand_path, exist_ok=True)

for sub in np.arange(1, 35):
    sub_file = 'CorrMatrix_Subject{0}.mat'.format(str(sub).zfill(3))
    print(sub_file)

    # load correlation
    fc_file = path + 'symetrical_corr_mat/' + sub_file
    fc = io.loadmat(fc_file)

    # 3. corr[corr>=0] = [0, 1]  # winning matrices !!!
    Xfc_thr = deepcopy(fc['CorrMatrix'])
    Xfc_thr[Xfc_thr <= 0] = 0
    np.fill_diagonal(Xfc_thr, 0)

    # generate 100 randomised versions
    random_samples = 100
    # node_strengths = []
    networks = []
    for r in range(random_samples):
        Xfc_rand = bct.null_model_und_sign(Xfc_thr, bin_swaps=np.shape(Xfc_thr)[0], wei_freq=1)
        networks.append(Xfc_rand[0])

        # Calculate the node strength
        # G = nx.from_numpy_array(Xfc_rand[0])  # to nx format
        # strength = np.array([v for k, v in G.degree(weight='weight')])  # strength
        # node_strengths.append(list(strength))

    adjacency_matrices = np.array([network for network in networks])

    # Save the concatenated matrix to a .mat file
    filename = os.path.join(rand_path, 'RandMatrices_Subject{0}.mat'.format(str(sub).zfill(3)))
    io.savemat(filename, {'RandMatrices': adjacency_matrices})

print(time.time() - t)
# Plot the node strength distribution for each network
# import matplotlib.pyplot as plt
# plt.figure(figsize=(10, 6))
# colors = plt.cm.get_cmap('rainbow', random_samples)
# bins = 25
# for i, strengths in enumerate(node_strengths):
#     plt.hist(strengths, bins=bins, range=(10, 40), alpha=0.5, color=colors(i))
#
# # Plot the mean node strength distribution
# mean_strengths = np.mean(node_strengths, axis=0)
# plt.hist(mean_strengths, bins=bins, range=(10, 40), alpha=0.9, color='black', label='Mean random')
#
# # Plot the strength distribution of the original network
# G = nx.from_numpy_array(Xfc_thr)  # to nx format
# strength = np.array([v for k, v in G.degree(weight='weight')])  # strength
# plt.hist(strength, bins=bins, range=(10, 40), alpha=0.9, color='white', label='Original', edgecolor='black', linewidth=1.2)
#
# plt.xlabel('Node Strength')
# plt.ylabel('Frequency')
# plt.title('Node Strength Distribution of Multiple Networks')
# plt.legend()
# plt.show()
