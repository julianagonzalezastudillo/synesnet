import os.path
import numpy as np
from scipy import io
import bct
from copy import deepcopy
import time
from net.randomize import rand_wu
import networkx as nx

# path = '/network/lustre/iss02/aramis/users/juliana.gonzalez/synesnet/'
path = '/Users/juliana.gonzalez/ownCloud/graph_analysis/'
num_sub = len(os.listdir(path + 'symetrical_corr_mat'))

strength_select = True
# Create a directory to save the adjacency matrices
rand_path = path + 'rand_mat'
os.makedirs(rand_path, exist_ok=True)

# Open strength t-val file
strength_stat_file = os.path.join(os.getcwd(), 'plots', 'glb', 'strength_thr_t-val.mat')
strength_stat = io.loadmat(strength_stat_file)
idx_select = strength_stat['names_idx'][0]

for sub in np.arange(1, 35):
    # t = time.time()
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
        # Xfc_rand = bct.null_model_und_sign(Xfc_thr, bin_swaps=np.shape(Xfc_thr)[0], wei_freq=1)

        if strength_select:
            Xfc_rand = rand_wu(Xfc_thr[idx_select][:, idx_select])
            filename = os.path.join(rand_path, 'RandMatrices_Subject{0}_strength_selection.mat'.format(str(sub).zfill(3)))
        else:
            Xfc_rand = rand_wu(Xfc_thr)
            filename = os.path.join(rand_path, 'RandMatrices_Subject{0}.mat'.format(str(sub).zfill(3)))

        # Check if the new rand network already exists
        if any(np.array_equal(Xfc_rand, network) for network in networks):
            r = r - 1
        else:
            networks.append(Xfc_rand)

        # Calculate the node strength
        # G = nx.from_numpy_array(Xfc_rand)  # to nx format
        # strength = np.array([v for k, v in G.degree(weight='weight')])  # strength
        # node_strengths.append(list(strength))

    adjacency_matrices = np.array([network for network in networks])

    # Save the concatenated matrix to a .mat file
    io.savemat(filename, {'RandMatrices': adjacency_matrices})

    # print(time.time() - t)

    #%% Plot the node strength distribution for each network
    # import matplotlib.pyplot as plt
    # plt.figure(figsize=(10, 6))
    # colors = plt.cm.get_cmap('rainbow', random_samples)
    # bins = 25
    # for i, strengths in enumerate(node_strengths):
    #     plt.hist(strengths, bins=bins, range=(10, 40), alpha=0.5, color=colors(i))
    #
    # # Plot the mean node strength distribution
    # mean_strengths = np.mean(node_strengths, axis=0)
    # # plt.hist(mean_strengths, bins=bins, range=(10, 40), alpha=0.9, color='black', label='Mean random')
    #
    # # Plot the strength distribution of the original network
    # G = nx.from_numpy_array(Xfc_thr)  # to nx format
    # strength = np.array([v for k, v in G.degree(weight='weight')])  # strength
    # plt.hist(strength, bins=bins, range=(10, 40), alpha=0.9, color='white', label='Original', edgecolor='black', linewidth=1.2)
    #
    # plt.xlabel('Node Strength')
    # plt.ylabel('Frequency')
    # plt.title('SUBJECT {0}: Node Strength Distribution of Multiple Networks'.format(sub))
    # plt.legend()
    # plt.show()
