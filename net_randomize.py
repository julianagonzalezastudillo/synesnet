"""
=================================
            SYNESNET
=================================
Generate randomized network conserving strength distribution, to normalize coreness in net_core_random.py
"""

import numpy as np
from scipy import io
import bct
from tools import load_fc
import time
import networkx as nx
from config import RAND_DIR


# Load connectivity matrices for all subjects
Xfc_syn, Xfc_ctr = load_fc()
Xfc = np.concatenate((Xfc_syn, Xfc_ctr), axis=0)

# generate 100 randomised versions
random_samples = 100
for sub in np.arange(np.shape(Xfc)[0]):
    X = Xfc[sub]

    # node_strengths = []
    networks = []
    for r in range(random_samples):
        # Randomize conserving strength distribution
        Xfc_rand = bct.null_model_und_sign(X, bin_swaps=np.shape(X)[0], wei_freq=1)

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
    filename = RAND_DIR / f"RandMatrices_Subject{str(sub).zfill(3)}.mat"
    # io.savemat(filename, {"RandMatrices": adjacency_matrices})

    # print(time.time() - t)

    # %% Plot the node strength distribution for each network
    # import matplotlib.pyplot as plt
    #
    # plt.figure(figsize=(10, 6))
    # colors = plt.cm.get_cmap("rainbow", random_samples)
    # bins = 25
    #
    # # Plot the strength distribution of each random matrix
    # for i, strengths in enumerate(node_strengths):
    #     plt.hist(strengths, bins=bins, range=(10, 40), alpha=0.5, color=colors(i))
    #
    # # Plot the strength distribution of the original network
    # G = nx.from_numpy_array(X)  # to nx format
    # strength = np.array([v for k, v in G.degree(weight="weight")])  # strength
    # plt.hist(
    #     strength,
    #     bins=bins,
    #     range=(10, 40),
    #     alpha=0.9,
    #     color="white",
    #     label="Original",
    #     edgecolor="black",
    #     linewidth=1.2,
    # )
    #
    # plt.xlabel("Node Strength")
    # plt.ylabel("Frequency")
    # plt.title(
    #     "SUBJECT {0}: Node Strength Distribution of Multiple Networks".format(sub)
    # )
    # plt.legend()
    # plt.show()
