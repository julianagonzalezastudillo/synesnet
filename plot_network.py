#%% plot core-periphery
import bct
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from net.core import coreness
import os.path
from scipy import io
from tools import load_fc, load_xyz, load_node_names


path = '/Users/juliana.gonzalez/ownCloud/graph_analysis/'

# nodes positions
xyz = load_xyz()

# node names
n_name, full_n_name = load_node_names()

# Load connectivity matrices
Xfc_syn, Xfc_ctr = load_fc()

# Reshape all matrices
Xfc_syn_mean = Xfc_syn.mean(axis=0)
Xfc_ctr_mean = Xfc_ctr.mean(axis=0)

thr_i = 11
for matrix in [Xfc_syn_mean, Xfc_ctr_mean]:

    # Get core & periphery
    C, isCore = coreness(matrix)

    # Apply selected threshold
    N = np.shape(matrix)[0]
    thr = np.arange(1, N) / (N - 1)
    X = bct.threshold_proportional(matrix, thr[thr_i])
    G = nx.Graph(X)

    # Set color for core and periphery nodes
    core_periph = isCore[:, thr_i]
    node_colors = ['b' if core_periph[i] == 0 else 'r' for i in range(N)]

    # Extract node positions
    node_positions = {node: xyz[node] for node in np.arange(len(G))}

    # Draw nodes in 3D
    fig = plt.figure(figsize=(8, 8), dpi=300)
    ax = fig.add_subplot(111, projection='3d')
    for node, pos in node_positions.items():
        x, y, z = pos
        ax.scatter(x, y, z, label=node, s=1, c=node_colors[node])

    # Draw edges in 3D
    for edge in G.edges():
        start, end = edge
        x1, y1, z1 = node_positions[start]
        x2, y2, z2 = node_positions[end]

        edge_weight = X[start, end]
        edge_width = ((edge_weight * 10) ** 2) / 100
        ax.plot([x1, x2], [y1, y2], [z1, z2], 'k-', linewidth=edge_width, alpha=edge_width)

    # Set axis labels
    ax.view_init(elev=90, azim=-90)
    ax.grid(False)
    ax.axis('off')
    ax.axis('equal')
    plt.show()
