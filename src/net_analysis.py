"""
=================================
            SYNESNET
=================================
This module is designed to extract network properties from FC matrices.
Diagonal is put to 0 and only positive FCC are retained.
"""

import os.path
import scipy.io as sio
import numpy as np
import networkx as nx
import bct
from net.core import coreness
from tools import load_fc
from config import NET_DIR


# Load connectivity matrices for all subjects
Xfc_syn, Xfc_ctr = load_fc()
Xfc = np.concatenate((Xfc_syn, Xfc_ctr), axis=0)

for sub in np.arange(np.shape(Xfc)[0]):
    net_file = NET_DIR / "_".join(
        ("net_metrics", "Subject{0}".format(str(sub + 1).zfill(3)), "thr.mat")
    )
    X = Xfc[sub]

    # network metrics
    G = nx.from_numpy_array(X)  # to nx format
    strength = np.array([v for k, v in G.degree(weight="weight")])  # strength
    geff = bct.efficiency_wei(X)
    leff = bct.efficiency_wei(X, local=True)
    C, isCore = coreness(X)  # coreness

    net = {
        "strength": strength,
        "global_efficiency": geff,
        "local_efficiency": leff,
        "coreness": C,
        "isCore": isCore,
    }

    # save
    if os.path.exists(net_file):
        user_input = input(
            f"The file {net_file} already exists. Do you want to continue? (y/n): "
        )
        if user_input.lower() != "y":
            print("Operation aborted.")
        else:
            Xnet.update(net)
    else:
        Xnet = net

    ###### sio.savemat(net_file, Xnet)
