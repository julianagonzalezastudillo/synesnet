"""
=================================
            SYNESNET
=================================
This script is designed to normalized coreness by randomized networks.
"""

import os.path
import scipy.io as sio
import numpy as np
from scipy import io
from net.core import coreness
from copy import deepcopy
import time
from tools import load_net_metrics
from config import NET_DIR, RAND_DIR


# Load original coreness computation from file generated in net_analysis.py
Xnet_syn, Xnet_ctr = load_net_metrics("coreness")
C_orig = np.concatenate((Xnet_syn, Xnet_ctr), axis=0)

for sub in np.arange(np.shape(C_orig)[0]):
    print(f"Subject: {sub + 1}")
    t = time.time()

    # Open random matrices (the randomized matrices are computed on Xfc>0)
    fc_rand_file = RAND_DIR / f"RandMatrices_Subject{str(sub + 1).zfill(3)}.mat"
    fc_rand = io.loadmat(fc_rand_file)
    Xfc_rand = deepcopy(fc_rand["RandMatrices"])

    C_rand = []  # random networks coreness
    for m in np.arange(np.shape(Xfc_rand)[0]):
        print("random network n#{0}".format(m))
        X = Xfc_rand[m]
        C, isCore = coreness(X)  # coreness
        C_rand.append(C)

    C_rand = np.array(C_rand)
    C_rand_mean = np.mean(C_rand, axis=0)

    C_norm = C_orig[sub] / C_rand_mean  # coreness ratio
    C_zscore = (C_orig - C_rand_mean) / np.std(C_rand, axis=0)  # coreness z-score

    # Replace NaN and infinite values with zeros
    C_norm[np.isnan(C_norm) | np.isinf(C_norm)] = 0
    C_zscore[np.isnan(C_zscore) | np.isinf(C_zscore)] = 0

    print(time.time() - t)

    # save
    net_file = NET_DIR / f"net_metrics_Subject{str(sub + 1).zfill(3)}_thr.mat"
    if os.path.exists(net_file):
        Xnet = sio.loadmat(net_file)
        Xnet.update(
            {
                "coreness_norm_by_rand": np.nan_to_num(C_norm),
                "coreness_zcore_by_rand": np.nan_to_num(C_zscore),
                "coreness_rand": np.nan_to_num(C_rand),
            }
        )
    # sio.savemat(net_file, Xnet)
