import os.path
import scipy.io as sio
import numpy as np
from scipy import io
from net.core import coreness
from copy import deepcopy
import time


class bcolors:
    WARNING = '\033[93m'
    ENDC = '\033[0m'


path = '/Users/juliana.gonzalez/ownCloud/graph_analysis/'
net_path = '/Users/juliana.gonzalez/ownCloud/graph_analysis/net_metrics/'
if not os.path.exists(net_path):
    os.makedirs(net_path)

n_name = []
with open(path + "BN_Atlas_246_LUT_reoriented.txt", "r") as filestream:
    for line in filestream:
        n_name.append(line.split(",")[0])

# for sub_file in os.listdir(path + 'symetrical_corr_mat'):
for sub in range(25, 35):
    sub_file = 'CorrMatrix_Subject{0}'.format(str(sub).zfill(3))
    print(sub_file)
    t = time.time()

    # load correlation
    fc_file = path + 'symetrical_corr_mat/' + sub_file
    fc = io.loadmat(fc_file)

    # connectivity matrix
    # Xfc = np.nan_to_num(fc['CorrMatrix'])
    Xfc_thr = deepcopy(fc['CorrMatrix'])
    Xfc_thr[Xfc_thr <= 0] = 0
    # net_file_thr = net_path + '_'.join(('net_metrics', 'Subject{0}'.format(str(sub).zfill(3)), 'thr.mat'))
    np.fill_diagonal(Xfc_thr, 0)

    # open random matrices
    # the randomized matrices (already) only contain Xfc>0
    fc_rand_file = path + 'rand_mat/' + 'RandMatrices_Subject{0}'.format(str(sub).zfill(3))
    fc_rand = io.loadmat(fc_rand_file)
    Xfc_rand = deepcopy(fc_rand['RandMatrices'])
    Xfc_rand[Xfc_rand <= 0] = 0
    C_rand = []  # random networks coreness
    isCore_rand = []
    for m in np.arange(np.shape(Xfc_rand)[0]):
        print('random network n#{0}'.format(m))
        X = Xfc_rand[m]
        np.fill_diagonal(X, 0)
        t = time.time()
        C, isCore = coreness(X)  # coreness
        C_rand.append(C)
        isCore_rand.append(isCore)

    isC_rand = np.array([isC for isC in isCore_rand])
    C_rand = np.array(C_rand)
    C_rand_mean = np.mean(C_rand, axis=0)
    C_orig, isCore_orig = coreness(Xfc_thr)
    C_norm = C_orig / C_rand_mean  # coreness ratio
    C_zscore = (C_orig - C_rand_mean) / np.std(C_rand, axis=0)  # coreness z-score
    print(time.time() - t)

    # save
    net_file = net_path + '_'.join(('net_metrics', sub_file.split("_")[-1], 'thr.mat'))
    if os.path.exists(net_file):
        Xnet = sio.loadmat(net_file)
        Xnet.update({'coreness_norm': C_norm,
                     'coreness_zcore': C_zscore})
    sio.savemat(net_file, Xnet)
