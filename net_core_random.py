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

# Constants
RAND_PRESERVE_STRENGTH = True

if not os.path.exists(net_path):
    os.makedirs(net_path)

n_name = []
with open(path + "BN_Atlas_246_LUT_reoriented.txt", "r") as filestream:
    for line in filestream:
        n_name.append(line.split(",")[0])

# Open strength t-val file
strength_stat_file = os.path.join(os.getcwd(), 'plots', 'glb', 'strength_thr_t-val.mat')
strength_stat = io.loadmat(strength_stat_file)
idx_select = strength_stat['names_idx'][0]

# for sub_file in os.listdir(path + 'symetrical_corr_mat'):
strength_select = False
if RAND_PRESERVE_STRENGTH:
    selec = '_conserve_strenght_distribution'
    rand_folder = 'rand_mat_strength(equal)'
else:
    selec = ''
    rand_folder = 'rand_mat'

for sub in range(6, 22):
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

    selec, Xfc_thr = ('_strength_selection', Xfc_thr[idx_select][:, idx_select]) if strength_select else ('', Xfc_thr)

    # Open random matrices
    # the randomized matrices (already) only contain Xfc>0
    fc_rand_file = os.path.join(path, rand_folder, 'RandMatrices_Subject{0}{1}'.format(str(sub).zfill(3), selec))
    fc_rand = io.loadmat(fc_rand_file)
    Xfc_rand = deepcopy(fc_rand['RandMatrices'])
    Xfc_rand[Xfc_rand <= 0] = 0
    C_rand = []  # random networks coreness
    # isCore_rand = []
    for m in np.arange(np.shape(Xfc_rand)[0]):
        print('random network n#{0}'.format(m))
        X = Xfc_rand[m]
        np.fill_diagonal(X, 0)
        C, isCore = coreness(X)  # coreness
        C_rand.append(C)
        # isCore_rand.append(isCore)

    # isC_rand = np.array([isC for isC in isCore_rand])
    C_rand = np.array(C_rand)
    C_rand_mean = np.mean(C_rand, axis=0)
    C_orig, isCore_orig = coreness(Xfc_thr)
    C_norm = C_orig / C_rand_mean  # coreness ratio

    # Replace NaN and infinite values with zeros
    C_norm[np.isnan(C_norm) | np.isinf(C_norm)] = 0

    C_zscore = (C_orig - C_rand_mean) / np.std(C_rand, axis=0)  # coreness z-score
    print(time.time() - t)

    # save
    net_file = net_path + '_'.join(('net_metrics', sub_file.split("_")[-1], 'thr.mat'))
    if os.path.exists(net_file):
        Xnet = sio.loadmat(net_file)
        Xnet.update({'coreness_norm_by_rand{0}'.format(selec): np.nan_to_num(C_norm),
                     'coreness_zcore_by_rand{0}'.format(selec): np.nan_to_num(C_zscore),
                     'coreness_rand{0}'.format(selec): np.nan_to_num(C_rand)})
    sio.savemat(net_file, Xnet)



#%%
# matrix = np.random.uniform(low=-1, high=1, size=(246, 246))
# np.fill_diagonal(matrix, 0)
# matrix[matrix <= 0] = 0
