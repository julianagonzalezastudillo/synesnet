import os.path
import scipy.io as sio
import numpy as np
from scipy import io
import networkx as nx
from net_metrics import local_laterality
import bct


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

for sub_file in os.listdir(path + 'symetrical_corr_mat'):
    print(sub_file)

    # load correlation
    fc_file = path + 'symetrical_corr_mat/' + sub_file
    fc = io.loadmat(fc_file)

    # connectivity matrix
    # Xfc = np.nan_to_num(fc['CorrMatrix'])
    Xfc = abs(fc['CorrMatrix'])
    np.fill_diagonal(Xfc, 0)

    # network metrics
    G = nx.from_numpy_matrix(Xfc)  # to nx format
    strength = np.array([v for k, v in G.degree(weight = 'weight')])  # strength
    lat = local_laterality(Xfc, n_name)  # laterality
    geff = bct.efficiency_wei(Xfc)
    leff = bct.efficiency_wei(Xfc, local=True)
    # sw = nx.sigma(Xfc)

    # save
    Xnet = {'strength': strength,
            'laterality': lat,
            'global_efficiency': geff,
            'local_efficiency': leff}
    net_file = net_path + '_'.join(('net_metrics', sub_file.split("_")[-1]))
    sio.savemat(net_file, Xnet)
    # net_file = net_path + '_'.join(('laterality', sub_file.split("_")[-1]))
    # sio.savemat(net_file, {'laterality': lat})


