import os.path
import scipy.io as sio
import numpy as np
from scipy import io
import networkx as nx


path = '/Users/juliana.gonzalez/ownCloud/graph_analysis/'
net_path = '/Users/juliana.gonzalez/ownCloud/graph_analysis/net_metrics/'
if not os.path.exists(net_path):
    os.makedirs(net_path)

for sub_file in os.listdir(path + 'symetrical_corr_mat'):
    print(sub_file)

    # load correlation
    fc_file = path + 'symetrical_corr_mat/' + sub_file
    fc = io.loadmat(fc_file)

    # connectivity matrix
    # Xfc = np.nan_to_num(fc['CorrMatrix'])
    Xfc = fc['CorrMatrix']
    np.fill_diagonal(Xfc, 0)

    # network metrics
    G = nx.from_numpy_matrix(Xfc)  # to nx format
    strength = [v for k, v in G.degree(weight = 'weight')]

    # save
    Xnet = {'strength': strength}
    net_file = net_path + '_'.join(('strength', sub_file.split("_")[-1]))
    sio.savemat(net_file, {'strength': strength})

