import numpy as np
from scipy import io
import os.path


def load_net_metrics(net_path, metric, corr_type, num_sub, idx_select):
    Xnet = np.array([
        io.loadmat(os.path.join(net_path, f'net_metrics_Subject{str(sub).zfill(3)}{corr_type}'))[metric][0]
        for sub in range(1, num_sub + 1)])
    X_ = np.zeros(np.shape(Xnet))
    X_[:, idx_select] = Xnet[:, idx_select]
    return X_