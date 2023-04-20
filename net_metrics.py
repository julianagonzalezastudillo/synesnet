import numpy as np


def channel_idx(ch_names):
    RH_idx = []
    LH_idx = []
    for i in range(len(ch_names)):
        if ch_names[i][-1] == 'R':
            RH_idx = np.append(RH_idx, i)
            LH_idx = np.append(LH_idx, ch_names.index(ch_names[i][:-1] + 'L'))  # find homologue
            # print('R: {0}, L: {1}'.format(ch_names[int(i)], ch_names[ch_names.index(ch_names[i][:-1] + 'L')]))
        elif ch_names[i][-2:] == 'R\n':
            RH_idx = np.append(RH_idx, i)
            LH_idx = np.append(LH_idx, ch_names.index(ch_names[i][:-2] + 'L\n'))  # find homologue
            # print('R: {0}, L: {1}'.format(ch_names[int(i)], ch_names[ch_names.index(ch_names[i][:-2] + 'L\n')]))

    return RH_idx.astype(int), LH_idx.astype(int)


def local_laterality(X, ch_names):
    """
    order LEFT-RIGHT
    :param X:
    :param ch_names:
    :return:
    """
    RH_idx, LH_idx = channel_idx(ch_names)
    # LH_names = [ch_names[index] for index in LH_idx]
    # print('LH: {0}'.format(LH_names))
    # RH_names = [ch_names[index] for index in RH_idx]
    # print('RH: {0}'.format(RH_names))

    LH = X[LH_idx, :][:, LH_idx]
    RH = X[RH_idx, :][:, RH_idx]

    lat_homol = np.zeros(len(LH) + len(RH))
    # for each node
    for j in range(len(RH_idx)):
        # HOMOLOGOUS NODES
        lat_homol[j] = (np.sum(LH[j]) - np.sum(RH[j]))
        lat_homol[j + len(RH_idx)] = (np.sum(RH[j]) - np.sum(LH[j]))

    return lat_homol

