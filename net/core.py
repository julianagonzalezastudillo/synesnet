"""
=================================
        Core - Periphery
=================================
This module is design to compute coreness in weighted networks, using strength as richness metric.

   References:

       [1] Battiston F, Guillon J, Chavez M, Latora V, De Vico Fallani F.
       Multiplex core-periphery organization of the human connectome.
       J R Soc Interface. 2018 Sep 12;15(146):20180514.
       doi: 10.1098/rsif.2018.0514. PMID: 30209045; PMCID: PMC6170773.
"""
import numpy as np
import matplotlib.pyplot as plt
import bct


def coreness(X):
    """
    Coreness C_i of each node i is defined as the normalized number of thresholds at
    which the corresponding node is present in the rich core.

    :param
        X: NxN np.ndarray
        undirected weighted connection matrix

    :return:
        C: corness
        isCore: core (0), periphery (0)

    """

    N = np.shape(X)[0]

    # Density max
    m = np.count_nonzero(X)
    if m == 0 or N <= 1:
        d_max = 0
    else:
        d_max = m / (N * (N - 1))

    # Density thresholds
    THRESHOLDS = np.arange(1, N) / (N - 1)

    # Check to use total number of links
    thr_idx = np.where(THRESHOLDS <= d_max)[0]
    L = np.count_nonzero(bct.threshold_proportional(X, THRESHOLDS[thr_idx[-1]]))
    if L < np.count_nonzero(X):
        thr_idx = np.append(thr_idx, thr_idx[-1] + 1)

    thr = THRESHOLDS[thr_idx]

    isCore = np.empty((N, len(thr)))
    isCore[:] = np.nan

    for iThresh, Thresh in enumerate(thr):
        # threshold the weighted network by the matching density network
        T = bct.threshold_proportional(X, Thresh)

        # normalize min-max
        T = normalize(T)

        # binarize
        # T[T > 0] = 1

        # find core
        isCore[:, iThresh] = coreperiphery(T)

    C = np.mean(isCore, axis=1)  # to check if it's the right axis
    return C, isCore


def normalize(X):
    # normalize between 0 and one
    # existing zeros to nan to keep them
    X[X == 0] = np.nan
    minX = np.nanmin(X)
    maxX = np.nanmax(X)

    # normalize min max
    normX = (X - minX) / (maxX - minX)

    # converted nan to zero
    normX[np.isnan(normX)] = 0
    return normX


def coreperiphery(X):
    """
    Compute Core - Periphery

    :param X: NxN np.ndarray
        undirected weighted connection matrix

    :return:
        isCore: core (1), periphery (0)
    """
    N = np.shape(X)[0]

    # get richness (k) and richness towards nodes with higher richness (kPlus)
    [k, kMinus, kPlus] = richness(X)

    # rank according to descending richness
    rankingInd = np.argsort(-k)

    # the value of the rank corresponding to the maximum of kPlus finally determines the
    # core–periphery structure
    rankOfMaxkPlus = np.argmax(kPlus[rankingInd])

    # nodes with rank lower than rankOfMaxkPlus are assigned to the core,
    # the remaining ones become part of the periphery
    isCore = np.zeros(N)
    isCore[rankingInd[0: rankOfMaxkPlus + 1]] = 1
    isCore = isCore.astype(bool)

    plot_rank = False
    if plot_rank:
        plot_core_periphery(X, kPlus, rankOfMaxkPlus, rankingInd)

    return isCore


def richness(A):
    """
    Computes the richness based on strength.

    :param A: NxN np.ndarray
        undirected weighted connection matrix

    :return:
        k: richness
        kMinus: richness of towards lower richer nodes
        kPlus: richness of towards richer nodes
    """
    N = np.shape(A)[0]

    # degree/strength vector --> richness
    k = bct.strengths_und(A)

    kMinus = np.zeros(len(k))
    kPlus = np.zeros(len(k))

    # divide the links of a node i in two groups,
    # those towards nodes with lower richness (lr) and
    # those towards nodes with higher richness (hr).
    for i in range(N):
        lrInd = k <= k[i]  # Indices of nodes with Lower Richness (lr)
        hrInd = k > k[i]  # Indices of nodes with Higher Richness (hr)

        lrA = np.copy(A)
        lrA[i, hrInd] = 0
        lrA[hrInd, i] = 0
        hrA = np.copy(A)
        hrA[i, lrInd] = 0
        hrA[lrInd, i] = 0

        kMinusForI = bct.strengths_und(lrA)
        kMinus[i] = kMinusForI[i]
        kPlusForI = bct.strengths_und(hrA)
        kPlus[i] = kPlusForI[i]  # richness of node i towards richer nodes

    return k, kMinus, kPlus


def plot_core_periphery(X, kPlus, rankOfMaxkPlus, rankingInd):
    """
    Plot k+ as a function of rank_i

    :param X:
    :param kPlus:
    :param rankOfMaxkPlus:
    :param rankingInd:
    :return:
    """

    N = np.shape(X)[0]
    density = (np.count_nonzero(X)) / (N*(N-1))

    k7 = 100/(N-1)
    if density < round(k7, 4):
        fig = plt.figure(figsize=(12, 6), dpi=100)
        plt.axvline(x=rankOfMaxkPlus, color='grey', linestyle='--', linewidth=4)
        plt.plot(np.arange(rankOfMaxkPlus, N), kPlus[rankingInd[rankOfMaxkPlus:]] / np.max(kPlus),
                 color='#203FB6', linewidth=4)
        plt.plot(np.arange(0, rankOfMaxkPlus + 1), kPlus[rankingInd[:rankOfMaxkPlus + 1]] / np.max(kPlus),
                 color='#F62336', linewidth=4)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.ylim(0, 1)
        plt.xlim(0, len(kPlus))

        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['bottom'].set_visible(True)
        plt.gca().spines['left'].set_visible(True)
        L = np.count_nonzero(X)
        plot_name = f'plots/core-periphery/rank_k+_{str(L).zfill(6)}.png'
        plt.savefig(plot_name, bbox_inches='tight', pad_inches=0, dpi=300, transparent=True)
        plt.show()
