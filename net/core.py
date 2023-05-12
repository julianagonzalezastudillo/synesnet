"""
=================================
        Core - Periphery
=================================
This module is design to compute coreness

   References:

       [1] Battiston F, Guillon J, Chavez M, Latora V, De Vico Fallani F.
       Multiplex core-periphery organization of the human connectome.
       J R Soc Interface. 2018 Sep 12;15(146):20180514.
       doi: 10.1098/rsif.2018.0514. PMID: 30209045; PMCID: PMC6170773.
"""
import numpy as np
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

    # Density thresholds
    THRESHOLDS = np.arange(1, N) / (N-1)
    isCore = np.empty((N, len(THRESHOLDS)))
    isCore[:] = np.nan

    for iThresh in np.arange(len(THRESHOLDS)):
        # threshold the weighted network by the matching density network
        T = bct.threshold_proportional(X, THRESHOLDS[iThresh])

        # normalize min-max
        T = normalize(T)

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
    isCore[rankingInd[0: rankOfMaxkPlus+1]] = 1
    isCore = isCore.astype(bool)
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
