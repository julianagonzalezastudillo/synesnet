import numpy as np
import bct


def coreness(X, f, c):
    N = np.shape(X)[0]
    # thresholds of
    THRESHOLDS = np.arange(1, N) / (N-1)
    isCore = np.empty((N, len(THRESHOLDS)))
    isCore[:] = np.nan

    for iThresh in np.arange(len(THRESHOLDS)):
        T = bct.threshold_proportional(X, THRESHOLDS[iThresh])  #threshold the weighted network by the matching density network
        T = normalize(T)

        # find core
        isCore[:, iThresh] = coreperiphery(T, f, c)

    C = np.mean(isCore, axis=1)  # to check if it's the right axis
    return C, isCore


def normalize(X):
    W = X
    # existing zeros to nan to keep them
    W[W==0] = np.nan
    minX = np.nanmin(W)
    maxX = np.nanmax(W)
    normX = (W - minX) / (maxX - minX)  # normalize min max
    normX[np.isnan(normX)] = 0
    return normX


def coreperiphery(X, f, c):
    N = np.shape(X)[0]
    c = 1  # Richness coefficients default values
    f = richness

    # Compute Core - Periphery
    [k, kMinus, kPlus] = richness(X)

    rankingInd = np.argsort(-k)  # rank descending order
    rankOfMaxkPlus = np.argmax(kPlus[rankingInd])

    isCore = np.zeros(N)
    isCore[rankingInd[0: rankOfMaxkPlus+1]] = 1
    isCore = isCore.astype(bool)
    return isCore


def richness(A):
    N = np.shape(A)[0]
    k = bct.strengths_und(A)
    kMinus = np.zeros(len(k))
    kPlus = np.zeros(len(k))

    for i in range(N):
        lrInd = k <= k[i]  # Indices of nodes with Lower Richness (LR)
        hrInd = k > k[i]  # Indices of nodes with Higher Richness (HR)

        lrA = np.copy(A)
        lrA[i, hrInd] = 0
        lrA[hrInd, i] = 0
        hrA = np.copy(A)
        hrA[i, lrInd] = 0
        hrA[lrInd, i] = 0

        kMinusForI = bct.strengths_und(lrA)
        kMinus[i] = kMinusForI[i]
        kPlusForI = bct.strengths_und(hrA)
        kPlus[i] = kPlusForI[i]

    return k, kMinus, kPlus
