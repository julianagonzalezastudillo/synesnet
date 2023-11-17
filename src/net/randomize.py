import numpy as np
import random


def check_equal_elements(matrix1, matrix2):
    row, col = np.triu_indices(len(matrix1), k=1)
    for i, j in zip(row, col):
        if matrix1[i][j] == matrix2[i][j]:
            print(i, j)
            return True
    return False


def rand_wu(W):
    """
    This function randomizes an undirected weighted network with positive and
    negative weights. The function does not preserve the strength distribution.

    Parameters
    ----------
    W : NxN np.ndarray
        undirected binary/weighted connection matrix

    Returns
    -------
    R : NxN np.ndarray
        randomized network
    """

    R = W.copy()
    R0 = W.copy()

    # while check_equal_elements(R, R0):
    row, col = np.triu_indices(len(R), k=1)
    while not (len(row) / 2).is_integer():
        add_n = random.sample(range(len(row)), 1)
        row = np.append(row, row[add_n])
        col = np.append(col, col[add_n])

    while len(row) > 0:
        a, c = random.sample(range(len(row)), 2)

        r0_a = R0[row[a], col[a]]
        r0_c = R0[row[c], col[c]]

        R[row[a], col[a]] = R[col[a], row[a]] = r0_c
        R[row[c], col[c]] = R[col[c], row[c]] = r0_a

        row = np.delete(row, [a, c])
        col = np.delete(col, [a, c])

    return R
