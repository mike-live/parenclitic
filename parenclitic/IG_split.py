import numpy as np

#from numba import float64, int32, int64, uint64, int8, jit
#import numba_scipy

class numba_config:
    cache = False
    nopython = True
    nogil = True
    parallel = True


import math

#@jit(nopython = numba_config.nopython, nogil = numba_config.nogil, cache = numba_config.cache)
def xlogy(x, y):
    if x == 0:
        return 0
    else:
        return x * math.log(y)

#@jit(nopython = numba_config.nopython, nogil = numba_config.nogil, cache = numba_config.cache)
def npxlogy(x, y):
    res = np.empty(x.shape)
    for i in range(len(res)):
        res[i] = xlogy(x[i], y[i])
    return res

#@jit(nopython = numba_config.nopython, nogil = numba_config.nogil, cache = numba_config.cache)
def npentropy(a, b, n):
    #from scipy.special import xlogy
    p0 = a / n
    p1 = b / n
    return -npxlogy(p0, p0) - npxlogy(p1, p1)

#@jit(nopython = numba_config.nopython, nogil = numba_config.nogil, cache = numba_config.cache)
def entropy(a, b, n):
    #from scipy.special import xlogy
    p0 = a / n
    p1 = b / n
    return -xlogy(p0, p0) - xlogy(p1, p1)

#@jit(nopython = numba_config.nopython, nogil = numba_config.nogil, cache = numba_config.cache)
def IG_split(p, y):
    ids = np.argsort(p)
    n = len(p)
    n1 = 0
    n2 = 0
    for tid in ids:
        n1 += y[tid] == -1
        n2 += y[tid] == +1
    
    sm = 0
    sp = 0
    en = entropy(n1, n2, n)
    max_gain_information = -np.inf
    for i in range(len(ids) - 1):
        tid = ids[i]
        sm += y[tid] == -1
        sp += y[tid] == +1
        n11 = sm
        n21 = sp
        n12 = n1 - n11
        n22 = n2 - n21
        gain_information = en - (entropy(n11, n21, n11 + n21) * (n11 + n21) + entropy(n12, n22, n12 + n22) * (n12 + n22)) / n
        if gain_information > max_gain_information:
            max_gain_information = gain_information
            id_imp = i

    thr = (p[ids[id_imp]] + p[ids[id_imp + 1]]) / 2
    best_gain_information = max_gain_information / (2 * np.log(2))
    acc = ((thr > p) == (y > 0)).mean()
    return thr, best_gain_information, acc
    
def IG_split_old(p, y):
    ids = np.argsort(p)
    p = p[ids]
    y = y[ids]
    n = np.array(len(y), dtype = p.dtype)

    n11 = np.cumsum(y == -1).astype(p.dtype)
    n21 = np.cumsum(y == +1).astype(p.dtype)
    n1 = n11[-1]
    n2 = n21[-1]
    n11 = n11[:-1]
    n21 = n21[:-1]
    n12 = n1 - n11
    n22 = n2 - n21

    gain_information = entropy(n1, n2, n) - (npentropy(n11, n21, n11 + n21) * (n11 + n21) + 
                                             npentropy(n12, n22, n12 + n22) * (n12 + n22)) / n
    id_imp = np.argmax(gain_information)
    thr = (p[id_imp] + p[id_imp + 1]) / 2
    best_gain_information = gain_information[id_imp] / (2 * np.log(2))
    acc = ((thr > p) == (y > 0)).mean()
    return thr, best_gain_information, acc

