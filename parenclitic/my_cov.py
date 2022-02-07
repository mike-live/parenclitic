import numpy as np
from numba import float64, int32, int64, uint64, int8, jit

class numba_config:
    cache = False
    nopython = True
    nogil = True
    parallel = True


@jit(nopython = numba_config.nopython, nogil = numba_config.nogil, cache = numba_config.cache)
def my_cov(data):
    m1 = 0.0
    for x in data[0]:
        m1 += x
    m2 = 0.0
    for x in data[1]:
        m2 += x
    m1 /= data.shape[1]
    m2 /= data.shape[1]
    c11 = 0.0
    c12 = 0.0
    c22 = 0.0
    for i in range(data.shape[1]):
        c12 += (data[0, i] - m1) * (data[1, i] - m2)
        c11 += (data[0, i] - m1) * (data[0, i] - m2)
        c22 += (data[1, i] - m1) * (data[1, i] - m2)
    c11 /= data.shape[1] - 1
    c12 /= data.shape[1] - 1
    c22 /= data.shape[1] - 1
    return np.array([[c11, c12], [c12, c22]])


@jit(nopython = numba_config.nopython, nogil = numba_config.nogil, cache = numba_config.cache)
def my_corrcoef(data):
    c = my_cov(data)
    d = np.diag(c)
    stddev = np.sqrt(d)

    c[0, 1] /= stddev[0] * stddev[1]
    c[1, 0] /= stddev[0] * stddev[1]
    c[0, 0] = 1
    c[1, 1] = 1

    np.clip(c, -1, 1, out=c)
    return c

"""
cov = my_cov
np.cov = my_cov
import numpy
numpy.cov = my_cov
"""
