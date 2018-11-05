import timeit
import numpy as np
from sklearn.model_selection import train_test_split
import sys

def load_data_down_GSE63347():
    from configurations.config_down_GSE63347 import config
    start = timeit.default_timer()
    X = np.genfromtxt(config.ifname("x"), dtype='float32', delimiter=' ')[:, 1:]

    genes_names = np.genfromtxt(config.ifname("x"), dtype='str', usecols = 0)
    config.params["num_genes"].value = min(genes_names.size, config.params["num_genes"].value)

    genes_dict = dict((v, i) for i, v in enumerate(genes_names))

    #ranged_genes = np.genfromtxt(config.ifname("ranged_genes"), dtype='str', usecols = 0)

    stop = timeit.default_timer()
    print 'Data loaded: ', stop - start
    print X.dtype, X.shape

    sys.stdout.flush()
    #X = np.random.rand(len(genes_names), 656)
    #genes_names = ranged_genes[:config.params["num_genes"].value]
    genes_names = genes_names[:config.params["num_genes"].value]

    indices = np.array([genes_dict[x] for x in genes_names])

    X = X[indices, :].T
    #X = X.T

    y = np.zeros((X.shape[0], 1), dtype = 'uint8')

    patients_info = np.genfromtxt(config.ifname("patients_info"), dtype='str', usecols = 1)
    y = (patients_info == "N").astype(np.uint8)

    config.params["normal_mask"].value = np.flatnonzero(y == 1)
    config.params["down_mask"].value = np.flatnonzero(y == 0)

    config.params["kde_mask"].value = "normal_mask"

    print X.shape, config.params["num_genes"].value
    sys.stdout.flush()

    mask = (y == 1)
    return X, y, X[mask, :], genes_names

def load_data_down_GSE63347_cpg_hannum():
    from configurations.config_down_GSE63347_cpg_hannum import config
    start = timeit.default_timer()
    X = np.genfromtxt(config.ifname("hannum_cpgs_beta"), dtype='float32', delimiter=' ')[:, 1:]

    cpgs_names  = np.genfromtxt(config.ifname("hannum_cpgs_beta"), dtype='str', usecols = 0)
    config.params["num_cpgs"].value = min(cpgs_names.size, config.params["num_cpgs"].value)

    stop = timeit.default_timer()
    print 'Data loaded: ', stop - start
    print X.dtype, X.shape

    sys.stdout.flush()

    X = X.T

    patients_info = np.genfromtxt(config.ifname("patients_info"), dtype='str', usecols = 1)
    y = (patients_info == "N").astype(np.uint8)

    config.params["normal_mask"].value = np.flatnonzero(y == 1)
    config.params["down_mask"].value = np.flatnonzero(y == 0)

    config.params["kde_mask"].value = "normal_mask"

    print X.shape, config.params["num_cpgs"].value
    sys.stdout.flush()

    mask = (y == 1)
    return X, y, X[mask, :], cpgs_names

def load_data_down_GSE63347_cpg_horvath():
    from configurations.config_down_GSE63347_cpg_horvath import config
    start = timeit.default_timer()
    X = np.genfromtxt(config.ifname("horvath_cpgs_beta"), dtype='float32', delimiter=' ')[:, 1:]

    cpgs_names  = np.genfromtxt(config.ifname("horvath_cpgs_beta"), dtype='str', usecols = 0)
    config.params["num_cpgs"].value = min(cpgs_names.size, config.params["num_cpgs"].value)

    stop = timeit.default_timer()
    print 'Data loaded: ', stop - start
    print X.dtype, X.shape

    sys.stdout.flush()

    X = X.T

    patients_info = np.genfromtxt(config.ifname("patients_info"), dtype='str', usecols = 1)
    y = (patients_info == "N").astype(np.uint8)

    config.params["normal_mask"].value = np.flatnonzero(y == 1)
    config.params["down_mask"].value = np.flatnonzero(y == 0)

    config.params["kde_mask"].value = "normal_mask"

    print X.shape, config.params["num_cpgs"].value
    sys.stdout.flush()

    mask = (y == 1)
    return X, y, X[mask, :], cpgs_names