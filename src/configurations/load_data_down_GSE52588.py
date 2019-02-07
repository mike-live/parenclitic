import timeit
import numpy as np
from sklearn.model_selection import train_test_split
import sys

def load_data_down_GSE52588():
    from configurations.config_down_GSE52588 import config
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
    y[config.params["mongoloids_mask"].value] = 0
    y[config.params["siblings_mask"].value] = 1
    y[config.params["mothers_mask"].value] = 2
    print X.shape, config.params["num_genes"].value
    sys.stdout.flush()

    mask = config.params[config.params["kde_mask"].value].value
    return X, y, X[mask, :], genes_names

def load_data_down_GSE52588_cpgs():
    from configurations.config_down_GSE52588_cpg import config
    start = timeit.default_timer()
    import pandas as pd
    X = np.genfromtxt(config.ifname("x"), dtype='float32', delimiter='\t', skip_header = 1)[:, 1:]

    cpgs_names  = np.genfromtxt(config.ifname("x"), dtype='str', usecols = 0)
    config.params["num_cpgs"].value = min(cpgs_names.size, config.params["num_cpgs"].value)

    stop = timeit.default_timer()
    print 'Data loaded: ', stop - start
    print X.dtype, X.shape

    sys.stdout.flush()
    
    import pandas as pd
    cpgs_info = pd.read_csv(config.ifname("cpgs"), delimiter='\t')
    cpgs_island = cpgs_info["ID_REF"][cpgs_info["RELATION_TO_UCSC_CPG_ISLAND"] == "Island"]
    
    cpgs_island = np.array(cpgs_island)
    cpgs_names = np.array(cpgs_names)
    _, indices, _ = np.intersect1d(cpgs_names, cpgs_island, return_indices = True)

    print indices, X.shape

    good_cpgs = ~np.isnan(X).any(axis=0)
    X = X[indices, :].T
    X = X[:, good_cpgs]
    cpgs_names = cpgs_names[good_cpgs]
    print X.shape
    #X = X.T

    y = np.zeros((X.shape[0], ), dtype = 'uint8')
    y[config.params["mongoloids_mask"].value] = 0
    y[config.params["siblings_mask"].value] = 1
    y[config.params["mothers_mask"].value] = 2
    
    config.params["num_cpgs"].value = min(X.shape[1], config.params["num_cpgs"].value)
    
    print X.shape, config.params["num_cpgs"].value, y.shape
    sys.stdout.flush()

    mask = config.params[config.params["kde_mask"].value].value
    return X, y, X[mask, :], cpgs_names
    
def load_data_mongoloids_cpgs():
    from mongoloids_cpgs_config import config
    start = timeit.default_timer()
    X = np.genfromtxt(config.ifname("cpgs"), dtype='float32', delimiter=' ')[:, 2:]

    genes_names = np.genfromtxt(config.ifname("cpgs"), dtype='str', usecols = 0)
    cpgs_names  = np.genfromtxt(config.ifname("cpgs"), dtype='str', usecols = 1)
    config.params["num_cpgs"].value = min(cpgs_names.size, config.params["num_cpgs"].value)

    genes_dict = dict((v, i) for i, v in enumerate(genes_names))

    #ranged_genes = np.genfromtxt(config.ifname("ranged_genes"), dtype='str', usecols = 0)

    stop = timeit.default_timer()
    print 'Data loaded: ', stop - start
    print X.dtype, X.shape

    sys.stdout.flush()
    #X = np.random.rand(len(genes_names), 656)
    #genes_names = ranged_genes[:config.params["num_genes"].value]
    #genes_names = genes_names[:config.params["num_genes"].value]

    #indices = np.array([genes_dict[x] for x in genes_names])

    #X = X[indices, :].T
    X = X.T
    #X = X.T

    y = np.zeros((X.shape[0], 1), dtype = 'uint8')
    y[config.params["mongoloids_mask"].value] = 0
    y[config.params["siblings_mask"].value] = 1
    y[config.params["mothers_mask"].value] = 2
    print X.shape, config.params["num_cpgs"].value
    sys.stdout.flush()

    mask = config.params[config.params["kde_mask"].value].value
    return X, y, X[mask, :], genes_names, cpgs_names

def load_data_mongoloids_horvath_cpgs():
    from config_mongoloids_cpg_horvath import config
    start = timeit.default_timer()
    X = np.genfromtxt(config.ifname("horvath_cpgs_beta"), dtype='float32', delimiter=' ')[:, 1:]

    cpgs_names  = np.genfromtxt(config.ifname("horvath_cpgs_beta"), dtype='str', usecols = 0)
    config.params["num_cpgs"].value = min(cpgs_names.size, config.params["num_cpgs"].value)

    stop = timeit.default_timer()
    print 'Data loaded: ', stop - start
    print X.dtype, X.shape

    sys.stdout.flush()

    X = X.T

    y = np.zeros((X.shape[0], 1), dtype = 'uint8')
    y[config.params["mongoloids_mask"].value] = 0
    y[config.params["siblings_mask"].value] = 1
    y[config.params["mothers_mask"].value] = 2
    print X.shape, config.params["num_cpgs"].value
    sys.stdout.flush()

    mask = config.params[config.params["kde_mask"].value].value
    return X, y, X[mask, :], cpgs_names

def load_data_mongoloids_hannum_cpgs():
    from mongoloids_cpg_hannum_config import config
    start = timeit.default_timer()
    X = np.genfromtxt(config.ifname("hannum_cpgs_beta"), dtype='float32', delimiter=' ')[:, 1:]

    cpgs_names  = np.genfromtxt(config.ifname("hannum_cpgs_beta"), dtype='str', usecols = 0)
    config.params["num_cpgs"].value = min(cpgs_names.size, config.params["num_cpgs"].value)

    stop = timeit.default_timer()
    print 'Data loaded: ', stop - start
    print X.dtype, X.shape

    sys.stdout.flush()

    X = X.T

    y = np.zeros((X.shape[0], 1), dtype = 'uint8')
    y[config.params["mongoloids_mask"].value] = 0
    y[config.params["siblings_mask"].value] = 1
    y[config.params["mothers_mask"].value] = 2
    print X.shape, config.params["num_cpgs"].value
    sys.stdout.flush()

    mask = config.params[config.params["kde_mask"].value].value
    return X, y, X[mask, :], cpgs_names