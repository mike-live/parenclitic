import timeit
import numpy as np
from sklearn.model_selection import train_test_split
import sys

def get_classes(config, X):
    y = np.zeros((X.shape[0], 1), dtype = 'uint8')

    patients_info = np.genfromtxt(config.ifname("patients_info"), dtype='str', usecols = (1))
    y = (patients_info == "N").astype(np.uint8)

    config.params["normal_mask"].value = np.flatnonzero(y == 1)
    config.params["down_mask"].value = np.flatnonzero(y == 0)

    config.params["kde_mask"].value = "normal_mask"

    return y, mask

def load_data_down_GSE63347():
    from configurations.config_down_GSE63347 import config
    start = timeit.default_timer()
    X = np.genfromtxt(config.ifname("x"), dtype='float32', delimiter=' ')[:, 1:]

    genes_names = np.genfromtxt(config.ifname("x"), dtype='str', usecols = 0)
    config.params["num_genes"].value = min(genes_names.size, config.params["num_genes"].value)

    genes_dict = dict((v, i) for i, v in enumerate(genes_names))

    #ranged_genes = np.genfromtxt(config.ifname("ranged_genes"), dtype='str', usecols = 0)

    stop = timeit.default_timer()
    print('Data loaded: ', stop - start)
    print(X.dtype, X.shape)

    sys.stdout.flush()
    #X = np.random.rand(len(genes_names), 656)
    #genes_names = ranged_genes[:config.params["num_genes"].value]
    genes_names = genes_names[:config.params["num_genes"].value]

    indices = np.array([genes_dict[x] for x in genes_names])

    X = X[indices, :].T
    #X = X.T

    y = np.zeros((X.shape[0], 1), dtype = 'uint8')

    patients_info = np.genfromtxt(config.ifname("patients_info"), dtype='str', usecols = (1))
    y = (patients_info == "N").astype(np.uint8)

    config.params["normal_mask"].value = np.flatnonzero(y == 1)
    config.params["down_mask"].value = np.flatnonzero(y == 0)

    config.params["kde_mask"].value = "normal_mask"

    print(X.shape, config.params["num_genes"].value)
    sys.stdout.flush()
   
    mask = (y == 1)
    return X, y, X[mask, :], genes_names

    
def load_data_down_GSE63347_cpgs():
    from configurations.config_down_GSE63347_cpg import config
    start = timeit.default_timer()
    import pandas as pd
    X = np.genfromtxt(config.ifname("x"), dtype='float32', delimiter='\t', skip_header = 1)[:, 1:]

    cpgs_names  = np.genfromtxt(config.ifname("x"), dtype='str', skip_header = 1, usecols = 0)
    config.params["num_cpgs"].value = min(cpgs_names.size, config.params["num_cpgs"].value)

    stop = timeit.default_timer()
    print('Data loaded: ', stop - start)
    print(X.dtype, X.shape)

    sys.stdout.flush()
    
    import pandas as pd
    cpgs_info = pd.read_csv(config.ifname("cpgs_annotations"), delimiter='\t')
    cpgs_island = cpgs_info["ID_REF"][cpgs_info["RELATION_TO_UCSC_CPG_ISLAND"] == "Island"]
    
    cpgs_all = np.array(cpgs_info["ID_REF"])
    
    cpgs_island = np.array(cpgs_island)
    cpgs_names = np.array(cpgs_names)
    _, indices, _ = np.intersect1d(cpgs_names, cpgs_all, return_indices = True)
    
    cpgs_names = cpgs_names[indices]
    X = X[indices, :]
    X = X.T
    
    good_cpgs = ~np.isnan(X).any(axis=0)
    X = X[:, good_cpgs]
    cpgs_names = cpgs_names[good_cpgs]
    print(X.shape)
    #X = X.T
    
    config.params["num_cpgs"].value = min(X.shape[1], config.params["num_cpgs"].value)
    y, mask = get_classes(config, X)
    #print y 
    #print mask
    print(X.shape, config.params["num_cpgs"].value, y.shape, cpgs_names.shape)
    sys.stdout.flush()
        
    return X, y, mask, cpgs_names

def load_data_down_GSE63347_cpg_hannum():
    from configurations.config_down_GSE63347_cpg_hannum import config
    start = timeit.default_timer()
    X = np.genfromtxt(config.ifname("hannum_cpgs_beta"), dtype='float32', delimiter=' ')[:, 1:]

    cpgs_names  = np.genfromtxt(config.ifname("hannum_cpgs_beta"), dtype='str', usecols = 0)
    config.params["num_cpgs"].value = min(cpgs_names.size, config.params["num_cpgs"].value)

    stop = timeit.default_timer()
    print('Data loaded: ', stop - start)
    print(X.dtype, X.shape)

    sys.stdout.flush()

    X = X.T

    patients_info = np.genfromtxt(config.ifname("patients_info"), dtype='str', usecols = 1)
    y = (patients_info == "N").astype(np.uint8)

    config.params["normal_mask"].value = np.flatnonzero(y == 1)
    config.params["down_mask"].value = np.flatnonzero(y == 0)

    config.params["kde_mask"].value = "normal_mask"

    print(X.shape, config.params["num_cpgs"].value)
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
    print('Data loaded: ', stop - start)
    print(X.dtype, X.shape)

    sys.stdout.flush()

    X = X.T

    patients_info = np.genfromtxt(config.ifname("patients_info"), dtype='str', usecols = 1)
    y = (patients_info == "N").astype(np.uint8)

    config.params["normal_mask"].value = np.flatnonzero(y == 1)
    config.params["down_mask"].value = np.flatnonzero(y == 0)

    config.params["kde_mask"].value = "normal_mask"

    print(X.shape, config.params["num_cpgs"].value)
    sys.stdout.flush()

    mask = (y == 1)
    return X, y, X[mask, :], cpgs_names