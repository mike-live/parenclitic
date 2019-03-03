import timeit
import numpy as np
from sklearn.model_selection import train_test_split
import sys
import pandas as pd

def get_classes(config, X):
    patients_info = pd.read_csv(config.ifname("patients_info"), sep = ' ')
    age = patients_info['age']
    num_groups = config.params["num_groups"].value
    #age_group = config.params["age_group"].value
    bins = np.percentile(age, np.linspace(0, 100, num_groups + 1))
    y = np.digitize(age, bins) - 1

    config.params["young_mask"].value = np.flatnonzero(y == 1)
    config.params["old_mask"].value = np.flatnonzero(y == 0)

    #config.params["kde_mask"].value = "young_mask"
    mask = y == 0
    return y, mask, age


def load_data_age_GSE87571():
    from configurations.config_age_GSE87571 import config
    start = timeit.default_timer()
    data = np.load(config.ifname("x"))
    X = data['X']
    genes_names = data['genes_names']
    
    config.params["num_genes"].value = min(genes_names.size, config.params["num_genes"].value)

    stop = timeit.default_timer()
    print 'Data loaded: ', stop - start
    print X.dtype, X.shape

    sys.stdout.flush()
    genes_names = genes_names[:config.params["num_genes"].value]

    y, mask, age = get_classes(config, X)

    print X.shape, config.params["num_genes"].value
    sys.stdout.flush()

    return X, y, mask, genes_names, age

def load_data_age_GSE87571_cpgs():
    from configurations.config_age_GSE87571_cpg import config
    start = timeit.default_timer()
    data = np.load(config.ifname("cpgs"))
    X = data['cpg_data']
    #cpgs_names  = np.genfromtxt(config.ifname("cpgs_annotations.txt"), dtype='str', usecols = 0)
    data = np.load(config.ifname("cpgs_names"))
    cpgs_names = data['cpgs_names']
    
    stop = timeit.default_timer()
    print 'Data loaded: ', stop - start
    print X.dtype, X.shape

    sys.stdout.flush()

    cpgs_info = pd.read_csv(config.ifname("cpgs_annotations"), delimiter='\t')
    cpgs_island = cpgs_info["ID_REF"][cpgs_info["RELATION_TO_UCSC_CPG_ISLAND"] == "Island"]
    
    cpgs_island = np.array(cpgs_island)
    cpgs_names = np.array(cpgs_names)
    _, indices, _ = np.intersect1d(cpgs_names, cpgs_names, return_indices = True)
    cpgs_names = cpgs_names[indices]

    print indices, X.shape

    X = X.T
    #cpgs_names = cpgs_names[indices]
    #X = X[:, indices]
    
    good_cpgs = ~np.isnan(X).any(axis=0)
    X = X[:, good_cpgs]
    cpgs_names = cpgs_names[good_cpgs]
    print X.shape
    
    config.params["num_cpgs"].value = min(X.shape[1], config.params["num_cpgs"].value)

    #patients_info = np.genfromtxt(config.ifname("patients_info"), dtype='str', usecols = 1)
    return X, y, X[mask, :], cpgs_names, age    
    
def load_data_age_GSE87571_cpg_horvath():
    from configurations.config_age_GSE87571_cpg_horvath import config
    start = timeit.default_timer()
    X = np.genfromtxt(config.ifname("horvath_cpgs_beta"), dtype='float32', delimiter=' ')[:, 1:]

    cpgs_names  = np.genfromtxt(config.ifname("horvath_cpgs_beta"), dtype='str', usecols = 0)
    config.params["num_cpgs"].value = min(cpgs_names.size, config.params["num_cpgs"].value)

    stop = timeit.default_timer()
    print 'Data loaded: ', stop - start
    print X.dtype, X.shape

    sys.stdout.flush()

    X = X.T

    #patients_info = np.genfromtxt(config.ifname("patients_info"), dtype='str', usecols = 1)
    patients_info = pd.read_csv(config.ifname("patients_info"), sep = ' ')
    y = (np.array(patients_info['age'] < config.params["delimiter_age"].value)).astype(np.uint8)

    config.params["young_mask"].value = np.flatnonzero(y == 1)
    config.params["old_mask"].value = np.flatnonzero(y == 0)

    config.params["kde_mask"].value = "young_mask"

    print X.shape, config.params["num_cpgs"].value
    sys.stdout.flush()

    age = patients_info['age']
    
    mask = (y == 1)
    return X, y, X[mask, :], cpgs_names, age