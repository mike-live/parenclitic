import timeit
import numpy as np
from sklearn.model_selection import train_test_split
import sys

def get_classes(config, X):
    y = np.zeros((X.shape[0], ), dtype = 'int8')
    mask = y.copy()
    '''
    y[config.params["mongoloids_mask"].value] = 0
    y[config.params["siblings_mask"].value] = 1
    y[config.params["mothers_mask"].value] = 2
    if config.params["kde_mask"].value == "mothers_mask":
        mask[config.params["mongoloids_mask"].value] = 1
        mask[config.params["mothers_mask"].value] = -1
        mask[config.params["siblings_mask"].value] = -2
    elif config.params["kde_mask"].value == "siblings_mask":
        mask[config.params["mongoloids_mask"].value] = 1
        mask[config.params["siblings_mask"].value] = -1
        mask[config.params["mothers_mask"].value] = -2
    elif config.params["kde_mask"].value == "healthy_mask":
        mask[config.params["mongoloids_mask"].value] = 1
        mask[config.params["siblings_mask"].value] = -1
        mask[config.params["mothers_mask"].value] = -1
    elif config.params["kde_mask"].value == "mongoloids_mask":
        mask[config.params["mongoloids_mask"].value] = -1
        mask[config.params["siblings_mask"].value] = 1
        mask[config.params["mothers_mask"].value] = 2
    elif config.params["kde_mask"].value == "nonhealthy_mask":
        mask[config.params["mongoloids_mask"].value] = -1
        mask[config.params["siblings_mask"].value] = 1
        mask[config.params["mothers_mask"].value] = 1
    '''
    return y, mask

def load_data_down_GSE100825():
    from configurations.config_down_GSE100825 import config
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

    patients_info = pd.read_csv(config.ifname("patients_info"), sep = '\t')
    age = patients_info['age']
        
    return X, y, mask, cpgs_names, age
    