import timeit
import numpy as np
from sklearn.model_selection import train_test_split
import sys
import pandas as pd

def get_classes(config, X):
    y = np.zeros((X.shape[0], ), dtype = 'int8')
    mask = y.copy()
    y[config.params["control_mask"].value] = 0
    y[config.params["schizophrenia_mask"].value] = 1
    if config.params["kde_mask"].value == "control_mask":
        mask[config.params["control_mask"].value] = -1
        mask[config.params["schizophrenia_mask"].value] = 1
    elif config.params["kde_mask"].value == "schizophrenia_mask":
        mask[config.params["control_mask"].value] = 1
        mask[config.params["schizophrenia_mask"].value] = -1
    return y, mask

def load_data_schizophrenia_cpgs(is_small=False, is_train=True):
    from configurations.config_schizophrenia_cpg import config
    
    #import pickle
    
    start = timeit.default_timer()
    #with open(config.ifname('beta_values'), 'rb') as fp:
    #    X = pickle.load(fp)
    if is_train:
        X = pd.read_pickle(config.ifname('beta_values'))
    else:
        X = pd.read_pickle(config.ifname('beta_values_full'))
    
    config.params["control_mask"].value = X['Status'].values == 'Control'
    config.params["schizophrenia_mask"].value = X['Status'].values == 'Schizophrenia'
    cpgs_names = X.columns
    X = X.iloc[:, 2:].values.astype('float32')
    print('3')
    #X = X[:, :100]
    #cpgs_names = cpgs_names[:100]

    stop = timeit.default_timer()
    print('Data loaded: ', stop - start)
    print(X.dtype, X.shape)

    sys.stdout.flush()
    
    if is_small:
        from annotations.cpgs import cpgs_annotation
        cpgs = cpgs_annotation(config.ifname('cpgs_annotations'))
        bad_cpgs = np.loadtxt(config.ifname('bad_cpgs'), dtype='str')

        subset_cpg_names, ids = cpgs.get_cpgs({'gene_out': [np.NaN], 
                                         'cpgs_in': cpgs_names, 
                                         'chr_out': ['X', 'Y'], 
                                         'geotype_in': ['Island'],
                                         'cpgs_out': bad_cpgs})
        cpgs_names, indices, _ = np.intersect1d(cpgs_names, subset_cpg_names, return_indices = True)
        X = X[:, indices]
    
    good_cpgs = ~np.isnan(X[:config.params["num_train"].value]).any(axis=0)
    X = X[:, good_cpgs]
    cpgs_names = cpgs_names[good_cpgs]
    print(X.shape)
    
    config.params["num_cpgs"].value = min(X.shape[1], config.params["num_cpgs"].value)
    y, mask = get_classes(config, X)
    if not is_train:
        mask[config.params["num_train"].value:] *= 2

    print(X.shape, config.params["num_cpgs"].value, y.shape, cpgs_names.shape)
    print('Fraction of zeros in y', (y == 0).mean())
    print('Fraction of ones  in y', (y == 1).mean())
    print('mask info:')
    print(' 0:', (mask == 0).sum())
    print('-1:', (mask == -1).sum(), '1:', (mask == 1).sum())
    print('-2:', (mask == -2).sum(), '2:', (mask == 2).sum())
    print()
    sys.stdout.flush()
    
    return X, y, mask, cpgs_names
    

def _load_data_schizophrenia_cpgs(is_small = False):
    from configurations.config_schizophrenia_cpg import config
    start = timeit.default_timer()
    import pandas as pd
    X = np.genfromtxt(config.ifname("x"), dtype='float32', delimiter='\t', skip_header = 0)[:, 1:]

    cpgs_names = np.genfromtxt(config.ifname("x"), dtype='str', skip_header = 0, usecols = 0)
    config.params["num_cpgs"].value = min(cpgs_names.size, config.params["num_cpgs"].value)

    stop = timeit.default_timer()
    print('Data loaded: ', stop - start)
    print(X.dtype, X.shape)

    sys.stdout.flush()
    
    if is_small:
        from annotations.cpgs import cpgs_annotation
        cpgs = cpgs_annotation(config.ifname('cpgs_annotations'))
        bad_cpgs = np.loadtxt(config.ifname('bad_cpgs'), dtype='str')

        subset_cpg_names, ids = cpgs.get_cpgs({'gene_out': [np.NaN], 
                                         'cpgs_in': cpgs_names, 
                                         'chr_out': ['X', 'Y'], 
                                         'geotype_in': ['Island'],
                                         'cpgs_out': bad_cpgs})
        cpgs_names, indices, _ = np.intersect1d(cpgs_names, subset_cpg_names, return_indices = True)

        #cpgs_island = np.array(cpgs_island)
        #cpgs_names = np.array(cpgs_names)
        #_, indices, _ = np.intersect1d(cpgs_names, cpgs_all, return_indices = True)

        #cpgs_names = cpgs_names[indices]
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