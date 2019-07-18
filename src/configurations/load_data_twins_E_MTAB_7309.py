import timeit
import numpy as np
import sys
import pandas as pd

def get_classes(config, X, ):
    y = np.zeros((X.shape[0], ), dtype = 'uint8')
    patients_info = pd.read_csv(config.ifname('patients_info'))
    age = patients_info['age'].values
    config.params["less_65_mask"].value = age < 65
    config.params["greater_80_mask"].value = age > 80
    config.params["between_65_80_mask"].value = \
        np.bitwise_not(np.bitwise_or(config.params["less_65_mask"].value, config.params["greater_80_mask"].value))
    
    y[config.params["less_65_mask"].value] = 0
    y[config.params["greater_80_mask"].value] = 1
    y[config.params["between_65_80_mask"].value] = 2
    mask = y.copy()
    return y, mask, age

def load_data_twins_E_MTAB_7309():
    from configurations.config_twins_E_MTAB_7309 import config
    start = timeit.default_timer()
    print(config.ifname("x"))
    data = np.load(config.ifname("x"))
    print(data.files)
    X = data['data']

    gene_names = np.genfromtxt(config.ifname("gene_names"), dtype='str', usecols = 0)
    config.params["num_genes"].value = gene_names.size

    stop = timeit.default_timer()
    print('Data loaded: ', stop - start)
    print(X.dtype, X.shape)

    sys.stdout.flush()
    
    X = X.T
    y, mask, age = get_classes(config, X)

    print(X.shape, config.params["num_genes"].value)
    sys.stdout.flush()

    return X, y, mask, gene_names, age
