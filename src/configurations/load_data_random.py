import timeit
import numpy as np
from sklearn.model_selection import train_test_split
import sys

def get_classes(config, X):
    y = np.zeros((X.shape[0], ), dtype = 'uint8')
    y[29:] = 1
    mask = y == 0
    return y, mask

def load_data_random():
    from configurations.config_random_data import config
    start = timeit.default_timer()
    X = np.zeros(shape = (config.params["num_genes"].value, config.params["num_samples"].value), dtype = np.float32)
    for i in range(config.params["num_genes"].value):
        loc = np.random.uniform(low = 0.01, high = 0.99)
        scale = np.exp(np.random.uniform(low = np.log(0.002), high = np.log(0.05)))
        X[i, :] = np.clip(np.random.normal(loc = loc, scale = scale, size = (1, config.params["num_samples"].value)), 0, 1).astype(np.float32)

    genes_names = list(map(str, range(config.params["num_genes"].value)))
    
    stop = timeit.default_timer()
    print('Data loaded: ', stop - start)
    print(X.dtype, X.shape)

    sys.stdout.flush()
    X = X.T
    y, mask = get_classes(config, X)

    print(X.shape, config.params["num_genes"].value)
    sys.stdout.flush()

    return X, y, mask, genes_names

