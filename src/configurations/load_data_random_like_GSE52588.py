import timeit
import numpy as np
from sklearn.model_selection import train_test_split
import sys

def get_classes(config, X):
    y = np.zeros((X.shape[0], ), dtype = 'uint8')
    y[config.params["mongoloids_mask"].value] = 0
    y[config.params["siblings_mask"].value] = 1
    y[config.params["mothers_mask"].value] = 2
    mask = y.copy()
    if config.params["kde_mask"].value == "mothers_mask":
        mask[config.params["mongoloids_mask"].value] = 0
        mask[config.params["mothers_mask"].value] = 1
        mask[config.params["siblings_mask"].value] = 2
    elif config.params["kde_mask"].value == "siblings_mask":
        mask[config.params["mongoloids_mask"].value] = 0
        mask[config.params["siblings_mask"].value] = 1
        mask[config.params["mothers_mask"].value] = 2
    elif config.params["kde_mask"].value == "healthy_mask":
        mask[config.params["mongoloids_mask"].value] = 0
        mask[config.params["siblings_mask"].value] = 1
        mask[config.params["mothers_mask"].value] = 1
    elif config.params["kde_mask"].value == "mongoloids_mask":
        mask[config.params["mongoloids_mask"].value] = 1
        mask[config.params["siblings_mask"].value] = 0
        mask[config.params["mothers_mask"].value] = 0
    return y, mask

def load_data_random_like_GSE52588():
    from configurations.load_data_down_GSE52588 import load_data_down_GSE52588
    Xw, _, _, genes_names = load_data_down_GSE52588()
    
    from configurations.config_random_data_like_GSE52588 import config
    start = timeit.default_timer()
    X = np.zeros(shape = (config.params["num_genes"].value, config.params["num_samples"].value), dtype = np.float32)
    masks = ["mongoloids_mask", "siblings_mask", "mothers_mask"]
    for j in range(3):
        mask = config.params[masks[j]].value
        for i in range(config.params["num_genes"].value):
            loc = np.mean(Xw[mask, i])
            scale = np.std(Xw[mask, i])
            x = np.random.normal(loc = loc, scale = scale, size = (1, len(mask)))
            X[i, mask] = np.clip(x, 0, 1).astype(np.float32)
    
    stop = timeit.default_timer()
    print('Data loaded: ', stop - start)
    print(X.dtype, X.shape)

    sys.stdout.flush()
    X = X.T
    y, mask = get_classes(config, X)

    print(X.shape, config.params["num_genes"].value)
    sys.stdout.flush()

    return X, y, mask, genes_names

