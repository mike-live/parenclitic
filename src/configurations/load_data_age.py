import timeit
import numpy as np
from ages_config import config
from sklearn.model_selection import train_test_split
import sys

def load_data_age():
    start = timeit.default_timer()
    X = np.genfromtxt(config.ifname("x"), dtype='float32', delimiter=' ')[:, 1:]

    genes_names = np.genfromtxt(config.ifname("x"), dtype='str', usecols = 0)
    config.params["num_genes"].value = min(genes_names.size, config.params["num_genes"].value)

    genes_dict = dict((v, i) for i, v in enumerate(genes_names))
    y = np.loadtxt(config.ifname("y"), dtype='float32')

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

    print(X.shape, config.params["num_genes"].value)

    median_age = np.median(y)
    print(median_age)
    min_age = np.min(y)
    max_age = np.max(y)

    ages_edges = [min_age, median_age, max_age]

    X_prob, _, y_prob, _ = train_test_split(
        X, y, test_size=0.9, random_state=42)

    mask = y_prob < ages_edges[1]
    sys.stdout.flush()
    return X, y, X_prob[mask, :], genes_names