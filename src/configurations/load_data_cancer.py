import timeit
import numpy as np
from sklearn.model_selection import train_test_split
import sys

def load_data_cancer():
    from cancer_config import config
    start = timeit.default_timer()
    X = np.genfromtxt(config.ifname("x"), dtype='float32', delimiter=' ')[0:, 1:]

    genes_names = np.genfromtxt(config.ifname("x"), dtype='str', usecols = 0)[1:]
    config.params["num_genes"].value = min(genes_names.size, config.params["num_genes"].value)
    genes_names[:] = [s.strip('"') for s in genes_names]

    genes_dict = dict((v, i) for i, v in enumerate(genes_names))

    patients_names = np.genfromtxt(config.ifname("patients_id"), dtype='str', usecols = 0)[1:]
    patients_names[:] = [s.strip('"') for s in patients_names]
    patients_dict = dict((v, i) for i, v in enumerate(patients_names))

    #ranged_genes = np.genfromtxt(config.ifname("ranged_genes"), dtype='str', usecols = 0)

    stop = timeit.default_timer()
    print('Data loaded: ', stop - start)
    print(X.dtype, X.shape)

    sys.stdout.flush()
    #X = np.random.rand(len(genes_names), 656)
    #genes_names = ranged_genes[:config.params["num_genes"].value]
    cur = np.genfromtxt(config.ifname("classes_blood"), dtype='str', skip_header = 1)

    patients_names = np.genfromtxt(config.ifname("classes_blood"), dtype='str', skip_header = 1).T[0]
    y_status = np.genfromtxt(config.ifname("classes_blood"), dtype='str', skip_header = 1).T[1]
    immune_frac = np.genfromtxt(config.ifname("classes_blood"), dtype='str', skip_header = 1).T[2]

    patients_names[:] = [s.strip('"') for s in patients_names]
    y_status[:] = [s.strip('"') for s in y_status]
    immune_frac = np.array([float(s.strip('"')) for s in immune_frac], dtype = "float")

    genes_names = genes_names[:config.params["num_genes"].value]

    #indices = np.array([genes_dict[x] for x in genes_names])
    patients_indices = np.array([patients_dict[x] for x in patients_names])

    #X = X[indices, :].T
    X = X.T
    y = ((y_status == "BRCAmut").astype(np.uint8))[patients_indices]
    mask = np.zeros(y.shape, dtype = 'uint8')

    ids = (np.flatnonzero(y == 0))
    np.random.shuffle(ids)
    cnt = int(len(ids) * 0.1)
    mask[ids[:cnt]] = 1

    config.params["control_mask"].value = ids[:cnt]
    config.params["health_mask"].value = ids
    config.params["cancer_mask"].value = np.flatnonzero(y == 1)

    print(X.shape, config.params["num_genes"].value, mask.sum(), mask)
    sys.stdout.flush()
    
    return X, y, X[mask, :], genes_names
