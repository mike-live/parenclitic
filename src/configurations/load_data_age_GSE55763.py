import timeit
import numpy as np
from sklearn.model_selection import train_test_split
import sys
import pandas as pd

def get_classes(config, X):
    patients_info = pd.read_csv(config.ifname("patients_info"), sep = ',')
    age = patients_info['age']
    if "num_groups" in config.params:
        num_groups = config.params["num_groups"].value
        bins = np.percentile(age, np.linspace(0, 100, num_groups + 1))
        
    elif "age_delimiter" in config.params:
        age_delimiter = config.params["age_delimiter"].value
        num_groups = 2
        bins = np.array([0, age_delimiter, 100])
        
    y = (np.minimum(np.digitize(age, bins), num_groups) - 1)
    config.params["young_mask"].value = np.flatnonzero(y == 0)
    config.params["old_mask"].value = np.flatnonzero(y == 1)
    print('young:', sum(y == 0), 'old:', sum(y == 1))
    #config.params["kde_mask"].value = "young_mask"
    mask = np.array(y)
    mask[y == 0] = -1
    mask[y == 1] = +1
    return y, mask, age


def load_data_age_GSE55763_cpgs(is_small = True):
    from configurations.config_age_GSE55763_cpg import config
    start = timeit.default_timer()
    data = np.load(config.ifname("cpgs"))
    X = data['data']
    cpgs_names = pd.read_csv(config.ifname("cpgs_names"))["nameCpG"].values
    print('cpg_names', cpgs_names.shape, cpgs_names[0])
    
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
        X = X[indices, :]

    X = X.T
    
    good_cpgs = ~np.isnan(X).any(axis=0)
    X = X[:, good_cpgs]
    cpgs_names = cpgs_names[good_cpgs]
    print(X.shape)
    
    y, mask, age = get_classes(config, X)
    
    config.params["num_cpgs"].value = min(X.shape[1], config.params["num_cpgs"].value)
    
    return X, y, mask, cpgs_names, age    
