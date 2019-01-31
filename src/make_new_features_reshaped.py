from threading import Event, Lock, Semaphore
import timeit
import numpy as np
from transform_data import *
import sys
from infrastructure.auxillary import *
import time
import fcntl
import os

#from ages_config import config
#from load_data_age import load_data_age

#from load_data_mongoloids import load_data_mongoloids
#from mongoloids_config import config

from configurations.load_data_down_GSE63347 import load_data_down_GSE63347
from configurations.config_down_GSE63347 import config

def calc_new_feautures(x, g, i):
    #disable_print()
    g = np.unpackbits(g, axis = 1)[:, :x.size]
    res = (i, parenclitic_transform(x, G = g))
    #enable_print()
    return res

def make_new_features(X):
    i = config.params["id_sample"].value
    thr_p = -1
    if "thr_p" in config.params:
        thr_p = config.params["thr_p"].value

    print 'Make new features for', i, '/', X.shape[0], 'and thr_p', thr_p
    sys.stdout.flush()
    start = timeit.default_timer()

    data = np.load(config.ofname(["graphs", "g"], ext = ".npz", include_set = config.params_sets["graph"]))
    g = data['G']
    print g.shape
    res = calc_new_feautures(X[i, :], g, i)

    stop = timeit.default_timer()
    print 'New features were calculated in', stop - start, 'for', i, '/', X.shape[0], 'and thr_p', thr_p
    sys.stdout.flush()

    return res
    

if __name__ == '__main__':
    config.params["id_part"].manual_ticks = True

    config.upd_ticks()

    i = config.params["id_sample"].value
    thr_p = -1
    if "thr_p" in config.params:
        thr_p = config.params["thr_p"].value

    print 'Start', i, '/', config.params['num_samples'].value, 'and thr_p', thr_p

    #if config.info["task_id"] in set([677, 676, 675]):
    #    sys.exit(0)

    degrees_path = config.ofname(["degrees"], ext = ".txt", include_set = config.params_sets["degrees_sample"])
    parenclitic_path = config.ofname(["parenclitic"], ext = ".txt", include_set = config.params_sets["parenclitic_sample"])
    if os.path.exists(degrees_path) and os.path.exists(parenclitic_path):
        sys.exit(0)

    X, y, _, genes_names = load_data_down_GSE63347()    
    #X, y, _, genes_names = load_data_age()
    config.params["num_genes"].value = min(genes_names.size, config.params["num_genes"].value)

    id_sample, parenclitic = make_new_features(X)

    config.save_params(include_set = config.params_sets["degrees"])
    
    '''
    def sync_save(filename, writer):
        f = os.open(filename,os.O_WRONLY | os.O_APPEND)
        fcntl.lockf(f,fcntl.LOCK_EX)

        writer(f)
        f.write('\n')

        fcntl.lockf(f,fcntl.LOCK_UN)
        os.close(f)
    '''
    #np.savetxt(degrees_path, degrees, delimiter = '\t', fmt = '%6d')
    #np.savetxt(parenclitic_path, parenclitic, delimiter = '\t')
    parenclitic.to_pickle(config.ofname(["parenclitic"], ext = ".pkl", include_set = config.params_sets["parenclitic_sample"]), compression = "gzip")
