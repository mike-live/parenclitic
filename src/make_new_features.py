from threading import Event, Lock, Semaphore
from multiprocessing import Pool
import timeit
import numpy as np
from .make_graphs_parts import *
import sys
from .infrastructure.auxillary import *
from .infrastructure.configuration import param
import pandas as pd
import time

#from configurations.load_data_down_GSE63347 import load_data_down_GSE63347_cpg_hannum
#from configurations.config_down_GSE63347_cpg_hannum import config

#from configurations.load_data_down_GSE63347 import load_data_down_GSE63347_cpg_horvath
#from configurations.config_down_GSE63347_cpg_horvath import config

from .configurations.load_data_age_GSE87571 import load_data_age_GSE87571_cpg_horvath
from .configurations.config_age_GSE87571_cpg_horvath import config

#from ages_config import config
#from load_data_age import load_data_age

#from configurations.load_data_mongoloids import load_data_mongoloids
#from configurations.config_mongoloids import config

def traverse_graphs(config, X, G, calc_func, upd_func):
    print('Traverse graphs started')
    sys.stdout.flush()
    global done_tasks
    done_tasks = 0
    global ready
    ready = Semaphore(config.params["num_workers"].value * 2)
    for i, x in enumerate(X):
        g = G[i, :, :]
        def upd_callback(res):
            global done_tasks, ready
            upd_func(res)
            done_tasks += 1
            ready.release()
            print('Done tasks:', done_tasks, '/', len(X))
            sys.stdout.flush()

        ready.acquire()
        pool.apply_async(calc_func, args = (x, g, i), callback = upd_callback)
    
    while done_tasks < X.shape[0]:
        ready.acquire()
    print('Traverse graphs finished')
    sys.stdout.flush()

def calc_new_feautures(x, g, i):
    disable_print()
    g = np.unpackbits(g, axis = 1)[:, :x.size]
    res = parenclitic_transform(x, G = g)
    enable_print()
    res['id_sample'] = i
    return res

def make_new_features_graph(X, id_thr = None):
    #disable_print()
    parenclitic = [pd.DataFrame()]
    G = read_graphs(config, X, id_thr)
    G = extract_graphs(G, config.params["num_features"].value, X.shape[0])
    print(G.shape, G.dtype)
    sys.stdout.flush()

    def upd_new_features(res, id_thr = id_thr):
        parenclitic[0] = parenclitic[0].append(res, ignore_index=True)

    traverse_graphs(config, X, G, calc_new_feautures, upd_new_features)
    #enable_print()
    parenclitic[0].set_index('id_sample', inplace=True)
    parenclitic[0].sort_index(inplace=True)
    return parenclitic[0]
    
    
def make_new_features(X):
    #degrees = np.zeros((config.params["thr_p"].num_ticks, X.shape[0], X.shape[1]), dtype = np.int32)
    #parenclitic = np.zeros((config.params["thr_p"].num_ticks, X.shape[0], parenclitic_num_features()), dtype = np.float32)
    print('Make new features')
    sys.stdout.flush()
    start = timeit.default_timer()

    if "thr_p" in config.params:
        parenclitic = [] * config.params["thr_p"].num_ticks
        for id_thr, thr in enumerate(config.params["thr_p"].get_values()):
            config.params["thr_p"].set_tick(id_thr)
            start_thr = timeit.default_timer()
            parenclitic[id_thr] = make_new_features_graph(X, id_thr)
            print('Threshold:', id_thr, thr, stop_thr - start_thr)
            stop_thr = timeit.default_timer()
            sys.stdout.flush()
    else:
        parenclitic = make_new_features_graph(X)

    stop = timeit.default_timer()
    print('New features were calculated in ', stop - start)
    sys.stdout.flush()
    return parenclitic

if __name__ == '__main__':
    global pool
    pool = Pool(config.params["num_workers"].value)

    X, y, _, features_names = load_data_age_GSE87571_cpg_horvath()

    config.params["num_features"] = param(features_names.size, name = 'num_features')

    if "thr_p" in config.params:
        config.params["thr_p"].whole_values = False
    parenclitic = make_new_features(X)

    if "thr_p" in config.params:
        for id_thr, thr in enumerate(config.params["thr_p"].get_values()):
            config.params["thr_p"].set_tick(id_thr)
            config.save_params(include_set = config.params_sets["parenclitic"])
            parenclitic[id_thr].to_pickle(config.ofname(["parenclitic"], ext = ".pkl", include_set = config.params_sets["parenclitic"]))
    else:
        parenclitic.to_pickle(config.ofname(["parenclitic"], ext = ".pkl", include_set = config.params_sets["parenclitic"]))
    