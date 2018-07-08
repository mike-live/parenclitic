from threading import Event, Lock, Semaphore
from multiprocessing import Pool
import timeit
import numpy as np
from make_graphs_parts import *
import sys
from auxillary import *
import time

from ages_config import config
from load_data_age import load_data_age

#from load_data_mongoloids import load_data_mongoloids
#from mongoloids_config import config

def traverse_graphs(config, X, G, calc_func, upd_func):
    print 'Traverse graphs started'
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
            print 'Done tasks:', done_tasks, '/', len(X)
            sys.stdout.flush()

        ready.acquire()
        pool.apply_async(calc_func, args = (x, g, i), callback = upd_callback)
    
    while done_tasks < X.shape[0]:
        ready.acquire()
    print 'Traverse graphs finished'
    sys.stdout.flush()

def calc_new_feautures(x, g, i):
    #disable_print()
    g = np.unpackbits(g, axis = 1)[:, :x.size]
    res = (i, get_degrees(g), parenclitic_transform(x, G = g))
    #enable_print()
    return res

def make_new_features(X):
    degrees = np.zeros((config.params["thr_p"].num_ticks, X.shape[0], X.shape[1]), dtype = np.int32)
    parenclitic = np.zeros((config.params["thr_p"].num_ticks, X.shape[0], parenclitic_num_features()), dtype = np.float32)
    print 'Make new features'
    sys.stdout.flush()
    start = timeit.default_timer()
    for id_thr, thr in enumerate(config.params["thr_p"].get_values()):
        config.params["thr_p"].set_tick(id_thr)
        start_thr = timeit.default_timer()

        #disable_print()
        G = read_graphs(config, X, id_thr)
        G = extract_graphs(G, config.params["num_genes"].value, X.shape[0])
        print G.shape, G.dtype
        sys.stdout.flush()

        def upd_new_features(res, id_thr = id_thr):
            degrees[id_thr, res[0]] = res[1]
            parenclitic[id_thr, res[0]] = res[2]
    
        traverse_graphs(config, X, G, calc_new_feautures, upd_new_features)
        #enable_print()
        stop_thr = timeit.default_timer()
        print 'Threshold:', id_thr, thr, stop_thr - start_thr
        sys.stdout.flush()

    stop = timeit.default_timer()
    print 'New features were calculated in ', stop - start
    sys.stdout.flush()
    return degrees, parenclitic

if __name__ == '__main__':
    global pool
    pool = Pool(config.params["num_workers"].value)

    X, y, _, genes_names = load_data_age()
    #X, y, _, genes_names = load_data_mongoloids()

    config.params["num_genes"].value = min(genes_names.size, config.params["num_genes"].value)
    config.params["thr_p"].whole_values = False
    degrees, parenclitic = make_new_features(X)

    for id_thr, thr in enumerate(config.params["thr_p"].get_values()):
        config.params["thr_p"].set_tick(id_thr)
    
        config.save_params(include_set = config.params_sets["degrees"])
    
        degrees_cur = degrees[id_thr, :, :]
        print degrees_cur.dtype
        np.savetxt(config.ofname(["degrees"], ext = ".txt", include_set = config.params_sets["degrees"]), degrees_cur, delimiter = '\t', fmt = '%6d')
    
        parenclitic_cur = parenclitic[id_thr, :, :]
        np.savetxt(config.ofname(["parenclitic"], ext = ".txt", include_set = config.params_sets["parenclitic"]), parenclitic_cur, delimiter = '\t')
