from threading import Event, Lock, Semaphore
from multiprocessing import Pool
import timeit
import numpy as np
from make_graphs_parts import *
import sys
from auxillary import *
import time

#from ages_config import config
#from load_data_age import load_data_age

#from load_data_mongoloids import load_data_mongoloids
#from mongoloids_config import config

from load_data_cancer import load_data_cancer
from cancer_config import config


def reshape_graphs(X):
    print 'Reshape graph'
    sys.stdout.flush()
    start = timeit.default_timer()
    for id_thr, thr in enumerate(config.params["thr_p"].get_values()):
        config.params["thr_p"].set_tick(id_thr)
        start_thr = timeit.default_timer()

        G = read_graphs(config, X, id_thr)
        G = extract_graphs(G, config.params["num_genes"].value, X.shape[0])
        for i, x in enumerate(X):
            config.params["id_sample"].set_tick(i)
            np.savez_compressed(config.ofname(["graphs", "g"], ext = ".npz", include_set = config.params_sets["graph"]), G = G[i, :, :])

    stop = timeit.default_timer()
    print 'Graph reshaped in ', stop - start
    sys.stdout.flush()

if __name__ == '__main__':
    #X, y, _, genes_names = load_data_age()
    X, y, _, genes_names = load_data_cancer()

    config.params["num_genes"].value = min(genes_names.size, config.params["num_genes"].value)
    config.params["thr_p"].whole_values = False
    config.params["id_sample"].manual_ticks = True 

    reshape_graphs(X)
