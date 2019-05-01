from threading import Event, Lock, Semaphore
from multiprocessing import Pool
import timeit
import numpy as np
from .make_graphs_parts import *
import sys
from .infrastructure.auxillary import *
import time

#from ages_config import config
#from load_data_age import load_data_age

#from configurations.load_data_mongoloids import load_data_mongoloids
#from configurations.config_mongoloids import config

#from load_data_mongoloids import load_data_mongoloids_hannum_cpgs
#from mongoloids_cpg_hannum_config import config

from .configurations.load_data_down_GSE52588 import load_data_down_GSE52588_cpgs
from .configurations.config_down_GSE52588_cpg import config

#from configurations.load_data_down_GSE63347 import load_data_down_GSE63347
#from configurations.config_down_GSE63347 import config

#from configurations.load_data_down_GSE63347 import load_data_down_GSE63347_cpg_hannum
#from configurations.config_down_GSE63347_cpg_hannum import config

#from configurations.load_data_down_GSE63347 import load_data_down_GSE63347_cpg_horvath
#from configurations.config_down_GSE63347_cpg_horvath import config

#from configurations.load_data_age_GSE87571 import load_data_age_GSE87571_cpg_horvath
#from configurations.config_age_GSE87571_cpg_horvath import config

#from load_data_cancer import load_data_cancer
#from cancer_config import config

#from configurations.load_data_down_GSE63347 import load_data_down_GSE63347
#from configurations.config_down_GSE63347 import config

#from configurations.load_data_mongoloids import load_data_mongoloids_horvath_cpgs
#from configurations.config_mongoloids_cpg_horvath import config

def reshape_graph(X, id_thr = None, need_G = True):
    G, D, IDS = read_graphs(config, X, need_G, id_thr)
    if need_G:
        G = extract_graphs(G, X.shape[1], X.shape[0])
    for i, x in enumerate(X):
        config.params["id_sample"].set_tick(i)
        fname = config.ofname(["graphs", "g"], ext = ".npz", include_set = config.params_sets["graph"])
        if need_G:
            np.savez_compressed(fname, G = G[i, :, :], D = D[i, :], IDS = IDS)
        else:
            np.savez_compressed(fname, D = D[i, :], IDS = IDS)

def reshape_graphs(X):
    print('Reshape graph')
    sys.stdout.flush()
    start = timeit.default_timer()
    if "thr_p" in config.params:
        for id_thr, thr in enumerate(config.params["thr_p"].get_values()):
            config.params["thr_p"].set_tick(id_thr)
            reshape_graph(X, id_thr)
    else:
        reshape_graph(X, need_G = False)
    stop = timeit.default_timer()
    print('Graph reshaped in ', stop - start)
    sys.stdout.flush()

if __name__ == '__main__':
    #X, y, _, genes_names = load_data_age()
    #X, y, _, genes_names = load_data_cancer()
    X, y, _, genes_names = load_data_down_GSE52588_cpgs()
    
    if "thr_p" in config.params:
        config.params["thr_p"].whole_values = False
    config.params["id_sample"].manual_ticks = True 

    reshape_graphs(X)
