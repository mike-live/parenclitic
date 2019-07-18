import numpy as np   
import pandas as pd
import os
import scipy 

def load_parenclitics(config, by_sample = False, id_thr = 0, from_mat = False):
    if from_mat:
        fname = config.ofname(["parenclitic"], ext = ".mat", include_set = config.params_sets["parenclitic"])
        print(fname)
        data = scipy.io.loadmat(fname)
        parenclitics = data['parenclitics']
        parenclitic_features = data['parenclitic_features']

        parenclitics = list(map(lambda x: tuple(map(lambda y: y.item() if y.size == 1 else y.flatten(), tuple(x[0]))), parenclitics))
        parenclitic_features = list(map(lambda x: tuple(x[0])[0], parenclitic_features))
        parenclitics = pd.DataFrame(parenclitics, columns = parenclitic_features)
        
    else:
        if by_sample:
            if "thr_p" in config.params:
                config.params["thr_p"].set_tick(id_thr)
    
            parenclitics = None
            for id_sample in config.params["id_sample"]:
                parenclitic = None
                pkl_file = config.ofname(["parenclitic"], ext = ".pkl", include_set = config.params_sets["parenclitic_sample"])
                if os.path.exists(pkl_file):
                    parenclitic = pd.read_pickle(pkl_file, compression = "gzip")
                
                txt_file = config.ofname(["parenclitic"], ext = ".txt", include_set = config.params_sets["parenclitic_sample"])
                if os.path.exists(txt_file):
                    parenclitic = pd.read_csv(txt_file, sep=" ", header=None).T
    
                if parenclitics is None:
                    parenclitics = parenclitic
                else:
                    parenclitics = pd.concat([parenclitics, parenclitic], ignore_index=True, sort=True)
        else:
            parenclitics = []
            if "thr_p" in config.params:
                for j, thr_p in enumerate(config.params["thr_p"]):
                    cur = pd.read_pickle(config.ofname(["parenclitic"], ext = ".pkl", include_set = config.params_sets["parenclitic"]))
                    parenclitics.append(cur)
            else:
                parenclitics = pd.read_pickle(config.ofname(["parenclitic"], ext = ".pkl", include_set = config.params_sets["parenclitic"]))
    return parenclitics