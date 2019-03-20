import pandas as pd
import os
def load_parenclitics(config, by_sample = False, id_thr = 0):
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