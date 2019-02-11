import pandas as pd
def load_parenclitics(config, by_sample = False):
    if by_sample:
        parenclitics = None
        for id_sample in config.params["id_sample"]:
            parenclitic = pd.read_pickle(config.ofname(["parenclitic"], ext = ".pkl", include_set = config.params_sets["parenclitic_sample"]), compression = "gzip")
            if parenclitics is None:
                parenclitics = parenclitic
            else:
                parenclitics = pd.concat([parenclitics, parenclitic], ignore_index=True)
    else:
        parenclitics = []
        if "thr_p" in config.params:
            for j, thr_p in enumerate(config.params["thr_p"]):
                cur = pd.read_pickle(config.ofname(["parenclitic"], ext = ".pkl", include_set = config.params_sets["parenclitic"]))
                parenclitics.append(cur)
        else:
            parenclitics = pd.read_pickle(config.ofname(["parenclitic"], ext = ".pkl", include_set = config.params_sets["parenclitic"]))
    return parenclitics