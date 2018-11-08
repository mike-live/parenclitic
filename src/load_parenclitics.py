import pandas as pd
def load_parenclitics(config):
    parenclitics = []
    for j, thr_p in enumerate(config.params["thr_p"]):
        cur = pd.read_pickle(config.ofname(["parenclitic"], ext = ".pkl", include_set = config.params_sets["parenclitic"]))
        parenclitics.append(cur)
    return parenclitics