from tqdm.notebook import tqdm

def sieve_graphs(config, gs, g_or, all_feature_names, get_4network):
    gss = []
    for g in tqdm(gs, desc = 'Sieve'):
        gc = g.copy()
        gc.vs['label'] = all_feature_names

        gn = g_or[get_4network(config)].copy()
        gn.vs['label'] = all_feature_names

        gc.delete_vertices([v.index for v in gn.vs if v.degree() < 1])
        gc.vs["label"]

        gss += [gc]
    return gss