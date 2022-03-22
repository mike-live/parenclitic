import pandas as pd
from tqdm.notebook import tqdm
import numpy as np
from graphs_aux import graph_to_dataframe

def strip_graph_cpgs(g, cpgs_erase, all_cpg_names):
    _, cpg_ids, _ = np.intersect1d(all_cpg_names, cpgs_erase, return_indices = True)
    erase_edges = [edge for edge in g.es if edge.tuple[0] in cpg_ids or edge.tuple[1] in cpg_ids]
    g.delete_edges(erase_edges)
    return g
    
def strip_graphs(gs, all_cpg_names, cpgs_erase):
    gs_stripped = []
    for g in tqdm(gs, desc = "Graph stripping"):
        gc = strip_graph_cpgs(g.copy(), cpgs_erase, all_cpg_names)
        gs_stripped.append(gc)
    return gs_stripped

def save_striped_graphs(config, gs_stripped):
    for id_sample in tqdm(config.params["id_sample"], desc = "Save stripped graphs"):
        g = gs_stripped[id_sample]
        graph_path = config.ofname(["stripped_graphs", "g"], ext = ".tsv", include_set = config.params_sets["graph"])
        graph_to_dataframe(g).to_csv(graph_path, sep = '\t')
    return

def make_stripped_graphs(config, cpgs, gs, all_feature_names):
    is_cpgs = "num_cpgs" in config.params

    path_erase_cpgs = config.ifname("singmann_sex_related")
    df_cpgs_erase = pd.read_csv(path_erase_cpgs)
    cpgs_erase = df_cpgs_erase["nameCpG"].values

    if is_cpgs:
        features_erase = cpgs_erase
    else:
        genes_erase = cpgs.get_crit_col_values("gene", criterions = {'cpgs_in': cpgs_erase})
        features_erase = genes_erase

    gs_stripped = strip_graphs(gs, all_feature_names, features_erase)
    save_striped_graphs(config, gs_stripped)
    return gs_stripped
