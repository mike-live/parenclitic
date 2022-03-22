import pandas as pd
from graphs_aux import graph_to_dataframe
from tqdm.notebook import tqdm
import numpy as np
import igraph

def export_best_cpg_pairs(config, X, gc, case_name, path = []):
    df_graph = graph_to_dataframe(gc)
    best_edges = df_graph.iloc[df_graph[case_name + "_percent"].argsort()[-3:]]
    for index, edge in best_edges.iterrows():
        feature_1, feature_2 = edge["feature_1"], edge["feature_2"]
        v1, v2 = edge["v1"], edge["v2"]
        x1, x2 = X[:, v1], X[:, v2]
        df = pd.DataFrame({feature_1: x1, feature_2: x2})
        feature_subset_path = config.ofname(path + [[feature_1, feature_2]], ext = ".tsv", include_set = config.params_sets["diff_graph"])
        #print(feature_subset_path)
        df.to_csv(feature_subset_path, sep = '\t', index = False)
        #break
    return

def save_plot_graph(gc, graph_name, X, group_masks, group_names, all_feature_names, config, case_name = None, 
                    is_plot_graph = True, is_print_name = True, remove_isolated = True):
    fig_graph_path = config.ofname(graph_name, ext = ".png", include_set = config.params_sets["diff_graph"])
    tsv_graph_path = config.ofname(graph_name, ext = ".tsv", include_set = config.params_sets["diff_graph"])
    from copy import deepcopy
    export_cpg_path = deepcopy(graph_name)
    def add_graph_name_suffix(graph_name, suffix):
        if type(graph_name) is str:
            graph_name += "_" + suffix
        elif type(graph_name[-1]) is str:
            graph_name[-1] += "_" + suffix
        else:
            graph_name[-1][-1] += "_" + suffix
        return graph_name
    
    tsv_vertices_path = config.ofname(add_graph_name_suffix(graph_name, 'vertices'), 
                                      ext = ".tsv", include_set = config.params_sets["diff_graph"])
    tsv_list_path = config.ofname(add_graph_name_suffix(graph_name, 'list'), 
                                  ext = ".tsv", include_set = config.params_sets["diff_graph"])

    g = gc.copy()
    g.vs["label"] = all_feature_names
    
    if not case_name is None:
        export_best_cpg_pairs(g, case_name, path = export_cpg_path)
    
    if remove_isolated:
        g.delete_vertices([v.index for v in g.vs if v.degree() < 1])
    
    if is_plot_graph:
        layout = g.layout("kk")
        #layout = g.layout("circle")
        #igraph.plot(g, layout = layout)

        visual_style = {}
        visual_style["font_size"] = 5
        visual_style["vertex_label_size"] = 10
        #visual_style["vertex_color"] = [color_dict[gender] for gender in g.vs["gender"]]
        visual_style["vertex_label"] = g.vs["label"]
        #visual_style["edge_width"] = [1 + abs(weight) for weight in g.es["weight"]]
        visual_style["layout"] = layout
        visual_style["bbox"] = (700, 700)
        visual_style["margin"] = 50
        visual_style["vertex_label_dist"] = 1.1
        visual_style["vertex_shape"] = "circle"

        fig_g = igraph.plot(g, **visual_style)

        fig_g.save(fig_graph_path)
    
    #os.startfile(fig_graph_path)
    name = "cpg" if "num_cpgs" in config.params else "gene"
    df_graph = graph_to_dataframe(g, name)
    #print(df_graph)
    df_graph.to_csv(tsv_graph_path, sep = '\t')
    feature_list = list(set(df_graph[name + "_1"].values.tolist() + df_graph[name + "_2"].values.tolist()))
    _, indices, _ = np.intersect1d(all_feature_names, feature_list, return_indices = True)
    feature_list = all_feature_names[indices]
        
    def add_statistics(d, X, group_masks, group_names):
        _, indices, _ = np.intersect1d(all_feature_names, d[name], return_indices = True)
        if len(indices) > 0:
            for group_mask, group_name in zip(group_masks, group_names):
                scur = X[group_mask]
                d[group_name + '_mean'] = scur[:, indices].mean(axis = 0)
            for group_mask, group_name in zip(group_masks, group_names):
                scur = X[group_mask]
                d[group_name + '_std'] = scur[:, indices].std(axis = 0)
        return d
    
    d = {name: feature_list}
    d = add_statistics(d, X, group_masks, group_names)
    df_vlist = pd.DataFrame(d)
    df_vlist.to_csv(tsv_list_path, sep = '\t')
    
    d = {name: g.vs["label"]}
    d = add_statistics(d, X, group_masks, group_names)
    df_vertices = pd.DataFrame(d)
    df_vertices.to_csv(tsv_vertices_path, sep = '\t')
    
    if is_print_name:
        print(graph_name)

def export_all(config, gs, gss, g_diff, g_and, g_or, X, group_masks, group_names, all_feature_names, path = [["final_graph"]], \
               is_plot_graph = False, is_plot_groups = False, is_plot_subjects = True):
    kwargs = {
        'config': config,
        'X': X, 
        'group_masks': group_masks, 
        'group_names': group_names,
        'all_feature_names': all_feature_names,
        'is_plot_graph': is_plot_graph
    }
    
    if is_plot_groups:
        save_plot_graph(g_diff, path + [["diff_graph"]], **kwargs)

        save_plot_graph(g_and[0], path + [["ds_intersection_graph"]], **kwargs)
        save_plot_graph(g_or[0], path + [["ds_union_graph"]], case_name = "ds",  **kwargs)
        save_plot_graph(g_and[1], path + [["siblings_intersection_graph"]], **kwargs)
        save_plot_graph(g_or[1], path + [["siblings_union_graph"]], case_name = "siblings", **kwargs)
        save_plot_graph(g_and[2], path + [["mothers_intersection_graph"]], **kwargs)
        save_plot_graph(g_or[2], path + [["mothers_union_graph"]], case_name = "mothers", **kwargs)
    
    if is_plot_subjects:
        for id_sample in tqdm(config.params["id_sample"], desc = "Graphs export"):
            save_plot_graph(gs[id_sample], path + [["graphs_export"],["graph_id_sample_" + str(id_sample)]], \
                            **kwargs, is_print_name = False)
            if gss:
                save_plot_graph(gss[id_sample], path + [["graphs_sieved_export"],["graph_id_sample_" + str(id_sample)]], \
                                **kwargs, is_print_name = False, remove_isolated = False)
        
    return
