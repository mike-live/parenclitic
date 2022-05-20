def load_graphs(config, all_feature_names):
    from parenclitic import parenclitic
    clf = parenclitic()
    paths = []
    for id_sample in config.params["id_sample"]:
        paths.append(config.ofname(["graphs", "g"], ext = ".npz", include_set = config.params_sets["graph"]))
    #print(len(paths), config.params["num_samples"].value)
    clf.set_num_samples(config.params["num_samples"].value)
    clf.set_graph_paths(paths)
    gs = clf.get_graphs(features_names = all_feature_names)
    return gs