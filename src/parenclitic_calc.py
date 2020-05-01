import parenclitic
import numpy as np
import sys

#from configurations.load_data_down_GSE52588 import load_data_down_GSE52588
#from configurations.config_down_GSE52588 import config

#from configurations.load_data_age_GSE87571 import load_data_age_GSE87571
#from configurations.config_age_GSE87571 import config

#from configurations.load_data_down_GSE52588 import load_data_down_GSE52588_cpgs
#from configurations.config_down_GSE52588_cpg import config

from configurations.load_data_age_GSE55763 import load_data_age_GSE55763_cpgs
from configurations.config_age_GSE55763_cpg import config


import multiprocessing as mp
import multiprocessing, logging
#mpl = multiprocessing.log_to_stderr()
#mpl.setLevel(logging.DEBUG)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        num_workers = int(sys.argv[1])
    else:
        num_workers = 2
    print('Start with', num_workers, 'workers')
    #X, y, mask, all_features_names = load_data_down_GSE52588()
    #X, y, mask, all_features_names, age = load_data_age_GSE87571()
    #X, y, mask, all_features_names = load_data_down_GSE52588_cpgs(True)
    X, y, mask, all_features_names, age = load_data_age_GSE55763_cpgs(True)

    print('Graph path', config.ofname(["graphs", "g"], ext = ".tsv", include_set = config.params_sets["graph"]))

    import time
    num_samples = X.shape[0]
    num_features = X.shape[1]
    
    #by_group = config.params["by_group"].value
    max_score_1d = config.params["max_score_1d"].value
    min_score = config.params["min_score"].value
    thr_type = config.params["thr_type"].value
    division_rule = config.params["division_rule"].value
    kernel = parenclitic.pdf_kernel(thr_type = thr_type, min_score = min_score, division_rule = division_rule)
    pair_filter = parenclitic.IG_filter(max_score = max_score_1d)
    clf = parenclitic.parenclitic(kernel = kernel, pair_filter = pair_filter, verbose = 0)
    # partition = parenclitic.graph_partition_subset()
    #clf = parenclitic.parenclitic(kernel = parenclitic.classifier_kernel(min_score = min_score, by_group = by_group)) # partition = parenclitic.graph_partition_subset()
    #clf = parenclitic.parenclitic(kernel = parenclitic.pdf_kernel(thr_p = 0.88))#, partition = parenclitic.graph_partition_subset())
    be = time.time()
    #clf.fit(X, y, mask, num_workers = 5)
    clf.fit(X[:, :], y[:], mask[:], num_workers = num_workers, chunk_size = 1000)
    en = time.time()
    print(en - be)
    #clf.calc_parenclitic()
    #gr = clf.get_graphs(features_names = features_names)
    paths = []
    for id_sample in config.params["id_sample"]:
        paths.append(config.ofname(["graphs", "g"], ext = ".tsv", include_set = config.params_sets["graph"]))
    clf.set_graph_paths(paths = paths[:])
    clf.save_graphs(gtype = 'csv')
    
    paths = []
    for id_sample in config.params["id_sample"]:
        paths.append(config.ofname(["graphs", "g"], ext = ".npz", include_set = config.params_sets["graph"]))
    clf.set_graph_paths(paths = paths[:])
    clf.save_graphs(gtype = 'npz')