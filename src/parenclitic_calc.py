import parenclitic
import numpy as np
import sys

#from configurations.load_data_down_GSE52588 import load_data_down_GSE52588
#from configurations.config_down_GSE52588 import config

#from configurations.load_data_age_GSE87571 import load_data_age_GSE87571
#from configurations.config_age_GSE87571 import config

from configurations.load_data_down_GSE52588 import load_data_down_GSE52588_cpgs
from configurations.config_down_GSE52588_cpg import config

#from configurations.load_data_age_GSE55763 import load_data_age_GSE55763_cpgs
#from configurations.config_age_GSE55763_cpg import config


import multiprocessing as mp
import multiprocessing, logging
from tqdm import tqdm
#mpl = multiprocessing.log_to_stderr()
#mpl.setLevel(logging.DEBUG)

def build_parenclitic(config, X, y, mask, all_features_names, num_workers, subset = None):
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
    #clf = parenclitic.parenclitic(kernel = kernel, pair_filter = pair_filter, verbose = 0)
    if subset is None:
        partition = parenclitic.graph_partition()
    else:    
        partition = parenclitic.graph_partition_subset()
    clf = parenclitic.parenclitic(kernel = kernel, pair_filter = pair_filter, verbose = 0, partition = partition)

    #clf = parenclitic.parenclitic(kernel = parenclitic.classifier_kernel(min_score = min_score, by_group = by_group)) # partition = parenclitic.graph_partition_subset()
    #clf = parenclitic.parenclitic(kernel = parenclitic.pdf_kernel(thr_p = 0.88))#, partition = parenclitic.graph_partition_subset())
    be = time.time()
    #clf.fit(X, y, mask, num_workers = 5)
    clf.fit(X[:, :], y[:], mask[:], num_workers = num_workers, chunk_size = 1000, subset = subset)
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

def get_edges_subset(config):
    file_name = config.ofname(["graphs", "g"], ext = ".npz", include_set = config.params_sets["graph"])    
    print(file_name)
    data = np.load(file_name)
    print(np.array(data['E']).shape)
    return data['E']

if __name__ == "__main__":
    if len(sys.argv) > 1:
        num_workers = int(sys.argv[1])
    else:
        num_workers = 2
    print('Start with', num_workers, 'workers')
    #X, y, mask, all_features_names = load_data_down_GSE52588()
    #X, y, mask, all_features_names, age = load_data_age_GSE87571()
    X, y, mask, all_features_names = load_data_down_GSE52588_cpgs(True)
    #X, y, mask, all_features_names, age = load_data_age_GSE55763_cpgs(True)
    if "LOO" in config.params:
        for id_leave in tqdm(config.params["LOO"]):
            if id_leave < 11:
                continue
            maskc = mask.copy()
            for i in [0, 29, 29 * 2]:
                value = mask[id_leave + i]
                maskc[id_leave + i] = {-2: -2, -1: -2, 0: 0, 1: 2, 2: 2}[value]
            print(id_leave)
            print(maskc)
            #subset = get_edges_subset(config)
            subset = None
            build_parenclitic(config, X, y, maskc, all_features_names, num_workers, subset = subset)
    else:
        build_parenclitic(config, X, y, mask, all_features_names, num_workers)
