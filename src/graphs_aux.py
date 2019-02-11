import numpy as np
import igraph
import scipy
import struct

def get_graph_file(config, id_thr = 0, id_sample = 0):
    config.params["id_sample"].set_tick(id_sample)
    if "thr_p" in config.params:
        config.params["thr_p"].set_tick(id_thr)
    fname = config.ofname(["graphs", "g"], ext = ".npz", include_set = config.params_sets["graph"])
    print fname
    data = np.load(fname)
    return data
    
def load_graph(config, X, features_names, id_thr = 0, id_sample = 0):
    data = get_graph_file(config, id_thr, id_sample)
    if 'IDS' in data:
        d = data['D']
        ids = data['IDS']
        
        ids = ids[d > 0]
        d = d[d > 0]
        
        g = igraph.Graph(n = X.shape[1], edges = zip(*ids.T))
        g.es["weight"] = d
    else:
        g = data['G']
        g = np.unpackbits(g, axis = 1)[:, :X.shape[1]].astype(np.bool)
        g = igraph.Graph.Weighted_Adjacency(g.tolist(), mode=igraph.ADJ_UNDIRECTED)
    
    g.vs["name"] = features_names
    g.vs["label"] = features_names
    return g

def graph_to_crs(g):
    m = g.get_adjacency()
    m = scipy.sparse.csr_matrix(m)
    return m

def save_crs(file_name, m):
    with open(file_name, 'wb') as f:
        n = int(m.shape[0])
        f.write(struct.pack('>i', n))
        f.write(m.nnz.to_bytes(4, byteorder='big', signed=True))
        f.write(m.indices.to_bytes(4, byteorder='big', signed=True))
        f.write(m.indptr.to_bytes(4, byteorder='big', signed=True))
        