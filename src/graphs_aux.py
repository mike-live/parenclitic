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

def make_graph(edges = None, weights = None, G = None, features_names = None, num_vertices = None):
    if not features_name is None:
        num_vertices = len(features_name)
    if not edges is None:
        if num_vertices is None:
            num_vertices = edges.max() + 1
        edges = edges[weights > 0]
        weights = weights[weights > 0]
        
        g = igraph.Graph(n = num_vertices, edges = zip(*edges.T))
        g.es["weight"] = weights
    elif not G is None:
        if num_vertices is None:
            num_vertices = G.shape[0]
        G = np.unpackbits(G, axis = 1)[:, :num_vertices].astype(np.bool)
        g = igraph.Graph.Weighted_Adjacency(G.tolist(), mode=igraph.ADJ_UNDIRECTED)
    else:
        return None
        
    if not features_names is None:
        g.vs["name"] = features_names
        g.vs["label"] = features_names
    return g


def load_graph(config, features_names = None, num_vertices = None, id_thr = 0, id_sample = 0):
    data = get_graph_file(config, id_thr, id_sample)
    if 'IDS' in data:
        g = make_graph(edges = data['IDS'], weights = data['D'], features_names = features_names, num_vertices = num_vertices)
    else:
        g = make_graph(G = data['G'], features_names = features_names, num_vertices = num_vertices)
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
        