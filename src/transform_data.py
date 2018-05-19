import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

import sys
import igraph
import scipy.stats as stats
import scipy.integrate as integrate
import numpy as np
import timeit

def calc_edge_weight(feature_i, feature_j, kde, p, I, min_x, max_x):
#    start = timeit.default_timer()
    
    q = kde([feature_i, feature_j])
    pos = np.searchsorted(p, q)
    w = (I[pos] + q * (p.size - pos)) / p.size
#    print p
#    w = np.random.rand()
#    stop = timeit.default_timer()
#    print 'Time', stop - start 
#    w2, err = integrate.dblquad(lambda y, x: min(kde([x, y]), q), min_x, max_x, lambda x: min_x, lambda x: max_x)
#    w2, err = integrate.nquad(lambda y, x: q, [[min_x, max_x], [min_x, max_x]])
#    print w, w2, abs(w - w2)
    return w


def make_weights(x, kdes, p, I, min_x, max_x):
    w = np.zeros((x.size, x.size), dtype = np.float32)
    for i, feature_i in enumerate(x):
        for j, feature_j in enumerate(x):
            if (i <= j): continue
            w[i, j] = calc_edge_weight(feature_i, feature_j, kdes[i, j], p[i, j], I[i, j], min_x, max_x)
            w[j, i] = w[i, j]
    return w


def make_graph(w, threshold_p):
    n = w.shape[0]
    g = igraph.Graph(n = n, directed = False)
    for i in range(n):
        for j in range(n):
            if i <= j: continue
            if w[i, j] < threshold_p:
                 g.add_edge(i, j, weight = w[i, j])
    return g


def robustness(g):
    cnt = 0
    while g.vcount() > 0:
        degrees = np.array(g.degree())
        max_degree = np.max(degrees)
        if max_degree == 0:
            break
        g.delete_vertices(np.flatnonzero(max_degree == degrees))
        cnt = cnt + 1
    return cnt


def calculate_metrics(g, w, need_weights = False):
    if need_weights:
        weights = 'weight'
    else:
        weights = None            
    parenclitic = np.zeros(parenclitic_num_features(), dtype = np.float32)

    start = timeit.default_timer()

    degrees = np.array(g.vs.degree())
    parenclitic[0] = np.max (degrees)
    parenclitic[1] = np.mean(degrees)
    parenclitic[2] = np.std (degrees)
    degrees = None

    stop = timeit.default_timer()
    print 'Parenclitic 1', stop - start
    
    start = timeit.default_timer()

    shortest_paths = np.array(g.shortest_paths(weights = weights))
    shortest_paths = shortest_paths[(shortest_paths > 0) & (shortest_paths != np.inf)]
    efficency = 0
    if len(shortest_paths) > 0:
        efficency = (1.0 / shortest_paths).sum() / (g.vcount() * (g.vcount() - 1))
    # In paper Latora V., Marchiori M.: Efficient behavior of small-world networks. Phys. Rev. Lett. 87 (Article No. 198701) (2001)
    # suggested to normalize by efficiency for threshold_p = 0 (cause graph has all edges when thr_p = 0)
    parenclitic[3] = efficency
    shortest_paths = None

    stop = timeit.default_timer()
    print 'Parenclitic 2', stop - start
    
    start = timeit.default_timer()


    betweenness = g.betweenness(weights = weights)
    parenclitic[4] = np.max (betweenness)
    parenclitic[5] = np.mean(betweenness)
    parenclitic[6] = np.std (betweenness)
    betweenness = None

    stop = timeit.default_timer()
    print 'Parenclitic 3', stop - start
    
    start = timeit.default_timer()


    closeness = g.closeness(weights = weights)
    parenclitic[7] = np.max (closeness)
    parenclitic[8] = np.mean(closeness)
    parenclitic[9] = np.std (closeness)
    closeness = None

    stop = timeit.default_timer()
    print 'Parenclitic 4', stop - start
    
    start = timeit.default_timer()


    pagerank = g.pagerank(weights = weights)
    parenclitic[10] = np.max (pagerank)
    parenclitic[11] = np.mean(pagerank)
    parenclitic[12] = np.std (pagerank)
    pagerank = None

    stop = timeit.default_timer()
    print 'Parenclitic 5', stop - start
    
    start = timeit.default_timer()

    # alpha centrality with alpha = 1
    eigenvector_centrality = g.eigenvector_centrality(weights = weights)
    parenclitic[13] = np.mean(eigenvector_centrality)
    eigenvector_centrality = None

    stop = timeit.default_timer()
    print 'Parenclitic 6', stop - start
    
    start = timeit.default_timer()

    parenclitic[14] = g.ecount()

    if g.ecount() > 0:
        weights = np.array(g.es["weight"])
        parenclitic[15] = np.max (weights)
        parenclitic[16] = np.mean(weights)
        parenclitic[17] = np.std (weights)
        weights = None

    stop = timeit.default_timer()
    print 'Parenclitic 7', stop - start
    
    start = timeit.default_timer()
    parenclitic[18] = g.community_edge_betweenness().optimal_count
    stop = timeit.default_timer()
    print 'Parenclitic 8', stop - start
    
    start = timeit.default_timer()
    parenclitic[19] = robustness(g)
    stop = timeit.default_timer()
    print 'Parenclitic 9', stop - start
    
    assert(parenclitic.size == parenclitic_num_features())
    return parenclitic

def parenclitic_num_features():
    return 20


def parenclitic_transform(x, kdes, p, I, threshold_p = 0.5, min_x = 0, max_x = 1, graph_path = '', id_patient = -1, genes_names = []):
    start = timeit.default_timer()
    w = make_weights(x, kdes, p, I, min_x, max_x)
    stop = timeit.default_timer()
    print 'Time weights', stop - start

    start = timeit.default_timer()
    g = make_graph(w, threshold_p)
    stop = timeit.default_timer()
    print 'Make graph', stop - start

    if False and graph_path != '':
        g.vs["label"] = genes_names
        if g.ecount() > 0:
            g.es["label"] = g.es["weight"]
        layout = g.layout("fr")
        igraph.plot(g, graph_path, bbox = (1024, 1024), layout = layout, vertex_size = 20)

    parenclitic = calculate_metrics(g, w)
    return parenclitic


def parenclitic_kdes(X, min_x = 0, max_x = 1):
    k = X.shape[1]
    kdes = np.empty((k, k), dtype=object)
    print kdes.shape
    num_points = 10000
    p = np.zeros((k, k, num_points), dtype=np.float32)
    I = np.zeros((k, k, num_points + 1), dtype=np.float32)
    for i in range(k):
        start = timeit.default_timer()
        for j in range(k):
            if (i == j): continue
            data = np.array([X[:, i], X[:, j]])
            kde = stats.gaussian_kde(data)
            points = kde.resample(num_points)
            p[i, j] = np.sort(np.array(kde(points)))
            I[i, j] = np.concatenate([[0], np.cumsum(p[i, j])])
            kdes[i, j] = kde
        stop = timeit.default_timer()
        print 'KDE for ', i, 'calculated in ', stop - start
        sys.stdout.flush()

    return kdes, p, I