import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

#from multiprocessing.pool import ThreadPool as Pool
from multiprocessing.pool import Pool
from threading import Event, Lock, Semaphore

import sys
import collections
import igraph
import scipy.stats as stats
import scipy.integrate as integrate
import numpy as np
import pandas as pd
import timeit
from sklearn.neighbors.kde import KernelDensity
from numpy import linalg as LA

feature_names = []

def calc_edge_weight(feature_i, feature_j, kde, p, I, min_x, max_x):
#    start = timeit.default_timer()
    
    q = kde(np.array([feature_i, feature_j])) # .reshape(1, -1)
    pos = np.searchsorted(p, q)
    #w = (I[pos] + q * (p.size - pos)) / p.size
    w = (I[-1] - I[pos]) / p.size
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


def make_graph(G, w):
    g = igraph.Graph.Weighted_Adjacency(w.tolist(), mode=igraph.ADJ_UNDIRECTED)

    print 'make_graph_finished'
    sys.stdout.flush()
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
    print 'Metrics'
    sys.stdout.flush()

    if need_weights:
        weights = 'weight'
    else:
        weights = None

    parenclitic = pd.DataFrame()
    
    start = timeit.default_timer()

    degrees = np.array(g.vs.degree())
    parenclitic['degrees'] = [degrees]
    parenclitic['max_degrees'] = np.max (degrees)
    parenclitic['mean_degrees'] = np.mean(degrees)
    parenclitic['std_degrees'] = np.std (degrees)
    degrees = None

    stop = timeit.default_timer()
    print 'Parenclitic 1', stop - start
    sys.stdout.flush()
        
    start = timeit.default_timer()

    shortest_paths = np.array(g.shortest_paths(weights = weights))
    shortest_paths = shortest_paths[(shortest_paths > 0) & (shortest_paths != np.inf)]
    efficiency = 0
    if len(shortest_paths) > 0:
        efficiency = (1.0 / shortest_paths).sum() / (g.vcount() * (g.vcount() - 1))
    # In paper Latora V., Marchiori M.: Efficient behavior of small-world networks. Phys. Rev. Lett. 87 (Article No. 198701) (2001)
    # suggested to normalize by efficiency for threshold_p = 0 (cause graph has all edges when thr_p = 0)
    parenclitic['efficiency'] = efficiency
    shortest_paths = None

    stop = timeit.default_timer()
    print 'Parenclitic 2', stop - start
    sys.stdout.flush()
        
    
    start = timeit.default_timer()

    betweenness = g.betweenness(weights = weights)
    parenclitic['betweenness'] = [betweenness]
    parenclitic['max_betweenness'] = np.max (betweenness)
    parenclitic['mean_betweenness'] = np.mean(betweenness)
    parenclitic['std_betweenness'] = np.std (betweenness)
    betweenness = None

    stop = timeit.default_timer()
    print 'Parenclitic 3', stop - start
    sys.stdout.flush()
        
    
    start = timeit.default_timer()

    closeness = g.closeness(weights = weights)
    parenclitic['closeness'] = [closeness]
    parenclitic['max_closeness'] = np.max (closeness)
    parenclitic['mean_closeness'] = np.mean(closeness)
    parenclitic['std_closeness'] = np.std (closeness)
    closeness = None

    stop = timeit.default_timer()
    print 'Parenclitic 4', stop - start
    sys.stdout.flush()

    start = timeit.default_timer()

    pagerank = g.pagerank(weights = weights)
    parenclitic['pagerank'] = [pagerank]
    parenclitic['max_pagerank'] = np.max (pagerank)
    parenclitic['mean_pagerank'] = np.mean(pagerank)
    parenclitic['std_pagerank'] = np.std (pagerank)
    pagerank = None

    stop = timeit.default_timer()
    print 'Parenclitic 5', stop - start
    sys.stdout.flush()
        
    start = timeit.default_timer()

    # alpha centrality with alpha = 1
    
    eigenvector_centrality = g.eigenvector_centrality(weights = weights)
    parenclitic['eigenvector_centrality'] = [eigenvector_centrality]
    parenclitic['mean_eigenvector_centrality'] = np.mean(eigenvector_centrality)
    eigenvector_centrality = None

    stop = timeit.default_timer()
    print 'Parenclitic centrality', stop - start
    sys.stdout.flush()
        
    start = timeit.default_timer()

    m = np.array(g.get_adjacency().data)
    sys.stdout.flush()

    eigenvalues, eigenvectors = LA.eig(m)
    print 'Eigenvectors', stop - start
    sys.stdout.flush()

    eigenvalues_intervals = np.diff(np.sort(eigenvalues)) 
    print 'intervals', stop - start
    sys.stdout.flush()

    eigenvalues_intervals_normalized = eigenvalues_intervals / np.mean(eigenvalues_intervals)

    print 'normalized', stop - start
    sys.stdout.flush()
    
    parenclitic['eigenvalues'] = [eigenvalues]
    parenclitic['eigenvalues_intervals'] = [eigenvalues_intervals]
    parenclitic['eigenvalues_intervals_normalized'] = [eigenvalues_intervals_normalized]

    stop = timeit.default_timer()
    print 'Parenclitic: eigenvalues', stop - start
    sys.stdout.flush()
    
    IPR = np.sum(np.power(eigenvectors, 4), axis=1) / np.power(np.sum(np.power(eigenvectors, 2), axis=1), 2)
    parenclitic['IPR'] = [IPR]
    parenclitic['max_IPR'] = np.max(IPR)
    parenclitic['mean_IPR'] = np.mean(IPR)

    stop = timeit.default_timer()
    print 'Parenclitic 7', stop - start
    sys.stdout.flush()
        
    start = timeit.default_timer()

    parenclitic['num_edges'] = g.ecount()

    if g.ecount() > 0:
        weights = np.array(g.es["weight"])
        parenclitic['weights'] = [weights]
        parenclitic['max_weights'] = np.max (weights)
        parenclitic['mean_weights'] = np.mean(weights)
        parenclitic['std_weights'] = np.std (weights)
        weights = None

    stop = timeit.default_timer()
    print 'Parenclitic 8', stop - start
    sys.stdout.flush()

    start = timeit.default_timer()
#    parenclitic['community_edge_betweenness_optimal'] = g.community_edge_betweenness().optimal_count
    stop = timeit.default_timer()
    print 'Parenclitic 9', stop - start
    sys.stdout.flush()
        
    start = timeit.default_timer()
    parenclitic['robustness'] = robustness(g)
    stop = timeit.default_timer()
    print 'Parenclitic 10', stop - start
    sys.stdout.flush()

    return parenclitic

def parenclitic_feature_names():
    feature_names = {}
    feature_names['degrees'] = 'Degrees'
    feature_names['max_degrees'] = 'Max degrees'
    feature_names['mean_degrees'] = 'Mean degrees'
    feature_names['std_degrees'] = 'Std degrees'
    
    feature_names['efficiency'] = 'Efficiency'
    feature_names['betweenness'] = 'Betweenness'
    feature_names['max_betweenness'] = 'Max betweenness'
    feature_names['mean_betweenness'] = 'Mean betweenness'
    feature_names['std_betweenness'] = 'Std betweenness'
    
    feature_names['closeness'] = 'Closeness'
    feature_names['max_closeness'] = 'Max closeness'
    feature_names['mean_closeness'] = 'Mean closeness'
    feature_names['std_closeness'] = 'Std closeness'
    
    feature_names['pagerank'] = 'Pagerank'
    feature_names['max_pagerank'] = 'Max pagerank'
    feature_names['mean_pagerank'] = 'Mean pagerank'
    feature_names['std_pagerank'] = 'Std pagerank'
    
    feature_names['eigenvalues'] = 'Eigenvalues'
    feature_names['mean_eigenvector_centrality'] = 'Mean eigenvector centrality'
    feature_names['num_edges'] = 'Number of edges'

    feature_names['eigenvalues_intervals'] = 'Eigenvalues intervals'
    feature_names['eigenvalues_intervals_normalized'] = 'Eigenvalues intervals normalized'

    feature_names['IPR'] = 'IPR'
    feature_names['max_IPR'] = 'Max IPR'
    feature_names['mean_IPR'] = 'Mean IPR'

    feature_names['weights'] = 'Weights'
    feature_names['max_weights'] = 'Max weights'
    feature_names['mean_weights'] = 'Mean weights'
    feature_names['std_weights'] = 'Std weights'
    
    feature_names['community_edge_betweenness_optimal'] = 'Community edge betweenness: optimal count'
        
    feature_names['robustness'] = 'Robustness'
    return feature_names    


def parenclitic_transform(x, kdes = None, p = None, I = None, G = None, threshold_p = 0.5, min_x = 0, max_x = 1, graph_path = '', id_patient = -1, genes_names = []):
    print 'parenclitic_transform', G.sum()
    sys.stdout.flush()
    if G is None:
        start = timeit.default_timer()
        w = make_weights(x, kdes, p, I, min_x, max_x)
        stop = timeit.default_timer()
        print 'Time weights', stop - start
    
        start = timeit.default_timer()
        g = make_graph(w > threshold_p, w)
        stop = timeit.default_timer()
        print 'Make graph', stop - start
    else:
        w = G
        start = timeit.default_timer()
        g = make_graph(G, w)
        stop = timeit.default_timer()
        print 'Make graph', stop - start
    sys.stdout.flush()

    '''
    if graph_path != '':
        g.vs["label"] = genes_names
        if g.ecount() > 0 and "weight" in g.es:
            g.es["label"] = g.es["weight"]
        layout = g.layout("fr")
        igraph.plot(g, graph_path, bbox = (1024, 1024), layout = layout, vertex_size = 20)
    '''

    parenclitic = calculate_metrics(g, w)
    sys.stdout.flush()
    return parenclitic

def parenclitic_kdes(X, min_x = 0, max_x = 1):
    k = X.shape[1]
    kdes = np.empty((k, k), dtype=object)
    print kdes.shape, X.shape
    num_points = 10000
    p = np.zeros((k, k, num_points), dtype=np.float32)
    I = np.zeros((k, k, num_points + 1), dtype=np.float32)
    for i in range(k):
        start = timeit.default_timer()
        for j in range(k):
            if (i == j): continue
            data = np.array([X[:, i], X[:, j]]).astype('float32')
#            kde = KernelDensity(kernel='gaussian', algorithm='kd_tree', rtol=1.0e-4, atol=1.0e-4).fit(data)
            kde = stats.gaussian_kde(data)
            points = kde.resample(num_points)
            pr = np.array(kde(points))
            #points = kde.sample(num_points)
            #pr = np.array(kde.score_samples(points))
            p[i, j] = np.sort(pr)
            I[i, j] = np.concatenate([[0], np.cumsum(p[i, j])])
            kdes[i, j] = kde
        stop = timeit.default_timer()
        print 'KDE for ', i, 'calculated in ', stop - start
        sys.stdout.flush()

    return kdes, p, I


def make_genes_edge(X_prob_i, X_prob_j, X_i, X_j, thresholds_p):
    num_points = 10000
    data = np.array([X_prob_i, X_prob_j])
    det = np.linalg.det(np.corrcoef(data))
    eps = 1.0e-9
    if abs(det) < eps:
        return np.zeros((len(thresholds_p), X_i.shape[0]), dtype=np.bool)
    kde = stats.gaussian_kde(data)

    data = np.array([X_i, X_j])
    p = np.array(kde(data))

    points = kde.resample(num_points)
    pr = np.array(kde(points))
    pr.sort()
    #I = np.concatenate([np.cumsum(pr), [1]])
    G = np.zeros((len(thresholds_p), X_i.shape[0]), dtype=np.bool)
    for i, thr_p in enumerate(thresholds_p):
        #ind = np.flatnonzero(I >= thr_p)[0]
        ind = int(thr_p * num_points)
        if ind < pr.size:
            q = pr[ind]
            G[i] = p < q
        else:
            G[i] = np.ones((X_i.shape[0]), dtype=np.bool)
    return G

def make_genes_edges(X_prob, X, threshold_p):
    G = np.zeros((X.shape[0], 1, k), dtype=np.bool)
    for j in range(k):
        G[:, 0, j] = make_genes_edge(X_prob[:, i], X_prob[:, j], X[:, i], X[:, j], threshold_p)
    return G


def parenclitic_graphs(X_prob, X, threshold_p = 0.5, num_workers = 1, skip_values = lambda i, j: i >= j):
    if not isinstance(threshold_p, collections.Iterable):
        threshold_p = [threshold_p]

    k = X.shape[1]
    print 'parenclitic_graphs'
    sys.stdout.flush()

    global num_done, num_pairs
    num_pairs = 0
    num_done = 0
    for i in range(k):
        for j in range(k):
            if skip_values(i, j): continue
            num_pairs += 1
    num_bytes = len(np.packbits(np.zeros((X.shape[0], 1),dtype = np.bool)))
    G = np.zeros((len(threshold_p), num_pairs, num_bytes), dtype = np.uint8)
    print 'G size:', G.nbytes, G.shape
    
    need_parallel = num_workers > 1
    if need_parallel:
        pool = Pool(num_workers)
        global done_tasks
        done_tasks = 0
        global ready
        ready = Semaphore(num_workers * 10)

    start = timeit.default_timer()
    each_progress = int(np.sqrt(num_pairs + 0.5))
    lid = 0
    for i in range(k):
        for j in range(k):
            if skip_values(i, j): continue

            def upd_graph(g, i = i, j = j, lid = lid):
                global num_done, done_tasks, ready
                gp = np.packbits(g, axis=1)
                for id_thr, gl in enumerate(gp):
                    G[id_thr][lid] = gl
                done_tasks += 1
                ready.release()

                num_done += 1
                if num_done % each_progress == 0 or num_done == num_pairs:
                    stop = timeit.default_timer()
                    print 'Graph for', num_done, 'pairs calculated in', stop - start
                    sys.stdout.flush()

            if need_parallel:
                ready.acquire()
                pool.apply_async(make_genes_edge, args = (X_prob[:, i], X_prob[:, j], X[:, i], X[:, j], threshold_p), callback = upd_graph)
            else:
                g = make_genes_edge(X_prob[:, i], X_prob[:, j], X[:, i], X[:, j], threshold_p)
                upd_graph(g)

            lid += 1

    if need_parallel:
        while done_tasks < num_pairs:
            ready.acquire()

    pool.close()
    pool.join()
    sys.stdout.flush()
    return G


def extract_graph(G, num_features, id_sample):
    g = np.zeros((num_features, num_features), dtype=np.bool)
    lid = 0
    for i in range(num_features):
        for j in range(num_features):
            if i >= j: continue
            lg = G[lid]
            g[i, j] = get_bit(lg, id_sample)
            lid += 1
    return g

def extract_graphs(G, num_features, num_samples):
    print 'Extract graphs'
    sys.stdout.flush()    
    start = timeit.default_timer()
    num_bytes = len(np.packbits(np.zeros((num_features, 1),dtype = np.bool)))    
    g = np.zeros((num_features, num_samples, num_bytes), dtype=np.uint8)
    lid = 0
    for i in range(num_features - 1):
        startf = timeit.default_timer()
        gp = np.zeros((num_samples, num_features), dtype=np.bool)
        cnt = num_features - i - 1
        #print i, lid, lid + cnt, G.shape[0]
        #sys.stdout.flush()

        gp[:, (i + 1):] = np.unpackbits(G[lid:(lid + cnt), :], axis = 1)[:, :num_samples].T
        lid += cnt
        '''
        for j in range(num_features):
            if i >= j: continue
            lg = G[lid]
            gp[:, j] = np.unpackbits(lg)[:num_samples]
            lid += 1
        '''
        #stopf = timeit.default_timer()
        #print 'Extraction of', i, 'feature:', stopf - startf
        #sys.stdout.flush()

        startf = timeit.default_timer()
        g[i, :, :] = np.packbits(gp, axis=1)
        #stopf = timeit.default_timer()
        #print 'Packing of', i, 'feature:', stopf - startf
        #sys.stdout.flush()

    assert (lid == G.shape[0])
    g = np.swapaxes(g, 0, 1)
    print g.shape
    stop = timeit.default_timer()
    print 'Graphs extracted in', stop - start
    sys.stdout.flush()
    return g

def get_bit(arr, pos):
    num_bits = 8 * arr.dtype.itemsize
    return (arr[pos // num_bits] >> (num_bits - pos % num_bits - 1)) & 1

def get_degrees(G):
    g = make_graph(G, G)
    degrees = np.array(g.vs.degree())
    return degrees

def read_graphs(config, X, id_thr):
    num_bytes = len(np.packbits(np.zeros((X.shape[0], 1), dtype = np.bool)))
    G = np.zeros((X.shape[1] * (X.shape[1] - 1) / 2, num_bytes), dtype = np.uint8)
    print 'Read graph'
    sys.stdout.flush()
    start = timeit.default_timer()
    lid = 0
    for i in range(config.params["num_parts"].value):
        config.params["id_part"].set_tick(i)
        data = np.load(config.ofname(["graphs", "g"], ext = ".npz", include_set = config.params_sets["graphs"]))
        stop = timeit.default_timer()
        print 'Part', i, 'of graphs was read in', stop - start
        sys.stdout.flush()
        #if G is None:
        cur = data['G'][id_thr]
        G[lid:(lid + cur.shape[0]), :] = cur
        lid += cur.shape[0]
        #else:
        #    G = np.concatenate([G, data['G'][id_thr]])
        stop = timeit.default_timer()
        print 'Part', i, 'of graphs was added in', stop - start
        sys.stdout.flush()
        data.close()
    assert (lid == G.shape[0])
    return G
