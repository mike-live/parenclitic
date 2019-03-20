#import matplotlib as mpl
#mpl.use('Agg')
#import matplotlib.pyplot as plt
#import matplotlib.pylab as pylab

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
import graphs_aux
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



def _robustness(g):
    cnt = 0
    while g.vcount() > 0:
        degrees = np.array(g.degree())
        max_degree = np.max(degrees)
        if max_degree == 0:
            break
        g.delete_vertices(np.flatnonzero(max_degree == degrees)[0])
        cnt = cnt + 1
    return cnt

def robustness(g, weights = None):
    cnt = 0
    while g.ecount() > 0:
        degrees = np.array(g.strength(weights = weights))
        g.delete_vertices(np.argmax(degrees))
        cnt = cnt + 1
    return cnt


def calculate_metrics(g, w, need_weights = True, get_big = True):
    print 'Metrics'
    sys.stdout.flush()

    if need_weights:
        weights = 'weight'
    else:
        weights = None

    parenclitic = pd.DataFrame(index=[0])
    start = timeit.default_timer()

    degrees = np.array(g.strength(weights = weights))
    if get_big: 
        parenclitic['degrees'] = [degrees]
    parenclitic['min_degrees'] = np.min (degrees)
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
    print 'here'
    sys.stdout.flush()

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
    if get_big: 
        parenclitic['betweenness'] = [np.array(betweenness)]
    parenclitic['min_betweenness'] = np.min (betweenness)
    parenclitic['max_betweenness'] = np.max (betweenness)
    parenclitic['mean_betweenness'] = np.mean(betweenness)
    parenclitic['std_betweenness'] = np.std (betweenness)
    betweenness = None

    stop = timeit.default_timer()
    print 'Parenclitic 3', stop - start
    sys.stdout.flush()
        
    
    start = timeit.default_timer()

    closeness = g.closeness(weights = weights)
    if get_big: 
        parenclitic['closeness'] = [np.array(closeness)]
    parenclitic['min_closeness'] = np.min (closeness)
    parenclitic['max_closeness'] = np.max (closeness)
    parenclitic['mean_closeness'] = np.mean(closeness)
    parenclitic['std_closeness'] = np.std (closeness)
    closeness = None

    stop = timeit.default_timer()
    print 'Parenclitic 4', stop - start
    sys.stdout.flush()

    start = timeit.default_timer()

    pagerank = g.pagerank(weights = weights)
    if get_big: 
        parenclitic['pagerank'] = [np.array(pagerank)]
    parenclitic['min_pagerank'] = np.min (pagerank)
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
    if get_big: 
        parenclitic['eigenvector_centrality'] = [np.array(eigenvector_centrality)]
    parenclitic['min_eigenvector_centrality'] = np.min(eigenvector_centrality)
    parenclitic['max_eigenvector_centrality'] = np.max(eigenvector_centrality)
    parenclitic['mean_eigenvector_centrality'] = np.mean(eigenvector_centrality)
    parenclitic['std_eigenvector_centrality'] = np.std(eigenvector_centrality)
    eigenvector_centrality = None

    stop = timeit.default_timer()
    print 'Parenclitic centrality', stop - start
    sys.stdout.flush()
        
    start = timeit.default_timer()
    
    largest = g.clusters().giant()    
    m = np.array(largest.get_adjacency().data)
    sys.stdout.flush()

    eigenvalues, eigenvectors = LA.eig(m)
    #Suppose symmetric matrix
    eigenvalues = np.real(eigenvalues)
    eigenvectors = np.real(eigenvectors)
    print 'Eigenvectors', stop - start
    sys.stdout.flush()

    eigenvalues_intervals = np.diff(np.sort(eigenvalues)) 
    print 'intervals', stop - start
    sys.stdout.flush()

    eigenvalues_intervals_normalized = eigenvalues_intervals / np.mean(eigenvalues_intervals)

    print 'normalized', stop - start
    sys.stdout.flush()
    
    if get_big: 
        parenclitic['eigenvalues'] = [np.array(eigenvalues)]
        parenclitic['eigenvalues_intervals'] = [np.array(eigenvalues_intervals)]
        parenclitic['eigenvalues_intervals_normalized'] = [np.array(eigenvalues_intervals_normalized)]

    stop = timeit.default_timer()
    print 'Parenclitic: eigenvalues', stop - start
    sys.stdout.flush()
    
    IPR = np.sum(np.power(eigenvectors, 4), axis=0) / np.power(np.sum(np.power(eigenvectors, 2), axis=0), 2)
    if get_big: 
        parenclitic['IPR'] = [np.array(IPR)]
    parenclitic['max_IPR'] = np.max(IPR)
    parenclitic['mean_IPR'] = np.mean(IPR)

    stop = timeit.default_timer()
    print 'Parenclitic 7', stop - start
    sys.stdout.flush()

    eigenvectors = None
    eigenvalues = None
    IPR = None
    eigenvalues_intervals = None
    eigenvalues_intervals_normalized = None
    
    start = timeit.default_timer()
    parenclitic['num_edges'] = g.ecount()

    if g.ecount() > 0:
        weights = np.array(g.es["weight"])
        if get_big: 
            parenclitic['weights'] = [np.array(weights)]
        parenclitic['sum_weights'] = np.sum (weights)
        parenclitic['min_weights'] = np.min (weights)
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
    parenclitic['robustness'] = robustness(g, weights)
    stop = timeit.default_timer()
    print 'Parenclitic 10', stop - start
    sys.stdout.flush()

    return parenclitic

def parenclitic_feature_logscale():
    feature_logscale = {
        'degrees': False,
        'max_degrees': True,
        'mean_degrees': False,
        'std_degrees': False,

        'efficiency': True,
        'betweenness': False,
        'max_betweenness': True,
        'mean_betweenness': True,
        'std_betweenness': True,

        'closeness': False,
        'max_closeness': False,
        'mean_closeness': False,
        'std_closeness': True,

        'pagerank': False,
        'max_pagerank': True,
        'mean_pagerank': True,
        'std_pagerank': True,

        'eigenvalues': False,
        'mean_eigenvector_centrality': False,
        'num_edges': False,

        'eigenvalues_intervals': False,
        'eigenvalues_intervals_normalized': False,

        'IPR': False,
        'max_IPR': False,
        'mean_IPR': False,

        'weights': False,
        'max_weights': False,
        'mean_weights': False,
        'std_weights': False,

        'community_edge_betweenness_optimal': False,

        'robustness': False,
    }
    return feature_logscale
    
def parenclitic_feature_names():
    feature_names = {}
    feature_names['degrees'] = 'Degrees'
    feature_names['min_degrees'] = 'Min degrees'
    feature_names['max_degrees'] = 'Max degrees'
    feature_names['mean_degrees'] = 'Mean degrees'
    feature_names['std_degrees'] = 'Std degrees'
    
    feature_names['efficiency'] = 'Efficiency'
    feature_names['betweenness'] = 'Betweenness'
    feature_names['min_betweenness'] = 'Min betweenness'
    feature_names['max_betweenness'] = 'Max betweenness'
    feature_names['mean_betweenness'] = 'Mean betweenness'
    feature_names['std_betweenness'] = 'Std betweenness'
    
    feature_names['closeness'] = 'Closeness'
    feature_names['min_closeness'] = 'Min closeness'
    feature_names['max_closeness'] = 'Max closeness'
    feature_names['mean_closeness'] = 'Mean closeness'
    feature_names['std_closeness'] = 'Std closeness'
    
    feature_names['pagerank'] = 'Pagerank'
    feature_names['min_pagerank'] = 'Min pagerank'
    feature_names['max_pagerank'] = 'Max pagerank'
    feature_names['mean_pagerank'] = 'Mean pagerank'
    feature_names['std_pagerank'] = 'Std pagerank'
    
    feature_names['eigenvalues'] = 'Eigenvalues'
    feature_names['min_eigenvector_centrality'] = 'Min eigenvector centrality'
    feature_names['max_eigenvector_centrality'] = 'Max eigenvector centrality'
    feature_names['mean_eigenvector_centrality'] = 'Mean eigenvector centrality'
    feature_names['std_eigenvector_centrality'] = 'Std eigenvector centrality'
    feature_names['num_edges'] = 'Number of edges'

    feature_names['eigenvalues_intervals'] = 'Eigenvalues intervals'
    feature_names['eigenvalues_intervals_normalized'] = 'Eigenvalues intervals normalized'

    feature_names['IPR'] = 'IPR'
    feature_names['max_IPR'] = 'Max IPR'
    feature_names['mean_IPR'] = 'Mean IPR'

    feature_names['weights'] = 'Weights'
    feature_names['min_weights'] = 'Min weights'
    feature_names['max_weights'] = 'Max weights'
    feature_names['mean_weights'] = 'Mean weights'
    feature_names['std_weights'] = 'Std weights'
    
    feature_names['community_edge_betweenness_optimal'] = 'Community edge betweenness: optimal count'
        
    feature_names['robustness'] = 'Robustness'
    return feature_names    


def parenclitic_transform(x = None, kdes = None, p = None, I = None, G = None, IDS = None, w = None, g = None, threshold_p = 0.5, min_x = 0, max_x = 1, graph_path = '', id_patient = -1, genes_names = []):
    print 'parenclitic_transform'
    sys.stdout.flush()
    if g is None:
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
            start = timeit.default_timer()
            g = graphs_aux.make_graph(G = G, weights = w, edges = IDS)
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
    g = None
    w = None
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

def make_genes_edge_svc(X_i, X_j, y, mask, min_score, by_group):
    from sklearn import svm, datasets
    from scipy import stats
    data = stats.zscore(np.array([X_i, X_j]).T)

    if by_group:
        classes = np.unique(y)
        G = np.zeros((len(y), ), np.bool)
        D = np.zeros((len(y), ), np.float32)
        for c in classes:
            fit_mask = (y == c) | (mask == 1)
            if len(np.unique(y[fit_mask])) == 1:
                continue
            clf = svm.LinearSVC(C = 1, class_weight = "balanced")
            clf.fit(data[fit_mask], y[fit_mask] == c)
            score = clf.score(data[fit_mask], y[fit_mask] == c)
            
            G[y == c] = clf.predict(data[y == c]) == 1
            D[y == c] = clf.decision_function(data[y == c])
            if score < min_score:
                G[fit_mask] = False
                
        G = G.reshape((1, len(G)))
    else:
        clf = svm.LinearSVC(C = 1, class_weight = "balanced")
    
        #clf = svm.SVC(kernel = 'linear', C = 1, class_weight = "balanced")
        fit_mask = (mask == 0) | (mask == 1)
        
        clf.fit(data[fit_mask], y[fit_mask] == 0)
        G = clf.predict(data) == 1
        score = clf.score(data[fit_mask], y[fit_mask] == 0)
            
        if score < min_score:
            G[:] = False
        G = G.reshape((1, len(G)))
        D = clf.decision_function(data)
    return [G, D]
    
def make_genes_edges(X_prob, X, y, threshold_p):
    G = np.zeros((X.shape[0], 1, k), dtype=np.bool)
    for j in range(k):
        #G[:, 0, j] = make_genes_edge(X_prob[:, i], X_prob[:, j], X[:, i], X[:, j], threshold_p)
        G[:, 0, j] = make_genes_edge_svc(X[:, i], X[:, j], y)
    return G

def parenclitic_graphs(mask, X, y, get_ids, min_score = 0.75, by_group = False, threshold_p = 0.5, num_workers = 1, algo = "svc"): # skip_values = lambda i, j: i >= j
    if not threshold_p is None and not isinstance(threshold_p, collections.Iterable):
        threshold_p = [threshold_p]

    k = X.shape[1]
    print 'parenclitic_graphs'
    sys.stdout.flush()

    global num_done, num_pairs
    num_pairs = 0
    num_done = 0
    for i, j in get_ids():
        num_pairs += 1
        
    '''
    for i in range(k):
        for j in range(k):
            if skip_values(i, j): continue
            num_pairs += 1
    '''
    num_bytes = len(np.packbits(np.zeros((X.shape[0], 1),dtype = np.bool)))
    D = []
    IDS = []
    if threshold_p is None:
        G = np.zeros((num_pairs, num_bytes), dtype = np.uint8)
    else:
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
    '''
    for i in range(k):
        for j in range(k):
            if skip_values(i, j): continue
    '''
    for i, j in get_ids():
        def upd_graph(ls, i = i, j = j, lid = lid):
            global num_done, done_tasks, ready
            g = ls[0]
            d = ls[1]
            gp = np.packbits(g, axis=1)
            if threshold_p is None:
                G[lid] = gp
                if g.any():
                    D.append(d)
                    IDS.append([i, j])
            else: 
                for id_thr, gl in enumerate(gp):
                    G[id_thr][lid] = gl
        
            if need_parallel:
                done_tasks += 1
                ready.release()

            num_done += 1
            if num_done % each_progress == 0 or num_done == num_pairs:
                stop = timeit.default_timer()
                print 'Graph for', num_done, 'pairs calculated in', stop - start
                sys.stdout.flush()

        if need_parallel:
            ready.acquire()
            if algo == "svc":
                sys.stdout.flush()
                pool.apply_async(make_genes_edge_svc, args = (X[:, i], X[:, j], y, mask, min_score, by_group), callback = upd_graph)
            else: 
                pool.apply_async(make_genes_edge, args = (X[mask == 1, i], X[mask == 1, j], X[:, i], X[:, j], threshold_p), callback = upd_graph)
        else:
            if algo == "svc":
                g = make_genes_edge_svc(X[:, i], X[:, j], y, mask, min_score, by_group)
            else:
                g = make_genes_edge(X[mask == 1, i], X[mask == 1, j], X[:, i], X[:, j], threshold_p)
            upd_graph(g)

        lid += 1

    if need_parallel:
        while done_tasks < num_pairs:
            ready.acquire()

    pool.close()
    pool.join()
    sys.stdout.flush()
    
    D = np.array(D)
    IDS = np.array(IDS)
    return G, D, IDS


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

def read_graphs_part(config, id_part, id_thr = None):
    config.params["id_part"].set_tick(id_part)
    data = np.load(config.ofname(["graphs", "g"], ext = ".npz", include_set = config.params_sets["graphs"]))
    if id_thr is None:
        cur = data['G']
        dcur = data['D']
        idscur = data['IDS']
    else:
        cur = data['G'][id_thr]
        dcur = None
        idscur = None
    data.close()
    return cur, dcur, idscur
    
def read_graphs(config, X, need_G = True, id_thr = None):
    num_bytes = len(np.packbits(np.zeros((X.shape[0], 1), dtype = np.bool)))
    if need_G:
        G = np.zeros((X.shape[1] * (X.shape[1] - 1) / 2, num_bytes), dtype = np.uint8)
        lid = 0
    else:
        G = None
    print 'Read graph'
    sys.stdout.flush()
    start = timeit.default_timer()
    D = []
    IDS = []
    for id_part in range(config.params["num_parts"].value):
        cur, dcur, idscur = read_graphs_part(config, id_part, id_thr)
        
        stop = timeit.default_timer()
        print 'Part', id_part, 'of graphs was read in', stop - start
        print 'Part', id_part, 'of graphs consist of', len(dcur), 'edges'
        sys.stdout.flush()
        
        if need_G:
            G[lid:(lid + cur.shape[0]), :] = cur
            lid += cur.shape[0]
            
        D.extend(dcur.tolist())
        IDS.extend(idscur.tolist())

        stop = timeit.default_timer()
        print 'Part', id_part, 'of graphs was added in', stop - start
        sys.stdout.flush()
    
    if need_G:
        print lid, G.shape[0]
        assert (lid == G.shape[0])
    D = np.array(D).T
    IDS = np.array(IDS)
    print D.shape, IDS.shape, IDS.dtype
    return G, D, IDS
