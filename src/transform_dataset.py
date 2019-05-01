from multiprocessing import Pool
import timeit
import numpy as np
from .transform_data import *
from sklearn.model_selection import train_test_split
from .gen_files.genes_mean_and_std import generate_genes_mean_and_std
import os.path
import sys
from .plot_kde import plot_parenclitic_kdes

num_genes = 5
threshold_p = 0.5
threshold_p_be = 0.5
threshold_p_en = 0.9
threshold_p_n = 5

import sys
if len(sys.argv) > 1:
    param_id = float(sys.argv[1])
    threshold_p = (threshold_p_en - threshold_p_be) / (threshold_p_n - 1) * param_id + threshold_p_be
    

start = timeit.default_timer()
project_path = '/home/krivonosov_m/Projects/Gerontology/'

data_name = 'GSE52588'
data_path = project_path + 'data/' + data_name + '/'

mongoloids_mask = np.arange(0, 28)
siblings_mask = np.arange(29, 57)
mothers_mask = np.arange(58, 86)

x_file_name = 'gene_mean_islands_shores.txt'
ranged_genes_file_name = 'linreg_genes_mean_islands_shores.txt'
#ranged_genes_file_name = 'all_gene_importance_RF_2class.csv'
name_genes = 'linreg_genes_mean_islands_shores'
#xp_file_name = 'parenclitic_features_' + name_genes + '_' + str(num_genes) + '_thr_p_' + "{:.1f}".format(threshold_p) + '.txt'
degrees_file_name = 'degrees_' + name_genes + '_num_genes_' + str(num_genes) + '_thr_p_' + "{:.1f}".format(threshold_p) + '.txt'
print(degrees_file_name)

#xp_file_name = 'parenclitic_features_RF.txt'

x_file_path = data_path + x_file_name
#xp_file_path = data_path + xp_file_name
degrees_file_path = data_path + degrees_file_name
ranged_genes_file_path = data_path + ranged_genes_file_name

kde_path = data_path + 'kdes/'
graph_path = data_path + 'graphs/'

if not os.path.exists(kde_path):
    os.makedirs(kde_path)

if not os.path.exists(graph_path):
    os.makedirs(graph_path)

#if not os.path.exists(x_file_path):
#    generate_genes_mean_and_std(data_path)

X = np.genfromtxt(x_file_path, dtype='float32', delimiter=' ')[:, 1:]

genes_names = np.genfromtxt(x_file_path, dtype='str', usecols = 0)
genes_dict = dict((v, i) for i, v in enumerate(genes_names))

ranged_genes = np.genfromtxt(ranged_genes_file_path, dtype='str', usecols = 0)

stop = timeit.default_timer()
print('Data loaded: ', stop - start) 
print(X.dtype, X.shape)

sys.stdout.flush()



#X = np.random.rand(len(genes_names), 656)
genes_names = ranged_genes[:num_genes]

indices = np.array([genes_dict[x] for x in genes_names])

X = X[indices, :].T

print(X.shape, num_genes)


X_prob = X[mothers_mask, :]

start = timeit.default_timer()
sys.stdout.flush()

G = parenclitic_graphs(X_prob, X, threshold_p, num_workers = 32)
stop = timeit.default_timer()
print('Parenclitic kdes elapsed: ', stop - start) 
sys.stdout.flush()

#plot_parenclitic_kdes(kde_path, X, y, mask, kdes, p, I, genes_names)


#parenclitic_transform(x, G = G[451, :, :])

start = timeit.default_timer()

#Xp = np.zeros((X.shape[0], parenclitic_num_features()), dtype = np.float32)
degrees = np.zeros((X.shape[0], X.shape[1]), dtype = np.float32)
pool = Pool(32)
results = np.empty((X.shape[0]), dtype=object)
for i, x in enumerate(X):
    g_path = graph_path + 'graph_' + str(i + 1) + '.png'
    g = extract_graph(G, num_genes, i)
    print(i)
    print(g)
    results[i] = pool.apply_async(get_degrees, kwds = {'G': g})
#    results[i] = pool.apply_async(parenclitic_transform, args = (x,), kwds = {'G': g, 'graph_path': g_path, 'genes_names': genes_names})
    #Xp[i, :] = parenclitic_transform(x, G = G[i, :, :], graph_path = g_path, genes_names = genes_names)
#    Xp[i, :] = parenclitic_transform(x, kdes, p, I, G, graph_path = g_path, genes_names = genes_names)
    sys.stdout.flush()

for i in range(results.shape[0]):
    degrees[i, :] = results[i].get()
    print('Results', i)

degrees_mongoloids = degrees[mongoloids_mask, :].mean(axis = 0)
degrees_siblings = degrees[siblings_mask, :].mean(axis = 0)

print(degrees_mongoloids.shape, degrees_siblings.shape, genes_names.shape)

T = np.zeros(genes_names.size, dtype = [('genes_names', 'S32'), ('degrees_mongoloids', np.float32), ('degrees_siblings', np.float32)])
T['genes_names'] = genes_names
T['degrees_mongoloids'] = degrees_mongoloids
T['degrees_siblings'] = degrees_siblings

np.savetxt(degrees_file_path, T, fmt = ('%s', '%.18e', '%.18e'), header = "\t".join(list(T.dtype.names)), delimiter = '\t', comments = '')

stop = timeit.default_timer()
#print 'Parenclitic transform elapsed: ', stop - start 
print('Degrees calc elapsed: ', stop - start) 

#np.savetxt(xp_file_path, Xp)

