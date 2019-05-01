from multiprocessing import Pool
import timeit
import numpy as np
from .transform_data import *
from sklearn.model_selection import train_test_split
from .gen_files.genes_mean_and_std import generate_genes_mean_and_std
import os.path
import sys
from .plot_kde import plot_parenclitic_kdes

num_genes = 10
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

data_name = 'GSE40279'
data_path = project_path + 'data/' + data_name + '/'

x_file_name = 'gene_mean.txt'
ranged_genes_file_name = 'pvals_mean_genes.txt'
#ranged_genes_file_name = 'all_gene_importance_RF_2class.csv'
y_file_name = 'ages.txt'
name_genes = 'anova'
xp_file_name = 'parenclitic_features_' + name_genes + '_' + str(num_genes) + '_thr_p_' + "{:.1f}".format(threshold_p) + '.txt'
print(xp_file_name)

#xp_file_name = 'parenclitic_features_RF.txt'

x_file_path = data_path + x_file_name
y_file_path = data_path + y_file_name
xp_file_path = data_path + xp_file_name
ranged_genes_file_path = data_path + ranged_genes_file_name

kde_path = data_path + 'kdes/'
graph_path = data_path + 'graphs/'

if not os.path.exists(kde_path):
    os.makedirs(kde_path)

if not os.path.exists(graph_path):
    os.makedirs(graph_path)

if not os.path.exists(x_file_path):
    generate_genes_mean_and_std(data_path)

X = np.genfromtxt(x_file_path, dtype='float32', delimiter=' ')[:, 1:]

genes_names = np.genfromtxt(x_file_path, dtype='str', usecols = 0)
genes_dict = dict((v, i) for i, v in enumerate(genes_names))
y = np.loadtxt(y_file_path, dtype='float32')

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

median_age = np.median(y)
print(median_age)
min_age = np.min(y)
max_age = np.max(y)

ages_edges = [min_age, median_age, max_age]

X_prob, _, y_prob, _ = train_test_split(
    X, y, test_size=0.9, random_state=42)

start = timeit.default_timer()
mask = y_prob < ages_edges[1]
print(mask.sum(), y_prob.shape, ages_edges[1])
sys.stdout.flush()
#kdes, p, I = parenclitic_kdes(X_prob[mask, :]) # Not old
G = parenclitic_graphs(X_prob, X, threshold_p)
stop = timeit.default_timer()
print('Parenclitic kdes elapsed: ', stop - start) 
sys.stdout.flush()

#plot_parenclitic_kdes(kde_path, X, y, mask, kdes, p, I, genes_names)


#parenclitic_transform(x, G = G[451, :, :])

start = timeit.default_timer()

Xp = np.zeros((X.shape[0], parenclitic_num_features()), dtype = np.float32)
pool = Pool(32)
results = np.empty((X.shape[0]), dtype=object)
for i, x in enumerate(X):
    g_path = graph_path + 'graph_' + str(i + 1) + '_' + str(y[i]) + '.png'
    g = extract_graph(G, num_genes, i)
    results[i] = pool.apply_async(parenclitic_transform, args = (x,), kwds = {'G': g, 'graph_path': g_path, 'genes_names': genes_names})
    #Xp[i, :] = parenclitic_transform(x, G = G[i, :, :], graph_path = g_path, genes_names = genes_names)
#    Xp[i, :] = parenclitic_transform(x, kdes, p, I, G, graph_path = g_path, genes_names = genes_names)
    sys.stdout.flush()

for i in range(results.shape[0]):
    Xp[i, :] = results[i].get()
    print('Results', i)

print((Xp[0, :]))
stop = timeit.default_timer()
print('Parenclitic transform elapsed: ', stop - start) 

np.savetxt(xp_file_path, Xp)
