import timeit
import numpy as np
from transform_data import *
from sklearn.model_selection import train_test_split
from gen_files.genes_mean_and_std import generate_genes_mean_and_std
import os.path
import sys
from plot_kde import plot_parenclitic_kdes

start = timeit.default_timer()
project_path = '/home/krivonosov_m/Projects/Gerontology/'

data_path = project_path + 'data/GSE40279/'

x_file_name = 'gene_mean.txt'
ranged_genes_file_name = 'pvals_mean_genes.txt'
#ranged_genes_file_name = 'all_gene_importance_RF_2class.csv'
y_file_name = 'ages.txt'
xp_file_name = 'parenclitic_features_anova.txt'
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
print 'Data loaded: ', stop - start 
print X.dtype, X.shape

sys.stdout.flush()

num_genes = 400

#X = np.random.rand(len(genes_names), 656)
genes_names = ranged_genes[:num_genes]

indices = np.array([genes_dict[x] for x in genes_names])

X = X[indices, :].T

print X.shape

median_age = np.median(y)
print median_age
min_age = np.min(y)
max_age = np.max(y)

ages_edges = [min_age, median_age, max_age]

X_prob, _, y_prob, _ = train_test_split(
    X, y, test_size=0.1, random_state=42)

start = timeit.default_timer()
mask = y_prob < ages_edges[1]
kdes, p, I = parenclitic_kdes(X_prob[mask]) # Not old
stop = timeit.default_timer()
print 'Parenclitic kdes elapsed: ', stop - start 
sys.stdout.flush()

#plot_parenclitic_kdes(kde_path, X, y, mask, kdes, p, I, genes_names)


start = timeit.default_timer()

Xp = np.zeros((X.shape[0], parenclitic_num_features()), dtype = np.float32)
for i, x in enumerate(X):
    g_path = graph_path + 'graph_' + str(i + 1) + '_' + str(y[i]) + '.png'
    Xp[i, :] = parenclitic_transform(x, kdes, p, I, graph_path = g_path, genes_names = genes_names)
    sys.stdout.flush()

print(Xp[0, :])
stop = timeit.default_timer()
print 'Parenclitic transform elapsed: ', stop - start 

np.savetxt(xp_file_path, Xp)

