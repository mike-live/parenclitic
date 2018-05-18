import timeit
import numpy as np
from transform_data import *
from sklearn.model_selection import train_test_split
from gen_files.genes_mean_and_std import generate_genes_mean_and_std
import os.path

project_path = '/home/krivonosov_m/Projects/Gerontology/'

data_path = project_path + 'data/GSE40279/'

x_file_name = 'gene_mean.txt'
y_file_name = 'ages.txt'
xp_file_name = 'parenclitic_features.txt'

x_file_path = data_path + x_file_name
y_file_path = data_path + y_file_name
xp_file_path = data_path + xp_file_name

if not os.path.exists(x_file_path):
    generate_genes_mean_and_std(data_path)

X = np.loadtxt(x_file_path)
y = np.loadtxt(y_file_path)

print X.shape
print y.shape


'''
ages_edges = [0, 58, 100]

X_prob, _, y_prob, _ = train_test_split(
    X, y, test_size=0.1, random_state=42)

start = timeit.default_timer()
kdes, p, I = parenclitic_kdes(X_prob[y_prob < ages_edges[1]]) # Not old
stop = timeit.default_timer()
print 'Parenclitic kdes elapsed: ', stop - start 

start = timeit.default_timer()

Xp = np.zeros((X.shape[0], parenclitic_num_features()), dtype = np.float32)
for i, x in enumerate(X):
    Xp[i, :] = parenclitic_transform(x, kdes, p, I)

print(Xp[0, :])
stop = timeit.default_timer()
print 'Parenclitic transform elapsed: ', stop - start 

np.savetxt(xp_file_path, Xp)
'''