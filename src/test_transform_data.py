import timeit
import numpy as np
from .transform_data import *
from sklearn.model_selection import train_test_split

np.random.seed(0)

n = 1000
k = 100
X = np.random.rand(n, k)
y = np.random.randint(2, size = n)

print(X)
print(y)

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.33, random_state=42)

start = timeit.default_timer()
kdes, p, I = parenclitic_kdes(X_train[y_train == 1]) # Healthy
stop = timeit.default_timer()
print('Parenclitic kdes elapsed: ', stop - start) 

par = parenclitic_transform(X_train[0, :], kdes, p, I)
print(par)

'''
start = timeit.default_timer()

Xp_train = np.zeros((X_train.shape[0], parenclitic_num_features()), dtype = np.float32)
for i, x in enumerate(X_train):
    Xp_train[i, :] = parenclitic_transform(x, kdes, p, I)

Xp_test = np.zeros((X_test.shape[0], parenclitic_num_features()), dtype = np.float32)
for i, x in enumerate(X_test):
    Xp_test[i, :] = parenclitic_transform(x, kdes, p, I)

print(Xp_train[0, :])
stop = timeit.default_timer()
print 'Parenclitic transform elapsed: ', stop - start 
'''