import os
print(os.environ)

import time
import numpy as np
from scipy.stats import gaussian_kde
from threading import Thread
from queue import Queue

def make_G(X_prob_i, X_prob_j, X_i, X_j, threshold_p = 0.9):
    num_points = 10000
    data = np.array([X_prob_i, X_prob_j])
    det = np.linalg.det(np.corrcoef(data))
    eps = 1.0e-9
    if abs(det) < eps:
        return np.zeros((X_i.shape[0]), dtype=np.bool)

    kde = gaussian_kde(data)
    points = kde.resample(num_points)
    pr = np.array(kde(points))
    pr.sort()
    I = np.cumsum(pr)
    I = I / I[-1]
    ind = np.flatnonzero(I >= threshold_p)[0]
    data = np.array([X_i, X_j])
    p = np.array(kde(data))
    if ind < pr.size:
        q = pr[ind]
        G = p < q
    else:
        G = np.ones((X_i.shape[0]), dtype=np.bool)
    
    return G

class CalcWorker(Thread):
    def __init__(self, queue):
        Thread.__init__(self)
        self.queue = queue
    
    def run(self):
        global x, G
        num = 65
        while True:
            i = self.queue.get()
            be = time.time()
            for j in range(G.shape[2]):
                G[:, i, j] = make_G(x[:num, i], x[:num, j], x[:, i], x[:, j])
            self.queue.task_done()
            en = time.time()
            print(i, 'Time: ', en - be)

def main():
    global x, G
    ts = time.time()
    n = 1000
    num = 65
    np.random.seed(42)
    x = np.random.rand(100, n)
    G = np.zeros((100, n, n), dtype = 'bool')
    
    #for i in range(n):
    '''
    be = time.time()
    i = 0
    for j in range(G.shape[2]):
        G[:, i, j] = make_G(x[:num, i], x[:num, j], x[:, i], x[:, j])
    en = time.time()
    print 'Time: ', en - be
    print G.mean()
    '''
    #return
    
    # Create a queue to communicate with the worker threads
    queue = Queue()
    # Create 8 worker threads
    for i in range(4):
        worker = CalcWorker(queue)
        # Setting daemon to True will let the main thread exit even though the workers are blocking
        worker.daemon = True
        worker.start()
    # Put the tasks into the queue as a tuple
    for i in range(4):
        queue.put(i)
    # Causes the main thread to wait for the queue to finish processing all the tasks
    queue.join()

    print('Took {}'.format(time.time() - ts))
    print(G.mean())

main()