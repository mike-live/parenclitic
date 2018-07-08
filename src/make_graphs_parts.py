from transform_data import *

def get_part(id_part, num_parts, num_features):
    cnt = num_features * (num_features - 1) / 2
    lens = np.tile(cnt // num_parts, (num_parts, 1))
    lens[:(cnt % num_parts)] += 1
    clens = np.concatenate([[0], np.cumsum(lens)])
    cur = 0
    cur_id = 0
    be = clens[id_part]
    en = clens[id_part + 1] - 1
    #print id_part, id_part + 1, clens[id_part + 1], clens[id_part + 1] - 1
    #print be, en, num_features
    #print lens, clens, cnt % num_parts, cnt, num_parts
    sys.stdout.flush()
    for i in range(num_features):
        for j in range(num_features):
            if i >= j: continue
            cur_id = i * num_features + j
            # print 'Cur: ', cur_id, cur
            if cur == be:
                be_part = cur_id
                print 'Begin part: ', i, j
            if cur == en:
                en_part = cur_id
                print 'End part: ', i, j
            cur += 1
    sys.stdout.flush()
    return be_part, en_part
    
def make_graphs_part(X_prob, X, threshold_p, id_part, num_parts, num_workers = 1):
    start = timeit.default_timer()
    be_part, en_part = get_part(id_part, num_parts, X.shape[1])
    G = parenclitic_graphs(X_prob, X, threshold_p, num_workers = num_workers, skip_values = lambda i, j: (i >= j) or (be_part > i * X.shape[1] + j or i * X.shape[1] + j > en_part))
    stop = timeit.default_timer()
    print 'Make graphs: ', stop - start 
    sys.stdout.flush()
    return G


