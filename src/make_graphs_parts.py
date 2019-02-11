from transform_data import *

def get_part(id_part, num_parts, num_features):
    cnt = num_features * (num_features - 1) / 2
    lens = np.tile(cnt // num_parts, (num_parts, 1))
    lens[:(cnt % num_parts)] += 1
    clens = np.concatenate([[0], np.cumsum(lens)])
    be = clens[id_part]
    en = clens[id_part + 1] - 1
    #print id_part, id_part + 1, clens[id_part + 1], clens[id_part + 1] - 1
    #print be, en, num_features
    #print lens, clens, cnt % num_parts, cnt, num_parts
    sys.stdout.flush()

    def get_ids(num_features = num_features, be = be, en = en):
        cur = 0
        for i in range(num_features):
            l = num_features - i - 1
            if (cur <= be and be < cur + l) or (cur <= en and en < cur + l) or (be <= cur and cur <= en):
                for j in range(num_features):
                    if i >= j: continue
                    cur_id = i * num_features + j
                    # print 'Cur: ', cur_id, cur
                    if be <= cur and cur <= en:
                        yield i, j
                    cur += 1
            else:        
                cur += l
    
    '''
    for i in range(num_features):
        l = num_features - i - 1
        if (cur <= be and be < cur + l) or (cur <= en and en < cur + l):
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
        else:        
            cur += l
    '''
    sys.stdout.flush()
    return get_ids
    
def make_graphs_part(mask, X, y, threshold_p, id_part, num_parts, num_workers = 1, algo = "svc"):
    start = timeit.default_timer()
    get_ids = get_part(id_part, num_parts, X.shape[1])
    G = parenclitic_graphs(mask, X, y, get_ids, threshold_p, num_workers = num_workers, algo = algo)
    #skip_values = lambda i, j: (i >= j) or (be_part > i * X.shape[1] + j or i * X.shape[1] + j > en_part)
    stop = timeit.default_timer()
    print 'Make graphs: ', stop - start 
    sys.stdout.flush()
    return G


