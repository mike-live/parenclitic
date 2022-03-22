import numpy as np
from tqdm.notebook import tqdm

def compute_percent_graph(gs, g, groups, group_masks, group_names):
    edges = [edge.tuple for edge in g.es]
    edge_id = dict(zip(edges, list(range(len(edges)))))
    for id_group, group in enumerate(groups):
        cnt = np.zeros(len(edges), dtype = np.int32)
        for id_sample in group_masks[id_group]:
            cur_edges = [edge.tuple for edge in gs[id_sample].es]
            for edge in cur_edges:
                if edge in edge_id:
                    cnt[edge_id[edge]] += 1
             
        g.es[group_names[id_group].lower() + "_percent"] = cnt / len(group_masks[id_group])
    return g

def make_diff_graph(config, gs, groups, group_masks, group_names, diff_graph_path):
    num_groups = len(groups)
    config.params["id_sample"].manual_ticks = True    
    g_or = [None] * num_groups
    g_and = [None] * num_groups
    print (diff_graph_path)
    
    for id_sample in tqdm(config.params["id_sample"], desc = "Graph merging"):
        #print (id_sample, '/', config.params["id_sample"].num_ticks, end = ' ')
        g = gs[id_sample].copy()
        
        #print (g.vcount(), g.ecount(), end = ' ')
        for id_group, group in enumerate(groups):
            if id_sample in config.params[group].value:
                #if id_group == 0:
                #    g.delete_edges(g.es.select(weight_lt = 0.0))
                if g_or[id_group] is None:
                    g_or[id_group] = g.copy()
                    g_and[id_group] = g.copy()
                else:
                    g_or[id_group] = g_or[id_group].union(g, byname=False)
                    g_and[id_group] = g_and[id_group].intersection(g, byname=False)
                #print(g_and[id_group].ecount())
        #print()
    #gn = g_or[get_4network(config)].copy()
    #print('set len', len(set([v.index for v in gn.vs if v.degree() < 1]) & set([v.index for v in gs[1].vs if v.degree() > 0])))
    #m_diff = np.logical_and(g_and[1], np.logical_not(g_or[0]))

    #edges_0 = [edge.tuple for edge in g_and[0].es]
    #print('edges_and')
    #print(edges_0)

    #edges_1 = [edge.tuple for edge in g_or[1].es]
    #print('edges_or')
    #print(edges_1)

    m_diff = g_and[0].copy()
    cur = g_and[0].intersection(g_or[1])
    #print(cur.vcount(), cur.ecount())
    edges = [edge.tuple for edge in cur.es]

    #print(m_diff.vcount(), m_diff.ecount())
    #print('edges_diff')
    #print(edges)

    
    m_diff.delete_edges(edges)
    for i in range(len(groups)):
        g_and[i] = compute_percent_graph(gs, g_and[i], groups, group_masks, group_names)
        g_or[i] = compute_percent_graph(gs, g_or[i], groups, group_masks, group_names)
    m_diff = compute_percent_graph(gs, m_diff, groups, group_masks, group_names)
    
    np.savez_compressed(diff_graph_path, 
                        m_diff = [m_diff], 
                        g_or = g_or, g_and = g_and)
    return

def load_diff_graphs(diff_graph_path):
    m_data = np.load(diff_graph_path, allow_pickle=True)
    m_diff = m_data['m_diff']
    g_or = m_data['g_or']
    g_and = m_data['g_and']
    m_data.close()
    print (g_or[0].ecount(), g_or[1].ecount(), g_and[0].ecount(), g_and[1].ecount())
    print (m_diff[0].vcount(), m_diff[0].ecount())
    #print(g_and[1].es["ds_percent"])
    return m_diff[0], g_or, g_and

