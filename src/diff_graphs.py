from configuration import param
import numpy as np

def make_diff_graph_mongoloids(config):
    diff_graph_path = config.ofname(["diff_graph"], ext = ".npz", include_set = config.params_sets["diff_graph"])
    
    config.params["id_sample"].manual_ticks = True
    m_mongoloids_or = None
    m_mongoloids_and = None
    m_siblings_or = None
    m_siblings_and = None
    m_mothers_or = None
    m_mothers_and = None
    print diff_graph_path
    for thr_p in config.params["thr_p"]:
        for id_sample in config.params["id_sample"]:
            #print id_sample, '/', config.params["id_sample"].num_ticks
            print config.ofname(["graphs", "g"], ext = ".npz", include_set = config.params_sets["graph"])
            data = np.load(config.ofname(["graphs", "g"], ext = ".npz", include_set = config.params_sets["graph"]))
            g = data['G']
            g = np.unpackbits(g, axis = 1)[:, :config.params["num_genes"].value].astype(np.bool)
            if id_sample in config.params["mongoloids_mask"].value:
                if m_mongoloids_or is None:
                    m_mongoloids_or = g
                    m_mongoloids_and = g
                else:
                    m_mongoloids_or = np.logical_or(g, m_mongoloids_or)
                    m_mongoloids_and = np.logical_and(g, m_mongoloids_and)
            
            if id_sample in config.params["siblings_mask"].value:
                if m_siblings_or is None:
                    m_siblings_or = g
                    m_siblings_and = g
                else:
                    m_siblings_or = np.logical_or(g, m_siblings_or)
                    m_siblings_and = np.logical_and(g, m_siblings_and)
                    
            if id_sample in config.params["mothers_mask"].value:
                if m_mothers_or is None:
                    m_mothers_or = g
                    m_mothers_and = g
                else:
                    m_mothers_or = np.logical_or(g, m_mothers_or)
                    m_mothers_and = np.logical_and(g, m_mothers_and)
                    
    #m_diff = np.logical_and(m_mongoloids_and, np.logical_not(m_siblings_or))
    m_diff = np.logical_and(m_mongoloids_and, np.logical_not(m_mothers_or))
    np.savez_compressed(diff_graph_path, 
                        m_diff = m_diff, 
                        m_mongoloids_and = m_mongoloids_and, 
                        m_mongoloids_or = m_mongoloids_or, 
                        m_siblings_and = m_siblings_and, 
                        m_siblings_or = m_siblings_or,
                        m_mothers_or = m_mothers_or,
                        m_mothers_and = m_mothers_and)

def make_diff_graph_age(config, y):
    diff_graph_path = config.ofname(["diff_graph"], ext = ".npz", include_set = config.params_sets["diff_graph"])
    median_age = np.median(y)
    config.params["id_sample"].manual_ticks = True
    num_genes = config.params["num_genes"].value
    num_ages = len(config.params["ages_diff"].value)
    m_ages_or  = np.zeros((num_ages, num_genes, num_genes), dtype = np.bool)
    m_ages_and = np.ones ((num_ages, num_genes, num_genes), dtype = np.bool)
    print diff_graph_path
    #for thr_p in config.params["thr_p"]:
    config.params["thr_p"].set_tick(8);
    for id_sample in config.params["id_sample"]:
        #print id_sample, '/', config.params["id_sample"].num_ticks
        print config.ofname(["graphs", "g"], ext = ".npz", include_set = config.params_sets["graph"])
        data = np.load(config.ofname(["graphs", "g"], ext = ".npz", include_set = config.params_sets["graph"]))
        g = data['G']
        g = np.unpackbits(g, axis = 1)[:, :config.params["num_genes"].value].astype(np.bool)
        y_age = y[id_sample]
        ages_diff = config.params["ages_diff"].value
        for id_age in range(len(ages_diff) - 1):
            if ages_diff[id_age] <= y_age and y_age < ages_diff[id_age + 1]:
                m_ages_and[id_age, :, :] = np.logical_and(m_ages_and[id_age, :, :], g)
                m_ages_or [id_age, :, :] = np.logical_or (m_ages_or [id_age, :, :], g)
                
    np.savez_compressed(diff_graph_path, 
                        m_ages_and = m_ages_and,
                        m_ages_or = m_ages_or)