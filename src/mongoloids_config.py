#from slurm_info import info
from ws_info import info
from configuration import configuration, param
import collections
import numpy as np
from pathlib2 import Path

params = collections.OrderedDict([
    ("num_genes", param(15024, name = 'num_genes')), # 15024, 20270
    ("kde_mask", param('mothers_mask', name = 'kde_mask')),
    ("id_part", param(value_be = 0, value_en = 29, num_ticks = 30, name = 'id_part')),
    ("thr_p", param(value_be = 0.1, value_en = 0.9, num_ticks = 9, name = 'threshold_p')),
    ("id_sample", param(value_be = 0, value_en = 86, num_ticks = 87, name = 'id_sample')),
    ("num_parts", param(30, name = 'num_parts')),
    ("num_workers", param(30, name = 'num_workers')),
    ("num_samples", param(87, name = 'num_samples')),
    ("mongoloids_mask", param(np.arange(0, 29), name = 'mongoloids_mask')),
    ("siblings_mask", param(np.arange(29, 58), name = 'siblings_mask')),
    ("mothers_mask", param(np.arange(58, 87), name = 'mothers_mask')),
])

files = {
    "x": 'gene_mean_islands_shores.txt',
    "name_genes": 'linreg_genes_mean_islands_shores',
    #"ranged_genes": 'linreg_genes_mean_islands_shores.txt',
    "g": 'graph',
    "kdes": 'kdes',
    "graphs": 'graphs',
    "degrees": Path('degrees') / 'degrees',
    "parenclitic": Path("parenclitics") / "parenclitic",
    "degrees_boxplots": "degrees_boxplots",
    "parenclitic_boxplots": "parenclitic_boxplots",
    "degrees_all": 'degrees_all',
    "parenclitic_all": "parenclitic_all",
    "diff_graph": 'diff_graph',
}

params_sets = {
    "graphs": set(['kde_mask', 'num_genes', 'id_part', 'num_parts']),
    "graph": set(['kde_mask', 'num_genes', 'thr_p', 'id_sample']),
    "degrees": set(['kde_mask', 'num_genes']),
    "parenclitic": set(['kde_mask', 'num_genes']),
    "degrees_sample": set(['kde_mask', 'num_genes', 'thr_p', 'id_sample']),
    "parenclitic_sample": set(['kde_mask', 'num_genes', 'thr_p', 'id_sample']),
    "degrees_boxplots": set(['kde_mask', 'num_genes']),
    "parenclitic_boxplots": set(['kde_mask', 'num_genes']),
    "diff_graph": set(['kde_mask', 'num_genes']),
}

config = configuration(params, info, files, data_name = 'GSE52588', project_name = 'Gerontology', config_name = 'mongoloids', params_sets = params_sets)
