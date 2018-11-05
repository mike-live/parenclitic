#from slurm_info import info
from ws_info import info
from configuration import configuration, param
import collections
from pathlib2 import Path
import numpy as np

params = collections.OrderedDict([
    ("num_genes", param(20270, name = 'num_genes')), # 20270
    ("id_part", param(value_be = 0, value_en = 29, num_ticks = 30, name = 'id_part')),
    ("thr_p", param(value_be = 0.1, value_en = 0.9, num_ticks = 9, name = 'threshold_p')), 
    ("id_sample", param(value_be = 0, value_en = 655, num_ticks = 656, name = 'id_sample')),
    ("num_parts", param(30, name = 'num_parts')),
    ("num_workers", param(1, name = 'num_workers')),
    ("num_samples", param(656, name = 'num_samples')),
    ("ages_diff", param(np.arange(0, 110, 5), name = 'ages_diff'))
])

files = {"x": 'gene_mean.txt',
         #"ranged_genes": 'pvals_mean_genes.txt',
         "y": 'ages.txt',
         "g": 'graph',
         "kdes": 'kdes',
         "graphs": 'graphs',
         "degrees": Path('degrees') / 'degrees',
         "parenclitic": Path("parenclitics") / "parenclitic",
         "diff_graph": 'diff_graph',
}

params_sets = {
    "graphs": set(['num_genes', 'id_part', 'num_parts']),
    "graph": set(['num_genes', 'thr_p', 'id_sample']),
    "degrees": set(['num_genes', 'thr_p']),
    "parenclitic": set(['num_genes', 'thr_p']),
    "degrees_sample": set(['num_genes', 'thr_p', 'id_sample']),
    "parenclitic_sample": set(['num_genes', 'thr_p', 'id_sample']),
    "diff_graph": set(['num_genes']),
}

config = configuration(params, info, files, data_name = 'GSE40279', project_name = 'Gerontology', config_name = 'ages', params_sets = params_sets)
