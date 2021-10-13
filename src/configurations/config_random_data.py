#from slurm_info import info
from .ws_info import info
from infrastructure.configuration import configuration, param
import collections
import numpy as np
from pathlib2 import Path

params = collections.OrderedDict([
    ("num_genes", param(14756, name = 'num_genes')), # 15024, 20270
    ("algorithm", param('kde', name = 'algorithm')), # svc, kde
    ("thr_p", param(0.88, name = 'thr_p')),
    #("min_score", param(0.9, name = 'min_score')),
    ("num_genes", param(14756, name = 'num_genes')), # 15024, 20270
    ("id_part", param(value_be = 0, value_en = 29, num_ticks = 30, name = 'id_part')),
    ("id_sample", param(value_be = 0, value_en = 28, num_ticks = 29, name = 'id_sample')),
    ("num_parts", param(30, name = 'num_parts')),
    ("num_workers", param(10, name = 'num_workers')),
    ("num_samples", param(29, name = 'num_samples')),
])

files = {
    "g": 'graph',
    "kdes": Path('kdes') / 'kdes',
    "graphs": 'graphs',
    "degrees": Path('degrees') / 'degrees',
    "parenclitic": Path("parenclitics") / "parenclitic",
    "degrees_boxplots": "degrees_boxplots",
    "parenclitic_boxplots": "parenclitic_boxplots",
    "degrees_all": 'degrees_all',
    "parenclitic_all": "parenclitic_all",
    "diff_graph": 'diff_graph',
    "pair_genes": Path('pair_genes') / 'pair_genes',
    "kdes_dist": Path('kdes_dist') / 'kdes_dist',
    "parenclitic_boxplot": Path("parenclitic_boxplots") / "parenclitic_boxplot",
    "down_phenotypes": "down_phenotypes",
}

params_sets = {
    "graphs": set(['num_genes', 'algorithm', 'thr_p', 'id_part', 'num_parts']),
    "graph": set(['num_genes', 'algorithm', 'thr_p', 'id_sample']),
    "degrees": set(['num_genes', 'algorithm', 'thr_p']),
    "parenclitic": set(['num_genes', 'algorithm', 'thr_p']),
    "degrees_sample": set(['num_genes', 'algorithm', 'thr_p', 'id_sample']),
    "parenclitic_sample": set(['num_genes', 'algorithm', 'thr_p', 'id_sample']),
    "degrees_boxplots": set(['num_genes', 'algorithm', 'thr_p']),
    "parenclitic_boxplots": set(['num_genes', 'algorithm', 'thr_p']),
    "diff_graph": set(['num_genes', 'algorithm', 'thr_p']),
    "pair_genes": set(['num_genes', 'algorithm', 'thr_p', 'id_pair']),
    "kdes": set(['num_genes', 'algorithm', 'thr_p', 'id_pair']),
    "parenclitic_boxplot": set(['num_genes', 'algorithm', 'thr_p', 'id_parenclitic']),
    "down_phenotypes": set(['num_genes', 'algorithm', 'thr_p']),
}

config = configuration(params, info, files, data_name = 'random_data', project_name = 'Gerontology', config_name = 'random_data', params_sets = params_sets)
