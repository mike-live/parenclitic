#from slurm_info import info
from .ws_info import info
from infrastructure.configuration import configuration, param
import collections
import numpy as np
from pathlib2 import Path

params = collections.OrderedDict([
    ("num_genes", param(14756, name = 'num_genes')), # 15024, 20270
    ("id_part", param(value_be = 0, value_en = 29, num_ticks = 30, name = 'id_part')),
    ("id_sample", param(value_be = 0, value_en = 86, num_ticks = 87, name = 'id_sample')),
    ("num_parts", param(30, name = 'num_parts')),
    ("num_workers", param(10, name = 'num_workers')),
    ("num_samples", param(87, name = 'num_samples')),
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
    "graphs": set(['num_genes', 'id_part', 'num_parts']),
    "graph": set(['num_genes', 'id_sample']),
    "degrees": set(['num_genes']),
    "parenclitic": set(['num_genes']),
    "degrees_sample": set(['num_genes', 'id_sample']),
    "parenclitic_sample": set(['num_genes', 'id_sample']),
    "degrees_boxplots": set(['num_genes']),
    "parenclitic_boxplots": set(['num_genes']),
    "diff_graph": set(['num_genes']),
    "pair_genes": set(['num_genes', 'id_pair']),
    "kdes": set(['num_genes', 'id_pair']),
    "parenclitic_boxplot": set(['num_genes', 'id_parenclitic']),
    "down_phenotypes": set(['num_genes']),
}

config = configuration(params, info, files, data_name = 'random_graph', project_name = 'Gerontology', config_name = 'random_graph', params_sets = params_sets)
