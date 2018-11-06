#from slurm_info import info
from ws_info import info
from configuration import configuration, param
import collections
import numpy as np
from pathlib2 import Path

params = collections.OrderedDict([
    ("num_cpgs", param(71, name = 'num_cpgs')), # 15024, 20270
    ("input", param('hannum_cpg', name = 'input')),
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
    "gene_chromosome": 'gene_chr.txt',
    "hannum_cpgs_beta": "hannum_cpgs_beta.txt",
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
}

params_sets = {
    "graphs": set(['input', 'kde_mask', 'num_cpgs', 'id_part', 'num_parts']),
    "graph": set(['input', 'kde_mask', 'num_cpgs', 'thr_p', 'id_sample']),
    "degrees": set(['input', 'kde_mask', 'num_cpgs']),
    "parenclitic": set(['input', 'kde_mask', 'num_cpgs', 'thr_p']),
    "degrees_sample": set(['input', 'kde_mask', 'num_cpgs', 'thr_p', 'id_sample']),
    "parenclitic_sample": set(['input', 'kde_mask', 'num_cpgs', 'thr_p', 'id_sample']),
    "degrees_boxplots": set(['input', 'kde_mask', 'num_cpgs']),
    "parenclitic_boxplots": set(['input', 'kde_mask', 'num_cpgs']),
    "diff_graph": set(['input', 'kde_mask', 'num_cpgs']),
    "pair_genes": set(['input', 'kde_mask', 'num_cpgs', 'id_pair']),
    "kdes": set(['input', 'kde_mask', 'num_cpgs', 'id_pair']),
    "parenclitic_boxplot": set(['input', 'kde_mask', 'num_cpgs', 'thr_p', 'id_parenclitic']),
}

config = configuration(params, info, files, data_name = 'GSE52588', project_name = 'Gerontology', config_name = 'mongoloids', params_sets = params_sets)
