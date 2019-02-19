from slurm_info import info
#from ws_info import info
from infrastructure.configuration import configuration, param
import collections
import numpy as np
from pathlib2 import Path

params = collections.OrderedDict([
    ("num_cpgs", param(150254, name = 'num_cpgs')), # 150254
    ("kde_mask", param('healthy_mask', name = 'kde_mask')),
    ("algorithm", param('svc', name = 'algorithm')),
    ("geotypes", param(['Island'], name = 'geotypes')),
    ("min_score", param(0.9, name = 'min_score')),
    ("id_part", param(value_be = 0, value_en = 899, num_ticks = 900, name = 'id_part')),
    ("id_sample", param(value_be = 0, value_en = 86, num_ticks = 87, name = 'id_sample')),
    ("num_parts", param(900, name = 'num_parts')),
    ("num_workers", param(10, name = 'num_workers')),
    ("num_samples", param(87, name = 'num_samples')),
    ("mongoloids_mask", param(np.arange(0, 29), name = 'mongoloids_mask')),
    ("siblings_mask", param(np.arange(29, 58), name = 'siblings_mask')),
    ("mothers_mask", param(np.arange(58, 87), name = 'mothers_mask')),
])

files = {
    "gene_chromosome": 'gene_chr.txt',
    "x": 'average_beta.txt',
    "cpgs": "cpgs_annotations.txt",
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
    "graphs": set(['kde_mask', 'num_cpgs', 'algorithm', 'id_part', 'num_parts']),
    "graph": set(['kde_mask', 'num_cpgs', 'algorithm', 'id_sample']),
    "degrees": set(['kde_mask', 'num_cpgs', 'algorithm']),
    "parenclitic": set(['kde_mask', 'num_cpgs', 'algorithm']),
    "degrees_sample": set(['kde_mask', 'num_cpgs', 'algorithm', 'id_sample']),
    "parenclitic_sample": set(['kde_mask', 'num_cpgs', 'algorithm', 'id_sample']),
    "degrees_boxplots": set(['kde_mask', 'num_cpgs', 'algorithm']),
    "parenclitic_boxplots": set(['kde_mask', 'num_cpgs', 'algorithm']),
    "diff_graph": set(['kde_mask', 'num_cpgs', 'algorithm']),
    "pair_genes": set(['kde_mask', 'num_cpgs', 'algorithm', 'id_pair']),
    "kdes": set(['kde_mask', 'num_cpgs', 'algorithm', 'id_pair']),
    "parenclitic_boxplot": set(['kde_mask', 'num_cpgs', 'algorithm', 'id_parenclitic']),
}

config = configuration(params, info, files, data_name = 'GSE52588', project_name = 'Gerontology', config_name = 'down_syndrome_cpgs', params_sets = params_sets)
