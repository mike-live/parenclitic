from configurations.slurm_info import info
#from configurations.ws_info import info
from infrastructure.configuration import configuration, param
import collections
import numpy as np
from pathlib2 import Path

params = collections.OrderedDict([
    ("num_cpgs", param(150254, name = 'num_cpgs')), # 150254
    ("kde_mask", param('normal_mask', name = 'kde_mask')),
    ("algorithm", param('svc', name = 'algorithm')),
    ("geotypes", param(['Island'], name = 'geotypes')),
    ("id_part", param(value_be = 0, value_en = 29, num_ticks = 30, name = 'id_part')),
    ("id_sample", param(value_be = 0, value_en = 70, num_ticks = 71, name = 'id_sample')),
    ("num_parts", param(900, name = 'num_parts')),
    ("num_workers", param(10, name = 'num_workers')),
    ("num_samples", param(71, name = 'num_samples')),
    ("down_mask", param([], name = 'down_mask')),
    ("normal_mask", param([], name = 'normal_mask')),
])

files = {
    "gene_chromosome": 'gene_chr.txt',
    "x": 'GSE63347_series_matrix.txt',
    "cpgs": "cpgs_annotations.txt",
    "patients_info": "patients_info.txt",
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
    "graphs": set(['input', 'kde_mask', 'num_cpgs', 'algorithm', 'id_part', 'num_parts']),
    "graph": set(['input', 'kde_mask', 'num_cpgs', 'algorithm', 'id_sample']),
    "degrees": set(['input', 'kde_mask', 'num_cpgs', 'algorithm']),
    "parenclitic": set(['input', 'kde_mask', 'num_cpgs', 'algorithm']),
    "degrees_sample": set(['input', 'kde_mask', 'num_cpgs', 'algorithm', 'id_sample']),
    "parenclitic_sample": set(['input', 'kde_mask', 'num_cpgs', 'algorithm', 'id_sample']),
    "degrees_boxplots": set(['input', 'kde_mask', 'num_cpgs', 'algorithm']),
    "parenclitic_boxplots": set(['input', 'kde_mask', 'num_cpgs', 'algorithm']),
    "diff_graph": set(['input', 'kde_mask', 'num_cpgs', 'algorithm']),
    "pair_genes": set(['input', 'kde_mask', 'num_cpgs', 'algorithm', 'id_pair']),
    "kdes": set(['input', 'kde_mask', 'num_cpgs', 'algorithm', 'id_pair']),
    "parenclitic_boxplot": set(['input', 'kde_mask', 'num_cpgs', 'algorithm', 'id_parenclitic']),
}

config = configuration(params, info, files, data_name = 'GSE63347', project_name = 'Gerontology', config_name = 'down_syndrome', params_sets = params_sets)
