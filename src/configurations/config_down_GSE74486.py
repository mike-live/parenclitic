#from configurations.slurm_info import info
from configurations.ws_info import info
from infrastructure.configuration import configuration, param
import collections
import numpy as np
from pathlib2 import Path

params = collections.OrderedDict([
    ("num_genes", param(15009, name = 'num_genes')), # 15009
    ("kde_mask", param('normal_mask', name = 'kde_mask')),
    ("algorithm", param('svc', name = 'algorithm')),
    ("id_part", param(value_be = 0, value_en = 29, num_ticks = 30, name = 'id_part')),
#    ("thr_p", param(value_be = 0.1, value_en = 0.9, num_ticks = 9, name = 'threshold_p')),
    ("id_sample", param(value_be = 0, value_en = 70, num_ticks = 71, name = 'id_sample')),
    ("num_parts", param(30, name = 'num_parts')),
    ("num_workers", param(30, name = 'num_workers')),
    ("num_samples", param(71, name = 'num_samples')),
    ("down_mask", param([], name = 'down_mask')),
    ("normal_mask", param([], name = 'normal_mask')),
])

files = {
    "x": 'gene_data.txt',
    "patients_info": "GSE74486_header.csv",
    #"ranged_genes": 'linreg_genes_mean_islands_shores.txt',
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
    "graphs": set(['kde_mask', 'num_genes', 'algorithm', 'id_part', 'num_parts']),
    "graph": set(['kde_mask', 'num_genes', 'algorithm', 'id_sample']),
    "degrees": set(['kde_mask', 'num_genes', 'algorithm']),
    "parenclitic": set(['kde_mask', 'num_genes', 'algorithm']),
    "degrees_sample": set(['kde_mask', 'num_genes', 'algorithm', 'id_sample']),
    "parenclitic_sample": set(['kde_mask', 'num_genes', 'algorithm', 'id_sample']),
    "degrees_boxplots": set(['kde_mask', 'num_genes', 'algorithm']),
    "parenclitic_boxplots": set(['kde_mask', 'num_genes', 'algorithm']),
    "diff_graph": set(['kde_mask', 'num_genes', 'algorithm']),
    "pair_genes": set(['kde_mask', 'num_genes', 'algorithm', 'id_pair']),
    "kdes": set(['kde_mask', 'num_genes', 'algorithm', 'id_pair']),
    "parenclitic_boxplot": set(['kde_mask', 'num_genes', 'algorithm', 'id_parenclitic']),
}

config = configuration(params, info, files, data_name = 'GSE74486', project_name = 'Gerontology', config_name = 'down_syndrome', params_sets = params_sets)
