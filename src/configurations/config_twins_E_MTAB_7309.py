#from slurm_info import info
from .ws_info import info
from infrastructure.configuration import configuration, param
import collections
import numpy as np
from pathlib2 import Path

params = collections.OrderedDict([
    ("num_genes", param(14743, name = 'num_genes')), # 15024, 20270
    ("kde_mask", param('less_65_mask', name = 'kde_mask')),
    ("algorithm", param('svc', name = 'algorithm')), # svc, kde
    #("thr_p", param(0.9, name = 'thr_p')),
    ("min_score", param(0.7, name = 'min_score')),
    ("id_part", param(value_be = 0, value_en = 29, num_ticks = 30, name = 'id_part')),
    ("id_sample", param(value_be = 0, value_en = 1058, num_ticks = 1059, name = 'id_sample')),
    ("num_parts", param(30, name = 'num_parts')),
    ("num_workers", param(10, name = 'num_workers')),
    ("num_samples", param(1059, name = 'num_samples')),
    ("less_65_mask", param([], name = 'less_65_mask')),
    ("greater_80_mask", param([], name = 'greater_80_mask')),
    ("between_65_80_mask", param([], name = 'between_65_80_mask')),
])

files = {
    "x": 'genes_norm(quantile)_source(betas).npz',
    "gene_names": 'genes_list.txt',
    "patients_info": 'observables.csv',
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
    "graphs": set(['kde_mask', 'num_genes', 'algorithm', 'thr_p', 'id_part', 'num_parts']),
    "graph": set(['kde_mask', 'num_genes', 'algorithm', 'thr_p', 'id_sample']),
    "degrees": set(['kde_mask', 'num_genes', 'algorithm']),
    "parenclitic": set(['kde_mask', 'num_genes', 'algorithm', 'thr_p']),
    "degrees_sample": set(['kde_mask', 'num_genes', 'algorithm', 'thr_p', 'id_sample']),
    "parenclitic_sample": set(['kde_mask', 'num_genes', 'algorithm', 'thr_p', 'id_sample']),
    "degrees_boxplots": set(['kde_mask', 'num_genes', 'algorithm']),
    "parenclitic_boxplots": set(['kde_mask', 'num_genes', 'algorithm', 'thr_p']),
    "diff_graph": set(['kde_mask', 'num_genes', 'algorithm']),
    "pair_genes": set(['kde_mask', 'num_genes', 'algorithm', 'thr_p', 'id_pair']),
    "kdes": set(['kde_mask', 'num_genes', 'algorithm', 'thr_p', 'id_pair']),
    "parenclitic_boxplot": set(['kde_mask', 'num_genes', 'algorithm', 'thr_p', 'id_parenclitic']),
    "down_phenotypes": set(['num_genes']),
}

config = configuration(params, info, files, data_name = 'E-MTAB-7309', project_name = 'Gerontology', config_name = 'old_twins', params_sets = params_sets)
