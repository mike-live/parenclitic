#from slurm_info import info
from ws_info import info
from configuration import configuration, param
import collections
import numpy as np
from pathlib2 import Path

params = collections.OrderedDict([
    ("num_genes", param(18484, name = 'num_genes')), 
    ("id_part", param(value_be = 0, value_en = 29, num_ticks = 30, name = 'id_part')),
    ("thr_p", param(value_be = 0.1, value_en = 0.9, num_ticks = 9, name = 'threshold_p')),
    ("id_sample", param(value_be = 0, value_en = 488, num_ticks = 489, name = 'id_sample')),
    ("num_parts", param(30, name = 'num_parts')),
    ("num_workers", param(30, name = 'num_workers')),
    ("num_samples", param(489, name = 'num_samples')),
    ("kde_mask", param('control_mask', name = 'kde_mask')),
    ("health_mask", param([], name = 'health_mask')),
    ("control_mask", param([], name = 'control_mask')),
    ("cancer_mask", param([], name = 'cancer_mask')),
])

files = {
    "patients_id": 'patients_id.csv',
    "classes_blood": 'classes_blood.csv',
    "x": 'gene_mean_blood.csv',

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
    "graphs": set(['kde_mask', 'num_genes', 'id_part', 'num_parts']),
    "graph": set(['kde_mask', 'num_genes', 'thr_p', 'id_sample']),
    "degrees": set(['kde_mask', 'num_genes']),
    "parenclitic": set(['kde_mask', 'num_genes', 'thr_p']),
    "degrees_sample": set(['kde_mask', 'num_genes', 'thr_p', 'id_sample']),
    "parenclitic_sample": set(['kde_mask', 'num_genes', 'thr_p', 'id_sample']),
    "degrees_boxplots": set(['kde_mask', 'num_genes']),
    "parenclitic_boxplots": set(['kde_mask', 'num_genes']),
    "diff_graph": set(['kde_mask', 'num_genes']),
    "pair_genes": set(['kde_mask', 'num_genes', 'id_pair']),
    "kdes": set(['kde_mask', 'num_genes', 'id_pair']),
    "parenclitic_boxplot": set(['kde_mask', 'num_genes', 'thr_p', 'id_parenclitic']),
}

config = configuration(params, info, files, data_name = '2018.09.01', project_name = 'Cancer', config_name = 'cancer', params_sets = params_sets)
