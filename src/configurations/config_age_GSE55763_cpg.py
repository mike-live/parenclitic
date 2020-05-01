#from configurations.slurm_info import info
from configurations.ws_info import info
from infrastructure.configuration import configuration, param
import collections
import numpy as np
from pathlib2 import Path

params = collections.OrderedDict([
    ("num_cpgs", param(431906, name = 'num_cpgs')), 
    ("algorithm", param('pdf', name = 'algorithm')), # svc, kde
    ("thr_type", param('best', name = 'thr_type')), # best, one
    ("division_rule", param('non_control', name = 'division_rule')), # non_control, atypical
    ("min_score", param(0.85, name = 'min_score')),
    ("max_score_1d", param(0.6, name = 'max_score_1d')),
    ("age_delimiter", param(38, name = 'age_delimiter')),
    #("num_groups", param(4, name = 'num_groups')),  
    ("id_part", param(value_be = 0, value_en = 29, num_ticks = 30, name = 'id_part')),
    ("id_sample", param(value_be = 0, value_en = 2710, num_ticks = 2711, name = 'id_sample')),
    ("num_parts", param(2, name = 'num_parts')),
    ("num_workers", param(2, name = 'num_workers')),
    ("num_samples", param(2711, name = 'num_samples')),
    ("young_mask", param([], name = 'young_mask')),
    ("old_mask", param([], name = 'old_mask')),
])

files = {
    "cpgs": "GSE55763_betas.npz",
    "cpgs_names": "cpgs_names.tsv",
    "patients_info": "GSE55763_samples.csv",
    "cpgs_annotations": Path("..") / "common" / "cpgs_annotations.txt",
    "bad_cpgs": Path("..") / "common" / "bad_cpgs.txt",
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
    "graphs": set(['input', 'kde_mask', 'num_cpgs', 'algorithm', 'thr_type', 'division_rule', 'age_delimiter', 'id_part', 'num_parts']),
    "graph": set(['input', 'kde_mask', 'num_cpgs', 'algorithm', 'thr_type', 'division_rule', 'age_delimiter', 'id_sample']),
    "degrees": set(['input', 'kde_mask', 'num_cpgs', 'algorithm', 'thr_type', 'division_rule', 'age_delimiter']),
    "parenclitic": set(['input', 'kde_mask', 'num_cpgs', 'algorithm', 'thr_type', 'division_rule', 'age_delimiter']),
    "degrees_sample": set(['input', 'kde_mask', 'num_cpgs', 'algorithm', 'thr_type', 'division_rule', 'age_delimiter', 'id_sample']),
    "parenclitic_sample": set(['input', 'kde_mask', 'num_cpgs', 'algorithm', 'thr_type', 'division_rule', 'age_delimiter', 'id_sample']),
    "degrees_boxplots": set(['input', 'kde_mask', 'num_cpgs', 'algorithm', 'thr_type', 'division_rule', 'age_delimiter']),
    "parenclitic_boxplots": set(['input', 'kde_mask', 'num_cpgs', 'algorithm', 'thr_type', 'division_rule', 'age_delimiter']),
    "diff_graph": set(['input', 'kde_mask', 'num_cpgs', 'algorithm', 'thr_type', 'division_rule', 'age_delimiter']),
    "pair_genes": set(['input', 'kde_mask', 'num_cpgs', 'algorithm', 'thr_type', 'division_rule', 'age_delimiter', 'id_pair']),
    "kdes": set(['input', 'kde_mask', 'num_cpgs', 'algorithm', 'thr_type', 'division_rule', 'age_delimiter', 'id_pair']),
    "parenclitic_boxplot": set(['input', 'kde_mask', 'num_cpgs', 'algorithm', 'thr_type', 'division_rule', 'age_delimiter', 'id_parenclitic']),
}

config = configuration(params, info, files, data_name = 'GSE55763', project_name = 'Gerontology', config_name = 'age', params_sets = params_sets)
