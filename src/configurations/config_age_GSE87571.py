#from slurm_info import info
from .ws_info import info
from infrastructure.configuration import configuration, param
import collections
import numpy as np
from pathlib2 import Path

params = collections.OrderedDict([
    ("num_genes", param(14768, name = 'num_genes')), # 15024, 20270
    ("algorithm", param('pdf', name = 'algorithm')), # svc, kde
    ("thr_type", param('best', name = 'thr_type')), # best, one
    ("division_rule", param('non_control', name = 'division_rule')), # non_control, atypical
    ("thr_p", param(0.88, name = 'thr_p')),
    #("by_group", param(True, name = 'by_group')),
    ("min_score", param(0.85, name = 'min_score')),
    ("age_delimiter", param(28, name = 'age_delimiter')),
#    ("num_groups", param(4, name = 'num_groups')),  
#    ("age_group", param(value_be = 1, value_en = 4, num_ticks = 4, name = 'age_group')),  
    ("by_group", param(True, name = 'by_group')),
    ("id_part", param(value_be = 0, value_en = 29, num_ticks = 30, name = 'id_part')),
    ("id_sample", param(value_be = 0, value_en = 728, num_ticks = 729, name = 'id_sample')),
    ("num_parts", param(30, name = 'num_parts')),
    ("num_workers", param(10, name = 'num_workers')),
    ("num_samples", param(729, name = 'num_samples')),
    ("young_mask", param([], name = 'young_mask')),
    ("old_mask", param([], name = 'old_mask'))
])

files = {
    "gene_chromosome": 'gene_chr.txt',
    "x": 'GSE87571_beta_qf.npz',
    "beta_gene_mean": 'GSE87571_beta_gene_mean.npz',
    "patients_info": "GSE87571_samples.txt",
    "cpgs_annotations": Path("..") / "common" / "cpgs_annotations.txt",
    "betas_dict": "betas_dict.pkl",
    "bad_cpgs": Path("..") / "common" / "bad_cpgs.txt",
    "g": 'graph',
    "kdes": Path('kdes') / 'kdes',
    "graphs": 'graphs',
    "degrees": Path('degrees') / 'degrees',
    "parenclitic": Path("parenclitics") / "parenclitic",
    "degrees_boxplots": "degrees_boxplots",
    "parenclitic_boxplots": "parenclitic_boxplots",
    "parenclitic_ages": "parenclitic_ages",
    "degrees_all": 'degrees_all',
    "parenclitic_all": "parenclitic_all",
    "diff_graph": 'diff_graph',
    "pair_genes": Path('pair_genes') / 'pair_genes',
    "kdes_dist": Path('kdes_dist') / 'kdes_dist',
    "parenclitic_boxplot": Path("parenclitic_boxplots") / "parenclitic_boxplot",
    "parenclitic_age": Path("parenclitic_ages") / "parenclitic_age",
}

params_sets = {
    "graphs": set(['num_genes', 'algorithm', 'age_delimiter', 'num_groups']),
    "graph": set(['num_genes', 'algorithm', 'age_delimiter', 'num_groups', 'id_sample']),
    "degrees": set(['num_genes', 'algorithm', 'age_delimiter', 'num_groups']),
    "parenclitic": set(['num_genes', 'algorithm', 'age_delimiter', 'num_groups']),
    "degrees_sample": set(['num_genes', 'algorithm', 'age_delimiter', 'num_groups', 'id_sample']),
    "parenclitic_sample": set(['num_genes', 'algorithm', 'age_delimiter', 'num_groups', 'id_sample']),
    "degrees_boxplots": set(['num_genes', 'algorithm', 'age_delimiter', 'num_groups']),
    "parenclitic_boxplots": set(['num_genes', 'algorithm', 'age_delimiter', 'num_groups']),
    "diff_graph": set(['num_genes', 'algorithm', 'age_delimiter', 'num_groups']),
    "pair_genes": set(['num_genes', 'algorithm', 'age_delimiter', 'num_groups', 'id_pair']),
    "kdes": set(['num_genes', 'algorithm', 'age_delimiter', 'num_groups', 'id_pair']),
    "parenclitic_boxplot": set(['num_genes', 'algorithm', 'age_delimiter', 'num_groups', 'id_parenclitic']),
}

config = configuration(params, info, files, data_name = 'GSE87571', project_name = 'Gerontology', config_name = 'age', params_sets = params_sets)
