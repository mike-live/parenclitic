#from slurm_info import info
from .ws_info import info
from infrastructure.configuration import configuration, param
import collections
import numpy as np
from pathlib2 import Path

params = collections.OrderedDict([
    ("num_genes", param(14756, name = 'num_genes')), # 15024, 20270
    ("algorithm", param('pdf', name = 'algorithm')), # svc, kde
    ("thr_type", param('best', name = 'thr_type')), # best, one
    ("division_rule", param('non_control', name = 'division_rule')), # non_control, atypical
    #("thr_p", param(0.88, name = 'thr_p')),
    #("by_group", param(True, name = 'by_group')),
    ("min_score", param(0.9, name = 'min_score')),
    ("id_sample", param(value_be = 0, value_en = 86, num_ticks = 87, name = 'id_sample')),
    ("num_workers", param(10, name = 'num_workers')),
    ("num_samples", param(87, name = 'num_samples')),
])

files = {
    "gene_chromosome": 'gene_chr.txt',
    #"x": 'gene_mean_islands_shores.txt',
    "x": 'GSE52588_average_beta.txt',
    "horvath_cpgs_beta": "horvath_cpgs_beta.txt",
    "name_genes": 'linreg_genes_mean_islands_shores',
    "good_pairs": "good_pairs.npz",
    "down_phenotypes_table": "DOWN_FENOTIPO_No4,8,12_PerCorrelazioni.tsv",
    "patients_info": "GSE52588_samples.txt",
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
    "down_phenotypes": "down_phenotypes",
    "cpgs_annotations": Path("..") / "common" / "cpgs_annotations.txt",
    "age_related": Path("..") / "common" / "age_related",
}

params_sets = {
    "graphs": set(['kde_mask', 'num_genes', 'algorithm', 'by_group', 'thr_p', 'thr_type', 'division_rule', 'id_part', 'num_parts']),
    "graph": set(['kde_mask', 'num_genes', 'algorithm', 'by_group', 'thr_p', 'thr_type', 'division_rule', 'id_sample']),
    "degrees": set(['kde_mask', 'num_genes', 'algorithm', 'by_group']),
    "parenclitic": set(['kde_mask', 'num_genes', 'algorithm', 'by_group', 'thr_p', 'thr_type', 'division_rule']),
    "degrees_sample": set(['kde_mask', 'num_genes', 'algorithm', 'by_group', 'thr_p', 'thr_type', 'division_rule', 'id_sample']),
    "parenclitic_sample": set(['kde_mask', 'num_genes', 'algorithm', 'by_group', 'thr_p', 'thr_type', 'division_rule', 'id_sample']),
    "degrees_boxplots": set(['kde_mask', 'num_genes', 'algorithm', 'by_group']),
    "parenclitic_boxplots": set(['kde_mask', 'num_genes', 'algorithm', 'by_group', 'thr_p', 'division_rule', 'thr_type']),
    "diff_graph": set(['kde_mask', 'num_genes', 'algorithm', 'by_group','thr_type', 'division_rule']),
    "pair_genes": set(['kde_mask', 'num_genes', 'algorithm', 'by_group', 'thr_p', 'thr_type', 'division_rule', 'id_pair']),
    "kdes": set(['kde_mask', 'num_genes', 'algorithm', 'by_group', 'thr_p', 'thr_type', 'division_rule', 'id_pair']),
    "parenclitic_boxplot": set(['kde_mask', 'num_genes', 'algorithm', 'by_group', 'thr_p', 'thr_type', 'division_rule', 'id_parenclitic']),
    "down_phenotypes": set(['num_genes']),
}

config = configuration(params, info, files, data_name = 'GSE52588', project_name = 'Gerontology', config_name = 'down_syndrome', params_sets = params_sets)
