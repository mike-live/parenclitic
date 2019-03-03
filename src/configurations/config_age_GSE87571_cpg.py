#from configurations.slurm_info import info
from configurations.ws_info import info
from infrastructure.configuration import configuration, param
import collections
import numpy as np
from pathlib2 import Path

params = collections.OrderedDict([
    ("num_cpgs", param(353, name = 'num_cpgs')), 
    ("input", param('horvath_cpg', name = 'input')),
    ("kde_mask", param('young_mask', name = 'kde_mask')),
    ("algorithm", param('svc', name = 'algorithm')),
    ("delimiter_age", param(65, name = 'delimiter_age')), 
    ("id_part", param(value_be = 0, value_en = 29, num_ticks = 30, name = 'id_part')),
#    ("thr_p", param(value_be = 0.1, value_en = 0.9, num_ticks = 9, name = 'threshold_p')),
    ("id_sample", param(value_be = 0, value_en = 728, num_ticks = 729, name = 'id_sample')),
    ("num_parts", param(2, name = 'num_parts')),
    ("num_workers", param(2, name = 'num_workers')),
    ("num_samples", param(729, name = 'num_samples')),
    ("young_mask", param([], name = 'young_mask')),
    ("old_mask", param([], name = 'old_mask')),
])

files = {
    "gene_chromosome": 'gene_chr.txt',
    "cpgs": "cpg_beta.npz",
    "cpgs_annotations": "cpgs_annotations.txt",
    "cpgs_names": "cpgs_names.npz",
    "patients_info": "GSE87571_samples.txt",
    "name_genes": 'linreg_genes_mean_islands_shores',
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
    "graphs": set(['input', 'kde_mask', 'num_cpgs', 'algorithm', 'id_part', 'num_parts']),
    "graph": set(['input', 'kde_mask', 'num_cpgs', 'algorithm', 'id_sample']),
    "degrees": set(['input', 'kde_mask', 'algorithm', 'num_cpgs']),
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

config = configuration(params, info, files, data_name = 'GSE87571', project_name = 'Gerontology', config_name = 'age', params_sets = params_sets)
