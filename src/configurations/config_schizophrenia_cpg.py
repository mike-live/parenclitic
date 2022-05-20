#from slurm_info import info
from .ws_info import info
from infrastructure.configuration import configuration, param
import collections
import numpy as np
from pathlib2 import Path

params = collections.OrderedDict([
    ("num_cpgs", param(422801, name = 'num_cpgs')), # 150254 114674 422801
    #("normalization", param('fn', name = 'normalization')), # qf fn
    #("geotypes", param(['Island'], name = 'geotypes')), # ONLY ISLANDS!!!
    ("kde_mask", param('control_mask', name = 'kde_mask')),
    ("algorithm", param('pdf', name = 'algorithm')), # svc, kde
    ("thr_type", param('best', name = 'thr_type')), # best, one
    ("division_rule", param('non_control', name = 'division_rule')), # non_control, atypical
    #("LOO", param(value_be = 0, value_en = 28, num_ticks = 29, name = 'LOO')),
    #("thr_p", param(0.88, name = 'thr_p')),
    #("by_group", param(True, name = 'by_group')),
    ("min_score", param(0.6, name = 'min_score')),
    ("max_score_1d", param(0.65, name = 'max_score_1d')),
    ("id_part", param(value_be = 0, value_en = 899, num_ticks = 900, name = 'id_part')),
    ("id_sample", param(value_be = 0, value_en = 1521, num_ticks = 1522, name = 'id_sample')),
    ("num_parts", param(900, name = 'num_parts')),
    ("num_workers", param(10, name = 'num_workers')),
    ("num_samples", param(1522, name = 'num_samples')),
    ("control_mask", param([], name = 'control_mask')),
    ("schizophrenia_mask", param([], name = 'schizoprenia_mask')),
])

files = {
    "gene_chromosome": 'gene_chr.txt',
    "beta_values": 'data_trn_val.pkl',
    "patients_info": "pheno_trn_val.xlsx",
    "g": 'graph',
    "cpgs_annotations": Path("..") / "common" / "cpgs_annotations.txt",
    "bad_cpgs": Path("..") / "common" / "bad_cpgs.txt",
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
    "age_related": Path("..") / "common" / "age_related",
    "sex_related": Path("..") / "common" / "sex_related",
    "singmann_sex_related": Path("..") / "common" / "sex_related" / "Singmann_metaanalysis_sex_Bonf_0.05.tsv",
    "stripped_graphs": "stripped_graphs",
    "stripped_parenclitic": Path("stripped_parenclitics") / "parenclitic",
    "parenclitics_sieved": Path("parenclitics_sieved") / "parenclitic",
}

base_params = ['kde_mask', 'num_cpgs'] # 'algorithm', 'by_group', 'thr_p', 'thr_type', 'division_rule'
params_sets = {
    "graphs": set(base_params + ['id_part', 'num_parts']),
    "graph": set(base_params + ['id_sample']),
    "degrees": set(base_params),
    "parenclitic": set(base_params),
    "degrees_sample": set(base_params + ['id_sample']),
    "parenclitic_sample": set(base_params + ['id_sample']),
    "degrees_boxplots": set(base_params),
    "parenclitic_boxplots": set(base_params),
    "diff_graph": set(base_params),
    "pair_genes": set(base_params + ['id_pair']),
    "kdes": set(base_params + ['id_pair']),
    "parenclitic_boxplot": set(base_params + ['id_parenclitic']),
    "all_cpg_names": set(['num_cpgs']),
    "down_phenotypes": set(['normalization', 'num_cpgs']),
    "down_epimutations": set(['normalization'] + base_params),
    "variance": set(['normalization', 'num_cpgs']),
    "age_related": set([]),
}

config = configuration(
    params, info, files, 
    data_name = 'Schizophrenia', 
    project_name = 'Gerontology', 
    config_name = 'schizophrenia_cpg', 
    params_sets = params_sets
)
