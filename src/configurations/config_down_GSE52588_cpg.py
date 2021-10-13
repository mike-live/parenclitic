#from slurm_info import info
from .ws_info import info
from infrastructure.configuration import configuration, param
import collections
import numpy as np
from pathlib2 import Path

params = collections.OrderedDict([
    ("num_cpgs", param(422801, name = 'num_cpgs')), # 150254 114674 422801
    ("normalization", param('fn', name = 'normalization')), # qf fn
    ("geotypes", param(['Island'], name = 'geotypes')), # ONLY ISLANDS!!!
    ("kde_mask", param('siblings_mask', name = 'kde_mask')),
    ("algorithm", param('pdf', name = 'algorithm')), # svc, kde
    ("thr_type", param('best', name = 'thr_type')), # best, one
    ("division_rule", param('non_control', name = 'division_rule')), # non_control, atypical
    #("thr_p", param(0.88, name = 'thr_p')),
    #("by_group", param(True, name = 'by_group')),
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
    #"x": 'average_beta.txt',
    "x": 'GSE52588_beta_fn.txt',
    "beta_values": 'GSE52588_beta_fn.npz',
    "patients_info": "GSE52588_samples.txt",
    "cells": "Down_beta_Funnorm_filtered_samples_MAge_Horvath_NEW_AgeCalculator.output.csv", 
        # Down_beta_Funnorm_filtered_samples_MAge_Horvath_NEW_AgeCalculator.output.csv 
        # Down_beta_Quantile_filtered_samples.output.csv 
        # Down_beta_Noob_filtered_samples.output.csv
    "epimutations": "epimutations.txt",
    "down_related_2015_cpgs": "down_related_2015_cpgs.txt",
    "var_diff_pv": "var_diff_pv.npz",
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
    "down_phenotypes": "down_phenotypes",
    "down_epimutations": "down_epimutations",
    "age_related": Path("..") / "common" / "age_related",
    "sex_related": Path("..") / "common" / "sex_related",
    "singmann_sex_related": Path("..") / "common" / "sex_related" / "Singmann_metaanalysis_sex_Bonf_0.05.tsv",
    "stripped_graphs": "stripped_graphs",
    "stripped_parenclitic": Path("stripped_parenclitics") / "parenclitic",
    "parenclitics_sieved": Path("parenclitics_sieved") / "parenclitic",
}

params_sets = {
    "graphs": set(['kde_mask', 'num_cpgs', 'algorithm', 'by_group', 'thr_p', 'thr_type', 'division_rule', 'id_part', 'num_parts']),
    "graph": set(['kde_mask', 'num_cpgs', 'algorithm', 'by_group', 'thr_p', 'thr_type', 'division_rule', 'id_sample']),
    "degrees": set(['kde_mask', 'num_cpgs', 'algorithm', 'by_group', 'thr_p', 'thr_type', 'division_rule']),
    "parenclitic": set(['kde_mask', 'num_cpgs', 'algorithm', 'by_group', 'thr_p', 'thr_type', 'division_rule']),
    "degrees_sample": set(['kde_mask', 'num_cpgs', 'algorithm', 'by_group', 'thr_p', 'thr_type', 'division_rule', 'id_sample']),
    "parenclitic_sample": set(['kde_mask', 'num_cpgs', 'algorithm', 'by_group', 'thr_p', 'thr_type', 'division_rule', 'id_sample']),
    "degrees_boxplots": set(['kde_mask', 'num_cpgs', 'algorithm', 'by_group', 'thr_p', 'thr_type', 'division_rule']),
    "parenclitic_boxplots": set(['kde_mask', 'num_cpgs', 'algorithm', 'by_group', 'thr_p', 'thr_type', 'division_rule']),
    "diff_graph": set(['kde_mask', 'num_cpgs', 'algorithm', 'by_group', 'thr_p', 'thr_type', 'division_rule']),
    "pair_genes": set(['kde_mask', 'num_cpgs', 'algorithm', 'by_group', 'thr_p', 'thr_type', 'division_rule', 'id_pair']),
    "kdes": set(['kde_mask', 'num_cpgs', 'algorithm', 'by_group', 'thr_p', 'thr_type', 'division_rule', 'id_pair']),
    "parenclitic_boxplot": set(['kde_mask', 'num_cpgs', 'algorithm', 'by_group', 'thr_p', 'thr_type', 'division_rule', 'id_parenclitic']),
    "all_cpg_names": set(['num_cpgs']),
    "down_phenotypes": set(['normalization', 'num_cpgs']),
    "down_epimutations": set(['normalization', 'kde_mask', 'num_cpgs', 'algorithm', 'by_group', 'thr_p', 'thr_type', 'division_rule']),
    "variance": set(['normalization', 'num_cpgs']),
    "age_related": set([]),
}

config = configuration(params, info, files, data_name = 'GSE52588', project_name = 'Gerontology', config_name = 'down_syndrome_cpgs', params_sets = params_sets)
