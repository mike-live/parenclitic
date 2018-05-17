import os.path
import math
import numpy as np
import pandas as pd
import scipy.stats as stats
from dicts import get_dict_cpg_gene

def generate_genes_mean_and_std(data_path):
    print_rate = 1000

    dict_cpg_gene = get_dict_cpg_gene(data_path)

    fn = 'ages.txt'
    ages = []
    full_path = data_path + fn
    with open(full_path) as f:
        for line in f:
            ages.append(int(line))

    fn = 'average_beta.txt'
    full_path = data_path + fn
    f = open(full_path)
    first_line = f.readline()
    col_names = first_line.split('\t')

    num_lines = 0
    gene_raw_dict = {}
    for line in f:

        col_vals = line.split('\t')
        CpG = col_vals[0]
        vals = list(map(float, col_vals[1::]))

        genes = dict_cpg_gene.get(CpG)

        if genes is not None:
            for gene in genes:
                if gene in gene_raw_dict:
                    for list_id in range(0, len(ages)):
                        gene_raw_dict[gene][list_id].append(vals[list_id])
                else:
                    gene_raw_dict[gene] = []
                    for list_id in range(0, len(ages)):
                        gene_raw_dict[gene].append([vals[list_id]])

        num_lines += 1
        if num_lines % print_rate == 0:
            print('num_lines: ' + str(num_lines))

    gene_mean_dict = {}
    gene_std_dict = {}
    gene_mean_str_list = []
    gene_std_str_list = []
    num_genes = 0
    for gene in gene_raw_dict:

        mean_list = []
        std_list = []
        for curr_list in gene_raw_dict[gene]:
            mean_list.append(np.mean(curr_list))
            std_list.append(np.std(curr_list))

        gene_mean_dict[gene] = mean_list
        gene_std_dict[gene] = std_list

        curr_mean_str = gene
        curr_std_str = gene
        for id in range(0, len(mean_list)):
            curr_mean_str += (' ' + str(format(mean_list[id], '0.8e')))
            curr_std_str += (' ' + str(format(std_list[id], '0.8e')))

        gene_mean_str_list.append(curr_mean_str)
        gene_std_str_list.append(curr_std_str)

        num_genes += 1
        if num_genes % print_rate == 0:
            print('num_genes: ' + str(num_genes))

    np.savetxt(data_path + 'gene_mean.txt', gene_mean_str_list, fmt="%s")
    np.savetxt(data_path + 'gene_std.txt', gene_std_str_list, fmt="%s")


