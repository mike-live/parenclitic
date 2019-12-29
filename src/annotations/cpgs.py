
import pandas as pd
import math
import numpy as np
from collections import defaultdict

def is_intersect(a, b):
    #if not hasattr(a, '__len__'):
    #    a = [a]
    #if not hasattr(b, '__len__'):
    #    b = [b]
    #   print a, b
    if len(a) > len(b):
        a, b = b, a
    if type(b) is set:
        return bool(b.intersection(a))
    return bool(set(a).intersection(b))
        

class cpgs_annotation:
    crit_cols = {
        'cpgs': 'ID_REF',
        'chr': 'CHR',
        'gene': 'UCSC_REFGENE_NAME',
        'geotype': 'RELATION_TO_UCSC_CPG_ISLAND',
        'crossr': 'CROSS_R',
        'class': 'Class',
        'genepart': 'UCSC_REFGENE_GROUP'
    }
    crit_types = ['in', 'ex', 'out']
    
    def __init__(self, fname, crit_cols = None):
        if not crit_cols is None:
            self.crit_cols = crit_cols
        self.load(fname)
        for crit_name, crit_col in self.crit_cols.items():
            if not crit_col in self.df_cpgs:
                raise ValueError(str(crit_col) + ' is not in DataFrame fname.')
        
        self.df_cpgs_set = self.df_cpgs.applymap(lambda x: set(x.split(';')) if type(x) is str else [x])
        self.df_cpgs_smart = self.df_cpgs_set.applymap(lambda x: list(x)[0] if type(x) is set and len(x) == 1 else x)

    def load(self, fname):
        self.df_cpgs = pd.read_csv(fname, delimiter='\t')
        return self.df_cpgs
    
    def get_crit_mask(self, crit_col, crit_type, crit_list):
        nan = any(type(x) is float and math.isnan(x) for x in crit_list)
        crit_list = set(crit_list)
        col = self.crit_cols[crit_col]
        if crit_type == 'in':
            crit_mask = self.df_cpgs_set[col].map(lambda x: is_intersect(x, crit_list))
            if nan:
                crit_mask = crit_mask | self.df_cpgs_set[col].map(lambda x: (type(x) is float and math.isnan(x)))
        elif crit_type == 'ex' or crit_type == 'out':
            crit_mask = self.df_cpgs_set[col].map(lambda x: not (is_intersect(x, crit_list)))
            if nan:
                crit_mask = crit_mask & self.df_cpgs_set[col].map(lambda x: not (type(x) is float and math.isnan(x)))
        return crit_mask

    def split_crit_name(self, crit_name):
        pos = crit_name.rfind('_')
        return crit_name[:pos], crit_name[pos + 1:]


    # O(n*k*log(m)*l)
    # n - number of overall cpgs (485578)
    # k - number of genes for one cpg (~2)
    # m - sum of length of all criterion lists (<3n)
    # l - number of criterions (<10)
    # O(2*n*log(n)) ~ 1e7 op. (per one criterion)
    # example:
    # cpgs.get_crits_mask({'cpgs_in': ['cg00000905']})
    # cpgs.get_crits_mask({'gene_in': ['TLR2', 'KLHL29']}) # 0.5 sec
    # cpgs.get_crits_mask({'gene_out': ['TLR2', 'KLHL29']})
    # cpgs.get_crits_mask({'gene_out': ['TLR2', 'KLHL29'], 'chr_in': ['1','2','22']})
    def get_crits_mask(self, criterions = None):
        crit_all = None
        if not criterions is None:
            for crit_name, crit_list in criterions.items():
                if type(crit_list) is str or not hasattr(type(crit_list), '__iter__'):
                    crit_list = [crit_list]
                crit_col, crit_type = self.split_crit_name(crit_name)
                if crit_col in self.crit_cols and crit_type in self.crit_types:
                    crit_mask = self.get_crit_mask(crit_col, crit_type, crit_list)
                    if crit_all is None:
                        crit_all = crit_mask
                    else:
                        crit_all = crit_all & crit_mask    
        if crit_all is None:
            crit_all = np.ones((len(self.df_cpgs_set),), dtype = np.bool)
        return crit_all
        
    # cur = cpgs.get_sub_frame({'chr_in': ['22'], 'geotype_ex': ['Island']})
    def get_sub_frame(self, criterions = None, original = False):
        df = self.df_cpgs if original else self.df_cpgs_smart
        return df.loc[self.get_crits_mask(criterions)]
                
    # cur, ids = cpgs.get_cpgs({'chr_in': ['22'], 'geotype_ex': ['Island'], 'gene_ex': [np.NaN]})
    def get_cpgs(self, criterions = None):
        cpgs = self.df_cpgs[self.crit_cols['cpgs']][self.get_crits_mask(criterions)]
        return cpgs.values, cpgs.index.values
    
    def get_crit_col_values(self, crit_col, criterions = None, original = False):
        df = self.get_sub_frame(criterions = criterions, original = original)
        values = list(df[self.crit_cols[crit_col]])
        values = list(map(lambda x: x if type(x) is list else (list(x) if type(x) is set else [x]), values))
        values = [y for x in values for y in x]
        return list(set(values))

    def aggregate(self, crit_col_key, crit_col_value, criterions = None, original = False):
        df = self.get_sub_frame(criterions = criterions, original = original)
        d = defaultdict(list)
        for index, x in df.iterrows():
            keys = x[self.crit_cols[crit_col_key]]
            keys = list(keys) if type(keys) is set else [keys]
            values = x[self.crit_cols[crit_col_value]]
            values = list(values) if type(values) is set else [values]
            for key in keys:
                for value in values:
                    d[key].append(value)
        for key in d:
            d[key] = list(set(d[key]))
        return d

