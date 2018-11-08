import timeit
import numpy as np
from make_graphs_parts import *
from infrastructure.configuration import param

#from load_data_age import load_data_age
#from ages_config import config

#from load_data_mongoloids import load_data_mongoloids
#from mongoloids_config import config

#from load_data_mongoloids import load_data_mongoloids_horvath_cpgs
#from mongoloids_cpg_horvath_config import config

#from load_data_mongoloids import load_data_mongoloids_hannum_cpgs
#from mongoloids_cpg_hannum_config import config

#from configurations.load_data_down_GSE63347 import load_data_down_GSE63347_cpg_hannum
#from configurations.config_down_GSE63347_cpg_hannum import config

#from configurations.load_data_down_GSE63347 import load_data_down_GSE63347_cpg_horvath
#from configurations.config_down_GSE63347_cpg_horvath import config

from configurations.load_data_age_GSE87571 import load_data_age_GSE87571_cpg_horvath
from configurations.config_age_GSE87571_cpg_horvath import config

#from load_data_cancer import load_data_cancer
#from cancer_config import config

import sys

config.params["thr_p"].manual_ticks = True 
config.params["thr_p"].whole_values = True
config.params["id_sample"].manual_ticks = True 
 
config.params["id_part"] = param(value_be = 0, value_en = config.info['run_num'] - 1, num_ticks = config.info['run_num'], name = 'id_part')
config.params["num_parts"] = param(config.info['run_num'], name = 'num_parts')
config.upd_ticks()


X, y, X_prob, _ = load_data_age_GSE87571_cpg_horvath()


config.save_params(include_set = config.params_sets["graphs"])
# config.params["thr_p"].get_values()

G = make_graphs_part(X_prob, X, y, None, config.params["id_part"].value, config.params["num_parts"].value, num_workers = config.params["num_workers"].value, algo = config.params["algorithm"].value)
np.savez_compressed(config.ofname(["graphs", "g"], ext = ".npz", include_set = config.params_sets["graphs"]), G = G)