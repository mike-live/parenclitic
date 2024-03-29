import timeit
import numpy as np
from .make_graphs_parts import *
from .infrastructure.configuration import param

#from configurations.load_data_age import load_data_age
#from configurations.ages_config import config

from configurations.load_data_down_GSE52588 import load_data_down_GSE52588_cpgs
from configurations.config_down_GSE52588_cpg import config

#from configurations.load_data_down_GSE52588 import load_data_down_GSE52588
#from configurations.config_down_GSE52588 import config

#from configurations.load_data_mongoloids import load_data_mongoloids_horvath_cpgs
#from configurations.config_mongoloids_cpg_horvath import config


#from configurations.load_data_mongoloids import load_data_mongoloids_hannum_cpgs
#from configurations.mongoloids_cpg_hannum_config import config

#from configurations.load_data_down_GSE63347 import load_data_down_GSE63347_cpg_hannum
#from configurations.config_down_GSE63347_cpg_hannum import config

#from configurations.load_data_down_GSE63347 import load_data_down_GSE63347_cpg_horvath
#from configurations.config_down_GSE63347_cpg_horvath import config

#from configurations.load_data_age_GSE87571 import load_data_age_GSE87571_cpg_horvath
#from configurations.config_age_GSE87571_cpg_horvath import config

#from .configurations.load_data_age_GSE87571 import load_data_age_GSE87571
#from .configurations.config_age_GSE87571 import config

#from configurations.load_data_cancer import load_data_cancer
#from configurations.cancer_config import config

#from configurations.load_data_down_GSE63347 import load_data_down_GSE63347
#from configurations.config_down_GSE63347 import config

#from configurations.load_data_down_GSE74486 import load_data_down_GSE74486
#from configurations.config_down_GSE74486 import config

import sys

if "thr_p" in config.params:
    config.params["thr_p"].manual_ticks = True 
    config.params["thr_p"].whole_values = True

config.params["id_sample"].manual_ticks = True 
 
config.params["id_part"] = param(value_be = 0, value_en = config.info['run_num'] - 1, num_ticks = config.info['run_num'], name = 'id_part')
config.params["num_parts"] = param(config.info['run_num'], name = 'num_parts')
config.upd_ticks()


#X, y, mask, _, _ = load_data_age_GSE87571()
X, y, mask, _ = load_data_down_GSE52588_cpgs()

config.save_params(include_set = config.params_sets["graphs"])
# config.params["thr_p"].get_values()

G, D, IDS = make_graphs_part(mask, X, y, None, config.params["id_part"].value, 
                                               config.params["num_parts"].value, 
                                               num_workers = config.params["num_workers"].value, 
                                               algo = config.params["algorithm"].value,
                                               min_score = config.params["min_score"].value,
                                               by_group = ("by_group" in config.params) & (config.params['by_group'].value))
np.savez_compressed(config.ofname(["graphs", "g"], ext = ".npz", include_set = config.params_sets["graphs"]), G = G, D = D, IDS = IDS)