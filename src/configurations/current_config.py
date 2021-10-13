from infrastructure.configuration import param

#from configurations.load_data_age_GSE87571 import load_data_age_GSE87571_cpg_horvath
#from configurations.config_age_GSE87571_cpg_horvath import config

#from configurations.load_data_down_GSE63347 import load_data_down_GSE63347
#from configurations.config_down_GSE63347 import config

#from configurations.load_data_mongoloids import load_data_mongoloids_horvath_cpgs
#from configurations.config_mongoloids_cpg_horvath import config

from configurations.load_data_down_GSE52588 import load_data_down_GSE52588_cpgs
from configurations.config_down_GSE52588_cpg import config

#from configurations.load_data_age_GSE87571 import load_data_age_GSE87571_cpgs
#from configurations.config_age_GSE87571_cpg import config

#from configurations.load_data_age_GSE87571 import load_data_age_GSE87571
#from configurations.config_age_GSE87571 import config

#from configurations.load_data_down_GSE52588 import load_data_down_GSE52588
#from configurations.config_down_GSE52588 import config

#from configurations.load_data_down_GSE100825 import load_data_down_GSE100825
#from configurations.config_down_GSE100825 import config

#from configurations.load_data_down_GSE131752 import load_data_down_GSE131752
#from configurations.config_down_GSE131752 import config

def get_data():
    #X, y, _, features_names, _ = load_data_age_GSE87571_cpg_horvath()
    #X, y, _, features_names = load_data_down_GSE63347()
    #X, y, _, features_names = load_data_mongoloids_horvath_cpgs()
    #X, y, _, features_names = load_data_down_GSE52588_cpgs()
    #X, y, _, features_names,  = load_data_age_GSE87571_cpgs()
    #return X, y, _, features_names, _
    return load_data_down_GSE52588_cpgs(True)
    #return load_data_down_GSE52588()
    
def get_config():
    return config
