
from transform_data import parenclitic_feature_names

def prepare_parenclitic(parenclitics, is_lists = False):
    numerics = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
    if is_lists:
        parenclitics = parenclitics.select_dtypes(exclude=numerics) # includ exclude
        parenclitics = parenclitics.rename(columns = {'weights': "weight"})
        parenclitics = parenclitics.drop(columns = ["id_components", "component_sizes", "weight"]) #   "weights" "component_sizes", 
    else:
        parenclitics = parenclitics.select_dtypes(include=numerics) # includ exclude
        renamed_columns = {"min_weights": "min_weight", 
                           "max_weights": "max_weight", 
                           "mean_weights": "mean_weight", 
                           "std_weights": "std_weight", 
                           "zeros_weights": "zeros_weight",
                           "nonzeros_weights": "nonzeros_weight"}
        parenclitics = parenclitics.rename(columns = renamed_columns)

    #print (parenclitics.columns.values, len(parenclitics.columns.values))
    parenclitic_names = parenclitic_feature_names()
    parenclitic_names = [parenclitic_names[name] for name in parenclitics.columns.values]
    return parenclitics, parenclitic_names