import sys
import os
#sys.path.append('C:\Workspace\part0075\GIT_REPOS')
from HydroModelBuilder.GWModelManager import GWModelManager

# Configuration Loader
from HydroModelBuilder.Utilities.Config.ConfigLoader import ConfigLoader

# MM is short for model manager
def run(model_folder, pest_folder):
    
    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder,"02_transient_flow_packaged.pkl"))
    
    name = list(MM.GW_build.keys())[0]

    MM.setupPEST(name, 
                 directory=pest_folder, 
                 csv_copy=True, models_ID=[name]) 
           
    MM.PEST.genParameters(method='csv')
    MM.PEST.genPESTpgp()
    MM.PEST.genPestfiles(models_ID=[name])       

if __name__ ==  "__main__":

# Get general model config information
    CONFIG = ConfigLoader('../../config/model_config.json')\
                    .set_environment("02_transient_flow")

    verbose=True
                    
    args = sys.argv
    if len(args) > 1:
        model_folder = sys.argv[1]
        data_folder = sys.argv[2]
    else:
        model_config = CONFIG.model_config
        model_folder = model_config['model_folder'] + model_config['grid_resolution'] + os.path.sep
        data_folder = model_config['data_folder']
    
    pest_folder = data_folder     
    
    run(model_folder, pest_folder)