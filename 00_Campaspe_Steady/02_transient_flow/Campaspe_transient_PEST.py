import sys
import os
sys.path.append('C:\Workspace\part0075\GIT_REPOS')
from HydroModelBuilder.GWModelManager import GWModelManager

# MM is short for model manager
def run(model_folder, pest_folder):
    
    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder,"02_transient_flow_packaged.pkl"))
    
    name = MM.GW_build.keys()[0]

    MM.setupPEST(name, 
                 directory=pest_folder, 
                 csv_copy=True, models_ID=[name]) 
           
    MM.PEST.genParameters(method='csv')
    MM.PEST.genPESTpgp()
    MM.PEST.genPestfiles(models_ID=[name])       

if __name__ ==  "__main__":

    grid_resolution = '5000'
    model_folder = r"C:\Workspace\part0075\MDB modelling\testbox\00_Campaspe_Cascade\02_transient_flow\structured_model_grid_" + grid_resolution + r"m\\" 
    pest_folder = r"C:\Workspace\part0075\MDB modelling\testbox\PEST5000\master"     
    
    run(model_folder, pest_folder)