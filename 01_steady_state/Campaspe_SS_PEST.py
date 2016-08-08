from HydroModelBuilder.HydroModelBuilder.GWModelManager import GWModelManager

# MM is short for model manager


MM = GWModelManager()
MM.load_GW_model(r"C:\Workspace\part0075\MDB modelling\testbox\01_steady_state\structured_model_grid_10000m\01_steady_state_packaged.pkl")

name = MM.GW_build.keys()[0]

data_folder = r"C:\Workspace\part0075\MDB modelling\testbox\\"    

MM.setupPEST(name, 
             directory=r"C:\Workspace\part0075\MDB modelling\testbox\\", 
             csv_copy=True) 
       
MM.PEST.genParameters(method='csv')
MM.PEST.genPESTpgp()
MM.PEST.genPestfiles(models_ID=[name])       

#print MM.GW_build[name].observations.obs['average head'].keys()
#print MM.GW_build[name].observations.obs['average head']['mapped_observations']


