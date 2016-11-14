from HydroModelBuilder.HydroModelBuilder.GWModelManager import GWModelManager

# MM is short for model manager

grid_resolution = '1000'

MM = GWModelManager()
MM.load_GW_model(r"C:\Workspace\part0075\MDB modelling\testbox\02_transient_flow\structured_model_grid_" + grid_resolution + r"m\02_transient_flow_packaged.pkl")

name = MM.GW_build.keys()[0]

pest_folder = r"C:\Workspace\part0075\MDB modelling\testbox\PEST"    

MM.setupPEST(name, 
             directory=pest_folder, 
             csv_copy=True, models_ID=[name]) 
       
MM.PEST.genParameters(method='csv')
MM.PEST.genPESTpgp()
MM.PEST.genPestfiles(models_ID=[name])       
