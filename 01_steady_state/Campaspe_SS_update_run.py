import sys
sys.path.append('C:\Workspace\part0075\GIT_REPOS')


from HydroModelBuilder.HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface
# MM is short for model manager


MM = GWModelManager()
MM.load_GW_model(r"C:\Workspace\part0075\MDB modelling\testbox\01_steady_state\structured_model_grid_10000m\01_steady_state_packaged.pkl")

name = MM.GW_build.keys()[0]
# modify data folder
#MM.GW_build['Campaspe'].data_folder
#modify output folder
#MM.GW_build['Campaspe'].out_data_folder

data_folder = r"C:\Workspace\part0075\MDB modelling\testbox\\"    

# Load in the new parameters based on parameters.txt or dictionary of new parameters
#MM.updateParameters('Campaspe', 'parameters.txt')

print "************************************************************************"
print " Updating HGU parameters "


print "************************************************************************"
print " Updating river parameters "

mapped_river = MM.GW_build[name].load_obj(r"C:\Workspace\part0075\MDB modelling\testbox\01_steady_state\structured_model_grid_10000m\Campaspe_Riv_model.shp_mapped.pkl")

simple_river = []
riv_width_avg = 10.0 #m
riv_bed_thickness = 0.10 #m
for riv_cell in mapped_river: #MM.GW_build[name].polyline_mapped['Campaspe_Riv_model.shp']:
    row = riv_cell[0][0]
    col = riv_cell[0][1]
    if MM.GW_build[name].model_mesh3D[1][0][row][col] == -1:
        continue
    #print test_model.model_mesh3D
    stage = MM.GW_build[name].model_mesh3D[0][0][row][col]
    bed = MM.GW_build[name].model_mesh3D[0][0][row][col] - MM.GW_build[name].parameters.param['bed_depress']
    cond = riv_cell[1] * riv_width_avg * MM.GW_build[name].parameters.param['Kv_riv'] / riv_bed_thickness
    simple_river += [[0, row, col, stage, cond, bed]]

riv = {}
riv[0] = simple_river
#MM.GW_build[name].boundaries.create_model_boundary_condition('Campaspe River', 'river', bc_static=True)
MM.GW_build[name].boundaries.assign_boundary_array('Campaspe River', riv)


print "************************************************************************"
print " Updating recharge boundary "

interp_rain = MM.GW_build[name].boundaries.bc['Rainfall']['bc_array'] * MM.GW_build[name].parameters.param['magic_rain']

rch = {}
rch[0] = interp_rain

#MM.GW_build[name].boundaries.create_model_boundary_condition('Rain_reduced', 'recharge', bc_static=True)
MM.GW_build[name].boundaries.assign_boundary_array('Rain_reduced', rch)

#print " Include irrigation in the recharge array"

#print MM.GW_build[name].observations.obs_group['head'].keys()
#print MM.GW_build[name].observations.obs_group['head']['mapped_observations']

print "************************************************************************"
print " Build and run MODFLOW model "

###########################################################################
###########################################################################
###########################################################################
## Currently using flopyInterface directly rather than running from the ModelManager ...
modflow_model = flopyInterface.ModflowModel(MM.GW_build[name], data_folder=data_folder)

modflow_model.executable = r"C:\Workspace\part0075\GIT_REPOS\CampaspeModel\MODFLOW-NWT_64.exe"

modflow_model.runMODFLOW()

#print " Return the stream-aquifer exchange for reaches as list "
#
#SWGWexchange = [1]
#
#print " Return the average depth to the GW table in areas as list "

#AvgDepthToGWTable = 1   
#DepthToGWTable = [1]

modflow_model.writeObservations()

#modflow_model.viewHeads()
#modflow_model.viewHeads2()    