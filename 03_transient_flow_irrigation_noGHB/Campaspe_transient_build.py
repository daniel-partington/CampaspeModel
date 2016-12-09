import os
import datetime

from osgeo import osr
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt

from HydroModelBuilder.GWModelBuilder import GWModelBuilder
from HydroModelBuilder.GISInterface.GDALInterface.GDALInterface import GDALInterface
from CampaspeModel.CustomScripts import processWeatherStations, getBoreData, get_GW_licence_info, processRiverGauges, readHydrogeologicalProperties

"""

"""
# Define basic model parameters:
Proj_CS = osr.SpatialReference()
Proj_CS.ImportFromEPSG(28355) # 28355 is the code for gda94 mga zone 55; http://spatialreference.org/ref/epsg/gda94-mga-zone-55/

Interface = GDALInterface()
Interface.projected_coordinate_system = Proj_CS 
Interface.pcs_EPSG = "EPSG:28355"

tr_model = GWModelBuilder(name="02_transient_flow_1966_2015", 
                          data_folder=r"C:\Workspace\part0075\MDB modelling\testbox\input_data\\",
                          out_data_folder=r"C:\Workspace\part0075\MDB modelling\testbox\02_transient_flow_1966_2015\\",
                          GISInterface=Interface,
                          model_type='Modflow',
                          mesh_type='structured')

# Cleanup
#tr_model.flush()

# Define the units for the project for consistency and to allow converions on input data
# tr_model.length = 'm'
# tr_model.time = 'd'

# Set the model boundary using a polygon shapefile:
print "************************************************************************"
print " Setting model boundary "

tr_model.set_model_boundary_from_polygon_shapefile("GW_model_area.shp", 
                                                      shapefile_path=tr_model.data_folder)

# Set data boundary for model data
print "************************************************************************"
print " Setting spatial data boundary "
tr_model.set_data_boundary_from_polygon_shapefile(tr_model.boundary_poly_file, 
                                                  buffer_dist=20000)

# Setup recharge:
# ... read in climate data using Custom_Scripts
weather_stations = ['Kyneton',  'Elmore', 'Rochester', 'Echuca']
print "************************************************************************"
print " Executing custom script: processWeatherStations "

rain_info_file = "rain_processed"
if os.path.exists(tr_model.out_data_folder + rain_info_file + '.h5'):
    long_term_historic_rainfall = tr_model.load_dataframe(tr_model.out_data_folder + rain_info_file + '.h5')
else:
    long_term_historic_rainfall = processWeatherStations.processWeatherStations(weather_stations, 
                                                                                path=r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Climate\\", 
                                                                                frequency='M')
    tr_model.save_dataframe(tr_model.out_data_folder + rain_info_file, long_term_historic_rainfall)

rain_gauges = tr_model.read_points_data(r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Climate\Rain_gauges.shp")
#points_dict = tr_model.getXYpairs(rain_gauges, feature_id='Name')

# $%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%
# INCLUDE NSW bores in this next part too for better head representation at the border, i.e. Murray River

# Read in bore data:
print "************************************************************************"
print " Executing custom script: getBoreData "

bore_levels_file = "bore_levels"
bore_info_file = "bore_info"
if os.path.exists(tr_model.out_data_folder + bore_levels_file + ".h5") & os.path.exists(tr_model.out_data_folder + bore_info_file + ".h5"):
    bore_data_levels = tr_model.load_dataframe(tr_model.out_data_folder + bore_levels_file + ".h5")
    bore_data_info = tr_model.load_dataframe(tr_model.out_data_folder + bore_info_file + ".h5")
else:
    bore_data_levels, bore_data_info = getBoreData.getBoreData()
    tr_model.save_dataframe(tr_model.out_data_folder + bore_levels_file, bore_data_levels)
    tr_model.save_dataframe(tr_model.out_data_folder + bore_info_file, bore_data_info)
# end if

# getBoreDepth ... assuming that midpoint of screen interval is representative location and assign to layer accordingly
bore_data_info['depth'] = (bore_data_info['TopElev'] + bore_data_info['BottomElev'])/2.0

bore_data_info["HydroCode"] = bore_data_info.index

# For steady state model, only use bore details containing average level, not 
#observation_bores = tr_model.read_points_data(r"C:\Workspace\part0075\MDB modelling\ngis_shp_VIC\ngis_shp_VIC\NGIS_Bores.shp")

print "************************************************************************"
print " Read in and filtering bore spatial data "

bores_shpfile = tr_model.read_points_data(r"C:\Workspace\part0075\MDB modelling\ngis_shp_VIC\ngis_shp_VIC\NGIS_Bores.shp")

bores_filtered_from_shpfile = tr_model.points_shapefile_obj2dataframe(bores_shpfile, feature_id="HydroCode")

# Get the intersection of bores_filtered_from_shpfile with bores_data_info

final_bores = pd.merge(bore_data_info, bores_filtered_from_shpfile, how='inner', on="HydroCode")

# Only consider bores whereby the measured values are above the bottom of screen
final_bores = final_bores[final_bores['mean level'] > final_bores['BottomElev']]

print 'Final number of bores within the data boundary that have level data and screen info: ', final_bores.shape[0]

#final_bores.plot(kind='scatter', x="Easting", y="Northing", c="mean level", cmap="Spectral") # , edgecolor='None'


# Load in the pumping wells data
filename = "Groundwater licence information for Dan Partington bc301115.xlsx"
path = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\GW\Bore data\\"    
out_path = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\GW\Bore data\\"
out_file = "pumping wells.shp"

print "************************************************************************"
print " Executing custom script: get_GW_licence_info "

pumping_data = get_GW_licence_info.get_GW_licence_info(filename, path=path, out_file=out_file, out_path=out_path)
pumps_points = tr_model.read_points_data(r"C:\Workspace\part0075\MDB modelling\Campaspe_data\GW\Bore data\pumping wells.shp")



print "************************************************************************"
print " Executing custom script: readHydrogeologicalProperties "

file_location = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\GW\Aquifer properties\Hydrogeologic_variables.xlsx"
HGU_props = readHydrogeologicalProperties.getHGUproperties(file_location)



#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************

print '########################################################################'
print '########################################################################'
print '## Model specific building '
print '########################################################################'
print '########################################################################'

print "************************************************************************"
print " Defining temporal aspects of the model"

start = datetime.date(1966, 01, 01)
#start = datetime.date(2014, 01, 01)
end = datetime.date(2015, 01, 01)
tr_model.model_time.set_temporal_components(steady_state=False, start_time=start, end_time=end, time_step='M')

# Define the grid width and grid height for the model mesh which is stored as a multipolygon shapefile GDAL object
print "************************************************************************"
print " Defining structured mesh"
tr_model.define_structured_mesh(5000, 5000) #10000,10000)

# Read in hydrostratigraphic raster info for layer elevations:
#hu_raster_path = r"C:\Workspace\part0075\MDB modelling\Campaspe_model\GIS\GIS_preprocessed\Hydrogeological_Unit_Layers\\"
hu_raster_path = r"C:\Workspace\part0075\MDB modelling\VAF_v2.0_ESRI_GRID\ESRI_GRID\\"
#hu_raster_path = r"C:\Workspace\part0075\MDB modelling\VAF_v2.0_ESRI_GRID\ESRI_GRID_raw\ESRI_GRID\\"

# Build basement file ... only need to do this once as it is time consuming so commented out for future runs
#tr_model.create_basement_bottom(hu_raster_path, "sur_1t", "bse_1t", "bse_2b", hu_raster_path)

hu_raster_files = ["qa_1t", "qa_2b", "utb_1t", "utb_2b", "utqa_1t", "utqa_2b", "utam_1t", "utam_2b", "utaf_1t", "utaf_2b", "lta_1t", "lta_2b", "bse_1t", "bse_2b.tif"]
#hu_raster_files = ["qa_1t", "qa_2b", "utb_1t", "utb_2b", "utqa_1t", "utqa_2b", "utam_1t", "utam_2b", "utaf_1t", "utaf_2b", "lta_1t", "lta_2b", "bse_1t", "bse_2b.tif"]
tr_model.read_rasters(hu_raster_files, path=hu_raster_path)
hu_raster_files_reproj = [x+"_reproj.bil" for x in hu_raster_files]

# Map HGU's to grid
print "************************************************************************"
print " Mapping rasters to grid "

hu_gridded_rasters = tr_model.map_rasters_to_grid(hu_raster_files, hu_raster_path)

# Build 3D grid
model_grid_raster_files = [x+"_model_grid.bil" for x in hu_raster_files]

# First two arguments of next function are arbitrary and not used ... need to rework module
print "************************************************************************"
print " Building 3D mesh "
tr_model.build_3D_mesh_from_rasters(model_grid_raster_files, tr_model.out_data_folder_grid, 1.0, 1000.0)

print "************************************************************************"
print " Assign properties to mesh based on zonal information"

# create list of HGU's from hu_raster_files
HGU = list(set([x.split('_')[0] for x in hu_raster_files]))

# NOTE *** utam is mapping to Shepparton Sands but it represents the 
# Loxton-Parilla Sand ... the HGU database needs updating to include this.
HGU_map = {'bse':'Bedrock', 'utb':'Newer Volcanics Basalts', 
           'utaf':'Calivil', 'lta':'Renmark', 
           'qa':'Coonambidgal Sands', 'utqa':'Shepparton Sands',
           'utam':'Shepparton Sands'}

for unit in HGU:
    tr_model.parameters.create_model_parameter('Kh_' + unit, value=HGU_props['Kh mean'][HGU_map[unit]])
    tr_model.parameters.parameter_options('Kh_' + unit, 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=HGU_props['Kh mean'][HGU_map[unit]] / 10., 
                                          PARUBND=HGU_props['Kh mean'][HGU_map[unit]] * 10., 
                                          PARGP='cond', 
                                          SCALE=1, 
                                          OFFSET=0)
    tr_model.parameters.create_model_parameter('Kv_' + unit, value=HGU_props['Kz mean'][HGU_map[unit]])
    tr_model.parameters.parameter_options('Kv_' + unit, 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=HGU_props['Kz mean'][HGU_map[unit]] / 10., 
                                          PARUBND=HGU_props['Kz mean'][HGU_map[unit]] * 10., 
                                          PARGP='cond', 
                                          SCALE=1, 
                                          OFFSET=0)
    tr_model.parameters.create_model_parameter('Sy_' + unit, value=HGU_props['Sy mean'][HGU_map[unit]])
    tr_model.parameters.parameter_options('Sy_' + unit, 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=0., 
                                          PARUBND=0.8, 
                                          PARGP='spec_yield', 
                                          SCALE=1, 
                                          OFFSET=0)
    tr_model.parameters.create_model_parameter('SS_' + unit, value=HGU_props['SS mean'][HGU_map[unit]])
    tr_model.parameters.parameter_options('SS_' + unit, 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=HGU_props['SS mean'][HGU_map[unit]] / 10., 
                                          PARUBND=HGU_props['SS mean'][HGU_map[unit]] * 10., 
                                          PARGP='spec_stor', 
                                          SCALE=1, 
                                          OFFSET=0)

# This needs to be automatically generated from with the map_raster2mesh routine ...
zone_map = {1:'qa', 2:'utb', 3:'utqa', 4:'utam', 5:'utaf', 6:'lta', 7:'bse'}

Kh = tr_model.model_mesh3D[1].astype(float)
Kv = tr_model.model_mesh3D[1].astype(float)
Sy = tr_model.model_mesh3D[1].astype(float)
SS = tr_model.model_mesh3D[1].astype(float)
for key in zone_map.keys():
    Kh[Kh == key] = tr_model.parameters.param['Kh_' + zone_map[key]]['PARVAL1']
    Kv[Kv == key] = tr_model.parameters.param['Kv_' + zone_map[key]]['PARVAL1']
    Sy[Sy == key] = tr_model.parameters.param['Sy_' + zone_map[key]]['PARVAL1']
    SS[SS == key] = tr_model.parameters.param['SS_' + zone_map[key]]['PARVAL1']

tr_model.properties.assign_model_properties('Kh', Kh)
tr_model.properties.assign_model_properties('Kv', Kv)
tr_model.properties.assign_model_properties('Sy', Sy)
tr_model.properties.assign_model_properties('SS', SS)

print "************************************************************************"
print " Interpolating rainfall data to grid "

long_term_historic_rainfall = long_term_historic_rainfall.ix[start:end]

interp_rain = {}
for step, month in enumerate(long_term_historic_rainfall.iterrows()):
    print step    
    interp_rain[step] = tr_model.interpolate_points2mesh(rain_gauges, month[1], feature_id='Name')
    # Adjust rainfall to m from mm 
    interp_rain[step] = interp_rain[step]/1000.0

tr_model.boundaries.create_model_boundary_condition('Rainfall', 'rainfall', bc_static=True)
tr_model.boundaries.assign_boundary_array('Rainfall', interp_rain)

# Adjust rainfall to recharge using rainfall reduction
for i in [1,2,3,7]:
    tr_model.parameters.create_model_parameter('rch_red_'+zone_map[i], value=0.1)
    tr_model.parameters.parameter_options('rch_red_'+zone_map[i], 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=0., 
                                          PARUBND=0.9, 
                                          PARGP='rech_mult', 
                                          SCALE=1, 
                                          OFFSET=0)
    for key in interp_rain.keys():
        interp_rain[key][tr_model.model_mesh3D[1][0]==i] = interp_rain[key][tr_model.model_mesh3D[1][0]==i] * tr_model.parameters.param['rch_red_'+zone_map[i]]['PARVAL1']

rch = {}
for key in interp_rain.keys():
    rch[key] = interp_rain[key]


print "************************************************************************"
print " Creating recharge boundary "

tr_model.boundaries.create_model_boundary_condition('Rain_reduced', 'recharge', bc_static=False)
tr_model.boundaries.assign_boundary_array('Rain_reduced', rch)

print "************************************************************************"
print " Mapping bores to grid "

tr_model.map_points_to_grid(bores_shpfile, feature_id = 'HydroCode')

print tr_model.points_mapped.keys()

bores_more_filter = []
for bores in tr_model.points_mapped["NGIS_Bores_clipped.shp"]:
    row = bores[0][0]
    col = bores[0][1]
    for bore in bores[1]: 
        try:
            bore_depth = bore_data_info.loc[bore, 'depth'] #[bore_data_info["HydroCode"] == HydroCode]['depth']        
        except:
            continue
        if bore_depth > tr_model.model_mesh3D[0][0][row][col]:
            #print 'Bore can't be above surface!!!        
            continue
        if bore_depth <= tr_model.model_mesh3D[0][-2][row][col]:
            #print 'Ignoring bores in bedrock!!!        
            continue
        bores_more_filter += [bore]        

print 'Final bores within aquifers: ', len(bores_more_filter)

final_bores = final_bores[final_bores["HydroCode"].isin(bores_more_filter)]

bore_points = [[final_bores.loc[x, "Easting"], final_bores.loc[x, "Northing"]] for x in final_bores.index]

bore_points3D = final_bores[["HydroCode", "Easting", "Northing", "depth"]] # [[final_bores.loc[x, "Easting"], final_bores.loc[x, "Northing"], final_bores.loc[x, "depth"]] for x in final_bores.index]
bore_points3D = bore_points3D.set_index("HydroCode")

bores_obs_time_series = bore_data_levels[bore_data_levels["HydroCode"].isin(final_bores["HydroCode"])]

# Modify into standard format for the GWModelBuilder class
bores_obs_time_series = bores_obs_time_series.rename(columns={'HydroCode':'name', 'bore_date':'datetime', 'result':'value'})

bores_obs_time_series['datetime'] = pd.to_datetime(bores_obs_time_series['datetime'])

tr_model.observations.set_as_observations('head', bores_obs_time_series, bore_points3D, domain='porous', obs_type='head', units='mAHD')

tr_model.map_obs_loc2mesh3D(method='nearest')
tr_model.map_obs2model_times()


bores_in_layers = tr_model.map_points_to_raster_layers(bore_points, final_bores["depth"].tolist(), hu_raster_files_reproj)

# Map bores to layers to create initial head maps for different hydrogeological units
#interp_heads = {}
#
#for i in range(len(hu_raster_files_reproj)/2):
#    bores_layer = np.array(bore_points)[np.array(bores_in_layers[i])]
#    print 'Creating head map for: ', hu_raster_files[2*i]
#    if bores_layer.shape[0] < 4: 
#        #interp_heads[hu_raster_files[2*i]] = (tr_model.model_mesh3D[0][i]+tr_model.model_mesh3D[0][i+1])/2
#        interp_heads[hu_raster_files[2*i]] = np.full(tr_model.model_mesh3D[1].shape[1:], np.NaN)
#    else:
#        bores_head_layer = np.array(final_bores["mean level"].tolist())[np.array(bores_in_layers[i])]
#        unique_bores = np.unique(bores_layer) 
#    
#        b = np.ascontiguousarray(bores_layer).view(np.dtype((np.void, bores_layer.dtype.itemsize * bores_layer.shape[1])))
#        _, idx = np.unique(b, return_index=True)
#    
#        unique_bores = bores_layer[idx]    
#    
#        interp_heads[hu_raster_files[2*i]] = tr_model.interpolate_points2mesh(bores_layer, bores_head_layer, use='griddata', method='linear')
        
##for key in interp_heads:
#    #bores_layer_df = pd.DataFrame()
#    #bores_layer_df["Easting"] = [x[0] for x in bores_layer] 
#    #bores_layer_df["Northing"] = [x[1] for x in bores_layer]
#    #bores_layer_df["mean level"] = bores_head_layer
#    #(XI, YI) = tr_model.model_mesh_centroids
#    #plt.figure()
#    #z_min = np.min(interp_heads[key])
#    #z_max = np.max(interp_heads[key])
#    #plt.pcolor(XI, YI, interp_heads[key], vmin=z_min, vmax=z_max)
#    #plt.scatter([x[0] for x in bores_layer], [x[1] for x in bores_layer], 50, bores_head_layer, vmin=z_min, vmax=z_max, cmap="jet")
#    #plt.colorbar()
#
#    #bores_layer_df.plot(kind='scatter', x="Easting", y="Northing", c="mean level", cmap="Spectral") # , edgecolor='None'
#    #plt.scatter(x=[x[0] for x in bores_layer], y=[x[1] for x in bores_layer], c=bores_head_layer)
#
# Initalise model with head from elevations
initial_heads_tr = np.full(tr_model.model_mesh3D[1].shape, 0.)

for i in range(len(hu_raster_files_reproj)/2):
    initial_heads_tr[i] = (tr_model.model_mesh3D[0][i]+tr_model.model_mesh3D[0][i+1])/2

tr_model.initial_conditions.set_as_initial_condition("Head", initial_heads_tr)#interp_heads[hu_raster_files[0]])
#
#initial_heads_tr = np.full(tr_model.model_mesh3D[1].shape, 0.)
#
#for i in range(len(hu_raster_files_reproj)/2):
#    initial_heads_tr[i] = (interp_heads[hu_raster_files[2*i]])
#
#tr_model.initial_conditions.set_as_initial_condition("OldHead", initial_heads_tr)#interp_heads[hu_raster_files[0]])
#

w

## Map river polyline feature to grid including length of river in cell
print "************************************************************************"
print " Mapping Campaspe river to grid"

river_poly = tr_model.read_polyline("Campaspe_Riv.shp", path=r"C:\Workspace\part0075\MDB modelling\Campaspe_model\GIS\GIS_preprocessed\Surface_Water\Streams\\") 
tr_model.map_polyline_to_grid(river_poly)
tr_model.parameters.create_model_parameter('bed_depress', value=0.01)
tr_model.parameters.parameter_options('bed_depress', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.001, 
                                      PARUBND=0.1, 
                                      PARGP='river', 
                                      SCALE=1, 
                                      OFFSET=0)
tr_model.parameters.create_model_parameter('Kv_riv', value=5E-3)
tr_model.parameters.parameter_options('Kv_riv', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=1E-8, 
                                      PARUBND=20, 
                                      PARGP='river', 
                                      SCALE=1, 
                                      OFFSET=0)

simple_river = []
riv_width_avg = 10.0 #m
riv_bed_thickness = 0.10 #m
for riv_cell in tr_model.polyline_mapped['Campaspe_Riv_model.shp']:
    row = riv_cell[0][0]
    col = riv_cell[0][1]
    if tr_model.model_mesh3D[1][0][row][col] == -1:
        continue
    #print tr_model.model_mesh3D
    stage = tr_model.model_mesh3D[0][0][row][col]
    bed = tr_model.model_mesh3D[0][0][row][col] - tr_model.parameters.param['bed_depress']['PARVAL1']
    cond = riv_cell[1] * riv_width_avg * tr_model.parameters.param['Kv_riv']['PARVAL1'] / riv_bed_thickness
    simple_river += [[0, row, col, stage, cond, bed]]

riv = {}
riv[0] = simple_river

#print tr_model.polyline_mapped
print "************************************************************************"
print " Creating Campaspe river boundary"

tr_model.boundaries.create_model_boundary_condition('Campaspe River', 'river', bc_static=True)
tr_model.boundaries.assign_boundary_array('Campaspe River', riv)

print "************************************************************************"
print " Mapping Murray River to grid"

river_poly = tr_model.read_polyline("River_Murray.shp", path=r"C:\Workspace\part0075\MDB modelling\Campaspe_model\GIS\GIS_preprocessed\Surface_Water\Streams\\") 
tr_model.map_polyline_to_grid(river_poly)

tr_model.parameters.create_model_parameter('RMstage', value=0.01)
tr_model.parameters.parameter_options('RMstage', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.001, 
                                      PARUBND=0.1, 
                                      PARGP='spec_stor', 
                                      SCALE=1, 
                                      OFFSET=0)
tr_model.parameters.create_model_parameter('Kv_RM', value=5E-3)
tr_model.parameters.parameter_options('Kv_RM', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=1E-8, 
                                      PARUBND=20, 
                                      PARGP='spec_stor', 
                                      SCALE=1, 
                                      OFFSET=0)


simple_river = []
riv_width_avg = 10.0 #m
riv_bed_thickness = 0.10 #m
for riv_cell in tr_model.polyline_mapped['River_Murray_model.shp']:
    row = riv_cell[0][0]
    col = riv_cell[0][1]
    if tr_model.model_mesh3D[1][0][row][col] == -1:
        continue
    #print tr_model.model_mesh3D
    stage = tr_model.model_mesh3D[0][0][row][col]
    bed = tr_model.model_mesh3D[0][0][row][col] - tr_model.parameters.param['RMstage']['PARVAL1']
    cond = riv_cell[1] * riv_width_avg * tr_model.parameters.param['Kv_RM']['PARVAL1'] / riv_bed_thickness
    simple_river += [[0, row, col, stage, cond, bed]]

riv = {}
riv[0] = simple_river

#print tr_model.polyline_mapped
print "************************************************************************"
print " Creating Murray River boundary"

tr_model.boundaries.create_model_boundary_condition('Murray River', 'river', bc_static=True)
tr_model.boundaries.assign_boundary_array('Murray River', riv)

print "************************************************************************"
print " Setting up Murray River GHB boundary"

tr_model.parameters.create_model_parameter('MGHB_stage', value=0.01)
tr_model.parameters.parameter_options('MGHB_stage', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.001, 
                                      PARUBND=0.1, 
                                      PARGP='ghb', 
                                      SCALE=1, 
                                      OFFSET=0)
tr_model.parameters.create_model_parameter('MGHBcond', value=5E-3)
tr_model.parameters.parameter_options('MGHBcond', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=1E-8, 
                                      PARUBND=20, 
                                      PARGP='ghb', 
                                      SCALE=1, 
                                      OFFSET=0)


MurrayGHB = []
MGHBcond = 5E-3 #m/day
for MurrayGHB_cell in tr_model.polyline_mapped['River_Murray_model.shp']:
    row = MurrayGHB_cell[0][0]
    col = MurrayGHB_cell[0][1]
    #print tr_model.model_mesh3D
    for lay in range(tr_model.model_mesh3D[1].shape[0]):    
        if tr_model.model_mesh3D[1][0][row][col] == -1:
            continue
        MurrayGHBstage = tr_model.model_mesh3D[0][0][row][col] + tr_model.parameters.param['MGHB_stage']['PARVAL1']
        dx = tr_model.model_mesh3D[0][0][0][1] - tr_model.model_mesh3D[0][0][0][0]
        MGHBconductance = dx * tr_model.parameters.param['MGHBcond']['PARVAL1']

ghb = {}
ghb[0] = MurrayGHB

print "************************************************************************"
print " Mapping Western GW boundary to grid"

WGWbound_poly = tr_model.read_polyline("western_head.shp", path=r"C:\Workspace\part0075\MDB modelling\testbox\input_data\\") 
tr_model.map_polyline_to_grid(WGWbound_poly)

print "************************************************************************"
print " Setting up Western GHB boundary"

tr_model.parameters.create_model_parameter('WGHB_stage', value=0.01)
tr_model.parameters.parameter_options('WGHB_stage', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.001, 
                                      PARUBND=0.1, 
                                      PARGP='ghb', 
                                      SCALE=1, 
                                      OFFSET=0)
tr_model.parameters.create_model_parameter('WGHBcond', value=5E-3)
tr_model.parameters.parameter_options('WGHBcond', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=1E-8, 
                                      PARUBND=20, 
                                      PARGP='ghb', 
                                      SCALE=1, 
                                      OFFSET=0)


WestGHB = []
for WestGHB_cell in tr_model.polyline_mapped['western_head_model.shp']:
    row = WestGHB_cell[0][0]
    col = WestGHB_cell[0][1]
    #print tr_model.model_mesh3D
    for lay in range(tr_model.model_mesh3D[1].shape[0]):    
        if tr_model.model_mesh3D[1][lay][row][col] == -1:
            continue
        WestGHBstage = tr_model.model_mesh3D[0][0][row][col] + tr_model.parameters.param['WGHB_stage']['PARVAL1']
        dx = tr_model.model_mesh3D[0][0][0][1] - tr_model.model_mesh3D[0][0][0][0]
        WGHBconductance = dx * tr_model.parameters.param['WGHBcond']['PARVAL1']
        WestGHB += [[lay, row, col, WestGHBstage, WGHBconductance]]


ghb[0] += WestGHB


print "************************************************************************"
print " Mapping Eastern GW boundary to grid"

EGWbound_poly = tr_model.read_polyline("eastern_head.shp", path=r"C:\Workspace\part0075\MDB modelling\testbox\input_data\\") 
tr_model.map_polyline_to_grid(EGWbound_poly)

print "************************************************************************"
print " Setting up Western GHB boundary"

tr_model.parameters.create_model_parameter('EGHB_stage', value=0.01)
tr_model.parameters.parameter_options('EGHB_stage', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.001, 
                                      PARUBND=0.1, 
                                      PARGP='ghb', 
                                      SCALE=1, 
                                      OFFSET=0)
tr_model.parameters.create_model_parameter('EGHBcond', value=5E-3)
tr_model.parameters.parameter_options('EGHBcond', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=1E-8, 
                                      PARUBND=20, 
                                      PARGP='ghb', 
                                      SCALE=1, 
                                      OFFSET=0)


EastGHB = []
for EastGHB_cell in tr_model.polyline_mapped['eastern_head_model.shp']:
    row = EastGHB_cell[0][0]
    col = EastGHB_cell[0][1]
    #print SS_model.model_mesh3D
    for lay in range(tr_model.model_mesh3D[1].shape[0]):    
        if tr_model.model_mesh3D[1][lay][row][col] == -1:
            continue
        #EastGHBstage = (SS_model.model_mesh3D[0][lay+1][row][col] + SS_model.model_mesh3D[0][lay][row][col])/2. + SS_model.parameters.param['EGHB_stage']['PARVAL1']
        EastGHBstage = tr_model.model_mesh3D[0][0][row][col] + tr_model.parameters.param['EGHB_stage']['PARVAL1']
        dx = tr_model.model_mesh3D[0][0][0][1] - tr_model.model_mesh3D[0][0][0][0]
        EGHBconductance = dx * tr_model.parameters.param['EGHBcond']['PARVAL1']
        EastGHB += [[lay, row, col, EastGHBstage, EGHBconductance]]

ghb[0] += EastGHB

print "************************************************************************"
print " Creating GHB boundary"

tr_model.boundaries.create_model_boundary_condition('GHB', 'general head', bc_static=True)
tr_model.boundaries.assign_boundary_array('GHB', ghb)

print "************************************************************************"
print " Mapping Drains to grid"

drain_poly = tr_model.read_polyline("Drain_Clip.shp", path=r"C:\Workspace\part0075\MDB modelling\testbox\input_data\SW\\") 
tr_model.map_polyline_to_grid(drain_poly)

tr_model.parameters.create_model_parameter('drain_drop', value=0.01)
tr_model.parameters.parameter_options('drain_drop', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.001, 
                                      PARUBND=0.1, 
                                      PARGP='drain', 
                                      SCALE=1, 
                                      OFFSET=0)
tr_model.parameters.create_model_parameter('Kv_drain', value=5E-3)
tr_model.parameters.parameter_options('Kv_drain', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=1E-8, 
                                      PARUBND=20, 
                                      PARGP='drain', 
                                      SCALE=1, 
                                      OFFSET=0)

simple_drain = []
Kv_drain = 5E-3 #m/day
drain_width_avg = 3.0 #m
drain_bed_thickness = 0.10 #m
for drain_cell in tr_model.polyline_mapped['Drain_Clip_model.shp']:
    row = drain_cell[0][0]
    col = drain_cell[0][1]
    if tr_model.model_mesh3D[1][0][row][col] == -1:
        continue
    #print tr_model.model_mesh3D
    drain_bed = tr_model.model_mesh3D[0][0][row][col] - tr_model.parameters.param['drain_drop']['PARVAL1']
    drain_cond = drain_cell[1] * drain_width_avg * tr_model.parameters.param['Kv_drain']['PARVAL1'] / drain_bed_thickness
    simple_drain += [[0, row, col, drain_bed, drain_cond]]

drain = {}
drain[0] = simple_drain

print "************************************************************************"
print " Creating Drains boundary"

tr_model.boundaries.create_model_boundary_condition('Drain', 'drain', bc_static=True)
tr_model.boundaries.assign_boundary_array('Drain', drain)


print "************************************************************************"
print " Mapping Channels to grid"

channel_poly = tr_model.read_polyline("Channel_Clip.shp", path=r"C:\Workspace\part0075\MDB modelling\testbox\input_data\SW\\") 
tr_model.map_polyline_to_grid(channel_poly)

tr_model.parameters.create_model_parameter('chan_drop', value=0.01)
tr_model.parameters.parameter_options('chan_drop', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.001, 
                                      PARUBND=0.1, 
                                      PARGP='channel', 
                                      SCALE=1, 
                                      OFFSET=0)
tr_model.parameters.create_model_parameter('Kv_chan', value=5E-3)
tr_model.parameters.parameter_options('Kv_chan', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=1E-8, 
                                      PARUBND=20, 
                                      PARGP='channel', 
                                      SCALE=1, 
                                      OFFSET=0)

simple_channel = []
channel_width_avg = 10.0 #m
channel_bed_thickness = 0.10 #m
for channel_cell in tr_model.polyline_mapped['Channel_Clip_model.shp']:
    row = channel_cell[0][0]
    col = channel_cell[0][1]
    if tr_model.model_mesh3D[1][0][row][col] == -1:
        continue
    #print tr_model.model_mesh3D
    channel_stage = tr_model.model_mesh3D[0][0][row][col]
    channel_bed = tr_model.model_mesh3D[0][0][row][col] - tr_model.parameters.param['chan_drop']['PARVAL1']
    channel_cond = channel_cell[1] * channel_width_avg * tr_model.parameters.param['Kv_chan']['PARVAL1'] / channel_bed_thickness
    simple_channel += [[0, row, col, channel_stage, channel_cond, channel_bed]]

channel = {}
channel[0] = simple_channel

print "************************************************************************"
print " Creating Channel boundary"

tr_model.boundaries.create_model_boundary_condition('Channel', 'channel', bc_static=True)
tr_model.boundaries.assign_boundary_array('Channel', channel)

print "************************************************************************"
print " Collate observations"

tr_model.observations.collate_observations()


print "************************************************************************"
print " Package up groundwater model builder object"

tr_model.package_model()
