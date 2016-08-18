import os

from osgeo import osr
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt

from HydroModelBuilder.HydroModelBuilder.GWModelBuilder import GWModelBuilder
from HydroModelBuilder.HydroModelBuilder.GISInterface.GDALInterface.GDALInterface import GDALInterface
from HydroModelBuilder.HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface
from CustomScripts import processWeatherStations, getBoreData, get_GW_licence_info, processRiverGauges, readHydrogeologicalProperties

"""
STEPS IN THE GROUNDWATER MODEL BUILDING PROCESS:

1. Define the projected coordinate system (PCS) and a GIS interface object to
   pass to the model builder. This will be used to do any necessary 
   geotransforms to put all the data in the same PCS. 

2. Define the attributes that make up the model:
    a. Model type
    b. Mesh type
    c. data folder for inputs
    d. output data folder for created files during building
    e. output data types (i.e. as binary or ascii files)
    f. pass a GISInterface object to use (from 1.), e.g. GDAL or ArcPy
    
3. Define the mesh either by simple corner coordinates or using a shapefile
   polygon.

4. Assuming the mesh doesn't change across layers, define the model array,
   i.e. the elevation of nodes/cells across each of the layers

5. If using geospatial data, then use the GISInterface to do any operations
   on the data and map it to the grid, e.g. define streams from polyline, 
   define pumping wells and observation bores from points.     

"""

# Define basic model parameters:
Proj_CS = osr.SpatialReference()
Proj_CS.ImportFromEPSG(28355) # 28355 is the code for gda94 mga zone 55; http://spatialreference.org/ref/epsg/gda94-mga-zone-55/

Interface = GDALInterface()
Interface.projected_coordinate_system = Proj_CS 
Interface.pcs_EPSG = "EPSG:28355"

test_model = GWModelBuilder(name="Campaspe", 
                            data_folder=r"C:\Workspace\part0075\MDB modelling\testbox\input_data\\",
                            out_data_folder=r"C:\Workspace\part0075\MDB modelling\testbox\model_files\\",
                            GISInterface=Interface,
                            model_type='Modflow',
                            mesh_type='structured')

# Cleanup
#test_model.flush()

# Define the units for the project for consistency and to allow converions on input data
test_model.length = 'm'
test_model.time = 'd'

# Set the model boundary using a polygon shapefile:
print "************************************************************************"
print " Setting model boundary "

test_model.set_model_boundary_from_polygon_shapefile("GW_model_area.shp", 
                                                      shapefile_path=test_model.data_folder)

# Set data boundary for model data
print "************************************************************************"
print " Setting spatial data boundary "
test_model.set_data_boundary_from_polygon_shapefile(test_model.boundary_poly_file, 
                                                    buffer_dist=20000)

# Setup recharge:
# ... read in climate data using Custom_Scripts
weather_stations = ['Kyneton',  'Elmore', 'Rochester', 'Echuca']
print "************************************************************************"
print " Executing custom script: processWeatherStations "

rain_info_file = "rain_processed.h5"
if os.path.exists(test_model.out_data_folder + rain_info_file):
    long_term_historic_rainfall = test_model.load_dataframe(test_model.out_data_folder + rain_info_file)
else:
    long_term_historic_rainfall = processWeatherStations.processWeatherStations(weather_stations, path=r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Climate\\")
    test_model.save_dataframe(test_model.out_data_folder + rain_info_file, long_term_historic_rainfall)

rain_gauges = test_model.read_points_data(r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Climate\Rain_gauges.shp")
#points_dict = test_model.getXYpairs(rain_gauges, feature_id='Name')

# $%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%
# INCLUDE NSW bores in this next part too for better head representation at the border, i.e. Murray River

# Read in bore data:
print "************************************************************************"
print " Executing custom script: getBoreData "

bore_levels_file = "bore_levels.h5"
bore_info_file = "bore_info.h5"
if os.path.exists(test_model.out_data_folder + bore_levels_file) & os.path.exists(test_model.out_data_folder + bore_info_file):
    bore_data_levels = test_model.load_dataframe(test_model.out_data_folder + bore_levels_file)
    bore_data_info = test_model.load_dataframe(test_model.out_data_folder + bore_info_file)
else:
    bore_data_levels, bore_data_info = getBoreData.getBoreData()
    test_model.save_dataframe(test_model.out_data_folder + bore_levels_file, bore_data_levels)
    test_model.save_dataframe(test_model.out_data_folder + bore_info_file, bore_data_info)
# end if

# getBoreDepth ... assuming that midpoint of screen interval is representative location and assign to layer accordingly
bore_data_info['depth'] = (bore_data_info['TopElev'] + bore_data_info['BottomElev'])/2.0

bore_data_info["HydroCode"] = bore_data_info.index

# For steady state model, only use bore details containing average level, not 
#observation_bores = test_model.read_points_data(r"C:\Workspace\part0075\MDB modelling\ngis_shp_VIC\ngis_shp_VIC\NGIS_Bores.shp")

print "************************************************************************"
print " Read in and filtering bore spatial data "

bores_shpfile = test_model.read_points_data(r"C:\Workspace\part0075\MDB modelling\ngis_shp_VIC\ngis_shp_VIC\NGIS_Bores.shp")

bores_filtered_from_shpfile = test_model.points_shapefile_obj2dataframe(bores_shpfile, feature_id="HydroCode")

# Get the intersection of bores_filtered_from_shpfile with bores_data_info

final_bores = pd.merge(bore_data_info, bores_filtered_from_shpfile, how='inner', on="HydroCode")

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
pumps_points = test_model.read_points_data(r"C:\Workspace\part0075\MDB modelling\Campaspe_data\GW\Bore data\pumping wells.shp")

print "************************************************************************"
print " Executing custom script: readHydrogeologicalProperties "

file_location = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\GW\Aquifer properties\Hydrogeologic_variables.xlsx"
HGU_props = readHydrogeologicalProperties.getHGUproperties(file_location)


print '########################################################################'
print '########################################################################'
print '## Mesh specific model building '
print '########################################################################'
print '########################################################################'

# Define the grid width and grid height for the model mesh which is stored as a multipolygon shapefile GDAL object
print "************************************************************************"
print " Defining structured mesh"
test_model.define_structured_mesh(10000, 10000) #10000,10000)

# Read in hydrostratigraphic raster info for layer elevations:
#hu_raster_path = r"C:\Workspace\part0075\MDB modelling\Campaspe_model\GIS\GIS_preprocessed\Hydrogeological_Unit_Layers\\"
hu_raster_path = r"C:\Workspace\part0075\MDB modelling\VAF_v2.0_ESRI_GRID\ESRI_GRID\\"
#hu_raster_path = r"C:\Workspace\part0075\MDB modelling\VAF_v2.0_ESRI_GRID\ESRI_GRID_raw\ESRI_GRID\\"

# Build basement file ... only need to do this once as it is time consuming so commented out for future runs
#test_model.create_basement_bottom(hu_raster_path, "sur_1t", "bse_1t", "bse_2b", hu_raster_path)

hu_raster_files = ["qa_1t", "qa_2b", "utb_1t", "utb_2b", "utqa_1t", "utqa_2b", "utam_1t", "utam_2b", "utaf_1t", "utaf_2b", "lta_1t", "lta_2b", "bse_1t", "bse_2b.tif"]
#hu_raster_files = ["qa_1t", "qa_2b", "utb_1t", "utb_2b", "utqa_1t", "utqa_2b", "utam_1t", "utam_2b", "utaf_1t", "utaf_2b", "lta_1t", "lta_2b", "bse_1t", "bse_2b.tif"]
test_model.read_rasters(hu_raster_files, path=hu_raster_path)
hu_raster_files_reproj = [x+"_reproj.bil" for x in hu_raster_files]


# Map HGU's to grid
print "************************************************************************"
print " Mapping rasters to grid "

hu_gridded_rasters = test_model.map_rasters_to_grid(hu_raster_files, hu_raster_path)

# Build 3D grid
model_grid_raster_files = [x+"_model_grid.bil" for x in hu_raster_files]

# First two arguments of next function are arbitrary and not used ... need to rework module
print "************************************************************************"
print " Building 3D mesh "
test_model.build_3D_mesh_from_rasters(model_grid_raster_files, test_model.out_data_folder_grid, 1.0, 1000.0)

print "************************************************************************"
print " Assign properties to mesh based on zonal information"

# create list of HGU's from hu_raster_files
HGU = [x.split('_')[0] for x in hu_raster_files]

# NOTE *** utam is mapping to Shepparton Sands but it represents the 
# Loxton-Parilla Sand ... the HGU database needs updating to include this.
HGU_map = {'bse':'Bedrock', 'utb':'Newer Volcanics Basalts', 
           'utaf':'Calivil', 'lta':'Renmark', 
           'qa':'Coonambidgal Sands', 'utqa':'Shepparton Sands',
           'utam':'Shepparton Sands'}

for unit in HGU:
    test_model.parameters.create_model_parameter('Kh_' + unit, value=HGU_props['Kh mean'][HGU_map[unit]])
    test_model.parameters.parameter_options('Kh_' + unit, 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=HGU_props['Kh mean'][HGU_map[unit]] / 10., 
                                          PARUBND=HGU_props['Kh mean'][HGU_map[unit]] * 10., 
                                          PARGP='cond', 
                                          SCALE=1, 
                                          OFFSET=0)
    test_model.parameters.create_model_parameter('Kv_' + unit, value=HGU_props['Kz mean'][HGU_map[unit]])
    test_model.parameters.parameter_options('Kv_' + unit, 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=HGU_props['Kz mean'][HGU_map[unit]] / 10., 
                                          PARUBND=HGU_props['Kz mean'][HGU_map[unit]] * 10., 
                                          PARGP='cond', 
                                          SCALE=1, 
                                          OFFSET=0)
    test_model.parameters.create_model_parameter('Sy_' + unit, value=HGU_props['Sy mean'][HGU_map[unit]])
    test_model.parameters.parameter_options('Sy_' + unit, 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=0., 
                                          PARUBND=0.8, 
                                          PARGP='spec_yield', 
                                          SCALE=1, 
                                          OFFSET=0)
    test_model.parameters.create_model_parameter('SS_' + unit, value=HGU_props['SS mean'][HGU_map[unit]])
    test_model.parameters.parameter_options('SS_' + unit, 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=HGU_props['SS mean'][HGU_map[unit]] / 10., 
                                          PARUBND=HGU_props['SS mean'][HGU_map[unit]] * 10., 
                                          PARGP='spec_stor', 
                                          SCALE=1, 
                                          OFFSET=0)

# This needs to be automatically generated from with the map_raster2mesh routine ...
zone_map = {1:'qa', 2:'utb', 3:'utqa', 4:'utam', 5:'utaf', 6:'lta', 7:'bse'}

Kh = test_model.model_mesh3D[1].astype(float)
Kv = test_model.model_mesh3D[1].astype(float)
Sy = test_model.model_mesh3D[1].astype(float)
SS = test_model.model_mesh3D[1].astype(float)
for key in zone_map.keys():
    Kh[Kh == key] = test_model.parameters.param['Kh_' + zone_map[key]]['PARVAL1']
    Kv[Kv == key] = test_model.parameters.param['Kv_' + zone_map[key]]['PARVAL1']
    Sy[Sy == key] = test_model.parameters.param['Sy_' + zone_map[key]]['PARVAL1']
    SS[SS == key] = test_model.parameters.param['SS_' + zone_map[key]]['PARVAL1']

test_model.properties.assign_model_properties('Kh', Kh)
test_model.properties.assign_model_properties('Kv', Kv)
test_model.properties.assign_model_properties('Sy', Sy)
test_model.properties.assign_model_properties('SS', SS)

print "************************************************************************"
print " Interpolating rainfall data to grid "

interp_rain = test_model.interpolate_points2mesh(rain_gauges, long_term_historic_rainfall, feature_id='Name')
# Adjust rainfall to m from mm and from year to days
interp_rain = interp_rain/1000.0/365.0
# Adjust rainfall to recharge using 10% magic number
test_model.parameters.create_model_parameter('magic_rain', 0.1)

interp_rain = interp_rain * test_model.parameters.param['magic_rain']['PARVAL1']

rch = {}
rch[0] = interp_rain

print "************************************************************************"
print " Creating recharge boundary "

test_model.boundaries.create_model_boundary_condition('Rain_reduced', 'recharge', bc_static=True)
test_model.boundaries.assign_boundary_array('Rain_reduced', rch)

print "************************************************************************"
print " Mapping bores to grid "

bore_points = [[final_bores.loc[x, "Easting"], final_bores.loc[x, "Northing"]] for x in final_bores.index]

bore_points3D = final_bores[["HydroCode", "Easting", "Northing", "depth"]] # [[final_bores.loc[x, "Easting"], final_bores.loc[x, "Northing"], final_bores.loc[x, "depth"]] for x in final_bores.index]
bore_points3D = bore_points3D.set_index("HydroCode")

bores_obs_time_series = final_bores[["HydroCode", "mean level"]]

# Modify into standard format for the GWModelBuilder class
bores_obs_time_series = bores_obs_time_series.rename(columns={'HydroCode':'name', 'mean level':'value'})

test_model.observations.set_as_observations('average head', bores_obs_time_series, bore_points3D, domain='porous', obs_type='head', units='mAHD')

test_model.map_obs_loc2mesh3D(method='nearest')

bores_in_layers = test_model.map_points_to_raster_layers(bore_points, final_bores["depth"].tolist(), hu_raster_files_reproj)

# Map bores to layers to create initial head maps for different hydrogeological units
interp_heads = {}

for i in range(len(hu_raster_files_reproj)/2):
    bores_layer = np.array(bore_points)[np.array(bores_in_layers[i])]
    print 'Creating head map for: ', hu_raster_files[2*i]
    if bores_layer.shape[0] < 4: 
        #interp_heads[hu_raster_files[2*i]] = (test_model.model_mesh3D[0][i]+test_model.model_mesh3D[0][i+1])/2
        interp_heads[hu_raster_files[2*i]] = np.full(test_model.model_mesh3D[1].shape[1:], np.NaN)
    else:
        bores_head_layer = np.array(final_bores["mean level"].tolist())[np.array(bores_in_layers[i])]
        unique_bores = np.unique(bores_layer) 
    
        b = np.ascontiguousarray(bores_layer).view(np.dtype((np.void, bores_layer.dtype.itemsize * bores_layer.shape[1])))
        _, idx = np.unique(b, return_index=True)
    
        unique_bores = bores_layer[idx]    
    
        interp_heads[hu_raster_files[2*i]] = test_model.interpolate_points2mesh(bores_layer, bores_head_layer, use='griddata', method='linear')
        

        
#for key in interp_heads:
    #bores_layer_df = pd.DataFrame()
    #bores_layer_df["Easting"] = [x[0] for x in bores_layer] 
    #bores_layer_df["Northing"] = [x[1] for x in bores_layer]
    #bores_layer_df["mean level"] = bores_head_layer
    #(XI, YI) = test_model.model_mesh_centroids
    #plt.figure()
    #z_min = np.min(interp_heads[key])
    #z_max = np.max(interp_heads[key])
    #plt.pcolor(XI, YI, interp_heads[key], vmin=z_min, vmax=z_max)
    #plt.scatter([x[0] for x in bores_layer], [x[1] for x in bores_layer], 50, bores_head_layer, vmin=z_min, vmax=z_max, cmap="jet")
    #plt.colorbar()

    #bores_layer_df.plot(kind='scatter', x="Easting", y="Northing", c="mean level", cmap="Spectral") # , edgecolor='None'
    #plt.scatter(x=[x[0] for x in bores_layer], y=[x[1] for x in bores_layer], c=bores_head_layer)

# Initalise model with head from elevations
initial_heads_SS = np.full(test_model.model_mesh3D[1].shape, 0.)

for i in range(len(hu_raster_files_reproj)/2):
    initial_heads_SS[i] = (test_model.model_mesh3D[0][i]+test_model.model_mesh3D[0][i+1])/2

test_model.initial_conditions.set_as_initial_condition("Head", initial_heads_SS)#interp_heads[hu_raster_files[0]])

initial_heads_SS = np.full(test_model.model_mesh3D[1].shape, 0.)

for i in range(len(hu_raster_files_reproj)/2):
    initial_heads_SS[i] = (interp_heads[hu_raster_files[2*i]])

test_model.initial_conditions.set_as_initial_condition("OldHead", initial_heads_SS)#interp_heads[hu_raster_files[0]])


print "************************************************************************"
print " Mapping pumping wells to grid "

test_model.map_points_to_grid(pumps_points, feature_id = 'OLD ID')

pump_array = []
for pump_cell in test_model.points_mapped['pumping wells_clipped.shp']:
    row = pump_cell[0][0]
    col = pump_cell[0][1]
    layers = [0]
    for pump in pump_cell[1]: 
        #HydroCode = pumping_data.loc[pump, 'Works ID']
        #if HydroCode not in bore_data_info.index:
        #    print HydroCode, ' not in index of bore_data_info'            
        #    continue
        #pump_depth = bore_data_info.loc[HydroCode, 'depth'] #[bore_data_info["HydroCode"] == HydroCode]['depth']        
        if pumping_data.loc[pump, 'Top screen depth (m)'] == 0.: 
            #print 'No data to place pump at depth ... ignoring ', pump            
            continue
        pump_depth = test_model.model_mesh3D[0][0][row][col] - pumping_data.loc[pump, 'Top screen depth (m)']        
        active = False
        for i in range(test_model.model_mesh3D[0].shape[0]-1):
            if pump_depth < test_model.model_mesh3D[0][i][row][col] and pump_depth > test_model.model_mesh3D[0][i+1][row][col]:
                active_layer = i
                active = True
                break
        if active == False: 
            #print 'Well not placed: ', pump            
            continue
        #Get top of screen layer and calculate length of screen in layer
        pump_rate = -pumping_data.loc[pump, 'Annual Volume'] / 365. * 1000. 
    
        #Get bottom of screen and calculate length of screen in layer    
        lay = 0 # Need to do a search of mesh to place pump in correct layers
        pump_array += [[lay, row, col, pump_rate]]

wel = {}

wel[0] = pump_array

print "************************************************************************"
print " Creating pumping boundary "

#test_model.boundaries.create_model_boundary_condition('licenced_wells', 'wells', bc_static=True)
#test_model.boundaries.assign_boundary_array('licenced_wells', wel)


# Map river polyline feature to grid including length of river in cell
print "************************************************************************"
print " Mapping Campaspe river to grid"

river_poly = test_model.read_polyline("Campaspe_Riv.shp", path=r"C:\Workspace\part0075\MDB modelling\Campaspe_model\GIS\GIS_preprocessed\Surface_Water\Streams\\") 
test_model.map_polyline_to_grid(river_poly)

test_model.parameters.create_model_parameter('bed_depression', 0.01)

simple_river = []
Kv_riv = 5E-3 #m/day
riv_width_avg = 10.0 #m
riv_bed_thickness = 0.10 #m
for riv_cell in test_model.polyline_mapped['Campaspe_Riv_model.shp']:
    row = riv_cell[0][0]
    col = riv_cell[0][1]
    if test_model.model_mesh3D[1][0][row][col] == -1:
        continue
    #print test_model.model_mesh3D
    stage = test_model.model_mesh3D[0][0][row][col]
    bed = test_model.model_mesh3D[0][0][row][col] - test_model.parameters.param['bed_depression']['PARVAL1']
    cond = riv_cell[1]*riv_width_avg*Kv_riv/riv_bed_thickness
    simple_river += [[0, row, col, stage, cond, bed]]

riv = {}
riv[0] = simple_river

#print test_model.polyline_mapped
print "************************************************************************"
print " Creating Campaspe river boundary"

test_model.boundaries.create_model_boundary_condition('Campaspe River', 'river', bc_static=True)
test_model.boundaries.assign_boundary_array('Campaspe River', riv)

print "************************************************************************"
print " Mapping Murray River to grid"

river_poly = test_model.read_polyline("River_Murray.shp", path=r"C:\Workspace\part0075\MDB modelling\Campaspe_model\GIS\GIS_preprocessed\Surface_Water\Streams\\") 
test_model.map_polyline_to_grid(river_poly)

test_model.parameters.create_model_parameter('bed_depression', 0.01)

simple_river = []
Kv_riv = 5E-3 #m/day
riv_width_avg = 10.0 #m
riv_bed_thickness = 0.10 #m
for riv_cell in test_model.polyline_mapped['River_Murray_model.shp']:
    row = riv_cell[0][0]
    col = riv_cell[0][1]
    if test_model.model_mesh3D[1][0][row][col] == -1:
        continue
    #print test_model.model_mesh3D
    stage = test_model.model_mesh3D[0][0][row][col]
    bed = test_model.model_mesh3D[0][0][row][col] - test_model.parameters.param['bed_depression']['PARVAL1']
    cond = riv_cell[1]*riv_width_avg*Kv_riv/riv_bed_thickness
    simple_river += [[0, row, col, stage, cond, bed]]

riv = {}
riv[0] = simple_river

#print test_model.polyline_mapped
print "************************************************************************"
print " Creating Murray River boundary"

test_model.boundaries.create_model_boundary_condition('Murray River', 'river', bc_static=True)
test_model.boundaries.assign_boundary_array('Murray River', riv)


print "************************************************************************"
print " Package up groundwater model builder object"

test_model.package_model()

#******************************************************************************
#
# AFTER data input pass data to model
#
#******************************************************************************
#print '########################################################################'
#print '########################################################################'
#print '## Model specific setup and running'
#print '########################################################################'
#print '########################################################################'
#
#print "************************************************************************"
#print " Set up Modflow and RUN"
#
# THIS CAN BE DONE FROM THE GWModelManager MODULE
#


#modflow_model = flopyInterface.ModflowModel(test_model, data_folder=r"C:\Workspace\part0075\MDB modelling\testbox\\")
#
#modflow_model.runMODFLOW()
#
#modflow_model.viewHeads()