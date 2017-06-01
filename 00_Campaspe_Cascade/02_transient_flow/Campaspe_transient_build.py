import os
import datetime

from osgeo import osr
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt

from HydroModelBuilder.GWModelBuilder import GWModelBuilder
from HydroModelBuilder.GISInterface.GDALInterface.GDALInterface import GDALInterface
from CampaspeModel.CustomScripts import processWeatherStations, getBoreData, get_GW_licence_info, processRiverStations, readHydrogeologicalProperties

"""

"""
# Define basic model parameters:
Proj_CS = osr.SpatialReference()
Proj_CS.ImportFromEPSG(28355) # 28355 is the code for gda94 mga zone 55; http://spatialreference.org/ref/epsg/gda94-mga-zone-55/

Interface = GDALInterface()
Interface.projected_coordinate_system = Proj_CS 
Interface.pcs_EPSG = "EPSG:28355"

tr_model = GWModelBuilder(name="02_transient_flow", 
                          data_folder=r"C:\Workspace\part0075\MDB modelling\testbox\input_data\\",
                          model_data_folder=r"C:\Workspace\part0075\MDB modelling\testbox\00_Campaspe_Cascade\02_transient_flow\\",
                          out_data_folder=r"C:\Workspace\part0075\MDB modelling\testbox\data_build\\",
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
                                                  buffer_dist=20000, shapefile_path=tr_model.data_folder)

# Setup recharge:
# ... read in climate data using Custom_Scripts
weather_stations = ['Kyneton',  'Elmore', 'Rochester', 'Echuca']
print "************************************************************************"
print " Executing custom script: processWeatherStations "

rain_info_file = "rain_processed_transient"
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
bore_salinity_file = "bore_salinity"
bore_info_file = "bore_info"
if os.path.exists(tr_model.out_data_folder + bore_levels_file + ".h5") & \
   os.path.exists(tr_model.out_data_folder + bore_info_file + ".h5") & \
   os.path.exists(tr_model.out_data_folder + bore_salinity_file + ".h5"):
    bore_data_levels = tr_model.load_dataframe(tr_model.out_data_folder + bore_levels_file + ".h5")
    bore_data_info = tr_model.load_dataframe(tr_model.out_data_folder + bore_info_file + ".h5")
    bore_data_salinity = tr_model.load_dataframe(tr_model.out_data_folder + bore_salinity_file + ".h5")
else:
    bore_data_levels, bore_data_info, bore_data_salinity = getBoreData.getBoreData()
    tr_model.save_dataframe(tr_model.out_data_folder + bore_levels_file, bore_data_levels)
    tr_model.save_dataframe(tr_model.out_data_folder + bore_info_file, bore_data_info)
    tr_model.save_dataframe(tr_model.out_data_folder + bore_salinity_file, bore_data_salinity)
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

print "************************************************************************"
print "Get the C14 data"

C14_points = tr_model.read_points_data(r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Chemistry\C14.shp")    

#C14_wells_info_file = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Chemistry\C14_bore_depth.csv"
#df_C14_info = pd.read_csv(C14_wells_info_file)    
#df_C14_info = df_C14_info.dropna()
#df_C14_info = df_C14_info.set_index('Bore_id')    

C14data = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Chemistry\C14_locs.xlsx"
df_C14 = pd.read_excel(C14data)
df_C14.drop_duplicates(subset=["Bore_id"], inplace=True)
df_C14.dropna(inplace=True)

print "************************************************************************"
print " Executing custom script: processRiverStations "

river_flow_file = "river_flow_processed"
# Check if this data has been processed and if not process it
if os.path.exists(tr_model.out_data_folder + river_flow_file + '.h5'):
    river_flow_data = tr_model.load_dataframe(tr_model.out_data_folder + river_flow_file + '.h5')
else:
    river_flow_data = processRiverStations.getFlow(path=r"C:\Workspace\part0075\MDB modelling\Campaspe_data\SW\All_streamflow_Campaspe_catchment\\")
    tr_model.save_dataframe(tr_model.out_data_folder + river_flow_file, river_flow_data)

river_stage_file = "river_stage_processed"
# Check if this data has been processed and if not process it
if os.path.exists(tr_model.out_data_folder + river_stage_file + '.h5'):
    river_stage_data = tr_model.load_dataframe(tr_model.out_data_folder + river_stage_file + '.h5')
else:
    river_stage_data = processRiverStations.getStage(path=r"C:\Workspace\part0075\MDB modelling\Campaspe_data\SW\All_streamflow_Campaspe_catchment\\")
    tr_model.save_dataframe(tr_model.out_data_folder + river_stage_file, river_stage_data)

river_gauges = tr_model.read_points_data(r"C:\Workspace\part0075\MDB modelling\Campaspe_data\SW\All_streamflow_Campaspe_catchment\processed_river_sites_stage.shp")

print "************************************************************************"
print "Load in the river shapefiles"

Campaspe_river_poly = tr_model.read_polyline("Campaspe_Riv.shp", path=r"C:\Workspace\part0075\MDB modelling\Campaspe_model\GIS\GIS_preprocessed\Surface_Water\Streams\\") 
Murray_river_poly = tr_model.read_polyline("River_Murray.shp", path=r"C:\Workspace\part0075\MDB modelling\Campaspe_model\GIS\GIS_preprocessed\Surface_Water\Streams\\") 

print "************************************************************************"
print "Load in the shapefiles defining groundwater boundaries"

WGWbound_poly = tr_model.read_polyline("western_head.shp", path=r"C:\Workspace\part0075\MDB modelling\testbox\input_data\\") 
EGWbound_poly = tr_model.read_polyline("eastern_head.shp", path=r"C:\Workspace\part0075\MDB modelling\testbox\input_data\\") 

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

start = datetime.date(1840, 01, 01)
#start = datetime.date(2014, 01, 01)
end_post_clearance = datetime.date(1880, 12, 31)
start_irrigation = datetime.date(1881, 01, 01)
before_pumping = datetime.date(1965, 12, 31)
start_pumping = datetime.date(1966, 01, 01)
end = datetime.date(2015, 12, 31)

date_index_post_clearance = pd.date_range(start=start, end=end_post_clearance, freq='10A')
date_index_post_irrigation = pd.date_range(start=start_irrigation, end=before_pumping, freq='4A')
date_index_post_pumping = pd.date_range(start=start_pumping, end=end, freq='M')

date_index_temp = date_index_post_clearance.append(date_index_post_irrigation)
date_index = date_index_temp.append(date_index_post_pumping)

#tr_model.model_time.set_temporal_components(steady_state=False, start_time=start, end_time=end, time_step='M')
tr_model.model_time.set_temporal_components(steady_state=False, start_time=start, end_time=end, date_index=date_index)

# Define the grid width and grid height for the model mesh which is stored as a multipolygon shapefile GDAL object
print "************************************************************************"
print " Defining structured mesh"
resolution = 5000
tr_model.define_structured_mesh(resolution, resolution) 

# Read in hydrostratigraphic raster info for layer elevations:
#hu_raster_path = r"C:\Workspace\part0075\MDB modelling\Campaspe_model\GIS\GIS_preprocessed\Hydrogeological_Unit_Layers\\"
hu_raster_path = r"C:\Workspace\part0075\MDB modelling\VAF_v2.0_ESRI_GRID\ESRI_GRID\\"
#hu_raster_path = r"C:\Workspace\part0075\MDB modelling\VAF_v2.0_ESRI_GRID\ESRI_GRID_raw\ESRI_GRID\\"

# Build basement file ... only need to do this once as it is time consuming so commented out for future runs
#tr_model.create_basement_bottom(hu_raster_path, "sur_1t", "bse_1t", "bse_2b", hu_raster_path)

hu_raster_files = ["qa_1t", "qa_2b", "utb_1t", "utb_2b", "utqa_1t", "utqa_2b", 
                   "utam_1t", "utam_2b", "utaf_1t", "utaf_2b", "lta_1t", 
                   "lta_2b", "bse_1t", "bse_2b.tif"]
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

# Cleanup any isolated cells:
tr_model.reclassIsolatedCells()

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

HGU_zone = {'bse':6, 'utb':1, 
            'utaf':4, 'lta':5, 
            'qa':0, 'utqa':2,
            'utam':3}

pilot_points = True
           
if pilot_points:
    pp_file = r"C:\Workspace\part0075\MDB modelling\testbox\00_Campaspe_Cascade\01_steady_state\structured_model_grid_{}m\pilot_points.pkl".format(resolution)       
    tr_model.load_pilot_points(pp_file)
    hk = tr_model.pilot_points['hk']
           
for unit in HGU:
    if pilot_points:
        tr_model.parameters.create_model_parameter_set('kh_' + unit, 
                                                       value=HGU_props['Kh mean'][HGU_map[unit]], 
                                                       num_parameters=hk.num_ppoints_by_zone[HGU_zone[unit]])
        tr_model.parameters.parameter_options_set('kh_' + unit, 
                                              PARTRANS='log', 
                                              PARCHGLIM='factor', 
                                              PARLBND=HGU_props['Kh mean'][HGU_map[unit]] / 10., 
                                              PARUBND=HGU_props['Kh mean'][HGU_map[unit]] * 10., 
                                              PARGP='cond_' + unit, 
                                              SCALE=1, 
                                              OFFSET=0)
    else:
        tr_model.parameters.create_model_parameter('kh_' + unit, 
                                                   value=HGU_props['Kh mean'][HGU_map[unit]])
        tr_model.parameters.parameter_options('kh_' + unit, 
                                              PARTRANS='log', 
                                              PARCHGLIM='factor', 
                                              PARLBND=HGU_props['Kh mean'][HGU_map[unit]] / 10., 
                                              PARUBND=HGU_props['Kh mean'][HGU_map[unit]] * 10., 
                                              PARGP='cond_' + unit, 
                                              SCALE=1, 
                                              OFFSET=0)

    tr_model.parameters.create_model_parameter('kv_' + unit, value=HGU_props['Kz mean'][HGU_map[unit]])
    tr_model.parameters.parameter_options('kv_' + unit, 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=HGU_props['Kz mean'][HGU_map[unit]] / 10., 
                                          PARUBND=HGU_props['Kz mean'][HGU_map[unit]] * 10., 
                                          PARGP='cond_' + unit, 
                                          SCALE=1, 
                                          OFFSET=0)
    tr_model.parameters.create_model_parameter('sy_' + unit, value=HGU_props['Sy mean'][HGU_map[unit]])
    tr_model.parameters.parameter_options('sy_' + unit, 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=0.01, 
                                          PARUBND=0.8, 
                                          PARGP='sy_' + unit, 
                                          SCALE=1, 
                                          OFFSET=0)
    tr_model.parameters.create_model_parameter('ss_' + unit, value=HGU_props['SS mean'][HGU_map[unit]])
    tr_model.parameters.parameter_options('ss_' + unit, 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=HGU_props['SS mean'][HGU_map[unit]] / 10., 
                                          PARUBND=HGU_props['SS mean'][HGU_map[unit]] * 10., 
                                          PARGP='ss_' + unit, 
                                          SCALE=1, 
                                          OFFSET=0)

# This needs to be automatically generated from with the map_raster2mesh routine ...
zone_map = {1:'qa', 2:'utb', 3:'utqa', 4:'utam', 5:'utaf', 6:'lta', 7:'bse'}

Kh = tr_model.model_mesh3D[1].astype(float)
Kv = tr_model.model_mesh3D[1].astype(float)
Sy = tr_model.model_mesh3D[1].astype(float)
SS = tr_model.model_mesh3D[1].astype(float)
for key in zone_map.keys():
    if not pilot_points:
        Kh[Kh == key] = tr_model.parameters.param['kh_' + zone_map[key]]['PARVAL1']
    Kv[Kv == key] = tr_model.parameters.param['kv_' + zone_map[key]]['PARVAL1']
    Sy[Sy == key] = tr_model.parameters.param['sy_' + zone_map[key]]['PARVAL1']
    SS[SS == key] = tr_model.parameters.param['ss_' + zone_map[key]]['PARVAL1']

if pilot_points:
    Kh = hk.val_array

tr_model.properties.assign_model_properties('Kh', Kh)
tr_model.properties.assign_model_properties('Kv', Kv)
tr_model.properties.assign_model_properties('Sy', Sy)
tr_model.properties.assign_model_properties('SS', SS)

print "************************************************************************"
print " Interpolating rainfall data to grid "


# To select a subset of the rainfall we can use:
long_term_historic_rainfall = long_term_historic_rainfall.ix[start:end]

#However if the time period for the model is longer we need to reindex the dataseries
#model_date_index = pd.date_range(start,end, freq='M')
long_term_historic_rainfall = long_term_historic_rainfall.reindex(date_index)
# Fill the nan values with the mean value of rainfall for recorded data
long_term_historic_rainfall = long_term_historic_rainfall.fillna(long_term_historic_rainfall.mean())


interp_rain = {}
for step, month in enumerate(long_term_historic_rainfall.iterrows()):
    #print step    
    interp_rain[step] = tr_model.interpolate_points2mesh(rain_gauges, month[1], feature_id='Name')
    # Adjust rainfall to m from mm 
    interp_rain[step] = interp_rain[step]/1000.0

tr_model.boundaries.create_model_boundary_condition('Rainfall', 'rainfall', bc_static=True)
tr_model.boundaries.assign_boundary_array('Rainfall', interp_rain)

# Adjust rainfall to recharge using rainfall reduction
for i in [1,2,3,7]:
    tr_model.parameters.create_model_parameter('ssrch_'+zone_map[i], value=0.01)
    tr_model.parameters.parameter_options('ssrch_'+zone_map[i], 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=0., 
                                          PARUBND=0.9, 
                                          PARGP='rech_mult', 
                                          SCALE=1, 
                                          OFFSET=0)


for i in [1,2,3,7]:
    tr_model.parameters.create_model_parameter('rch_red_'+zone_map[i], value=0.05)
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

tr_model.map_points_to_grid(bores_shpfile, feature_id='HydroCode')

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

print('Final bores within aquifers: '.format(len(bores_more_filter)))

final_bores = final_bores[final_bores["HydroCode"].isin(bores_more_filter)]

bore_points = [[final_bores.loc[x, "Easting"], final_bores.loc[x, "Northing"]] for x in final_bores.index]

bore_points3D = final_bores[["HydroCode", "Easting", "Northing", "depth"]] # [[final_bores.loc[x, "Easting"], final_bores.loc[x, "Northing"], final_bores.loc[x, "depth"]] for x in final_bores.index]
bore_points3D = bore_points3D.set_index("HydroCode")

bores_obs_time_series = bore_data_levels[bore_data_levels["HydroCode"].isin(final_bores["HydroCode"])]

# Modify into standard format for the GWModelBuilder class
bores_obs_time_series = bores_obs_time_series.rename(columns={'HydroCode':'name', 'bore_date':'datetime', 'result':'value'})

bores_obs_time_series['datetime'] = pd.to_datetime(bores_obs_time_series['datetime'])

# Perhaps some extra filtering to reduce the number of obs
#bores_obs_time_series[bores_obs_time_series['datetime'] > datetime.datetime(2000,01,01)]

# Kill all values where head observation is 0.0
bores_obs_time_series = bores_obs_time_series[bores_obs_time_series['value'] > 5.]

bores_obs_time_series = bores_obs_time_series[bores_obs_time_series['datetime'] != datetime.datetime(1899,12,30)]

bores_obs_time_series = bores_obs_time_series[bores_obs_time_series['datetime'] > datetime.datetime(2010,12,30)]
                                              
bores_obs_outliers_multiplier = 4.0

bores_obs_group = bores_obs_time_series.groupby('name')
bores_obs_means = bores_obs_group['value'].mean()
bores_obs_stds = bores_obs_group['value'].std()

bores_obs_cull = []
for bore in bores_obs_time_series['name'].unique():
    for row in bores_obs_time_series[bores_obs_time_series['name'] == bore].iterrows():                                              
        upper = bores_obs_means.loc[bore] + bores_obs_outliers_multiplier * bores_obs_stds.loc[bore]
        lower = bores_obs_means.loc[bore] - bores_obs_outliers_multiplier * bores_obs_stds.loc[bore]
        if row[1]['value'] < lower or row[1]['value'] > upper:
            bores_obs_cull += [int(row[0])]

# Now to remove those bores obs that have a reading greater or less than 5 * times sigma either side of the mean
bores_obs_time_series = bores_obs_time_series[~bores_obs_time_series.index.isin(bores_obs_cull)]
                               
# For the weigts of observations we need to specify them as 1/sigma, where sigma is the standard deviation of measurement error

tr_model.observations.set_as_observations('head', bores_obs_time_series, bore_points3D, 
                                          domain='porous', obs_type='head', units='mAHD', 
                                          weights=1.0/0.2, by_zone=True)

bores_in_layers = tr_model.map_points_to_raster_layers(bore_points, final_bores["depth"].tolist(), hu_raster_files_reproj)

# Map bores to layers to create initial head maps for different hydrogeological units
interp_heads = {}

for i in range(len(hu_raster_files_reproj)/2):
    bores_layer = np.array(bore_points)[np.array(bores_in_layers[i])]
    print 'Creating head map for: ', hu_raster_files[2*i]
    if bores_layer.shape[0] < 4: 
        #interp_heads[hu_raster_files[2*i]] = (tr_model.model_mesh3D[0][i]+tr_model.model_mesh3D[0][i+1])/2
        interp_heads[hu_raster_files[2*i]] = np.full(tr_model.model_mesh3D[1].shape[1:], np.NaN)
    else:
        bores_head_layer = np.array(final_bores["mean level"].tolist())[np.array(bores_in_layers[i])]
        unique_bores = np.unique(bores_layer) 
    
        b = np.ascontiguousarray(bores_layer).view(np.dtype((np.void, bores_layer.dtype.itemsize * bores_layer.shape[1])))
        _, idx = np.unique(b, return_index=True)
    
        unique_bores = bores_layer[idx]    
    
        interp_heads[hu_raster_files[2*i]] = tr_model.interpolate_points2mesh(bores_layer, bores_head_layer, use='griddata', method='linear')
        
#for key in interp_heads:
    #bores_layer_df = pd.DataFrame()
    #bores_layer_df["Easting"] = [x[0] for x in bores_layer] 
    #bores_layer_df["Northing"] = [x[1] for x in bores_layer]
    #bores_layer_df["mean level"] = bores_head_layer
    #(XI, YI) = tr_model.model_mesh_centroids
    #plt.figure()
    #z_min = np.min(interp_heads[key])
    #z_max = np.max(interp_heads[key])
    #plt.pcolor(XI, YI, interp_heads[key], vmin=z_min, vmax=z_max)
    #plt.scatter([x[0] for x in bores_layer], [x[1] for x in bores_layer], 50, bores_head_layer, vmin=z_min, vmax=z_max, cmap="jet")
    #plt.colorbar()

    #bores_layer_df.plot(kind='scatter', x="Easting", y="Northing", c="mean level", cmap="Spectral") # , edgecolor='None'
    #plt.scatter(x=[x[0] for x in bores_layer], y=[x[1] for x in bores_layer], c=bores_head_layer)

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

print "************************************************************************"
print "Create observation wells for C14"

tr_model.map_points_to_grid(C14_points, feature_id='Bore_id')

wel = {}

# Create another column in the pandas dataframe for the C14 data for the depth
# at which the sample was taken in mAHD ... which will be calculated in the next
# for loop

df_C14['z'] = 'null'

i = 0
well_name = {}
for C14wells in tr_model.points_mapped['C14_clipped.shp']:
    row = C14wells[0][0]
    col = C14wells[0][1]
    for well in C14wells[1]: 
        try:
            well_depth = df_C14.loc[df_C14[df_C14['Bore_id'] == int(well)].index.tolist()[0], 'avg_screen(m)']
            #well_depth = df_C14.loc[df_C14['Bore_id'] == int(well), 'avg_screen(m)']
        except:
            print 'Well was excluded due to lack of information: ', int(well)            
            continue
        
        well_depth = tr_model.model_mesh3D[0][0][row][col] - well_depth

        df_C14.set_value(df_C14['Bore_id'] == int(well), 'z', well_depth)
                
        active = False
        for i in range(tr_model.model_mesh3D[1].shape[0]):
            if well_depth < tr_model.model_mesh3D[0][i][row][col] and well_depth > tr_model.model_mesh3D[0][i+1][row][col]:
                active_layer = i
                active = True
                break
        if active == False: 
            #print 'Well not placed: ', pump            
            continue

        if tr_model.model_mesh3D[1][active_layer][row][col] == -1:
            continue

        # Well sits in the mesh, so assign to well boundary condition
        well_name[i] = well
        i=i+1            
        try:
            wel[0] += [[active_layer, row, col, 0.]]
        except:
            wel[0] = [[active_layer, row, col, 0.]]

tr_model.boundaries.create_model_boundary_condition('C14_wells', 'wells', bc_static=True)
tr_model.boundaries.assign_boundary_array('C14_wells', wel)

C14_obs_time_series = df_C14.copy() 
C14_obs_time_series = C14_obs_time_series[['Bore_id', 'a14C(pMC)']]
C14_obs_time_series['datetime'] = pd.to_datetime(datetime.date(2015,12,30))
C14_obs_time_series.rename(columns={'Bore_id':'name', 'a14C(pMC)':'value'}, inplace=True)
C14_bore_points3D = df_C14[['Bore_id', 'zone55_easting', 'zone55_northing', 'z']]
C14_bore_points3D = C14_bore_points3D.set_index("Bore_id")
C14_bore_points3D.rename(columns={'zone55_easting':'Easting', 'zone55_northing':'Northing'}, inplace=True)

tr_model.observations.set_as_observations('C14', C14_obs_time_series, C14_bore_points3D, domain='porous', obs_type='concentration', units='pMC', weights=1.0/5.0)

print "************************************************************************"
print " Mapping pumping wells to grid "

tr_model.map_points_to_grid(pumps_points, feature_id = 'OLD ID')


#tr_model.parameters.create_model_parameter('pump_use', value=0.6)
#tr_model.parameters.parameter_options('pump_use', 
#                                      PARTRANS='log', 
#                                      PARCHGLIM='factor', 
#                                      PARLBND=0.2, 
#                                      PARUBND=1., 
#                                      PARGP='pumping', 
#                                      SCALE=1, 
#                                      OFFSET=0)

# Convert pumping_data to time series

#pumping_data_ts = pd.DataFrame(cols=['Works ID', 'datetime'])

# Existing data is only for 10 years from 2005 to 2015
pump_date_index = pd.date_range(start=datetime.datetime(2005,07,01), end=datetime.datetime(2015,06,30), freq='AS-JUL')

wel = {}

# Need to establish what time step the pumping starts as an integer
def findInterval(row, times):
    key_time = pd.to_datetime(row)
    lower_time = times[0]
    for period, time in enumerate(times):
        if period > 0:
            if lower_time <= key_time < time:
                return period - 1
        lower_time = time
    return np.nan

wells_start = findInterval(start_pumping, date_index) 

pump_shallow = [] # Shallow (if <25m) or Deep (>= 25m)
for pump_cell in tr_model.points_mapped['pumping wells_clipped.shp']:
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
        pump_depth = tr_model.model_mesh3D[0][0][row][col] - pumping_data.loc[pump, 'Top screen depth (m)']        
        active = False
        for i in range(tr_model.model_mesh3D[0].shape[0]-1):
            if pump_depth < tr_model.model_mesh3D[0][i][row][col] and pump_depth > tr_model.model_mesh3D[0][i+1][row][col]:
                active_layer = i
                active = True
                break
        if active == False: 
            #print 'Well not placed: ', pump            
            continue

        # Specify if pump is shallow
        if pump_depth < 25:
            pump_shallow += [True]
        else:
            pump_shallow += [False]

        p05_06 = pumping_data.loc[pump, 'Use 2005/06'] / 365. * 1000.
        p06_07 = pumping_data.loc[pump, 'Use 2006/07'] / 365. * 1000.
        p07_08 = pumping_data.loc[pump, 'Use 2007/08'] / 365. * 1000.
        p08_09 = pumping_data.loc[pump, 'Use 2008/09'] / 365. * 1000.
        p09_10 = pumping_data.loc[pump, 'Use 2009/10'] / 365. * 1000.
        p10_11 = pumping_data.loc[pump, 'Use 2010/11'] / 365. * 1000.
        p11_12 = pumping_data.loc[pump, 'Use 2011/12'] / 365. * 1000.
        p12_13 = pumping_data.loc[pump, 'Use 2012/13'] / 365. * 1000.
        p13_14 = pumping_data.loc[pump, 'Use 2013/14'] / 365. * 1000.
        p14_15 = pumping_data.loc[pump, 'Use 2014/15'] / 365. * 1000.
        pump_rates = [p05_06, p06_07, p07_08, p08_09, p09_10, p10_11, p11_12, p12_13, p13_14, p14_15]        
        pumping_data_ts = pd.DataFrame(pump_rates, columns=[pump], index=pump_date_index)
        pump_install = pumping_data.loc[pump, 'Construction date']
        if isinstance(pump_install, datetime.time):
            pump_install = datetime.date(1950,01,01)    
        pump_date_index2 = pd.date_range(start=pump_install, end=datetime.datetime(2005,06,30), freq='AS-JUL')

        #pump_allocation = pumping_data.loc[pump, 'Annual Volume'] / 365. * 1000.

        # Assume historical pumping is a percentage of lowest non-zero use for well        
        non_zero_pumping = [x for x in pump_rates if x > 0.]         
        if non_zero_pumping == []:
            pumping_rate_old = 0.
        else:
            pumping_rate_old = np.min(non_zero_pumping)

        old_pumping_ts = pd.DataFrame(index=pump_date_index2)
        old_pumping_ts[pump] = pumping_rate_old * 0.6 #tr_model.parameters.param['pump_use']['PARVAL1']

        # Merge the old and measured data

        pumping_data_ts = pd.concat([pumping_data_ts, old_pumping_ts])

        # Now let's resample the data on a monthly basis, and we will take the mean    
        #pumping_data_ts = pumping_data_ts.resample(tr_model.model_time.t['time_step'], how='mean')
        pumping_data_ts = pumping_data_ts.resample('M').mean() #, how='mean')

        # Let's also get rid of NaN data and replace with backfilling
        pumping_data_ts = pumping_data_ts.fillna(method='bfill')

        # Let's only consider times in our date range though
        #date_index = pd.date_range(start=tr_model.model_time.t['start_time'], end=tr_model.model_time.t['end_time'], freq=tr_model.model_time.t['time_step'])
        date_index = pd.date_range(start=start_pumping, end=end, freq='M')
        pumping_data_ts = pumping_data_ts.reindex(date_index)    
        pumping_data_ts = pumping_data_ts.ix[start_pumping:end]
        pumping_data_ts = pumping_data_ts.fillna(0.0)

        # Now fill in the well dictionary with the values of pumping at relevant stress periods where Q is not 0.0
        for index, time in enumerate(pumping_data_ts.iterrows()):
            if index >= tr_model.model_time.t['steps']: 
                print index
                continue
            #if time[1]['m3/day used'] != 0.0 :
            try:
                wel[wells_start - 1 + index] += [[active_layer, row, col, -time[1][pump]]]
            except:
                wel[wells_start - 1 + index] = [[active_layer, row, col, -time[1][pump]]]
                
print "************************************************************************"
print " Creating pumping boundary "

tr_model.boundaries.create_model_boundary_condition('licenced_wells', 'wells', bc_static=True)
tr_model.boundaries.assign_boundary_array('licenced_wells', wel)

## Map river polyline feature to grid including length of river in cell
print "************************************************************************"
print " Mapping Campaspe river to grid"

print " Mapping Campaspe river to grid"

use_gauges = ['CAMPASPE RIVER @ EPPALOCK',
              'CAMPASPE RIVER @ DOAKS RESERVE',
              'CAMPASPE RIVER @ AXEDALE',
              'CAMPASPE RIVER @ BACKHAUS ROAD',
              'CAMPASPE RIVER @ BARNADOWN',
              'CAMPASPE RIVER @ ELMORE',
              'CAMPASPE RIVER @ CAMPASPE WEIR',
              'CAMPASPE RIVER @ CAMPASPE WEIR (HEAD GAUGE)',
              'CAMPASPE RIVER @ BURNEWANG-BONN ROAD',
              'CAMPASPE RIVER @ ROCHESTER D/S WARANGA WESTERN CH SYPHN',
              #'CAMPASPE RIVER @ FEHRINGS LANE',
              'CAMPASPE RIVER @ ECHUCA']

inflow_gauges = ['MILLEWA CREEK @ NORTHERN HIGHWAY ECHUCA',
                 'CAMPASPE DR NO 5 @ OUTFALL',
                 'CAMPASPE DR NO 4 U/S NORTHERN HIGHWAY',
                 'AXE CREEK @ LONGLEA',
                 'AXE CREEK @ STRATHFIELDSAYE']

tr_model.map_points_to_grid(river_gauges, feature_id='Site_Name')
#tr_model.map_points_to_grid(river_gauges, feature_id='Site_ID')

Campaspe_river_gauges = tr_model.points_mapped['processed_river_sites_stage_clipped.shp']

filter_gauges = []
for riv_gauge in Campaspe_river_gauges:
    #if riv_gauge[1][0] in use_gauges:
    if str(riv_gauge[1][0]) in use_gauges:
        filter_gauges += [riv_gauge]

tr_model.map_polyline_to_grid(Campaspe_river_poly)
tr_model.parameters.create_model_parameter('bed_depress', value=0.01)
tr_model.parameters.parameter_options('bed_depress', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.001, 
                                      PARUBND=0.1, 
                                      PARGP='camp_riv', 
                                      SCALE=1, 
                                      OFFSET=0)
tr_model.parameters.create_model_parameter('kv_riv', value=5E-3)
tr_model.parameters.parameter_options('kv_riv', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=1E-8, 
                                      PARUBND=20, 
                                      PARGP='camp_riv', 
                                      SCALE=1, 
                                      OFFSET=0)

simple_river = []
riv_width_avg = 10.0 #m
riv_bed_thickness = 0.10 #m

# Map river from high to low
new_riv = tr_model.polyline_mapped['Campaspe_Riv_model.shp']
for index, riv_cell in enumerate(tr_model.polyline_mapped['Campaspe_Riv_model.shp']):
    row = riv_cell[0][0]
    col = riv_cell[0][1]
    new_riv[index] += [tr_model.model_mesh3D[0][0][row][col]]

new_riv = sorted(new_riv, key=lambda x: (x[0][1]), reverse=False)    
new_riv = sorted(new_riv, key=lambda x: (x[0][0]), reverse=True)    

stages = np.full((len(new_riv)), np.nan, dtype=np.float64)
beds = np.full((len(new_riv)), np.nan, dtype=np.float64)

# Identify cells that correspond to river gauges
riv_gauge_logical = np.full((len(new_riv)), False, dtype=np.bool)


# To account for fact that river shapefile and gauges shapefile are not perfect
# we get the closest river cell to the gauge cell

def closest_node(node, nodes):
    nodes = np.asarray(nodes)
    dist_2 = np.sum((nodes - node)**2, axis=1)
    return np.argmin(dist_2)

# Define river gauges at start of river cell
new_riv_cells = [x[0] for x in new_riv]
filter_gauge_loc = [new_riv_cells[x] for x in [closest_node(x[0], new_riv_cells) for x in filter_gauges]]
    
# Set up the gauges as observations 

#SS_model.observations.set_as_observations('Campase_riv_gauges', CampaspeRiv_obs_time_series, CampaspeRiv_points3D, domain='surface', obs_type='stage', units='m') 
                    
for index, riv in enumerate(new_riv):
    # Create logical array to identify those which are gauges and those which are not
    if riv[0] in filter_gauge_loc:
        riv_gauge_logical[index] = True
        gauge_ind = [i for i, x in enumerate(filter_gauge_loc) if x == riv[0]]
        print filter_gauges[gauge_ind[0]][1][0]                     
        #stages[index] = river_stage_data["Mean stage (m)"].loc[river_stage_data["Site ID"] == filter_gauges[gauge_ind[0]][1][0]]
        stages[index] = river_stage_data["Mean stage (m)"].loc[river_stage_data["Site Name"] == filter_gauges[gauge_ind[0]][1][0]]
        beds[index] = stages[index] - 1.0 #river_stage_data["Mean stage (m)"].loc[river_stage_data["Site ID"]== ??]

    # Add chainage to new_riv array:
    if index == 0:
        new_riv[index] += [0.0]
    else:
        new_riv[index] += [new_riv[index-1][3] + new_riv[index-1][1]]        

# River x in terms of chainage:
river_x = np.array([x[3] for x in new_riv])
river_x_unknown = river_x[~riv_gauge_logical]
river_x_known = river_x[riv_gauge_logical]

# Now interpolate know values of stage and bed to unknown river locations:
stages[~riv_gauge_logical] = np.interp(river_x_unknown, river_x_known, stages[riv_gauge_logical])
beds[~riv_gauge_logical] = np.interp(river_x_unknown, river_x_known, beds[riv_gauge_logical])

# Create observations for stage or discharge at those locations



for index, riv_cell in enumerate(tr_model.polyline_mapped['Campaspe_Riv_model.shp']):
    row = riv_cell[0][0]
    col = riv_cell[0][1]
    if tr_model.model_mesh3D[1][0][row][col] == -1:
        continue
    #print SS_model.model_mesh3D
    stage_temp = stages[index]
    if stage_temp < tr_model.model_mesh3D[0][1][row][col]:
        stage = tr_model.model_mesh3D[0][1][row][col] + 0.01
    else:
        stage = stage_temp
    bed_temp = beds[index]
    if bed_temp < tr_model.model_mesh3D[0][1][row][col]:
        bed = tr_model.model_mesh3D[0][1][row][col]
    else:
        bed = bed_temp

    cond = riv_cell[1] * riv_width_avg * \
        tr_model.parameters.param['kv_riv']['PARVAL1'] / riv_bed_thickness
    simple_river += [[0, row, col, stage, cond, bed]]
#    stage = stages[index] 
#    bed = beds[index]
#    cond = riv_cell[1] * riv_width_avg * tr_model.parameters.param['Kv_riv']['PARVAL1'] / riv_bed_thickness
#    simple_river += [[0, row, col, stage, cond, bed]]

riv = {}
riv[0] = simple_river

fy_start = end - datetime.timedelta(days=363)
fy_end = end 

# Annual average net flux along river for final year of simulation
a_time_series = pd.DataFrame([{'name':'a_swgw', 'value':0.0, 'datetime':pd.to_datetime(datetime.date(2015,12,30))},])
tr_model.observations.set_as_observations('nrf_a', a_time_series, new_riv_cells, domain='surface', 
                            obs_type='swgw_a', units='m^3/d', weights=0.0, real=False)
# Seasonal average net flux along river for final year of simulation
entries = 4
s_time_series = pd.DataFrame({'name':['s_swgw{}'.format(x) for x in range(0, entries)], 
                            'value':[0.0] * entries, 
                            'datetime':pd.date_range(start=fy_start, end=fy_end, freq='Q-NOV')})
tr_model.observations.set_as_observations('nrf_s', s_time_series, new_riv_cells, domain='surface', 
                            obs_type='swgw_s', units='m^3/d', weights=0.0, real=False)
# Monthly average net flux along river for final year of simulation
entries = 12
m_time_series = pd.DataFrame({'name':['m_swgw{}'.format(x) for x in range(0, entries)],
                             'value':[0.0] * entries, 
                             'datetime':pd.date_range(start=fy_start, end=fy_end, freq='M')})
tr_model.observations.set_as_observations('nrf_m', m_time_series, new_riv_cells, domain='surface', 
                            obs_type='swgw_m', units='m^3/d', weights=0.0, real=False)

#print tr_model.polyline_mapped
print "************************************************************************"
print " Creating Campaspe river boundary"

tr_model.boundaries.create_model_boundary_condition('Campaspe River', 'river', bc_static=True)
tr_model.boundaries.assign_boundary_array('Campaspe River', riv)

print "************************************************************************"
print " Mapping Murray River to grid"

tr_model.map_polyline_to_grid(Murray_river_poly)

tr_model.parameters.create_model_parameter('rmstage', value=0.01)
tr_model.parameters.parameter_options('rmstage', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.001, 
                                      PARUBND=0.1, 
                                      PARGP='murr_riv', 
                                      SCALE=1, 
                                      OFFSET=0)
tr_model.parameters.create_model_parameter('kv_rm', value=5E-3)
tr_model.parameters.parameter_options('kv_rm', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=1E-8, 
                                      PARUBND=20, 
                                      PARGP='murr_riv', 
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
    bed = tr_model.model_mesh3D[0][0][row][col] - tr_model.parameters.param['rmstage']['PARVAL1']
    cond = riv_cell[1] * riv_width_avg * tr_model.parameters.param['kv_rm']['PARVAL1'] / riv_bed_thickness
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

tr_model.parameters.create_model_parameter('mghb_stage', value=0.01)
tr_model.parameters.parameter_options('mghb_stage', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.001, 
                                      PARUBND=0.1, 
                                      PARGP='ghb', 
                                      SCALE=1, 
                                      OFFSET=0)
tr_model.parameters.create_model_parameter('mghbcond', value=5E-3)
tr_model.parameters.parameter_options('mghbcond', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=1E-8, 
                                      PARUBND=20, 
                                      PARGP='ghb', 
                                      SCALE=1, 
                                      OFFSET=0)


MurrayGHB = []
Active_MurrayGHB_cells = []
for MurrayGHB_cell in tr_model.polyline_mapped['River_Murray_model.shp']:
    row = MurrayGHB_cell[0][0]
    col = MurrayGHB_cell[0][1]
    #print tr_model.model_mesh3D
    for lay in range(tr_model.model_mesh3D[1].shape[0]):    
        if tr_model.model_mesh3D[1][0][row][col] == -1:
            continue
        #MurrayGHBstage = (tr_model.model_mesh3D[0][lay+1][row][col] + tr_model.model_mesh3D[0][lay][row][col])/2. + tr_model.parameters.param['MGHB_stage']['PARVAL1']
        MurrayGHBstage = tr_model.model_mesh3D[0][0][row][col] + tr_model.parameters.param['mghb_stage']['PARVAL1']
        if MurrayGHBstage < tr_model.model_mesh3D[0][0][row][col]:
            continue
        if lay == 0:
            continue
        
        Active_MurrayGHB_cells += [[lay, row, col]]

# Now make sure that no cells are being caught surrounded by other GHB cells
Final_MurrayGHB_cells = []
for active_cell in Active_MurrayGHB_cells:
    # check if active GHB cell has any active non GHB cells N,E,S,W, above or below         
    lay, row, col = 0, 1, 2
    shape = tr_model.model_mesh3D[1].shape
    active_non_GHB = False
    zone = tr_model.model_mesh3D[1]
    
    ac = active_cell
    # Check above:
    ref_cell = [ac[lay] - 1, ac[row], ac[col]]        
#    # Make sure not at top boundary    
#    if active_cell[lay] != 0:
#        # Make sure above cell is active
#        if zone[ac[lay] - 1, ac[row], ac[col]] != -1:
#            # Make sure if active that is is not another GHB cell
#            if ref_cell not in Active_MurrayGHB_cells:
#                active_non_GHB = True
#
#    # Check below:
#    ref_cell = [ac[lay] + 1, ac[row], ac[col]]        
#    if active_cell[lay] != shape[lay] - 1:
#        if zone[ac[lay] + 1, ac[row], ac[col]] != -1:
#            if ref_cell not in Active_MurrayGHB_cells:
#                active_non_GHB = True
                
    # Check north:
    ref_cell = [ac[lay], ac[row] + 1, ac[col]]        
    if active_cell[row] != 0:
        if zone[ac[lay], ac[row] + 1, ac[col]] != -1:
            if ref_cell not in Active_MurrayGHB_cells:
                active_non_GHB = True

    # Check east:
    ref_cell = [ac[lay], ac[row], ac[col] + 1]        
    if active_cell[col] != shape[col] - 1:
        if zone[ac[lay], ac[row], ac[col] + 1] != -1:
            if ref_cell not in Active_MurrayGHB_cells:
                active_non_GHB = True

    # Check south:
    ref_cell = [ac[lay], ac[row] - 1, ac[col]]        
    if active_cell[row] != shape[row] - 1:
        if zone[ac[lay], ac[row] - 1, ac[col]] != -1:
            if ref_cell not in Active_MurrayGHB_cells:
                active_non_GHB = True

    # Check west:
    ref_cell = [ac[lay], ac[row], ac[col] - 1]        
    if active_cell[col] != 0:
        if zone[ac[lay], ac[row], ac[col] - 1] != -1:
            if ref_cell not in Active_MurrayGHB_cells:
                active_non_GHB = True

    if active_non_GHB:
        Final_MurrayGHB_cells += [active_cell]                
                
#for MurrayGHB_cell in tr_model.polyline_mapped['River_Murray_model.shp']:
for MurrayGHB_cell in Final_MurrayGHB_cells:
    #row = MurrayGHB_cell[0][0]
    #col = MurrayGHB_cell[0][1]

    lay = MurrayGHB_cell[0]
    row = MurrayGHB_cell[1]
    col = MurrayGHB_cell[2]
    #print tr_model.model_mesh3D
#    for lay in range(tr_model.model_mesh3D[1].shape[0]):    
#        if tr_model.model_mesh3D[1][0][row][col] == -1:
#            continue
        #MurrayGHBstage = (tr_model.model_mesh3D[0][lay+1][row][col] + tr_model.model_mesh3D[0][lay][row][col])/2. + tr_model.parameters.param['MGHB_stage']['PARVAL1']
#        if lay == 0:
            # To avoid having river cells in the same cells as GHB cells.
#            continue
        
    MurrayGHBstage = tr_model.model_mesh3D[0][0][row][col] + tr_model.parameters.param['mghb_stage']['PARVAL1']
    if MurrayGHBstage < tr_model.model_mesh3D[0][0][row][col]:
        continue
    dx = tr_model.gridHeight
    dz = tr_model.model_mesh3D[0][lay][row][col] - tr_model.model_mesh3D[0][lay+1][row][col]
    MGHBconductance = dx * dz * tr_model.parameters.param['mghbcond']['PARVAL1']
    MurrayGHB += [[lay, row, col, MurrayGHBstage, MGHBconductance]]


ghb = {}
ghb[0] = MurrayGHB

print "************************************************************************"
#print " Mapping Western GW boundary to grid"
#
#WGWbound_poly = tr_model.read_polyline("western_head.shp", path=r"C:\Workspace\part0075\MDB modelling\testbox\input_data\\") 
#tr_model.map_polyline_to_grid(WGWbound_poly)
#
#print "************************************************************************"
#print " Setting up Western GHB boundary"
#
#tr_model.parameters.create_model_parameter('WGHB_stage', value=0.01)
#tr_model.parameters.parameter_options('WGHB_stage', 
#                                      PARTRANS='log', 
#                                      PARCHGLIM='factor', 
#                                      PARLBND=0.001, 
#                                      PARUBND=0.1, 
#                                      PARGP='ghb', 
#                                      SCALE=1, 
#                                      OFFSET=0)
#tr_model.parameters.create_model_parameter('WGHBcond', value=5E-3)
#tr_model.parameters.parameter_options('WGHBcond', 
#                                      PARTRANS='log', 
#                                      PARCHGLIM='factor', 
#                                      PARLBND=1E-8, 
#                                      PARUBND=20, 
#                                      PARGP='ghb', 
#                                      SCALE=1, 
#                                      OFFSET=0)
#
#
#WestGHB = []
#for WestGHB_cell in tr_model.polyline_mapped['western_head_model.shp']:
#    row = WestGHB_cell[0][0]
#    col = WestGHB_cell[0][1]
#    #print tr_model.model_mesh3D
#    for lay in range(tr_model.model_mesh3D[1].shape[0]):    
#        if tr_model.model_mesh3D[1][lay][row][col] == -1:
#            continue
#        WestGHBstage = tr_model.model_mesh3D[0][0][row][col] + tr_model.parameters.param['WGHB_stage']['PARVAL1']
#        dx = tr_model.model_mesh3D[0][0][0][1] - tr_model.model_mesh3D[0][0][0][0]
#        WGHBconductance = dx * tr_model.parameters.param['WGHBcond']['PARVAL1']
#        WestGHB += [[lay, row, col, WestGHBstage, WGHBconductance]]
#
#
#ghb[0] += WestGHB
#
#
#print "************************************************************************"
#print " Mapping Eastern GW boundary to grid"
#
#EGWbound_poly = tr_model.read_polyline("eastern_head.shp", path=r"C:\Workspace\part0075\MDB modelling\testbox\input_data\\") 
#tr_model.map_polyline_to_grid(EGWbound_poly)
#
#print "************************************************************************"
#print " Setting up Western GHB boundary"
#
#tr_model.parameters.create_model_parameter('EGHB_stage', value=0.01)
#tr_model.parameters.parameter_options('EGHB_stage', 
#                                      PARTRANS='log', 
#                                      PARCHGLIM='factor', 
#                                      PARLBND=0.001, 
#                                      PARUBND=0.1, 
#                                      PARGP='ghb', 
#                                      SCALE=1, 
#                                      OFFSET=0)
#tr_model.parameters.create_model_parameter('EGHBcond', value=5E-3)
#tr_model.parameters.parameter_options('EGHBcond', 
#                                      PARTRANS='log', 
#                                      PARCHGLIM='factor', 
#                                      PARLBND=1E-8, 
#                                      PARUBND=20, 
#                                      PARGP='ghb', 
#                                      SCALE=1, 
#                                      OFFSET=0)
#
#
#EastGHB = []
#for EastGHB_cell in tr_model.polyline_mapped['eastern_head_model.shp']:
#    row = EastGHB_cell[0][0]
#    col = EastGHB_cell[0][1]
#    #print tr_model.model_mesh3D
#    for lay in range(tr_model.model_mesh3D[1].shape[0]):    
#        if tr_model.model_mesh3D[1][lay][row][col] == -1:
#            continue
#        #EastGHBstage = (tr_model.model_mesh3D[0][lay+1][row][col] + tr_model.model_mesh3D[0][lay][row][col])/2. + tr_model.parameters.param['EGHB_stage']['PARVAL1']
#        EastGHBstage = tr_model.model_mesh3D[0][0][row][col] + tr_model.parameters.param['EGHB_stage']['PARVAL1']
#        dx = tr_model.model_mesh3D[0][0][0][1] - tr_model.model_mesh3D[0][0][0][0]
#        EGHBconductance = dx * tr_model.parameters.param['EGHBcond']['PARVAL1']
#        EastGHB += [[lay, row, col, EastGHBstage, EGHBconductance]]
#
#ghb[0] += EastGHB

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
tr_model.parameters.create_model_parameter('kv_drain', value=5E-3)
tr_model.parameters.parameter_options('kv_drain', 
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
    drain_cond = drain_cell[1] * drain_width_avg * tr_model.parameters.param['kv_drain']['PARVAL1'] / drain_bed_thickness
    simple_drain += [[0, row, col, drain_bed, drain_cond]]

drain = {}
drain[0] = simple_drain

print "************************************************************************"
print " Creating Drains boundary"

tr_model.boundaries.create_model_boundary_condition('Drain', 'drain', bc_static=True)
tr_model.boundaries.assign_boundary_array('Drain', drain)
#
#
#print "************************************************************************"
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
tr_model.parameters.create_model_parameter('kv_chan', value=5E-3)
tr_model.parameters.parameter_options('kv_chan', 
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
    channel_cond = channel_cell[1] * channel_width_avg * tr_model.parameters.param['kv_chan']['PARVAL1'] / channel_bed_thickness
    simple_channel += [[0, row, col, channel_stage, channel_cond, channel_bed]]

channel = {}
channel[0] = simple_channel

print "************************************************************************"
print " Creating Channel boundary"

tr_model.boundaries.create_model_boundary_condition('Channel', 'channel', bc_static=True)
tr_model.boundaries.assign_boundary_array('Channel', channel)

print "************************************************************************"
print " Creating parameters for transport "

tr_model.parameters.create_model_parameter('porosity', value=0.25)
tr_model.parameters.parameter_options('porosity', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.05, 
                                      PARUBND=0.4, 
                                      PARGP='transport', 
                                      SCALE=1, 
                                      OFFSET=0)

tr_model.parameters.create_model_parameter('disp', value=0.01)
tr_model.parameters.parameter_options('disp', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=1E-5, 
                                      PARUBND=100., 
                                      PARGP='transport', 
                                      SCALE=1, 
                                      OFFSET=0)

print "************************************************************************"
print " Collate observations"
#
tr_model.map_obs_loc2mesh3D(method='nearest', ignore=[-1, 7])
tr_model.map_obs2model_times()
tr_model.observations.collate_observations()

# Looking at bores with large standard deviations for errors:
#import matplotlib.pyplot as plt
#b = bores_obs_time_series[bores_obs_time_series['value'] != 0.0]
#b = b[b['active'] == True]
#c = b[['name', 'value']].groupby('name').std()
#c.hist()
#d = c[c['value'] > 4]
#for bore in d.index.tolist():
#    b_df = bores_obs_time_series[bores_obs_time_series['name'] == bore]
#    b_df_min = b_df['datetime'].min()
#    b_df_max = b_df['datetime'].max()
#    b_df_sd = 4.0 * b_df['value'].std()
#    b_df_mean = b_df['value'].mean()
#    ax = b_df.plot(x='datetime', y='value', marker='o', label=bore)
#    plt.plot([b_df_min, b_df_max], [b_df_mean, b_df_mean], label='mean')
#    plt.plot([b_df_min, b_df_max], [b_df_mean + b_df_sd, b_df_mean + b_df_sd], label='2.5$\sigma$')
#    plt.plot([b_df_min, b_df_max], [b_df_mean - b_df_sd, b_df_mean - b_df_sd], label='2.5$\sigma$')       
    
print "************************************************************************"
print " Package up groundwater model builder object"

tr_model.package_model()
