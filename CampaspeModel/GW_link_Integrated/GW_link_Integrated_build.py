import datetime
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from osgeo import osr

from CampaspeModel.CustomScripts import (get_GW_licence_info, getBoreData,
                                         processRiverStations,
                                         processWeatherStations,
                                         readHydrogeologicalProperties)
from HydroModelBuilder.GISInterface.GDALInterface.GDALInterface import \
    GDALInterface
from HydroModelBuilder.GWModelBuilder import GWModelBuilder
from HydroModelBuilder.Utilities.Config.ConfigLoader import ConfigLoader

p_j = os.path.join
dir_name = os.path.dirname
CONFIG = ConfigLoader(p_j(dir_name(dir_name(os.path.realpath(__file__))),
                          "config", "model_config.json"))\
    .set_environment("GW_link_Integrated")

VERBOSE = True

# Define basic model parameters:
Proj_CS = osr.SpatialReference()

# 28355 is the code for gda94 mga zone 55; http://spatialreference.org/ref/epsg/gda94-mga-zone-55/
epsg_code = 28355
Proj_CS.ImportFromEPSG(epsg_code)

Interface = GDALInterface()
Interface.projected_coordinate_system = Proj_CS
Interface.pcs_EPSG = "EPSG:{}".format(epsg_code)

get_conf_set = CONFIG.get_setting
model_config = CONFIG.model_config
data_folder = model_config['data_folder']
mf_exe_folder = model_config['mf_exe_folder']
param_file = model_config['param_file']
climate_path = p_j(get_conf_set(['model_build', 'campaspe_data']), 'Climate')

temp_data_path = get_conf_set(['model_build', 'temp_data'])
input_data_path = get_conf_set(['model_build', 'input_data'])
river_path = p_j(input_data_path, "Waterways")
sw_data_path = p_j(temp_data_path, "Campaspe_data/SW/All_streamflow_Campaspe_catchment/Updated")

bore_levels_file = "bore_levels"
bore_info_file = "bore_info"
model_build_input_path = get_conf_set(['model_build', 'input_data'])

model_params = {
    "name": "GW_link_Integrated", 
    "data_folder": model_build_input_path,
    "campaspe_data": get_conf_set(['model_build', 'campaspe_data']), 
    "model_data_folder": model_config['data_folder'],
    "out_data_folder": get_conf_set(['model_build', 'data_build']),
    "GISInterface": Interface,
    "model_type": "Modflow",
    "mesh_type": "structured"
}
SS_model = GWModelBuilder(**model_params)

# Define the units for the project for consistency and to allow converions on input data

# ******************************************************************************
# Complimentary models requirements, i.e. bore and gauge data that should be
# referenceable to this model for parsing specific outputs and receiving inputs:

if VERBOSE:
    print "Attempting to map these bores..."
# End if
model_linking = r"../testbox/integrated/data/model_bores.csv"
with open(model_linking, 'r') as f:
    lines = f.readlines()

    def process_line(line):
        processed = [x.strip() for x in line.split(':')[1].strip().split(',')]
        return processed

    for line in lines:
        if line.split(':')[0] == 'Ecology':
            Ecology_bores = process_line(line)
            print('Ecology: {}'.format(Ecology_bores))
        elif line.split(':')[0] == 'Policy':
            Policy_bores = process_line(line)
            print('Policy: {}'.format(Policy_bores))
        elif line.split(':')[0] == 'SW_stream_gauges':
            Stream_gauges = process_line(line)
            print('SW: {}'.format(Stream_gauges))

# ******************************************************************************
# ******************************************************************************

# Set the model boundary using a polygon shapefile:
if VERBOSE:
    print "************************************************************************"
    print " Setting model boundary "

SS_model.set_model_boundary_from_polygon_shapefile("GW_model_area.shp",
                                                   shapefile_path=SS_model.data_folder)

# Set data boundary for model data
if VERBOSE:
    print "************************************************************************"
    print " Setting spatial data boundary "

SS_model.set_data_boundary_from_polygon_shapefile(SS_model.boundary_poly_file,
                                                  buffer_dist=20000)

# Setup recharge:
# ... read in climate data using Custom_Scripts
weather_stations = ['Kyneton',  'Elmore', 'Rochester', 'Echuca']

if VERBOSE:
    print "************************************************************************"
    print " Executing custom script: processWeatherStations "

rain_info_file = "rain_processed"
# Check if this data has been processed and if not process i
rain_info_h5 = p_j(SS_model.out_data_folder, rain_info_file + ".h5")

if os.path.exists(rain_info_h5):
    long_term_historic_rainfall = SS_model.load_dataframe(rain_info_h5)
else:
    long_term_historic_rainfall = processWeatherStations.processWeatherStations(weather_stations,
                                                                                path=climate_path)
    SS_model.save_dataframe(p_j(SS_model.out_data_folder, rain_info_file),
                            long_term_historic_rainfall)

rain_gauges = SS_model.read_points_data(p_j(climate_path, "Rain_gauges.shp"))
# INCLUDE NSW bores in this next part too for better head representation at the border,
# i.e. Murray River

# Read in bore data:
if VERBOSE:
    print "************************************************************************"
    print " Executing custom script: getBoreData "

bore_lf = p_j(SS_model.out_data_folder, bore_levels_file + ".h5")
bore_if = p_j(SS_model.out_data_folder, bore_info_file + ".h5")

if os.path.exists(bore_lf) & os.path.exists(bore_if):
    bore_data_levels = SS_model.load_dataframe(bore_lf)
    bore_data_info = SS_model.load_dataframe(bore_if)
else:
    bore_data_levels, bore_data_info, bore_data_salinity = \
        getBoreData.getBoreData(path=os.path.join(SS_model.campaspe_data, 
                                                  "ngis_shp_VIC"))
    SS_model.save_dataframe(p_j(SS_model.out_data_folder, bore_levels_file),
                            bore_data_levels)
    SS_model.save_dataframe(p_j(SS_model.out_data_folder, bore_info_file),
                            bore_data_info)
# End if

# getBoreDepth ... assuming that midpoint of screen interval
# is representative location and assign to layer accordingly
bore_data_info['depth'] = (bore_data_info['TopElev'] + bore_data_info['BottomElev']) / 2.0
bore_data_info["HydroCode"] = bore_data_info.index

# For steady state model, only use bore details containing average level, not
ngis_bore_shp = p_j(
    temp_data_path, "Campaspe_data", "ngis_shp_VIC", "ngis_shp_VIC", "NGIS_Bores.shp")

if VERBOSE:
    print "************************************************************************"
    print " Read in and filtering bore spatial data "

bores_shpfile = SS_model.read_points_data(ngis_bore_shp)
bores_filtered_from_shpfile = SS_model.points_shapefile_obj2dataframe(bores_shpfile,
                                                                      feature_id="HydroCode")

# Get the intersection of bores_filtered_from_shpfile with bores_data_info
final_bores = pd.merge(bore_data_info, bores_filtered_from_shpfile, how='inner', on="HydroCode")

# Only consider bores whereby the measured values are above the bottom of screen
final_bores = final_bores[final_bores['mean level'] > final_bores['BottomElev']]

if VERBOSE:
    print 'Final number of bores within the data boundary that have level data and screen info: ', \
        final_bores.shape[0]

# Load in the pumping wells data
filename = "Groundwater licence information for Dan Partington bc301115.xlsx"

path = out_path = p_j(temp_data_path, "Campaspe_data", "GW", "Bore data")
out_file = "pumping wells.shp"

if VERBOSE:
    print "************************************************************************"
    print " Executing custom script: get_GW_licence_info "

pumping_data = get_GW_licence_info.get_GW_licence_info(filename, path=path,
                                                       out_file=out_file, out_path=out_path)
pump_pnt_shp = p_j(path, "pumping wells.shp")
pumps_points = SS_model.read_points_data(pump_pnt_shp)

if VERBOSE:
    print "************************************************************************"
    print " Executing custom script: readHydrogeologicalProperties "

file_location = p_j(
    temp_data_path, "Campaspe_data/GW/Aquifer properties/Hydrogeologic_variables.xlsx")
HGU_props = readHydrogeologicalProperties.getHGUproperties(file_location)

if VERBOSE:
    print "************************************************************************"
    print " Executing custom script: processRiverStations "

river_flow_file = "river_flow_processed"
riv_flow_path = p_j(SS_model.out_data_folder, river_flow_file)

# Check if this data has been processed and if not process it
if os.path.exists(riv_flow_path + '.pkl'):
    river_flow_data = SS_model.load_obj(riv_flow_path + '.pkl')
else:
    river_flow_data = processRiverStations.getFlow(path=sw_data_path, summary=True)
    SS_model.save_obj(river_flow_data, riv_flow_path)
# End if

river_stage_file = "river_stage_processed"
riv_stage_path = p_j(SS_model.out_data_folder, river_stage_file)

# Check if this data has been processed and if not process it
if os.path.exists(riv_stage_path + '.pkl'):
    river_stage_data = SS_model.load_obj(riv_stage_path + '.pkl')
else:
    river_stage_data = processRiverStations.getStage(path=sw_data_path, summary=True)
    SS_model.save_obj(river_stage_data, riv_stage_path)

river_gauges = SS_model.read_points_data(p_j(sw_data_path, "Updated",
                                             r"processed_river_sites_stage.shp"))


river_data_folder = p_j(sw_data_path, "Updated")
site_details_file = "Site Details.csv"
site_details = pd.read_csv(os.path.join(river_data_folder, site_details_file))
# As all of the stream data for the whole of the Camaspe catchment is in the folder
# to be processed, we can prefilter sites to examine by specifying sites.
Campaspe_relevant = site_details[site_details['Site Name'].str.contains("CAMPASPE RIVER") | \
                        site_details['Site Name'].str.contains("MURRAY RIVER") | \
                        site_details['Site Name'].str.contains("AXE CREEK") | \
                        site_details['Site Name'].str.contains("MOUNT PLEASANT")]

Campaspe = site_details[site_details['Site Name'].str.contains("CAMPASPE RIVER")]
Campaspe = Campaspe[Campaspe['Northing'] >= \
                    Campaspe.loc[6]['Northing']]

if VERBOSE:
    print "************************************************************************"
    print "Load in the river shapefiles"

Campaspe_river_poly = SS_model.read_poly("Campaspe_Riv.shp", path=river_path)
Murray_river_poly = SS_model.read_poly("River_Murray.shp", path=river_path)
Campaspe_river_poly_file = p_j(river_path, "Campaspe_Riv.shp")
Murray_river_poly_file = p_j(river_path, "River_Murray.shp") 


if VERBOSE:
    print "************************************************************************"
    print "Load in the farms shapefile"

farms_path = p_j(temp_data_path, "Campaspe_data", "SW", "Farm")
farms_poly = SS_model.read_poly("farm_v1_prj.shp", path=farms_path, poly_type='polygon')


if VERBOSE:
    print "************************************************************************"
    print "Define the raster defining high res DEM"

surface_raster_high_res = p_j(temp_data_path, "ESRI_GRID_raw", "ESRI_GRID", "sur_1t")


if VERBOSE:
    print '########################################################################'
    print '## Mesh specific model building '
    print '########################################################################\n'

    print "************************************************************************"
    print " Defining temporal aspects of the model"

start, end = datetime.date(2014, 01, 01), datetime.date(2015, 01, 01)
SS_model.model_time.set_temporal_components(steady_state=False, start_time=start,
                                            end_time=end, time_step='A')

res = model_config['grid_resolution']
res = res[:-1] if res.endswith('m') else res

# Define the grid width and grid height for the model mesh which is stored as a multipolygon
# shapefile GDAL object
if VERBOSE:
    print "************************************************************************"
    print " Defining structured mesh at {res}x{res}".format(res=res)

SS_model.define_structured_mesh(int(res), int(res))

# Read in hydrostratigraphic raster info for layer elevations:
hu_raster_path = p_j(temp_data_path, "ESRI_GRID_raw", "ESRI_GRID")

# TODO RUN ON FLAG
SS_model.create_basement_bottom(hu_raster_path, "sur_1t", "bse_1t", "bse_2b", hu_raster_path)

hu_raster_files = ["qa_1t", "qa_2b", "utb_1t", "utb_2b", "utqa_1t", "utqa_2b", "utam_1t", "utam_2b",
                   "utaf_1t", "utaf_2b", "lta_1t", "lta_2b", "bse_1t", "bse_2b.tif"]
SS_model.read_rasters(hu_raster_files, path=hu_raster_path)
hu_raster_files_reproj = [x + "_reproj.tif" for x in hu_raster_files]

# Map HGU's to grid
if VERBOSE:
    print "************************************************************************"
    print " Mapping rasters to grid "

hu_gridded_rasters = SS_model.map_rasters_to_grid(hu_raster_files, hu_raster_path)

# Build 3D grid
model_grid_raster_files = [x + "_model_grid.tif" for x in hu_raster_files]

# First two arguments of next function are arbitrary and not used ... need to rework module
if VERBOSE:
    print "************************************************************************"
    print " Building 3D mesh "

SS_model.build_3D_mesh_from_rasters(model_grid_raster_files,
                                    SS_model.out_data_folder_grid,
                                    0.010,
                                    1000.0)
# Cleanup any isolated cells:
SS_model.reclassIsolatedCells()

if VERBOSE:
    print "************************************************************************"
    print " Assign properties to mesh based on zonal information"

# create list of HGU's from hu_raster_files
HGU = list(set([x.split('_')[0] for x in hu_raster_files]))

# NOTE *** utam is mapping to Shepparton Sands but it represents the
# Loxton-Parilla Sand ... the HGU database needs updating to include this.
HGU_map = {'bse': 'Bedrock', 'utb': 'Newer Volcanics Basalts',
           'utaf': 'Calivil', 'lta': 'Renmark',
           'qa': 'Coonambidgal Sands', 'utqa': 'Shepparton Sands',
           'utam': 'Shepparton Sands'}

SSparams = SS_model.parameters
for unit in HGU:
    SSparams.create_model_parameter('Kh_' + unit,
                                    value=HGU_props['Kh mean'][HGU_map[unit]])

    SSparams.parameter_options('Kh_' + unit,
                               PARTRANS='log',
                               PARCHGLIM='factor',
                               PARLBND=HGU_props['Kh mean'][HGU_map[unit]] / 10.,
                               PARUBND=HGU_props['Kh mean'][HGU_map[unit]] * 10.,
                               PARGP='cond',
                               SCALE=1,
                               OFFSET=0)

    SSparams.create_model_parameter('Kv_' + unit,
                                    value=HGU_props['Kz mean'][HGU_map[unit]])

    SSparams.parameter_options('Kv_' + unit,
                               PARTRANS='log',
                               PARCHGLIM='factor',
                               PARLBND=HGU_props['Kz mean'][HGU_map[unit]] / 10.,
                               PARUBND=HGU_props['Kz mean'][HGU_map[unit]] * 10.,
                               PARGP='cond',
                               SCALE=1,
                               OFFSET=0)

    SSparams.create_model_parameter('Sy_' + unit,
                                    value=HGU_props['Sy mean'][HGU_map[unit]])

    SSparams.parameter_options('Sy_' + unit,
                               PARTRANS='log',
                               PARCHGLIM='factor',
                               PARLBND=0.,
                               PARUBND=0.8,
                               PARGP='spec_yield',
                               SCALE=1,
                               OFFSET=0)

    SSparams.create_model_parameter('SS_' + unit,
                                    value=HGU_props['SS mean'][HGU_map[unit]])

    SSparams.parameter_options('SS_' + unit,
                               PARTRANS='log',
                               PARCHGLIM='factor',
                               PARLBND=HGU_props['SS mean'][HGU_map[unit]] / 10.,
                               PARUBND=HGU_props['SS mean'][HGU_map[unit]] * 10.,
                               PARGP='spec_stor',
                               SCALE=1,
                               OFFSET=0)

# This needs to be automatically generated from with the map_raster2mesh routine ...
zone_map = {1: 'qa', 2: 'utb', 3: 'utqa', 4: 'utam', 5: 'utaf', 6: 'lta', 7: 'bse'}

mesh3D_0 = SS_model.model_mesh3D[0]
mesh3D_1 = SS_model.model_mesh3D[1]

Kh, Kv, Sy, SS = [mesh3D_1.astype(float) for x in range(4)]
for key in zone_map.keys():
    z_map_key = zone_map[key]
    Kh[Kh == key] = SSparams.param['Kh_' + z_map_key]['PARVAL1']
    Kv[Kv == key] = SSparams.param['Kv_' + z_map_key]['PARVAL1']
    Sy[Sy == key] = SSparams.param['Sy_' + z_map_key]['PARVAL1']
    SS[SS == key] = SSparams.param['SS_' + z_map_key]['PARVAL1']

SSprop_assign = SS_model.properties.assign_model_properties
SSprop_assign('Kh', Kh)
SSprop_assign('Kv', Kv)
SSprop_assign('Sy', Sy)
SSprop_assign('SS', SS)

if VERBOSE:
    print "************************************************************************"
    print " Interpolating rainfall data to grid "

interp_rain = SS_model.interpolate_points2mesh(rain_gauges, long_term_historic_rainfall,
                                               feature_id='Name')
# Adjust rainfall to m from mm and from year to days
interp_rain = interp_rain / 1000.0 / 365.0
# Adjust rainfall to recharge using 10% magic number

SSbounds = SS_model.boundaries
SSbounds.create_model_boundary_condition('Rainfall', 'rainfall', bc_static=True)
SSbounds.assign_boundary_array('Rainfall', interp_rain)

for i in [1, 2, 3, 7]:
    SSparams.create_model_parameter('rch_red_' + zone_map[i], value=0.09)
    SSparams.parameter_options('rch_red_' + zone_map[i],
                               PARTRANS='log',
                               PARCHGLIM='factor',
                               PARLBND=0.,
                               PARUBND=0.9,
                               PARGP='rech_mult',
                               SCALE=1,
                               OFFSET=0)

    match = interp_rain[mesh3D_1[0] == i]
    interp_rain[mesh3D_1[0] == i] = match * SSparams.param['rch_red_' + zone_map[i]]['PARVAL1']
# End for

rch = {0: interp_rain}

if VERBOSE:
    print "************************************************************************"
    print " Creating recharge boundary "

SSbounds.create_model_boundary_condition('Rain_reduced', 'recharge', bc_static=True)
SSbounds.assign_boundary_array('Rain_reduced', rch)

if VERBOSE:
    print "************************************************************************"
    print " Mapping Campaspe river to grid"

# SW gauges being used in the SW model:
# '406214' AXE CREEK @ LONGLEA
# '406219' CAMPASPE RIVER @ LAKE EPPALOCK (HEAD GAUGE)
# '406201' CAMPASPE RIVER @ BARNADOWN
# '406224' MOUNT PLEASANT CREEK @ RUNNYMEDE
# '406218' CAMPASPE RIVER @ CAMPASPE WEIR (HEAD GAUGE)
# '406202' CAMPASPE RIVER @ ROCHESTER D/S WARANGA WESTERN CH SYPHN
# '406265' CAMPASPE RIVER @ ECHUCA

# sw_stream_gauges = ['406214', '406219', '406201', '406224', '406218', '406202', '406265']
# Removed 406214 as it represents Axe Creek which is not included in this model
# Remvoed 406224 as it represents Mount Pleasant Creek which is not included in the model.
# Remvoed 406219 as it represents the head in the reservoir which is not modelled.
# Added in 406203 as it represents stage just downstream of the weir,
#   NB not currently included in list of gauges online


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
              # 'CAMPASPE RIVER @ FEHRINGS LANE',
              'CAMPASPE RIVER @ ECHUCA']

inflow_gauges = ['MILLEWA CREEK @ NORTHERN HIGHWAY ECHUCA',
                 'CAMPASPE DR NO 5 @ OUTFALL',
                 'CAMPASPE DR NO 4 U/S NORTHERN HIGHWAY',
                 'AXE CREEK @ LONGLEA',
                 'AXE CREEK @ STRATHFIELDSAYE']

def find_layer(elev, col_vals):
    '''
    Function to find what layer a point lies within
    param elev: float, elevation of point of interest
    param col_vals: list of floats, layer elevations from top to bottom
    returns int
    '''
    for index, val in enumerate(col_vals):
        if elev > val:
            if index == 0:
                return index
            else:
                return index - 1
        #end if


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#******************************************************************************

SS_model.GISInterface.raster_reproject_by_grid(surface_raster_high_res,
                                               surface_raster_high_res[:-4] + '_reproj.tif',
                                               resample_method='min')

surface_raster_high_res = surface_raster_high_res[:-4] + '_reproj.tif'


SS_model.map_points_to_grid(river_gauges, feature_id='Site Name')
#SS_model.map_points_to_grid(river_gauges, feature_id='Site_ID')

Campaspe_river_gauges = SS_model.points_mapped['processed_river_sites_stage_clipped.shp']

filter_gauges = []
for riv_gauge in Campaspe_river_gauges:
    #if riv_gauge[1][0] in use_gauges:
    if str(riv_gauge[1][0]) in use_gauges:
        filter_gauges += [riv_gauge]


SS_model.create_river_dataframe('Campaspe', Campaspe_river_poly_file, surface_raster_high_res)

# Create reach data
river_seg = SS_model.river_mapping['Campaspe']
# Parameters are ordered from upstream to downstream
num_reaches = 4

SS_model.create_pilot_points('Campaspe', linear=True)
camp_pp = SS_model.pilot_points['Campaspe']
camp_pp.set_uniform_points(river_seg['rchlen'].sum(), num_reaches)

known_points = camp_pp.points

# Define split on river for which unique values will be given to props at 
# those points which will then be interpolated along the length of the river
#for reach in range(num_reaches):
# Setting up river bed hydraulic conductivity values
SS_model.parameters.create_model_parameter_set('kv_riv', 
                                           value=[1., 1., 0.01, 0.001], 
                                           num_parameters=num_reaches)
SS_model.parameters.parameter_options_set('kv_riv', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.01, 
                                      PARUBND=10.0, 
                                      PARGP='kv_riv', 
                                      SCALE=1, 
                                      OFFSET=0)
# Setting up river bed elevation correction parameter to account for 
# uncertainty in where bed lies relative to zero gauge
SS_model.parameters.create_model_parameter_set('beddep', 
                                           value=0.01, 
                                           num_parameters=num_reaches)
SS_model.parameters.parameter_options_set('beddep', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.001, 
                                      PARUBND=1.0, 
                                      PARGP='rivbed', 
                                      SCALE=1, 
                                      OFFSET=0)
# Setting up river bed roughness values
SS_model.parameters.create_model_parameter_set('mn_riv', 
                                           value=0.001, 
                                           num_parameters=num_reaches)
SS_model.parameters.parameter_options_set('mn_riv', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.001, 
                                      PARUBND=0.1, 
                                      PARGP='rough', 
                                      SCALE=1, 
                                      OFFSET=0)
# Setting up river width values
SS_model.parameters.create_model_parameter_set('rivwdth', 
                                           value=20.0, 
                                           num_parameters=num_reaches)
SS_model.parameters.parameter_options_set('rivwdth', 
                                      PARTRANS='fixed', 
                                      PARCHGLIM='factor', 
                                      PARLBND=4., 
                                      PARUBND=40., 
                                      PARGP='rivwdt', 
                                      SCALE=1, 
                                      OFFSET=0)
# Setting up riverbed thickness values
SS_model.parameters.create_model_parameter_set('bedthck', 
                                           value=0.10, 
                                           num_parameters=num_reaches)
SS_model.parameters.parameter_options_set('bedthck', 
                                      PARTRANS='fixed', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.01, 
                                      PARUBND=1., 
                                      PARGP='bedthk', 
                                      SCALE=1, 
                                      OFFSET=0)


strcond_val = [SS_model.parameters.param['kv_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
river_seg['strhc1'] = np.interp(river_seg['Cumulative Length'].tolist(), 
                                known_points, strcond_val)

strthick_val = [SS_model.parameters.param['bedthck{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
river_seg['strthick'] = np.interp(river_seg['Cumulative Length'].tolist(), 
                                  known_points, strthick_val)

strwidth_val = [SS_model.parameters.param['rivwdth{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
river_seg['width1'] = np.interp(river_seg['Cumulative Length'].tolist(), 
                                known_points, strwidth_val)

amalg_riv_points = []
for row in river_seg[['i', 'j']].iterrows():
    amalg_riv_points += [[row[1]['j'], row[1]['i']]]

# Sort out collocated stream reaches to avoid short circuiting:
river_seg['amalg_riv_points_tuple'] = river_seg['amalg_riv_points'].apply(lambda x: (x[0], x[1]))    
river_seg_group = river_seg.groupby(by='amalg_riv_points_tuple').count()
river_seg_group = river_seg_group[river_seg_group['amalg_riv_points'] > 1]

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# THINGS THAT NEED TO BE DONE FOR EACH COLUMN OF river_seg WHEN MERGING ROW
# BASED ON THE merge_group:
# For some testing
#
#  strtop = average weighted by rchlen
#  rchlen = sum()
#  amalg_riv_points = first entry
#  Cumulative length = last entry
#  strtop_raw = average weighted by rchlen
#  slope = average weighted by rchlen
#  k = first entry  
#  i = first entry
#  j = first entry
#  amalg_riv_points_collection = join lists of tuples into one
#  strhc1 = average weighted by rchlen
#  strthick = average weighted by rchlen
#  amalg_riv_points_tuple = first entry
#
#  1. Create new row based on above rules
#  2. Replace first row indexed in merge_group with new row
#  3. Delete all other rows indexed in merge_group
#
river_seg2 = river_seg.copy()

max_length = 3000.
merge_row = []
for ind in range(river_seg.shape[0]):
    if ind == 0:
        continue
    elif ind == river_seg.shape[0] - 1:
        prev = river_seg.iloc[ind - 1]    
        curr = river_seg.iloc[ind]    
    else:
        prev = river_seg.iloc[ind - 1]    
        curr = river_seg.iloc[ind]    
        nexx = river_seg.iloc[ind + 1]
        #def loc_tup(row):
        #    return (row['i'], row['j'])
        if prev['amalg_riv_points_tuple'] == nexx['amalg_riv_points_tuple']:
            if curr['rchlen'] < max_length:
                merge_row += [ind]
            
from operator import itemgetter
from itertools import groupby
merge_row_consec = []
for k, g in groupby(enumerate(merge_row), lambda (i,x):i-x):
    merge_row_consec.append(map(itemgetter(1), g))    

first_entry = lambda x: x[0]
last_entry = lambda x: x[-1]
flatten = lambda l: [item for sublist in l for item in sublist]

for merge_group in merge_row_consec:
    index_list = river_seg2.index.tolist()
    index_dict = {x:index for index, x in enumerate(index_list)}
    merge_group = [index_list[index_dict[merge_group[0]] - 1]] + merge_group
    merge_group = merge_group + [merge_group[-1] + 1] 
    river_seg_temp = river_seg2.loc[merge_group]
    rchlen_temp = river_seg_temp['rchlen']
    rchlen_sum = rchlen_temp.sum()
    rchlen_weights = rchlen_temp / rchlen_sum
    def weighted(col):
        return (col * rchlen_weights).sum()
    
    river_seg2.loc[merge_group[0], 'strtop'] = weighted(river_seg_temp['strtop']) 
    river_seg2.loc[merge_group[0], 'rchlen'] = rchlen_sum 
    river_seg2.set_value(merge_group[0], 'amalg_riv_points', first_entry(river_seg_temp['amalg_riv_points'].tolist()))
    river_seg2.loc[merge_group[0], 'Cumulative Length'] = last_entry(river_seg_temp['Cumulative Length'].tolist())     
    river_seg2.loc[merge_group[0], 'strtop_raw'] = weighted(river_seg_temp['strtop_raw']) 
    river_seg2.loc[merge_group[0], 'slope'] = weighted(river_seg_temp['slope']) 
    river_seg2.loc[merge_group[0], 'k'] = first_entry(river_seg_temp['k'].tolist()) 
    river_seg2.loc[merge_group[0], 'i'] = first_entry(river_seg_temp['i'].tolist()) 
    river_seg2.loc[merge_group[0], 'j'] = first_entry(river_seg_temp['j'].tolist()) 
    river_seg2.set_value(merge_group[0], 'amalg_riv_points_collection', flatten(river_seg_temp['amalg_riv_points_collection'])) 
    river_seg2.loc[merge_group[0], 'strhc1'] = weighted(river_seg_temp['strhc1']) 
    river_seg2.loc[merge_group[0], 'strthick'] = weighted(river_seg_temp['strthick']) 
    river_seg2.set_value(merge_group[0], 'amalg_riv_points_tuple', first_entry(river_seg_temp['amalg_riv_points_tuple'].tolist()))
    
    river_seg2.drop(merge_group[1:], inplace=True)

river_seg2.index = range(river_seg2.shape[0])

max_length = 500.
merge_row_too_short = []
for ind in range(river_seg2.shape[0]):
    if ind == 0:
        continue
    elif ind == river_seg2.shape[0] - 1:
        prev = river_seg2.iloc[ind - 1]    
        curr = river_seg2.iloc[ind]    
    else:
        prev = river_seg2.iloc[ind - 1]    
        curr = river_seg2.iloc[ind]    
        nexx = river_seg2.iloc[ind + 1]
        if prev['amalg_riv_points_tuple'] == nexx['amalg_riv_points_tuple']:
            pass
        else:
            if curr['rchlen'] < max_length:
                merge_row_too_short += [ind]

merge_row_too_short_consec = []
for k, g in groupby(enumerate(merge_row_too_short), lambda (i,x):i-x):
    merge_row_too_short_consec.append(map(itemgetter(1), g))    

for merge_group in merge_row_too_short_consec:
    index_list = river_seg2.index.tolist()
    index_dict = {x:index for index, x in enumerate(index_list)}
    merge_group = [index_list[index_dict[merge_group[0]] - 1]] + merge_group
    #merge_group = merge_group + [merge_group[-1] + 1] 
    river_seg_temp = river_seg2.loc[merge_group]
    rchlen_temp = river_seg_temp['rchlen']
    rchlen_sum = rchlen_temp.sum()
    rchlen_weights = rchlen_temp / rchlen_sum
    def weighted(col):
        return (col * rchlen_weights).sum()
    
    river_seg2.loc[merge_group[0], 'strtop'] = weighted(river_seg_temp['strtop']) 
    river_seg2.loc[merge_group[0], 'rchlen'] = rchlen_sum 
    river_seg2.set_value(merge_group[0], 'amalg_riv_points', first_entry(river_seg_temp['amalg_riv_points'].tolist()))
    river_seg2.loc[merge_group[0], 'Cumulative Length'] = last_entry(river_seg_temp['Cumulative Length'].tolist())     
    river_seg2.loc[merge_group[0], 'strtop_raw'] = weighted(river_seg_temp['strtop_raw']) 
    river_seg2.loc[merge_group[0], 'slope'] = weighted(river_seg_temp['slope']) 
    river_seg2.loc[merge_group[0], 'k'] = first_entry(river_seg_temp['k'].tolist()) 
    river_seg2.loc[merge_group[0], 'i'] = first_entry(river_seg_temp['i'].tolist()) 
    river_seg2.loc[merge_group[0], 'j'] = first_entry(river_seg_temp['j'].tolist()) 
    river_seg2.set_value(merge_group[0], 'amalg_riv_points_collection', flatten(river_seg_temp['amalg_riv_points_collection'])) 
    river_seg2.loc[merge_group[0], 'strhc1'] = weighted(river_seg_temp['strhc1']) 
    river_seg2.loc[merge_group[0], 'strthick'] = weighted(river_seg_temp['strthick']) 
    river_seg2.set_value(merge_group[0], 'amalg_riv_points_tuple', first_entry(river_seg_temp['amalg_riv_points_tuple'].tolist()))
    
    river_seg2.drop(merge_group[1], inplace=True)

river_seg = river_seg2
SS_model.river_mapping['Campaspe'] = river_seg
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

already_defined = []
old = []
for row in river_seg.iterrows():
    ind = row[0]
    row = row[1]
    new = row['amalg_riv_points']
    if new in old:
        already_defined += [ind]
    old += [new]
# end for

river_seg.loc[already_defined, 'strhc1'] = 0.0

new_k = []

for row in river_seg.iterrows():
    j_mesh = row[1]['i'] 
    i_mesh = row[1]['j']
    strtop = row[1]['strtop']
    strbot = row[1]['strtop'] - row[1]['strthick'] 
    new_k += [find_layer(strbot, SS_model.model_mesh3D[0][:, j_mesh, i_mesh])]
# end for

river_seg['k'] = new_k
       
# Remove any stream segments for which the elevation could not be mapped to a layer
river_seg['ireach'] = 1
river_seg['iseg'] = [x + 1 for x in range(river_seg.shape[0])]

                      
# Set up bed elevations based on the gauge zero levels:
gauge_points = [x for x in zip(Campaspe.Easting, Campaspe.Northing)]
river_gauge_seg = SS_model.get_closest_riv_segments('Campaspe', gauge_points)
river_seg.loc[:, 'bed_from_gauge'] = np.nan
river_seg.loc[:, 'stage_from_gauge'] = np.nan
river_seg.loc[:, 'gauge_id'] = 'none'

Campaspe['new_gauge'] = Campaspe[['Gauge Zero (Ahd)', 'Cease to flow level', 
                                  'Min value']].max(axis=1) 
Campaspe['seg_loc'] = river_gauge_seg         
Campaspe_gauge_zero = Campaspe[Campaspe['new_gauge'] > 10.]
# There are two values at the Campaspe weir, while it would be ideal to split the
# reach here it will cause problems for the segment
Campaspe_gauge_zero2 = Campaspe_gauge_zero[Campaspe_gauge_zero['Site Id'] != 406218]

Campaspe_stage = pd.merge(Campaspe, river_stage_data[1], on='Site Name', how='inner', suffixes=('','_r'))
Campaspe_stage = Campaspe_stage[[x for x in Campaspe_stage.columns if '_r' not in x]]

river_seg.loc[river_seg['iseg'].isin(Campaspe_gauge_zero2['seg_loc'].tolist()), 
              'bed_from_gauge'] = \
                  sorted(Campaspe_gauge_zero2['new_gauge'].tolist(), 
                         reverse=True)

river_seg.loc[river_seg['iseg'].isin(Campaspe_stage['seg_loc'].tolist()), 
              'gauge_id'] = \
                  Campaspe_stage.sort_values('Mean stage (m)', ascending=False)['Site ID'].tolist() 
                         

river_seg.loc[river_seg['iseg'].isin(Campaspe_stage['seg_loc'].tolist()), 
              'stage_from_gauge'] = \
                  sorted(Campaspe_stage['Mean stage (m)'].tolist(), 
                         reverse=True)

river_seg['bed_from_gauge'] = \
    river_seg.set_index(river_seg['Cumulative Length'])['bed_from_gauge']. \
        interpolate(method='values', limit_direction='both').tolist()

river_seg['stage_from_gauge'] = \
    river_seg.set_index(river_seg['Cumulative Length'])['stage_from_gauge']. \
        interpolate(method='values', limit_direction='both').tolist()


new_k = []
surface_layers = {}
bottom_layer = []

for row in river_seg.iterrows():
    j_mesh = row[1]['i'] 
    i_mesh = row[1]['j']
    strbot = row[1]['bed_from_gauge'] - row[1]['strthick']
    new_k += [find_layer(strbot, SS_model.model_mesh3D[0][:, j_mesh, i_mesh])]
    k = find_layer(strbot, SS_model.model_mesh3D[0][:, j_mesh, i_mesh])
    bottom_layer += [SS_model.model_mesh3D[0][k + 1, j_mesh, i_mesh]] 
    for layer in range(7):
        try:
            surface_layers[layer] += [SS_model.model_mesh3D[0][layer, j_mesh, 
                                                               i_mesh]]
        except:
            surface_layers[layer] = [SS_model.model_mesh3D[0][layer, j_mesh, 
                                                              i_mesh]]
        # end try
    # end for
    
for layer in range(7):
    river_seg["surf{}".format(layer)] = surface_layers[layer]
# end for

river_seg['k'] = new_k
river_seg['strtop'] = river_seg['bed_from_gauge']
river_seg['bottom_layer'] = bottom_layer
river_seg['stage'] = river_seg['stage_from_gauge']

river_seg.plot(x='Cumulative Length', y=['bed_from_gauge'] + \
                   ["surf{}".format(x) for x in range(7)])
river_seg.plot(x='Cumulative Length', y=['bed_from_gauge', 'bottom_layer'])

river_seg.plot(x='Cumulative Length', y=['bed_from_gauge', 'stage_from_gauge'], style='o')

# For stream reaches that didn't map properly to the mesh for z elevation we 
# can still include by setting to layer 0 with a bed hydraulic conductivity of 0
inds = np.where(river_seg['k'].isnull())[0]
river_seg.dropna(inplace=True)
          
river_seg['iseg'] = [x + 1 for x in range(river_seg.shape[0])]


SS_model.river_mapping['Campaspe'] = river_seg

###############################################################################

simple_river = []
for row in river_seg.iterrows():
    row = row[1]
    simple_river += [[row['k'], row['i'], row['j'], row['stage'], \
                      row['strhc1'] * row['rchlen'] * row['width1'], \
                      row['strtop']]]
# end for

riv = {0: simple_river}

if VERBOSE:
    print "************************************************************************"
    print " Creating Campaspe river boundary"

SSbounds.create_model_boundary_condition('Campaspe River', 'river', bc_static=True)
SSbounds.assign_boundary_array('Campaspe River', riv)




if VERBOSE:
    print "************************************************************************"
    print " Mapping Murray River to grid"

SS_model.parameters.create_model_parameter('rmstage', value=0.01)
SS_model.parameters.parameter_options('rmstage', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.001, 
                                      PARUBND=0.1, 
                                      PARGP='murr_riv', 
                                      SCALE=1, 
                                      OFFSET=0)
# Parameter to all shifting the location of the bed which is only estimated based on assumed depth below zero gauge
SS_model.parameters.create_model_parameter('rmbed', value=0.01)
SS_model.parameters.parameter_options('rmbed', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.001, 
                                      PARUBND=0.1, 
                                      PARGP='murr_riv', 
                                      SCALE=1, 
                                      OFFSET=0)
# Parameter for River Murray bed thickness
SS_model.parameters.create_model_parameter('rmbedthk', value=0.01)
SS_model.parameters.parameter_options('rmbedthk', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.001, 
                                      PARUBND=0.5, 
                                      PARGP='murr_riv', 
                                      SCALE=1, 
                                      OFFSET=0)
# Parameter for the vertical hydraulic conductivity of the River Murray
SS_model.parameters.create_model_parameter('kv_rm', value=5E-3)
SS_model.parameters.parameter_options('kv_rm', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=1E-8, 
                                      PARUBND=20, 
                                      PARGP='murr_riv', 
                                      SCALE=1, 
                                      OFFSET=0)
# Parameter for the width of the River Murray
SS_model.parameters.create_model_parameter('rmwdth', value=30)
SS_model.parameters.parameter_options('rmwdth', 
                                      PARTRANS='fixed', 
                                      PARCHGLIM='factor', 
                                      PARLBND=20, 
                                      PARUBND=50, 
                                      PARGP='murr_riv', 
                                      SCALE=1, 
                                      OFFSET=0)

SS_model.create_river_dataframe('Murray', Murray_river_poly_file, surface_raster_high_res)

mriver_seg = SS_model.river_mapping['Murray']
mriver_seg['strthick'] = SS_model.parameters.param['rmbedthk']['PARVAL1']

# Set up bed elevations based on the gauge zero levels:
Murray = Campaspe_relevant[Campaspe_relevant['Site Name'].str.contains("MURRAY RIVER")]
gauge_points = [x for x in zip(Murray.Easting, Murray.Northing)]
mriver_gauge_seg = SS_model.get_closest_riv_segments('Murray', gauge_points)
mriver_seg.loc[:, 'bed_from_gauge'] = np.nan
mriver_seg.loc[:, 'stage_from_gauge'] = np.nan

Murray.loc[:, 'new_gauge'] = Murray[['Gauge Zero (Ahd)', 'Cease to flow level', 'Min value']].max(axis=1) 
Murray.loc[:, 'seg_loc'] = mriver_gauge_seg         
Murray_gauge_zero = Murray[Murray['new_gauge'] > 10.]
mriver_seg['iseg'] = [x + 1 for x in range(mriver_seg.shape[0])]
#Murray_gauge_zero['Cumulative Length'] = mriver_seg.loc[Murray_gauge_zero['seg_loc'].tolist(), 'Cumulative Length'].tolist()

#Murray = pd.merge(Murray, river_stage_data[1], on='Site Name', how='inner', suffixes=('','_r'))
Murray_gauge_zero = pd.merge(Murray_gauge_zero, river_stage_data[1], on='Site Name', how='inner', suffixes=('','_r'))

Murray_gauge_zero = Murray_gauge_zero[[x for x in Murray_gauge_zero.columns if '_r' not in x]]

def values_from_gauge(column, reference):
    mriver_seg.loc[mriver_seg['iseg'].isin( \
        Murray_gauge_zero['seg_loc'].tolist()), column] = sorted( \
        Murray_gauge_zero[reference].tolist(), reverse=True)
    mriver_seg[column] = mriver_seg.set_index( \
        mriver_seg['Cumulative Length'])[column].interpolate( \
        method='values', limit_direction='both').tolist()
    mriver_seg[column].fillna(method='bfill', inplace=True)
# end values_from_gauge

values_to_edit = ['bed_from_gauge', 'stage_from_gauge']
references = ['new_gauge', 'Mean stage (m)']

for value, reference in zip(values_to_edit, references):
    values_from_gauge(value, reference)
# end for

new_k = []
active = []
surface_layers = {}
bottom_layer = []

for row in mriver_seg.iterrows():
    j_mesh = int(row[1]['i'])
    i_mesh = int(row[1]['j'])
    strtop = row[1]['bed_from_gauge']
    k = find_layer(strtop, SS_model.model_mesh3D[0][:, j_mesh, i_mesh])
    new_k += [k]
    bottom_layer += [SS_model.model_mesh3D[0][k+1, j_mesh, i_mesh]] 
    active += [SS_model.model_mesh3D[1][k, j_mesh, i_mesh]]
    for layer in range(7):
        try:
            surface_layers[layer] += [SS_model.model_mesh3D[0][layer, j_mesh, i_mesh]]
        except:
            surface_layers[layer] = [SS_model.model_mesh3D[0][layer, j_mesh, i_mesh]]
        # end try
    # end for
# end for

for layer in range(7):
    mriver_seg["surf{}".format(layer)] = surface_layers[layer]
# end for

mriver_seg['bottom_layer'] = bottom_layer

mriver_seg.plot(x='Cumulative Length', y=['bed_from_gauge'] + ["surf{}".format(x) for x in range(7)])

mriver_seg.plot(x='Cumulative Length', y=['bed_from_gauge', 'stage_from_gauge'])

mriver_seg['k'] = new_k
mriver_seg['active'] = active
      
# Remove any stream segments for which the elevation could not be mapped to a layer
mriver_seg[mriver_seg['active'] == -1] = np.nan
mriver_seg.dropna(inplace=True)
SS_model.river_mapping['Murray'] = mriver_seg

mriver_seg['strtop'] = mriver_seg['bed_from_gauge']                      
mriver_seg['strhc1'] = SS_model.parameters.param['kv_rm']['PARVAL1']                      
mriver_seg['width1'] = SS_model.parameters.param['rmwdth']['PARVAL1']
mriver_seg['stage'] = mriver_seg['stage_from_gauge'] + SS_model.parameters.param['rmstage']['PARVAL1']

# Avoid collisions with Campaspe River ...
def is_in_other_river(riv_df_testing, riv_df_other):
    riv_df_other_locs = riv_df_other['amalg_riv_points'].tolist()
    cell_used = []
    for row in riv_df_testing.iterrows():
        if row[1]['amalg_riv_points'] in riv_df_other_locs:
            cell_used += [0]
        else:
            cell_used += [1] 
        # end if
    # end for
    
    return cell_used
# end is_in_other_river0.

cells_overlapping = is_in_other_river(mriver_seg, river_seg)
mriver_seg['cell_used'] = cells_overlapping
mriver_seg[mriver_seg['cell_used'] == 0] = np.nan
mriver_seg.dropna(inplace=True)
SS_model.river_mapping['Murray'] = mriver_seg

mriver_seg['cell_loc_tuple'] = [(x[1]['k'], x[1]['i'], x[1]['j']) for x in mriver_seg.iterrows()]
mriver_seg = mriver_seg.groupby(by='cell_loc_tuple').mean()
mriver_seg.index = range(mriver_seg.shape[0])

simple_river = []

for row in mriver_seg.iterrows():
    row = row[1]
    simple_river += [[row['k'], row['i'], row['j'], row['stage'], \
                      row['strhc1'] * row['rchlen'] * row['width1'], \
                      row['strtop']]]
# end for

riv = {0: simple_river}

if VERBOSE:
    print "************************************************************************"
    print " Creating Murray River boundary"

SSbounds.create_model_boundary_condition('Murray River', 'river', bc_static=True)
SSbounds.assign_boundary_array('Murray River', riv)

if VERBOSE:
    print "************************************************************************"
    print " Setting up Murray River GHB boundary"

SS_model.parameters.create_model_parameter('mghb_stage', value=0.01)
SS_model.parameters.parameter_options('mghb_stage', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=-20.0, 
                                      PARUBND=50, 
                                      PARGP='ghb', 
                                      SCALE=1, 
                                      OFFSET=0)
SS_model.parameters.create_model_parameter('mghbk', value=10)
SS_model.parameters.parameter_options('mghbk', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=1E-8, 
                                      PARUBND=50, 
                                      PARGP='ghb', 
                                      SCALE=1, 
                                      OFFSET=0)

# First find which cells should make up the boundary based on the mapping 
# from the Murray river polyline to the grid

MurrayGHB = []
Active_MurrayGHB_cells = []
Murray_df_ind = []
checked = []

for mrow in mriver_seg.iterrows():
    ind = mrow[0]
    mrow = mrow[1]
    row = int(mrow['i'])
    col = int(mrow['j'])

    for lay in range(SS_model.model_mesh3D[1].shape[0]):    
        if [lay, row, col] in checked:
            continue
        # end if
        checked += [lay, row, col]
        if SS_model.model_mesh3D[1][0][row][col] == -1:
            continue
        # end if
        
        MurrayGHBstage = mrow['stage'] + \
                             SS_model.parameters.param['mghb_stage']['PARVAL1']

        if lay <= mrow['k']:
            continue
        # end if
        
        Murray_df_ind += [ind]        
        Active_MurrayGHB_cells += [[lay, row, col]]
    
    # end for
#end for

# Now make sure that no cells are being caught surrounded by other GHB cells to prevent short circuiting

def active_check(check_param, check_val, ref_cell, zone_cell, Active_MurrayGHB_cells):
    """
    Check target cell if it meets criteria to mark it as active

    :param check_param: int, target parameter to check
    :param check_val: int, value to check against, target parameter should not equal this
    :param ref_cell: list, reference to neighbouring cell ([layer, row, column])
    :param zone_cell: list, reference to HGU of cell with -1 indicating inactive cell 
    :param Active_MurrayGHB_cells: list, list of cell locations ([layer, row, column])

    :returns: bool, cell is active or not active
    """
    if (check_param != check_val) and \
            (zone_cell != -1) and \
            (ref_cell not in Active_MurrayGHB_cells):
        return True
    # end if
    return False
# end active_check()

Murray_df_ind2 = []
Final_MurrayGHB_cells = []
zone = SS_model.model_mesh3D[1]
shape = zone.shape
lay, row, col = 0, 1, 2

for index, ac in enumerate(Active_MurrayGHB_cells):
    # check if active GHB cell has any active non GHB cells N,E,S,W
    active_non_GHB = False
    acl, acr, acc = int(ac[lay]), int(ac[row]), int(ac[col])

    # Check north:
    ref_cell = [acl, acr + 1, acc]
    active_non_GHB = active_check(acr, 0, ref_cell, 
                                  zone[acl, acr + 1, acc], 
                                  Active_MurrayGHB_cells)

    # Check east:
    if not active_non_GHB:
        ref_cell = [acl, acr, acc + 1]
        active_non_GHB = active_check(acc, shape[col] - 1, ref_cell, 
                                      zone[acl, acr, acc + 1], 
                                      Active_MurrayGHB_cells)
    # end if

    # Check south:
    if not active_non_GHB:
        ref_cell = [acl, acr - 1, acc]
        active_non_GHB = active_check(acr, shape[row] - 1, ref_cell, 
                                      zone[acl, acr - 1, acc], 
                                      Active_MurrayGHB_cells)
    # end if

    # Check west:
    if not active_non_GHB:
        ref_cell = [acl, acr, acc - 1]
        active_non_GHB = active_check(acc, 0, ref_cell, 
                                      zone[acl, acr, acc - 1], 
                                      Active_MurrayGHB_cells)
    # end if
    
    if active_non_GHB:
        Final_MurrayGHB_cells += [ac]
        Murray_df_ind2 += [Murray_df_ind[index]]
    # end if
# end for

for index, MurrayGHB_cell in enumerate(Final_MurrayGHB_cells):

    lay = MurrayGHB_cell[0]
    row = MurrayGHB_cell[1]
    col = MurrayGHB_cell[2]
        
    MurrayGHBstage = mriver_seg['stage'].loc[Murray_df_ind2[index]] + \
                         SS_model.parameters.param['mghb_stage']['PARVAL1']
    dx = SS_model.gridHeight
    dz = SS_model.model_mesh3D[0][lay][row][col] - \
             SS_model.model_mesh3D[0][lay + 1][row][col]
    MGHBconductance = dx * dz * SS_model.parameters.param['mghbk']['PARVAL1']
    MurrayGHB += [[lay, row, col, MurrayGHBstage, MGHBconductance]]
# end for

ghb = {0: MurrayGHB}

if VERBOSE:
    print "************************************************************************"
    print " Creating GHB boundary"

SSbounds.create_model_boundary_condition('GHB', 'general head', bc_static=True)
SSbounds.assign_boundary_array('GHB', ghb)







if VERBOSE:
    print "************************************************************************"
    print " Mapping bores to grid "

SS_model.map_points_to_grid(bores_shpfile, feature_id='HydroCode')

bores_more_filter = []
bores_in_top_layer = []
for index, bores in enumerate(SS_model.points_mapped["NGIS_Bores_clipped.shp"]):
    row, col = bores[0][0], bores[0][1]
    for bore in bores[1]:
        try:
            # [bore_data_info["HydroCode"] == HydroCode]['depth']
            bore_depth = bore_data_info.loc[bore, 'depth']
        except Exception as e:
            # if bore in Ecology_bores:
            #    print 'Ecology bore not in info: ', bore
            #    sys.exit('Halting model build due to bore not being found')
            if bore in Policy_bores:
                print e
                print bore_data_info.head()
                print 'Policy bore not in info: ', bore
                # sys.exit('Halting model build due to bore not being found')
            continue
        # End try

        if bore_depth > mesh3D_0[0][row][col]:
            # print 'Bore can't be above surface!!!
            # if bore in Ecology_bores:
            #    print 'Ecology bore above surf: ', bore
            #    sys.exit('Halting model build due to bore not being mapped')
            if bore in Policy_bores:
                print 'Policy bore above surf: ', bore
                # sys.exit('Halting model build due to bore not being mapped')
            continue
        if bore_depth > mesh3D_0[1][row][col]:
            bores_in_top_layer += [bore]
            # if [row, col] in filter_gauge_loc:
            #    eco_gauge = [x[1] for x in filter_gauges if x[0] == [row, col]][0]
            #    print('candidate bore for ecology @ {0}: {1}'.format(eco_gauge, bore))
        if bore_depth <= mesh3D_0[-2][row][col]:
            # print 'Ignoring bores in bedrock!!!
            # if bore in Ecology_bores:
            #    print 'Ecology bore in  bedrock: ', bore
                # sys.exit('Halting model build due to bore not being mapped')
            #    continue
            if bore in Policy_bores:
                print('Policy bore in bedrock: {}'.format(bore))
                print('Bore depth is at: {}'.format(bore_depth))
                print('Bedrock top is at: {}'.format(mesh3D_0[-2][row][col]))
                print('Using above cell in Deep Lead for head by redefining bore_depth')
                final_bores.set_value(final_bores[final_bores['HydroCode'] == bore].index,
                                      'depth', mesh3D_0[-2][row][col] + 1.0)
                # sys.exit('Halting model build due to bore not being mapped')
        bores_more_filter += [bore]
    # End for
# End for

if VERBOSE:
    print 'Final bores within aquifers: ', len(bores_more_filter)

final_bores = final_bores[final_bores["HydroCode"].isin(bores_more_filter)]

# Now find ecology bores by distance to stream gauges of interest and that have
# mapped to the model domain and have sufficient data

bore_points = [[final_bores.loc[x, "Easting"], final_bores.loc[x, "Northing"]]
               for x in final_bores.index]

bore_points3D = final_bores[["HydroCode", "Easting", "Northing", "depth"]]
bore_points3D = bore_points3D.set_index("HydroCode")

bores_obs_time_series = final_bores[["HydroCode", "mean level"]]

# Modify into standard format for the GWModelBuilder class
bores_obs_time_series = bores_obs_time_series.rename(columns={'HydroCode': 'name',
                                                              'mean level': 'value'})

SS_model.observations.set_as_observations('head', bores_obs_time_series, bore_points3D,
                                          domain='porous', obs_type='head', units='mAHD')

bores_in_layers = SS_model.map_points_to_raster_layers(bore_points, final_bores["depth"].tolist(),
                                                       hu_raster_files_reproj)

# Map bores to layers to create initial head maps for different hydrogeological units
interp_heads = {}

for i in xrange(len(hu_raster_files_reproj) / 2):
    bores_layer = np.array(bore_points)[np.array(bores_in_layers[i])]
    print 'Creating head map for: ', hu_raster_files[2 * i]
    if bores_layer.shape[0] < 4:
        interp_heads[hu_raster_files[2 * i]] = np.full(mesh3D_1.shape[1:], np.NaN)
    else:
        bores_head_layer = np.array(
            final_bores["mean level"].tolist()
        )[np.array(bores_in_layers[i])]
        unique_bores = np.unique(bores_layer)

        b = np.ascontiguousarray(bores_layer).view(np.dtype(
            (np.void, bores_layer.dtype.itemsize *
             bores_layer.shape[1])))
        _, idx = np.unique(b, return_index=True)

        unique_bores = bores_layer[idx]

        interp_heads[hu_raster_files[2 * i]] = SS_model.interpolate_points2mesh(
            bores_layer,
            bores_head_layer,
            use='griddata',
            method='linear')

# initial_heads_SS = np.full(mesh3D_1.shape, 0.)
#
# for i in xrange(len(hu_raster_files_reproj) / 2):
#     initial_heads_SS[i] = mesh3D_1[0] - 10.0

initial_heads_SS = np.full(mesh3D_1.shape, 400.)

# interp_heads[hu_raster_files[0]])
SS_model.initial_conditions.set_as_initial_condition("Head", initial_heads_SS)

if VERBOSE:
    print "************************************************************************"
    print " Mapping pumping wells to grid "

SS_model.map_points_to_grid(pumps_points, feature_id='OLD ID')

SSparams.create_model_parameter('pump_use', value=0.6)
SSparams.parameter_options('pump_use',
                           PARTRANS='log',
                           PARCHGLIM='factor',
                           PARLBND=0.2,
                           PARUBND=1.,
                           PARGP='pumping',
                           SCALE=1,
                           OFFSET=0)

# Convert pumping_data to time series

# Existing data is only for 10 years from 2005 to 2015
pump_date_index = pd.date_range(start=datetime.datetime(2014, 07, 01),
                                end=datetime.datetime(2015, 06, 30), freq='AS-JUL')

# Normalise pumping data, i.e. find individual pumps contribution to
# total volume of water that is being pumped
total_pumping_rate = pumping_data['Use 2014/15'].mean()

wel = {}
pump_shallow = []  # Shallow (if <25m) or Deep (>= 25m)
for pump_cell in SS_model.points_mapped['pumping wells_clipped.shp']:
    row = pump_cell[0][0]
    col = pump_cell[0][1]
    layers = [0]
    for pump in pump_cell[1]:
        if pumping_data.loc[pump, 'Top screen depth (m)'] == 0.:
            # print 'No data to place pump at depth ... ignoring ', pump
            continue
        # End if

        pump_depth = mesh3D_0[0][row][col] - pumping_data.loc[
            pump,
            'Top screen depth (m)']
        active = False
        for i in xrange(mesh3D_0.shape[0] - 1):
            if pump_depth < mesh3D_0[i][row][col] and \
               pump_depth > mesh3D_0[i + 1][row][col]:
                active_layer = i
                active = True
                break
            # End if
        # End for

        if active is False:
            # print 'Well not placed: ', pump
            continue
        # Get top of screen layer and calculate length of screen in layer

        # Specify if pump is shallow
        if pump_depth < 25:
            pump_shallow += [True]
        else:
            pump_shallow += [False]
        # End if

        p14_15 = pumping_data.loc[pump, 'Use 2014/15'] / 365. * 1000.
        pump_rates = [p14_15]
        pumping_data_ts = pd.DataFrame(pump_rates, columns=[pump], index=pump_date_index)
        pump_install = pumping_data.loc[pump, 'Construction date']

        if isinstance(pump_install, datetime.time):
            pump_install = datetime.date(1950, 01, 01)
        # End if

        pump_date_index2 = pd.date_range(start=pump_install,
                                         end=datetime.datetime(2004, 06, 30),
                                         freq='AS-JUL')

        # Assume historical pumping is a percentage of lowest non-zero use for well
        non_zero_pumping = [x for x in pump_rates if x > 0.]
        if non_zero_pumping == []:
            pumping_rate_old = 0.
        else:
            pumping_rate_old = np.min(non_zero_pumping)
        # End if

        old_pumping_ts = pd.DataFrame(index=pump_date_index2)
        old_pumping_ts[pump] = pumping_rate_old * SS_model.parameters.param['pump_use']['PARVAL1']

        # Merge the old and measured data
        pumping_data_ts = pd.concat([pumping_data_ts, old_pumping_ts])

        # Now let's resample the data on a monthly basis, and we will take the mean
        pumping_data_ts = pumping_data_ts.resample(SS_model.model_time.t['time_step']).mean()

        # Let's also get rid of NaN data and replace with backfilling
        pumping_data_ts = pumping_data_ts.fillna(method='bfill')

        # Let's only consider times in our date range though
        s_time, e_time = SS_model.model_time.t['start_time'], SS_model.model_time.t['end_time']
        date_index = pd.date_range(start=s_time,
                                   end=e_time,
                                   freq=SS_model.model_time.t['time_step'])
        pumping_data_ts = pumping_data_ts.reindex(date_index)
        pumping_data_ts = pumping_data_ts.ix[s_time:e_time]
        pumping_data_ts = pumping_data_ts.fillna(0.0)

        # Now fill in the well dictionary with the values of pumping at relevant stress periods
        # where Q is not 0.0
        for index, timestep in enumerate(pumping_data_ts.iterrows()):
            if timestep[1][pump] == 0:
                continue
            # End try

            try:
                wel[index] += [[active_layer, row, col, -timestep[1][pump] / total_pumping_rate]]
            except:
                wel[index] = [[active_layer, row, col, -timestep[1][pump] / total_pumping_rate]]
            # End try
        # End for
    # End for
# End for

if VERBOSE:
    print "************************************************************************"
    print " Creating pumping boundary "

SSbounds.create_model_boundary_condition('licenced_wells', 'wells', bc_static=True)
SSbounds.assign_boundary_array('licenced_wells', wel)


if VERBOSE:
    print "************************************************************************"
    print " Collate observations"

SS_model.map_obs_loc2mesh3D(method='nearest')

obs_active_bores = bores_obs_time_series[bores_obs_time_series['zone'] != 'null']['name']
obs_active_bores = obs_active_bores[obs_active_bores.isin(bores_in_top_layer)].tolist()
obs_filter_bores = bore_points3D[bore_points3D.index.isin(obs_active_bores)]
obs_bores_list = zip(obs_filter_bores['Easting'], obs_filter_bores['Northing'])

stream_active = river_flow_data[0][river_flow_data[0]['Site ID'].isin([int(x) for x in Stream_gauges])]
stream_gauges_list = zip(stream_active['Easting'], stream_active['Northing'])

closest_bores_active = SS_model.find_closest_points_between_two_lists(obs_bores_list, stream_gauges_list)

ecol_bores = []
for ind in closest_bores_active:
    ecol_bores += [obs_filter_bores.index.tolist()[ind]]

ecol_bores_df = obs_filter_bores[obs_filter_bores.index.isin(ecol_bores)]

# ecology_found = [x for x in final_bores["HydroCode"] if x in Ecology_bores]
policy_found = [x for x in final_bores["HydroCode"] if x in Policy_bores]

# read in model link
link_file = os.path.join(data_folder, "model_linking.csv")
with open(link_file, 'r') as model_link:
    tmp = model_link.readlines()
    for i, line in enumerate(tmp):
        if "Ecology" in line:
            tmp[i] = "Ecology: {}\n".format(", ".join(ecol_bores))
        # End if

        if "Policy" in line:
            tmp[i] = "Policy: {}\n".format(", ".join(policy_found))
    # End for
# End with

# Write out ecology bore hydrocodes
with open(os.path.join(data_folder, "model_linking.csv"), 'w') as model_link:
    model_link.writelines(tmp)
# End with


# Setup the outputs for head based on location of stream gauges
# SS_model.observations.set_as_observations('head_stream_gauge', str_time_series,
#                                          new_riv_cells, domain='surface',
#                                          obs_type='head', units='m^3/d',
#                                          weights=0.0, real=False)


# Visuals checks on getting nearest mapped bore from top layer for the ecology part:
if VERBOSE:
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.scatter(obs_filter_bores['Easting'], obs_filter_bores['Northing'], label='bores')
    ax.scatter(stream_active['Easting'], stream_active[
               'Northing'], color='r', label='stream gauges')
    for idx in stream_active['Site ID']:
        # Annotate the closest gauge
        x_y = stream_active.loc[stream_active['Site ID'] == idx, ["Easting", "Northing"]]
        ax.annotate(idx, xy=x_y.values[0].tolist(), xycoords='data', xytext=(1.0, 0.8), textcoords='offset points')

    ax.scatter(ecol_bores_df['Easting'], ecol_bores_df['Northing'],
               marker='+', color='orange', label='closest bores')
    for idx in ecol_bores_df.index:
        # Annotate closest bore
        x_y = ecol_bores_df.loc[ecol_bores_df.index == idx, ["Easting", "Northing"]]
        ax.annotate(idx, xy=x_y.values[0].tolist(), xycoords='data', xytext=(-120.0, 50.0),
                    textcoords='offset points',
                    arrowprops=dict(facecolor='orange', shrink=0.05))

    # plt.legend()
    # fig.suptitle('Finding nearest bores to stream gauges')
    # plt.xlabel('Easting')
    # plt.ylabel('Northing')
    # plt.axis('equal')
    #
    # plt.show()

    # set_ylabel('Northing')

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # c = 'rgb'
    # for index, bore in enumerate(ecol_bores):
    #     bore_data_levels[bore_data_levels['HydroCode'] == bore][
    #         ['bore_date', 'result']].plot(x='bore_date', ax=ax, color=c[index])
    # # NOTE that the bores data is not going all the way to 2015, although bore filtering could include only those bores
    # # which have data that is recent .... can do this later!
    # # It is interesting to note that the distance can be quite far from gauge to bore
    # # Perhaps the restriction to top layer bores could be relaxed somewhat.
    #
    # plt.show()

if VERBOSE:
    print "************************************************************************"
    print " Mapping farm areas to grid"

SS_model.map_polygon_to_grid(farms_poly, feature_name="ZoneID")


print "************************************************************************"
print " Package up groundwater model builder object"

SS_model.package_model()

print "Packaged into {}".format(SS_model.out_data_folder_grid)
