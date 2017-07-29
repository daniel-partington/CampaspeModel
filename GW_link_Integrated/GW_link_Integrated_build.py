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
sw_data_path = p_j(temp_data_path, "Campaspe_data/SW/All_streamflow_Campaspe_catchment")

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
    bore_data_levels, bore_data_info = getBoreData.getBoreData(base_path=SS_model.campaspe_data)
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
if os.path.exists(riv_flow_path + '.h5'):
    river_flow_data = SS_model.load_dataframe(riv_flow_path + '.h5')
else:
    river_flow_data = processRiverStations.getFlow(path=sw_data_path)
    SS_model.save_dataframe(riv_flow_path, river_flow_data)
# End if

river_stage_file = "river_stage_processed"
riv_stage_path = p_j(SS_model.out_data_folder, river_stage_file)

# Check if this data has been processed and if not process it
if os.path.exists(riv_stage_path + '.h5'):
    river_stage_data = SS_model.load_dataframe(riv_stage_path + '.h5')
else:
    river_stage_data = processRiverStations.getStage(path=sw_data_path)
    SS_model.save_dataframe(riv_stage_path, river_stage_data)

river_gauges = SS_model.read_points_data(p_j(sw_data_path,
                                             r"processed_river_sites_stage.shp"))

if VERBOSE:
    print "************************************************************************"
    print "Load in the river shapefiles"

Campaspe_river_poly = SS_model.read_poly("Campaspe_Riv.shp", path=river_path)
Murray_river_poly = SS_model.read_poly("River_Murray.shp", path=river_path)


if VERBOSE:
    print "************************************************************************"
    print "Load in the farms shapefile"

farms_path = p_j(temp_data_path, "Campaspe_data", "SW", "Farm")
farms_poly = SS_model.read_poly("farm_v1_prj.shp", path=farms_path, poly_type='polygon')


if VERBOSE:
    print "************************************************************************"
    print "Load in the shapefiles defining groundwater boundaries"

WGWbound_poly = SS_model.read_poly("western_head.shp", path=model_build_input_path)
EGWbound_poly = SS_model.read_poly("eastern_head.shp", path=model_build_input_path)


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

# SS_model.map_points_to_grid(river_gauges, feature_id='Site_Name')
SS_model.map_points_to_grid(river_gauges, feature_id='Site_ID')

Campaspe_river_gauges = SS_model.points_mapped['processed_river_sites_stage_clipped.shp']

filter_gauges = []
for riv_gauge in Campaspe_river_gauges:
    if str(riv_gauge[1][0]) in Stream_gauges:
        filter_gauges += [riv_gauge]
    # End if
# End for

SS_model.map_polyline_to_grid(Campaspe_river_poly)
SSparams.create_model_parameter('bed_depress', value=0.01)
SSparams.parameter_options('bed_depress',
                           PARTRANS='log',
                           PARCHGLIM='factor',
                           PARLBND=0.001,
                           PARUBND=0.1,
                           PARGP='spec_stor',
                           SCALE=1,
                           OFFSET=0)
SSparams.create_model_parameter('Kv_riv', value=5E-3)
SSparams.parameter_options('Kv_riv',
                           PARTRANS='log',
                           PARCHGLIM='factor',
                           PARLBND=1E-8,
                           PARUBND=20,
                           PARGP='spec_stor',
                           SCALE=1,
                           OFFSET=0)

simple_river = []
riv_width_avg = 10.0  # m
riv_bed_thickness = 0.10  # m

# Map river from high to low
new_riv = SS_model.polyline_mapped['Campaspe_Riv_model.shp']
for index, riv_cell in enumerate(SS_model.polyline_mapped['Campaspe_Riv_model.shp']):
    row = riv_cell[0][0]
    col = riv_cell[0][1]
    new_riv[index] += [mesh3D_0[0][row][col]]

new_riv = sorted(new_riv, key=lambda x: (x[0][1]), reverse=False)
new_riv = sorted(new_riv, key=lambda x: (x[0][0]), reverse=True)

new_riv_len = (len(new_riv))
stages = np.full(new_riv_len, np.nan, dtype=np.float64)
beds = np.full(new_riv_len, np.nan, dtype=np.float64)

# Identify cells that correspond to river gauges
riv_gauge_logical = np.full(new_riv_len, False, dtype=np.bool)

# To account for fact that river shapefile and gauges shapefile are not perfect
# we get the closest river cell to the gauge cell


def closest_node(node, nodes):
    nodes = np.asarray(nodes)
    dist_2 = np.sum((nodes - node)**2, axis=1)
    return np.argmin(dist_2)


# Define river gauges at start of river cell
new_riv_cells = [x[0] for x in new_riv]
filter_gauge_loc = [new_riv_cells[x] for x in
                    [closest_node(x[0], new_riv_cells) for x in filter_gauges]]

# Set up the gauges as observations

for index, riv in enumerate(new_riv):
    # Create logical array to identify those which are gauges and those which are not
    if riv[0] in filter_gauge_loc:
        riv_gauge_logical[index] = True
        gauge_ind = [i for i, x in enumerate(filter_gauge_loc) if x == riv[0]]
        stages[index] = river_stage_data["Mean stage (m)"].loc[river_stage_data["Site ID"] ==
                                                               filter_gauges[gauge_ind[0]][1][0]]
        beds[index] = stages[index] - 1.0

    # Add chainage to new_riv array:
    if index == 0:
        new_riv[index] += [0.0]
    else:
        new_riv[index] += [new_riv[index - 1][3] + new_riv[index - 1][1]]
    # End if
# End for

# River x in terms of chainage:
river_x = np.array([x[3] for x in new_riv])
river_x_unknown = river_x[~riv_gauge_logical]
river_x_known = river_x[riv_gauge_logical]

# Now interpolate know values of stage and bed to unknown river locations:
stages[~riv_gauge_logical] = np.interp(river_x_unknown, river_x_known, stages[riv_gauge_logical])
beds[~riv_gauge_logical] = np.interp(river_x_unknown, river_x_known, beds[riv_gauge_logical])

# Create observations for stage or discharge at those locations

# TODO: GENERALISE THE TWO LOOPS BELOW THAT DOES THE SAME THING FOR CAMPASPE AND MURRAY
# locations = ['Campaspe_Riv_model.shp', 'River_Murray_model.shp']

kv_riv_parval1 = SS_model.parameters.param['Kv_riv']['PARVAL1']
for index, riv_cell in enumerate(SS_model.polyline_mapped['Campaspe_Riv_model.shp']):
    row = riv_cell[0][0]
    col = riv_cell[0][1]

    if mesh3D_1[0][row][col] == -1:
        continue
    # End if

    stage = stages[index]
    bed = beds[index]
    cond = (riv_cell[1] * riv_width_avg * kv_riv_parval1 / riv_bed_thickness)
    simple_river += [[0, row, col, stage, cond, bed]]

riv = {0: simple_river}

if VERBOSE:
    print "************************************************************************"
    print " Creating Campaspe river boundary"

SSbounds.create_model_boundary_condition('Campaspe River', 'river', bc_static=True)
SSbounds.assign_boundary_array('Campaspe River', riv)


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
    print " Mapping Murray River to grid"

SS_model.map_polyline_to_grid(Murray_river_poly)

SSparams.create_model_parameter('RMstage', value=0.01)
SSparams.parameter_options('RMstage',
                           PARTRANS='log',
                           PARCHGLIM='factor',
                           PARLBND=0.001,
                           PARUBND=0.1,
                           PARGP='spec_stor',
                           SCALE=1,
                           OFFSET=0)
SSparams.create_model_parameter('Kv_RM', value=5E-3)
SSparams.parameter_options('Kv_RM',
                           PARTRANS='log',
                           PARCHGLIM='factor',
                           PARLBND=1E-8,
                           PARUBND=20,
                           PARGP='spec_stor',
                           SCALE=1,
                           OFFSET=0)

simple_river = []
riv_width_avg = 10.0  # m
riv_bed_thickness = 0.10  # m

poly_mapped_murray_model = SS_model.polyline_mapped['River_Murray_model.shp']
Kv_RM_parval1 = SS_model.parameters.param['Kv_RM']['PARVAL1']
RMstage_parval1 = SS_model.parameters.param['RMstage']['PARVAL1']
for riv_cell in poly_mapped_murray_model:
    row = riv_cell[0][0]
    col = riv_cell[0][1]

    if mesh3D_1[0][row][col] == -1:
        continue
    # End if

    stage = mesh3D_0[0][row][col]
    bed = mesh3D_0[0][row][col] - RMstage_parval1
    cond = (riv_cell[1] * riv_width_avg * Kv_RM_parval1 / riv_bed_thickness)
    simple_river += [[0, row, col, stage, cond, bed]]

riv = {0: simple_river}

if VERBOSE:
    print "************************************************************************"
    print " Creating Murray River boundary"

SSbounds.create_model_boundary_condition('Murray River', 'river', bc_static=True)
SSbounds.assign_boundary_array('Murray River', riv)

if VERBOSE:
    print "************************************************************************"
    print " Setting up Murray River GHB boundary"

SSparams.create_model_parameter('MGHB_stage', value=0.01)
SSparams.parameter_options('MGHB_stage',
                           PARTRANS='log',
                           PARCHGLIM='factor',
                           PARLBND=-20.0,
                           PARUBND=50,
                           PARGP='ghb',
                           SCALE=1,
                           OFFSET=0)
SSparams.create_model_parameter('MGHBcond', value=5E-3)
SSparams.parameter_options('MGHBcond',
                           PARTRANS='log',
                           PARCHGLIM='factor',
                           PARLBND=1E-8,
                           PARUBND=50,
                           PARGP='ghb',
                           SCALE=1,
                           OFFSET=0)


MurrayGHB = []
dx = SS_model.gridHeight
MGHB_stage_parval1 = SSparams.param['MGHB_stage']['PARVAL1']
MGHBcond_parval1 = SSparams.param['MGHBcond']['PARVAL1']
for MurrayGHB_cell in poly_mapped_murray_model:
    row = MurrayGHB_cell[0][0]
    col = MurrayGHB_cell[0][1]

    for lay in xrange(mesh3D_1.shape[0]):
        if mesh3D_1[0][row][col] == -1:
            continue

        MurrayGHBstage = (mesh3D_0[0][row][col] +
                          MGHB_stage_parval1)
        if MurrayGHBstage < mesh3D_0[0][row][col]:
            continue
        # End if

        dz = mesh3D_0[lay][row][col] - mesh3D_0[lay + 1][row][col]
        MGHBconductance = dx * dz * MGHBcond_parval1
        MurrayGHB += [[lay, row, col, MurrayGHBstage, MGHBconductance]]
    # End for
# End for

ghb = {0: MurrayGHB}

if VERBOSE:
    print "************************************************************************"
    print " Creating GHB boundary"

SSbounds.create_model_boundary_condition('GHB', 'general head', bc_static=True)
SSbounds.assign_boundary_array('GHB', ghb)

if VERBOSE:
    print "************************************************************************"
    print " Collate observations"

SS_model.map_obs_loc2mesh3D(method='nearest')

obs_active_bores = bores_obs_time_series[bores_obs_time_series['zone'] != 'null']['name']
obs_active_bores = obs_active_bores[obs_active_bores.isin(bores_in_top_layer)].tolist()
obs_filter_bores = bore_points3D[bore_points3D.index.isin(obs_active_bores)]
obs_bores_list = zip(obs_filter_bores['Easting'], obs_filter_bores['Northing'])

stream_active = river_flow_data[river_flow_data['Site ID'].isin([int(x) for x in Stream_gauges])]
stream_gauges_list = zip(stream_active['Easting'], stream_active['Northing'])

closest_bores_active = SS_model.find_closest_points_between_two_lists(
    obs_bores_list, stream_gauges_list)

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
