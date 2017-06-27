import os

from osgeo import osr
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt


from HydroModelBuilder.GWModelBuilder import GWModelBuilder
from HydroModelBuilder.GISInterface.GDALInterface.GDALInterface import GDALInterface
from CampaspeModel.CustomScripts import processWeatherStations 
from CampaspeModel.CustomScripts import processRiverStations
from CampaspeModel.CustomScripts import readHydrogeologicalProperties

# Define basic model parameters:
Proj_CS = osr.SpatialReference()
Proj_CS.ImportFromEPSG(28355) # 28355 is the code for gda94 mga zone 55; http://spatialreference.org/ref/epsg/gda94-mga-zone-55/

Interface = GDALInterface()
Interface.projected_coordinate_system = Proj_CS 
Interface.pcs_EPSG = "EPSG:28355"

SS_model = GWModelBuilder(name="01_steady_state", 
                          data_folder=r"C:\Workspace\part0075\MDB modelling\testbox\input_data\\",
                          model_data_folder=r"C:\Workspace\part0075\MDB modelling\testbox\00_Campaspe_Cascade\01_steady_state\\",
                          out_data_folder=r"C:\Workspace\part0075\MDB modelling\testbox\data_build\\",
                          GISInterface=Interface,
                          model_type='Modflow',
                          mesh_type='structured')

# Cleanup
#SS_model.flush()

# Set the model boundary using a polygon shapefile:
print "************************************************************************"
print " Setting model boundary "

SS_model.set_model_boundary_from_polygon_shapefile("GW_model_area.shp", 
                                                   shapefile_path=SS_model.data_folder)

# Set data boundary for model data
print "************************************************************************"
print " Setting spatial data boundary "
SS_model.set_data_boundary_from_polygon_shapefile(SS_model.boundary_poly_file, 
                                                  shapefile_path=SS_model.data_folder,
                                                  buffer_dist=20000)

# Setup recharge:
# ... read in climate data using Custom_Scripts
weather_stations = ['Kyneton',  'Elmore', 'Rochester', 'Echuca']
print "************************************************************************"
print " Executing custom script: processWeatherStations "

rain_info_file = "rain_processed"
# Check if this data has been processed and if not process it
if os.path.exists(SS_model.out_data_folder + rain_info_file + '.h5'):
    long_term_historic_rainfall = SS_model.load_dataframe(SS_model.out_data_folder + rain_info_file + '.h5')
else:
    long_term_historic_rainfall = processWeatherStations.processWeatherStations(weather_stations, path=r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Climate\\")
    SS_model.save_dataframe(SS_model.out_data_folder + rain_info_file, long_term_historic_rainfall)

rain_gauges = SS_model.read_points_data(r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Climate\Rain_gauges.shp")

print "************************************************************************"
print " Executing custom script: readHydrogeologicalProperties "

file_location = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\GW\Aquifer properties\Hydrogeologic_variables.xlsx"
HGU_props = readHydrogeologicalProperties.getHGUproperties(file_location)

print "************************************************************************"
print " Executing custom script: processRiverStations "


river_flow_file = "river_flow_processed"
river_stage_file = "river_stage_processed"
river_ec_file = "river_ec_processed"
river_data_folder = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\SW\All_streamflow_Campaspe_catchment\Updated"
river_flow_file = os.path.join(SS_model.out_data_folder, river_flow_file)
river_stage_file = os.path.join(SS_model.out_data_folder, river_stage_file)
river_ec_file = os.path.join(SS_model.out_data_folder, river_ec_file)

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

from geopandas import GeoDataFrame
from shapely.geometry import Point

geometry = [Point(xy) for xy in zip(Campaspe.Easting, Campaspe.Northing)]
Campaspe_geo = pd.DataFrame.copy(Campaspe)
Campaspe_geo.drop(['Northing', 'Easting'], axis=1)
crs = {'init': Interface.pcs_EPSG}
Geo_df = GeoDataFrame(Campaspe_geo, crs=crs, geometry=geometry)
site_shapefile = os.path.join(river_data_folder, 'processed_river_sites_stage.shp') 
Geo_df.to_file(site_shapefile)
sites=Campaspe_relevant['Site Id'].tolist()

Campaspe_gauge_zero = Campaspe[(Campaspe['Gauge Zero (Ahd)'] > 0.0) | \
                               (Campaspe['Cease to flow level'] > 10.0) | \
                               (Campaspe['Min value'] > 10.0)]

# Check if this data has been processed and if not process it
if os.path.exists(river_flow_file + '.pkl'):
    river_flow_data = SS_model.load_obj(river_flow_file + '.pkl')
else:
    river_flow_data = processRiverStations.getFlow(path=river_data_folder, sites=sites)
    SS_model.save_obj(river_flow_data, river_flow_file)

# Check if this data has been processed and if not process it
if os.path.exists(river_stage_file + '.pkl'):
    river_stage_data = SS_model.load_obj(river_stage_file + '.pkl')
else:
    river_stage_data = processRiverStations.getStage(path=river_data_folder, sites=sites)
    SS_model.save_obj(river_stage_data, river_stage_file)

# Check if this data has been processed and if not process it
if os.path.exists(river_ec_file + '.pkl'):
    river_ec_data = SS_model.load_obj(river_ec_file + '.pkl')
else:
    river_ec_data = processRiverStations.getEC(path=river_data_folder, sites=sites)
    SS_model.save_obj(river_ec_data, river_ec_file)

river_gauges = SS_model.read_points_data(site_shapefile)

print "************************************************************************"
print "Load in the river shapefiles"
Campaspe_river_poly_file = r"C:\Workspace\part0075\MDB modelling\testbox\input_data\Waterways\Campaspe_Riv.shp"
Campaspe_river_poly = SS_model.read_poly("Campaspe_Riv.shp", path=r"C:\Workspace\part0075\MDB modelling\testbox\input_data\Waterways") 
Murray_river_poly_file = r"C:\Workspace\part0075\MDB modelling\testbox\input_data\Waterways\River_Murray.shp" 
Murray_river_poly = SS_model.read_poly("River_Murray.shp", path=r"C:\Workspace\part0075\MDB modelling\testbox\input_data\Waterways") 


# This is disgusting, but works for now ... needs cleaning up through first testing
# if raster is in the right projection and if not returning the name of the new
# reprojected file
#surface_raster_high_res = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Surface_DEM_Geoscience_Australia\Camp_1sHE_2026754\Camp_1sHE.tif"
surface_raster_high_res = r"C:\Workspace\part0075\MDB modelling\ESRI_GRID_raw\ESRI_GRID\sur_1t"

recharge_zones = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Linking_recharge\Zones_24.tif"

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
print '## Mesh specific model building '
print '########################################################################'
print '########################################################################'

# Define the grid width and grid height for the model mesh which is stored as a multipolygon shapefile GDAL object
print "************************************************************************"
print " Defining structured mesh"
resolution = 5000
SS_model.define_structured_mesh(resolution, resolution)

# Read in hydrostratigraphic raster info for layer elevations:
hu_raster_path = r"C:\Workspace\part0075\MDB modelling\ESRI_GRID_raw\ESRI_GRID"

# Build basement file ... only need to do this once as it is time consuming so commented out for future runs
# SS_model.create_basement_bottom(hu_raster_path, "sur_1t", "bse_1t", "bse_2b", hu_raster_path)

hu_raster_files = ["qa_1t", "qa_2b", "utb_1t", "utb_2b", "utqa_1t", "utqa_2b", 
                   "utam_1t", "utam_2b", "utaf_1t", "utaf_2b", "lta_1t", 
                   "lta_2b", "bse_1t", "bse_2b.tif"]

# This loads in the raster files and transforms them into the correct coordinate
# sytstem.
SS_model.read_rasters(hu_raster_files, path=hu_raster_path)


#hu_raster_files_reproj = [x + "_reproj.bil" for x in hu_raster_files]
hu_raster_files_reproj = [x + "_reproj.tif" for x in hu_raster_files]

# Map HGU's to grid
print "************************************************************************"
print " Mapping rasters to grid "

hu_gridded_rasters = SS_model.map_rasters_to_grid(hu_raster_files, 
                                                  hu_raster_path)

# Build 3D grid
model_grid_raster_files = [x + "_model_grid.tif" for x in hu_raster_files]

# First two arguments of next function are arbitrary and not used ... need to rework module
print "************************************************************************"
print " Building 3D mesh "
SS_model.build_3D_mesh_from_rasters(model_grid_raster_files, 
                                    SS_model.out_data_folder_grid, 1.0, 1000.0)
# Cleanup any isolated cells:
SS_model.reclassIsolatedCells()

print "************************************************************************"
print " Assign properties to mesh based on pilot points and zonal information"

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
           
zone_HGU = {HGU_zone[x]:x for x in HGU_zone.keys()}
           
pilot_points = True
           
if pilot_points:       
    # Set up pilot points:
    pilot_point_groups = ['hk', 'ss', 'sy']
    pp_group_dict = {}
    for pilot_points_group in pilot_point_groups:
        SS_model.create_pilot_points(pilot_points_group)           
        pp_group_dict[pilot_points_group] = SS_model.pilot_points[pilot_points_group]
        # create alias for brevity ...
        pp_grp = pp_group_dict[pilot_points_group]
        
        # Create some references to data inside the model builder object
        mesh_array = SS_model.model_mesh3D
        cell_centers = SS_model.model_mesh_centroids
        model_boundary = SS_model.model_boundary
        zones = len(np.unique(SS_model.model_mesh3D[1])) - 1
        
        # Create dict of zones and properties
        zone_prop_dict={zone: HGU_props['Kh mean'][HGU_map[HGU[zone]]] for zone in range(zones)}
        # Define some parameters for pilot point distribution
        if resolution == 1000:
            skip=[0, 0, 6, 0, 6, 6, 6] 
            skip_active=[49, 20, 0, 34, 0, 0, 0]
        elif resolution == 500:
            skip=[0, 0, 12, 0, 12, 12, 12] 
            skip_active=[100, 40, 0, 70, 0, 0, 0]
        else:
            skip=[0,  0, 3, 0, 2, 3, 3] 
            skip_active=[3, 20, 0, 4, 0, 0, 0]
        
        
        # Generate the pilot points 
        pp_grp.generate_points_from_mesh(mesh_array, cell_centers, 
            skip=skip, 
            skip_active=skip_active,
            zone_prop_dict=zone_prop_dict)
        
        # Create some necessary files fro pilot points utilities
        pp_grp.write_settings_fig()
        pp_grp.write_grid_spec(mesh_array, model_boundary, delc=resolution, delr=resolution)
        pp_grp.write_struct_file(mesh_array, nugget=0.0, 
                                 transform='log',numvariogram=1, variogram=0.15, 
                                 vartype=2, bearing=0.0, a=20000.0, anisotropy=1.0)
        
        # These search_radius values have been tested on the 1000m grid, would be good
        # to add in other resolution lists as they are developed
        if resolution == 1000:
            search_radius = [30000, 20000, 20000, 20000, 20000, 20000, 20000]
        else:
            search_radius = [30000, 20000, 20000, 20000, 20000, 20000, 20000]
            
        prefixes=['{}_{}'.format(pilot_points_group, zone_HGU[x]) for x in range(zones)]
        pp_grp.setup_pilot_points_by_zones(mesh_array, zones, search_radius, prefixes=prefixes)    
    
        pp_grp.generate_cov_mat_by_zones(zones)    
    
        #print("Running pyfac2real")
        pp_grp.run_pyfac2real_by_zones(zones)
    
    SS_model.save_pilot_points()
    hk = pp_group_dict['hk']     
    ss = pp_group_dict['ss']
    sy = pp_group_dict['sy']
    
for unit in HGU:
    if pilot_points:
        SS_model.parameters.create_model_parameter_set('kh_' + unit, 
                                                       value=HGU_props['Kh mean'][HGU_map[unit]], 
                                                       num_parameters=hk.num_ppoints_by_zone[HGU_zone[unit]])
        SS_model.parameters.parameter_options_set('kh_' + unit, 
                                              PARTRANS='log', 
                                              PARCHGLIM='factor', 
                                              PARLBND=HGU_props['Kh mean'][HGU_map[unit]] / 10., 
                                              PARUBND=HGU_props['Kh mean'][HGU_map[unit]] * 10., 
                                              PARGP='cond_' + unit, 
                                              SCALE=1, 
                                              OFFSET=0)
        SS_model.parameters.create_model_parameter_set('sy_' + unit, 
                                                       value=HGU_props['Sy mean'][HGU_map[unit]],
                                                       num_parameters=ss.num_ppoints_by_zone[HGU_zone[unit]])
        SS_model.parameters.parameter_options_set('sy_' + unit, 
                                              PARTRANS='log', 
                                              PARCHGLIM='factor', 
                                              PARLBND=1.0E-3, 
                                              PARUBND=0.8, 
                                              PARGP='sy_' + unit, 
                                              SCALE=1, 
                                              OFFSET=0)
        SS_model.parameters.create_model_parameter_set('ss_' + unit, 
                                                       value=HGU_props['SS mean'][HGU_map[unit]],
                                                       num_parameters=sy.num_ppoints_by_zone[HGU_zone[unit]])
        SS_model.parameters.parameter_options_set('ss_' + unit, 
                                              PARTRANS='log', 
                                              PARCHGLIM='factor', 
                                              PARLBND=HGU_props['SS mean'][HGU_map[unit]] / 10., 
                                              PARUBND=HGU_props['SS mean'][HGU_map[unit]] * 10., 
                                              PARGP='ss_' + unit, 
                                              SCALE=1, 
                                              OFFSET=0)


    else:
        SS_model.parameters.create_model_parameter('kh_' + unit, 
                                                   value=HGU_props['Kh mean'][HGU_map[unit]])
        SS_model.parameters.parameter_options('kh_' + unit, 
                                              PARTRANS='log', 
                                              PARCHGLIM='factor', 
                                              PARLBND=HGU_props['Kh mean'][HGU_map[unit]] / 10., 
                                              PARUBND=HGU_props['Kh mean'][HGU_map[unit]] * 10., 
                                              PARGP='cond_' + unit, 
                                              SCALE=1, 
                                              OFFSET=0)
        SS_model.parameters.create_model_parameter('sy_' + unit, value=HGU_props['Sy mean'][HGU_map[unit]])
        SS_model.parameters.parameter_options('sy_' + unit, 
                                              PARTRANS='log', 
                                              PARCHGLIM='factor', 
                                              PARLBND=1.0E-3, 
                                              PARUBND=0.8, 
                                              PARGP='sy_' + unit, 
                                              SCALE=1, 
                                              OFFSET=0)
        SS_model.parameters.create_model_parameter('ss_' + unit, value=HGU_props['SS mean'][HGU_map[unit]])
        SS_model.parameters.parameter_options('ss_' + unit, 
                                              PARTRANS='log', 
                                              PARCHGLIM='factor', 
                                              PARLBND=HGU_props['SS mean'][HGU_map[unit]] / 10., 
                                              PARUBND=HGU_props['SS mean'][HGU_map[unit]] * 10., 
                                              PARGP='ss_' + unit, 
                                              SCALE=1, 
                                              OFFSET=0)

        
    SS_model.parameters.create_model_parameter('kv_' + unit, value=HGU_props['Kz mean'][HGU_map[unit]])
    SS_model.parameters.parameter_options('kv_' + unit, 
                                          PARTRANS='fixed', 
                                          PARCHGLIM='factor', 
                                          PARLBND=HGU_props['Kz mean'][HGU_map[unit]] / 10., 
                                          PARUBND=HGU_props['Kz mean'][HGU_map[unit]] * 10., 
                                          PARGP='cond_' + unit, 
                                          SCALE=1, 
                                          OFFSET=0)

# This needs to be automatically generated from with the map_raster2mesh routine ...
zone_map = {1:'qa', 2:'utb', 3:'utqa', 4:'utam', 5:'utaf', 6:'lta', 7:'bse'}

Kh = SS_model.model_mesh3D[1].astype(float)
Kv = SS_model.model_mesh3D[1].astype(float)
Sy = SS_model.model_mesh3D[1].astype(float)
SS = SS_model.model_mesh3D[1].astype(float)
for key in zone_map.keys():
    if not pilot_points:
        Kh[Kh == key] = SS_model.parameters.param['kh_' + zone_map[key]]['PARVAL1']
        Sy[Sy == key] = SS_model.parameters.param['sy_' + zone_map[key]]['PARVAL1']
        SS[SS == key] = SS_model.parameters.param['ss_' + zone_map[key]]['PARVAL1']
    Kv[Kv == key] = SS_model.parameters.param['kv_' + zone_map[key]]['PARVAL1']

if pilot_points:
    Kh = hk.val_array
    # We are assuming an anisotropy parameter here where kh is 10 times kv ...
    Kv = hk.val_array * 0.1
    Sy = sy.val_array
    SS = ss.val_array

SS_model.properties.assign_model_properties('Kh', Kh)
SS_model.properties.assign_model_properties('Kv', Kv)
SS_model.properties.assign_model_properties('Sy', Sy)
SS_model.properties.assign_model_properties('SS', SS)

print "************************************************************************"
print " Interpolating rainfall data to grid "

interp_rain = SS_model.interpolate_points2mesh(rain_gauges, long_term_historic_rainfall, feature_id='Name', method='linear')
# Adjust rainfall to m from mm and from year to day
interp_rain = interp_rain/1000.0/365.0
# Adjust rainfall to recharge using 10% magic number

SS_model.boundaries.create_model_boundary_condition('Rainfall', 'rainfall', bc_static=True)
SS_model.boundaries.assign_boundary_array('Rainfall', interp_rain)

# Replace interp_rain with a copy to prevent alteration of assigned boundary array
interp_rain = np.copy(interp_rain)

recharge_zone_array = SS_model.map_raster_to_regular_grid_return_array(recharge_zones)

rch_zone_dict = {i:x for i, x in enumerate(np.unique(recharge_zone_array))}
rch_zones = len(rch_zone_dict.keys())

SS_model.parameters.create_model_parameter_set('ssrch', 
                                               value=0.01,
                                               num_parameters=rch_zones)
SS_model.parameters.parameter_options_set('ssrch', 
                                      PARTRANS='fixed', 
                                      PARCHGLIM='factor', 
                                      PARLBND=1.0E-3, 
                                      PARUBND=0.5, 
                                      PARGP='ssrch', 
                                      SCALE=1, 
                                      OFFSET=0)

for i in range(rch_zones - 1):
    interp_rain[recharge_zone_array == rch_zone_dict[i+1]] = \
        interp_rain[recharge_zone_array == rch_zone_dict[i+1]] * \
        SS_model.parameters.param['ssrch{}'.format(i)]['PARVAL1']

interp_rain[recharge_zone_array==rch_zone_dict[0]] = interp_rain[recharge_zone_array == rch_zone_dict[0]] * 0.0
interp_rain[SS_model.model_mesh3D[1][0] == -1] = 0.
    
rch = {}
rch[0] = interp_rain

print "************************************************************************"
print " Creating recharge boundary "

SS_model.boundaries.create_model_boundary_condition('Rain_reduced', 'recharge', bc_static=True)
SS_model.boundaries.associate_zonal_array_and_dict('Rain_reduced', recharge_zone_array, rch_zone_dict)
SS_model.boundaries.assign_boundary_array('Rain_reduced', rch)

# Initial heads using uniform head over entire model making whole domain saturated:
initial_heads_SS = np.full(SS_model.model_mesh3D[1].shape, 500.)

# Initial heads using ground surface elevation:
#for i in range(len(hu_raster_files_reproj)/2):
#    initial_heads_SS[i] = SS_model.model_mesh3D[1][0] #(interp_heads[hu_raster_files[2]]) # 2*i

SS_model.initial_conditions.set_as_initial_condition("Head", initial_heads_SS)#interp_heads[hu_raster_files[0]])

print "************************************************************************"
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


###############
##############
########  This needs to live elsewhere but will do for now ... TODO: Fix this
##############
###############
SS_model.GISInterface.raster_reproject_by_grid(surface_raster_high_res,
                                               surface_raster_high_res[:-4] + '_reproj.tif',
                                               resample_method='min')

surface_raster_high_res = surface_raster_high_res[:-4] + '_reproj.tif'


SS_model.map_points_to_grid(river_gauges, feature_id='Site_Name')
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
                                           value=[1., 1., 0.001, 0.0], 
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
                                           value=0.01, 
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
river_seg['strthick'] = np.interp(river_seg['Cumulative Length'].tolist(), known_points, strthick_val)

amalg_riv_points = []
for row in river_seg[['i', 'j']].iterrows():
    amalg_riv_points += [[row[1]['j'], row[1]['i']]]

# The depths in the column at row j and col i can be obtained using:
# SS_model.model_mesh3D[0][:,0,1]
def find_layer(elev, col_vals):
    for index, val in enumerate(col_vals):
        if elev > val:
            if index == 0:
                return index
            else:
                return index - 1
        #end if

# Sort out collocated stream reaches to avoid short circuiting:
river_seg['amalg_riv_points_tuple'] = river_seg['amalg_riv_points'].apply(lambda x: (x[0], x[1]))    
river_seg_group = river_seg.groupby(by='amalg_riv_points_tuple').count()
river_seg_group = river_seg_group[river_seg_group['amalg_riv_points'] > 1]

already_defined = []
old = []
for row in river_seg.iterrows():
    ind = row[0]
    row = row[1]
    new = row['amalg_riv_points']
    if new in old:
        already_defined += [ind]
    old += [new]

river_seg.loc[already_defined, 'strhc1'] = 0.0


new_k = []

for row in river_seg.iterrows():
    j_mesh = row[1]['i'] 
    i_mesh = row[1]['j']
    strtop = row[1]['strtop']
    strbot = row[1]['strtop'] - row[1]['strthick'] 
    new_k += [find_layer(strbot, SS_model.model_mesh3D[0][:, j_mesh, i_mesh])]

river_seg['k'] = new_k
       
# Remove any stream segments for which the elevation could not be mapped to a layer
#river_seg.dropna(inplace=True)

river_seg['ireach'] = 1
#river_seg['iseg'] = river_seg.index + 1
river_seg['iseg'] = [x + 1 for x in range(river_seg.shape[0])]

                      
# Set up bed elevations based on the gauge zero levels:
gauge_points = [x for x in zip(Campaspe.Easting, Campaspe.Northing)]
river_gauge_seg = SS_model.get_closest_riv_segments('Campaspe', gauge_points)
river_seg.loc[:, 'bed_from_gauge'] = np.nan

Campaspe['new_gauge'] = Campaspe[['Gauge Zero (Ahd)', 'Cease to flow level', 'Min value']].max(axis=1) 
Campaspe['seg_loc'] = river_gauge_seg         
Campaspe_gauge_zero = Campaspe[Campaspe['new_gauge'] > 10.]
# There are two values at the Campaspe weir, while it would be ideal to split the
# reach here it will cause problems for the segment
Campaspe_gauge_zero2 = Campaspe_gauge_zero[Campaspe_gauge_zero['Site Id'] != 406218]

#river_seg['bed_from_gauge'][river_seg['iseg'].isin(Campaspe_gauge_zero2['seg_loc'].tolist())] = sorted(Campaspe_gauge_zero2['new_gauge'].tolist(), reverse=True)
river_seg.loc[river_seg['iseg'].isin(Campaspe_gauge_zero2['seg_loc'].tolist()), 'bed_from_gauge'] = sorted(Campaspe_gauge_zero2['new_gauge'].tolist(), reverse=True)
river_seg['bed_from_gauge'] = river_seg.set_index(river_seg['Cumulative Length'])['bed_from_gauge'].interpolate(method='values', limit_direction='both').tolist()
#river_seg['bed_from_gauge'] = river_seg['bed_from_gauge'].interpolate(limit_direction='both')


new_k = []
for row in river_seg.iterrows():
    j_mesh = row[1]['i'] 
    i_mesh = row[1]['j']
    strbot = row[1]['bed_from_gauge'] - row[1]['strthick']
    new_k += [find_layer(strbot, SS_model.model_mesh3D[0][:, j_mesh, i_mesh])]

river_seg['k'] = new_k
river_seg['strtop'] = river_seg['bed_from_gauge']

# For stream reaches that didn't map properly to the mesh for z elevation we 
# can still include by setting to layer 0 with a bed hydraulic conductivity of 0
inds = np.where(river_seg['k'].isnull())[0]
#river_seg['strhc1'].loc[inds] = 0.0
#river_seg['k'].loc[inds] = 0
river_seg.dropna(inplace=True)
          
river_seg['iseg'] = [x + 1 for x in range(river_seg.shape[0])]
         
def slope_corrector(x):
    if  x  < 0.0001:
        return 0.0001
    else:
        return x
    # end if
    
river_seg['slope'] = river_seg['slope'].apply(lambda x: slope_corrector(x))

reach_df = river_seg[['k','i','j','iseg','ireach','rchlen','strtop','slope','strthick','strhc1']]
reach_data = reach_df.to_records(index=False)

# Create segment data
nseg = river_seg['iseg'].tolist()
icalc = [1] * len(nseg)
outseg = river_seg['iseg'] + 1
outseg = outseg.tolist()
outseg[-1] = 0
iupseg = [0] * len(nseg)
iprior = [0] * len(nseg)
nstrpts = [0] * len(nseg)
flow = [0] * len(nseg)
flow[0] = river_flow_data[406207]['Mean'].describe().loc['25%'] * 1000.
runoff = [0] * len(nseg)
etsw = [0] * len(nseg)
pptsw = [0] * len(nseg)

# Set the roughness for the channel
roughch_val = [SS_model.parameters.param['mn_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
roughch = np.interp(river_seg['Cumulative Length'].tolist(), 
                                known_points, roughch_val)
# Set the roughness for the banks
roughbk_val = [SS_model.parameters.param['mn_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
roughbk = np.interp(river_seg['Cumulative Length'].tolist(), 
                                known_points, roughbk_val)
river_seg['roughch'] = roughch
river_seg['roughbk'] = roughbk

cdpth = [0] * len(nseg)
fdpth = [0] * len(nseg)
awdth = [0] * len(nseg)
bwdth = [0] * len(nseg)

width1_val = [SS_model.parameters.param['rivwdth{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
width1 = width2 = np.interp(river_seg['Cumulative Length'].tolist(), 
                                known_points, width1_val)
river_seg['width2'] = river_seg['width1'] = width1

segment_data = pd.DataFrame({'nseg':nseg, 'icalc':icalc, 'outseg':outseg, 'iupseg':iupseg, 'iprior':iprior, 'nstrpts':nstrpts, \
                             'flow':flow, 'runoff':runoff, 'etsw':etsw, 'pptsw':pptsw, 'roughch':roughch, 'roughbk':roughbk, \
                             'cdpth':cdpth, 'fdpth':fdpth, 'awdth':awdth, 'bwdth':bwdth, 'width1':width1, 'width2':width2})
cols_ordered = ['nseg', 'icalc', 'outseg', 'iupseg', 'iprior', 'nstrpts', \
                'flow', 'runoff', 'etsw', 'pptsw', 'roughch', 'roughbk', \
                'cdpth', 'fdpth', 'awdth', 'bwdth', 'width1', 'width2']
segment_data = segment_data[cols_ordered]

SS_model.save_MODFLOW_SFR_dataframes('Campaspe', reach_df, segment_data)
SS_model.river_mapping['Campaspe'] = river_seg
                      
segment_data1 = segment_data.to_records(index=False)
seg_dict = {0: segment_data1}


###############################################################################

# Set up the gauges as observations 
#SS_model.observations.set_as_observations('Campase_riv_gauges', CampaspeRiv_obs_time_series, CampaspeRiv_points3D, domain='surface', obs_type='stage', units='m') 
                    
# Create observations for stage or discharge at those locations


print "************************************************************************"
print " Creating Campaspe river boundary"

SS_model.boundaries.create_model_boundary_condition('Campaspe River', 'river_flow', bc_static=True)
SS_model.boundaries.assign_boundary_array('Campaspe River', [reach_data, seg_dict])


print "************************************************************************"
print " Mapping Murray River to grid"

SS_model.map_polyline_to_grid(Murray_river_poly)

# Parameter to modify the stage, thus accounting for errors in values specified for stage
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

Murray = pd.merge(Murray, river_stage_data[1], on='Site Name', how='inner', suffixes=('','_r'))
Murray = Murray[[x for x in Murray.columns if '_r' not in x]]
def values_from_gauge(column):
    mriver_seg.loc[mriver_seg['iseg'].isin( \
        Murray_gauge_zero['seg_loc'].tolist()), column] = sorted( \
        Murray_gauge_zero['new_gauge'].tolist(), reverse=True)
    mriver_seg[column] = mriver_seg.set_index( \
        mriver_seg['Cumulative Length'])[column].interpolate( \
        method='values', limit_direction='both').tolist()
    mriver_seg[column].fillna(method='bfill', inplace=True)

values_to_edit = ['bed_from_gauge', 'stage_from_gauge']
for value in values_to_edit:
    values_from_gauge(value)

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

for layer in range(7):
    mriver_seg["surf{}".format(layer)] = surface_layers[layer]

mriver_seg['bottom_layer'] = bottom_layer

mriver_seg.plot(x='Cumulative Length', y=['bed_from_gauge'] + ["surf{}".format(x) for x in range(7)])
#mriver_seg.plot(x='Cumulative Length', y=['bed_from_gauge', 'bottom_layer'])

mriver_seg['k'] = new_k
mriver_seg['active'] = active
      
# Remove any stream segments for which the elevation could not be mapped to a layer
mriver_seg[mriver_seg['active'] == -1] = np.nan
mriver_seg.dropna(inplace=True)
SS_model.river_mapping['Murray'] = mriver_seg

mriver_seg['strtop'] = mriver_seg['stage_from_gauge']                      
                      
mriver_seg['strhc1'] = SS_model.parameters.param['kv_rm']['PARVAL1']                      

mriver_seg['width1'] = SS_model.parameters.param['rmwdth']['PARVAL1']

mriver_seg['stage'] = mriver_seg['strtop'] + SS_model.parameters.param['rmstage']['PARVAL1']

# Avoid collisions with Campaspe River ...
def is_in_other_river(riv_df_testing, riv_df_other):
    riv_df_other_locs = riv_df_other['amalg_riv_points'].tolist()
    cell_used = []
    for row in riv_df_testing.iterrows():
        if row[1]['amalg_riv_points'] in riv_df_other_locs:
            cell_used += [0]
        else:
            cell_used += [1]            
    #riv_df_testing['cell_used'] = cell_used
    return cell_used

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

riv = {}
riv[0] = simple_river

#print SS_model.polyline_mapped
print "************************************************************************"
print " Creating Murray River boundary"

SS_model.boundaries.create_model_boundary_condition('Murray River', 'river', bc_static=True)
SS_model.boundaries.assign_boundary_array('Murray River', riv)

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
        checked += [lay, row, col]
        if SS_model.model_mesh3D[1][0][row][col] == -1:
            continue
        MurrayGHBstage = mrow['stage'] + SS_model.parameters.param['mghb_stage']['PARVAL1']
        #if MurrayGHBstage < SS_model.model_mesh3D[0][lay][row][col]:
        #    continue
        if lay <= mrow['k']:
            continue

        Murray_df_ind += [ind]        
        Active_MurrayGHB_cells += [[lay, row, col]]

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

    return False
# End active_check()

Murray_df_ind2 = []
Final_MurrayGHB_cells = []
zone = SS_model.model_mesh3D[1]
shape = zone.shape
lay, row, col = 0, 1, 2
for index, ac in enumerate(Active_MurrayGHB_cells):
    # check if active GHB cell has any active non GHB cells N,E,S,W, above or below
    active_non_GHB = False
    acl, acr, acc = int(ac[lay]), int(ac[row]), int(ac[col])

    # Check north:
    ref_cell = [acl, acr + 1, acc]
    active_non_GHB = active_check(acr, 0, ref_cell, zone[acl, acr + 1, acc], Active_MurrayGHB_cells)

    # Check east:
    if not active_non_GHB:
        ref_cell = [acl, acr, acc + 1]
        active_non_GHB = active_check(acc, shape[col] - 1, ref_cell, zone[acl, acr, acc + 1], Active_MurrayGHB_cells)

    # Check south:
    if not active_non_GHB:
        ref_cell = [acl, acr - 1, acc]
        active_non_GHB = active_check(acr, shape[row] - 1, ref_cell, zone[acl, acr - 1, acc], Active_MurrayGHB_cells)

    # Check west:
    if not active_non_GHB:
        ref_cell = [acl, acr, acc - 1]
        active_non_GHB = active_check(acc, 0, ref_cell, zone[acl, acr, acc - 1], Active_MurrayGHB_cells)

    if active_non_GHB:
        Final_MurrayGHB_cells += [ac]
        Murray_df_ind2 += [Murray_df_ind[index]]

#for MurrayGHB_cell in SS_model.polyline_mapped['River_Murray_model.shp']:
for index, MurrayGHB_cell in enumerate(Final_MurrayGHB_cells):

    lay = MurrayGHB_cell[0]
    row = MurrayGHB_cell[1]
    col = MurrayGHB_cell[2]
        
    MurrayGHBstage = mriver_seg['stage'].loc[Murray_df_ind2[index]] + SS_model.parameters.param['mghb_stage']['PARVAL1']
    #if MurrayGHBstage < SS_model.model_mesh3D[0][0][row][col]:
    #    continue
    dx = SS_model.gridHeight
    dz = SS_model.model_mesh3D[0][lay][row][col] - SS_model.model_mesh3D[0][lay + 1][row][col]
    MGHBconductance = dx * dz * SS_model.parameters.param['mghbk']['PARVAL1']
    MurrayGHB += [[lay, row, col, MurrayGHBstage, MGHBconductance]]

ghb = {}
ghb[0] = MurrayGHB

print "************************************************************************"
print " Creating GHB boundary"

SS_model.boundaries.create_model_boundary_condition('GHB', 'general head', bc_static=True)
SS_model.boundaries.assign_boundary_array('GHB', ghb)

print "************************************************************************"
print " Creating parameters for transport "

# General parameters for transport
for unit in HGU:
    SS_model.parameters.create_model_parameter('por_{}'.format(unit), value=0.25)
    SS_model.parameters.parameter_options('por_{}'.format(unit), 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=0.05, 
                                          PARUBND=0.4, 
                                          PARGP='transport', 
                                          SCALE=1, 
                                          OFFSET=0)

SS_model.parameters.create_model_parameter('disp', value=0.01)
SS_model.parameters.parameter_options('disp', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=1E-4, 
                                      PARUBND=100., 
                                      PARGP='transport', 
                                      SCALE=1, 
                                      OFFSET=0)

print "************************************************************************"
print " Package up groundwater model builder object"

SS_model.package_model()
