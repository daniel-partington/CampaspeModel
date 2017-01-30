import sys
import os
import datetime

from osgeo import osr
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt


from HydroModelBuilder.GWModelBuilder import GWModelBuilder
from HydroModelBuilder.GISInterface.GDALInterface.GDALInterface import GDALInterface
from CampaspeModel.CustomScripts import processWeatherStations, getBoreData, get_GW_licence_info, processRiverStations, readHydrogeologicalProperties

# Define basic model parameters:
Proj_CS = osr.SpatialReference()
Proj_CS.ImportFromEPSG(28355) # 28355 is the code for gda94 mga zone 55; http://spatialreference.org/ref/epsg/gda94-mga-zone-55/

Interface = GDALInterface()
Interface.projected_coordinate_system = Proj_CS 
Interface.pcs_EPSG = "EPSG:28355"

SS_model = GWModelBuilder(name="GW_link_Integrated", 
                          data_folder=r"C:\Workspace\part0075\MDB modelling\testbox\input_data\\",
                          model_data_folder=r"C:\Workspace\part0075\MDB modelling\testbox\GW_link_Integrated\\",
                          out_data_folder=r"C:\Workspace\part0075\MDB modelling\testbox\data_build\\",
                          GISInterface=Interface,
                          model_type='Modflow',
                          mesh_type='structured')

# Cleanup
#SS_model.flush()

# Define the units for the project for consistency and to allow converions on input data
# SS_model.length = 'm'
# SS_model.time = 'd'

#******************************************************************************
#******************************************************************************
#
# Complimentary models requirements, i.e. bore and gauge data that should be 
# referenceable to this model for parsing specific outputs and receiving inputs:


model_linking = r"../testbox/integrated/data/model_linking.csv"
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

#******************************************************************************
#******************************************************************************


# Set the model boundary using a polygon shapefile:
print "************************************************************************"
print " Setting model boundary "
SS_model.set_model_boundary_from_polygon_shapefile("GW_model_area.shp", 
                                                      shapefile_path=SS_model.data_folder)

# Set data boundary for model data
print "************************************************************************"
print " Setting spatial data boundary "
SS_model.set_data_boundary_from_polygon_shapefile(SS_model.boundary_poly_file, 
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

# $%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%
# INCLUDE NSW bores in this next part too for better head representation at the border, i.e. Murray River

# Read in bore data:
print "************************************************************************"
print " Executing custom script: getBoreData "

bore_levels_file = "bore_levels"
bore_info_file = "bore_info"
if os.path.exists(SS_model.out_data_folder + bore_levels_file + ".h5") & os.path.exists(SS_model.out_data_folder + bore_info_file + ".h5"):
    bore_data_levels = SS_model.load_dataframe(SS_model.out_data_folder + bore_levels_file + ".h5")
    bore_data_info = SS_model.load_dataframe(SS_model.out_data_folder + bore_info_file + ".h5")
else:
    bore_data_levels, bore_data_info = getBoreData.getBoreData()
    SS_model.save_dataframe(SS_model.out_data_folder + bore_levels_file, bore_data_levels)
    SS_model.save_dataframe(SS_model.out_data_folder + bore_info_file, bore_data_info)
# end if

# getBoreDepth ... assuming that midpoint of screen interval is representative location and assign to layer accordingly
bore_data_info['depth'] = (bore_data_info['TopElev'] + bore_data_info['BottomElev'])/2.0

bore_data_info["HydroCode"] = bore_data_info.index

# For steady state model, only use bore details containing average level, not 
#observation_bores = SS_model.read_points_data(r"C:\Workspace\part0075\MDB modelling\ngis_shp_VIC\ngis_shp_VIC\NGIS_Bores.shp")

print "************************************************************************"
print " Read in and filtering bore spatial data "

bores_shpfile = SS_model.read_points_data(r"C:\Workspace\part0075\MDB modelling\ngis_shp_VIC\ngis_shp_VIC\NGIS_Bores.shp")

bores_filtered_from_shpfile = SS_model.points_shapefile_obj2dataframe(bores_shpfile, feature_id="HydroCode")

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
pumps_points = SS_model.read_points_data(r"C:\Workspace\part0075\MDB modelling\Campaspe_data\GW\Bore data\pumping wells.shp")

print "************************************************************************"
print " Executing custom script: readHydrogeologicalProperties "

file_location = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\GW\Aquifer properties\Hydrogeologic_variables.xlsx"
HGU_props = readHydrogeologicalProperties.getHGUproperties(file_location)

print "************************************************************************"
print "Get the C14 data"

C14_points = SS_model.read_points_data(r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Chemistry\C14.shp")    

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
if os.path.exists(SS_model.out_data_folder + river_flow_file + '.h5'):
    river_flow_data = SS_model.load_dataframe(SS_model.out_data_folder + river_flow_file + '.h5')
else:
    river_flow_data = processRiverStations.getFlow(path=r"C:\Workspace\part0075\MDB modelling\Campaspe_data\SW\All_streamflow_Campaspe_catchment\\")
    SS_model.save_dataframe(SS_model.out_data_folder + river_flow_file, river_flow_data)

river_stage_file = "river_stage_processed"
# Check if this data has been processed and if not process it
if os.path.exists(SS_model.out_data_folder + river_stage_file + '.h5'):
    river_stage_data = SS_model.load_dataframe(SS_model.out_data_folder + river_stage_file + '.h5')
else:
    river_stage_data = processRiverStations.getStage(path=r"C:\Workspace\part0075\MDB modelling\Campaspe_data\SW\All_streamflow_Campaspe_catchment\\")
    SS_model.save_dataframe(SS_model.out_data_folder + river_stage_file, river_stage_data)

#river_gauges = SS_model.read_points_data(r"C:\Workspace\part0075\MDB modelling\Campaspe_data\20_gauges\Site_info.shp")
river_gauges = SS_model.read_points_data(r"C:\Workspace\part0075\MDB modelling\Campaspe_data\SW\All_streamflow_Campaspe_catchment\processed_river_sites_stage.shp")

print "************************************************************************"
print "Load in the river shapefiles"

Campaspe_river_poly = SS_model.read_polyline("Campaspe_Riv.shp", path=r"C:\Workspace\part0075\MDB modelling\Campaspe_model\GIS\GIS_preprocessed\Surface_Water\Streams\\") 
Murray_river_poly = SS_model.read_polyline("River_Murray.shp", path=r"C:\Workspace\part0075\MDB modelling\Campaspe_model\GIS\GIS_preprocessed\Surface_Water\Streams\\") 

print "************************************************************************"
print "Load in the shapefiles defining groundwater boundaries"

WGWbound_poly = SS_model.read_polyline("western_head.shp", path=r"C:\Workspace\part0075\MDB modelling\testbox\input_data\\") 
EGWbound_poly = SS_model.read_polyline("eastern_head.shp", path=r"C:\Workspace\part0075\MDB modelling\testbox\input_data\\") 


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


print "************************************************************************"
print " Defining temporal aspects of the model"

start = datetime.date(2014, 01, 01)
#start = datetime.date(2014, 01, 01)
end = datetime.date(2015, 01, 01)
SS_model.model_time.set_temporal_components(steady_state=False, start_time=start, end_time=end, time_step='A')

# Define the grid width and grid height for the model mesh which is stored as a multipolygon shapefile GDAL object
print "************************************************************************"
print " Defining structured mesh"
SS_model.define_structured_mesh(1000, 1000)

# Read in hydrostratigraphic raster info for layer elevations:
hu_raster_path = r"C:\Workspace\part0075\MDB modelling\VAF_v2.0_ESRI_GRID\ESRI_GRID\\"

# Build basement file ... only need to do this once as it is time consuming so commented out for future runs
#SS_model.create_basement_bottom(hu_raster_path, "sur_1t", "bse_1t", "bse_2b", hu_raster_path)

hu_raster_files = ["qa_1t", "qa_2b", "utb_1t", "utb_2b", "utqa_1t", "utqa_2b", "utam_1t", "utam_2b", "utaf_1t", "utaf_2b", "lta_1t", "lta_2b", "bse_1t", "bse_2b.tif"]
#hu_raster_files = ["qa_1t", "qa_2b", "utb_1t", "utb_2b", "utqa_1t", "utqa_2b", "utam_1t", "utam_2b", "utaf_1t", "utaf_2b", "lta_1t", "lta_2b", "bse_1t", "bse_2b.tif"]
SS_model.read_rasters(hu_raster_files, path=hu_raster_path)
hu_raster_files_reproj = [x+"_reproj.bil" for x in hu_raster_files]

# Map HGU's to grid
print "************************************************************************"
print " Mapping rasters to grid "

hu_gridded_rasters = SS_model.map_rasters_to_grid(hu_raster_files, hu_raster_path)

# Build 3D grid
model_grid_raster_files = [x+"_model_grid.bil" for x in hu_raster_files]

# First two arguments of next function are arbitrary and not used ... need to rework module
print "************************************************************************"
print " Building 3D mesh "
SS_model.build_3D_mesh_from_rasters(model_grid_raster_files, SS_model.out_data_folder_grid, 1.0, 1000.0)
# Cleanup any isolated cells:
SS_model.reclassIsolatedCells()

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
    SS_model.parameters.create_model_parameter('Kh_' + unit, value=HGU_props['Kh mean'][HGU_map[unit]])
    SS_model.parameters.parameter_options('Kh_' + unit, 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=HGU_props['Kh mean'][HGU_map[unit]] / 10., 
                                          PARUBND=HGU_props['Kh mean'][HGU_map[unit]] * 10., 
                                          PARGP='cond', 
                                          SCALE=1, 
                                          OFFSET=0)
    SS_model.parameters.create_model_parameter('Kv_' + unit, value=HGU_props['Kz mean'][HGU_map[unit]])
    SS_model.parameters.parameter_options('Kv_' + unit, 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=HGU_props['Kz mean'][HGU_map[unit]] / 10., 
                                          PARUBND=HGU_props['Kz mean'][HGU_map[unit]] * 10., 
                                          PARGP='cond', 
                                          SCALE=1, 
                                          OFFSET=0)
    SS_model.parameters.create_model_parameter('Sy_' + unit, value=HGU_props['Sy mean'][HGU_map[unit]])
    SS_model.parameters.parameter_options('Sy_' + unit, 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=0., 
                                          PARUBND=0.8, 
                                          PARGP='spec_yield', 
                                          SCALE=1, 
                                          OFFSET=0)
    SS_model.parameters.create_model_parameter('SS_' + unit, value=HGU_props['SS mean'][HGU_map[unit]])
    SS_model.parameters.parameter_options('SS_' + unit, 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=HGU_props['SS mean'][HGU_map[unit]] / 10., 
                                          PARUBND=HGU_props['SS mean'][HGU_map[unit]] * 10., 
                                          PARGP='spec_stor', 
                                          SCALE=1, 
                                          OFFSET=0)

# This needs to be automatically generated from with the map_raster2mesh routine ...
zone_map = {1:'qa', 2:'utb', 3:'utqa', 4:'utam', 5:'utaf', 6:'lta', 7:'bse'}

Kh = SS_model.model_mesh3D[1].astype(float)
Kv = SS_model.model_mesh3D[1].astype(float)
Sy = SS_model.model_mesh3D[1].astype(float)
SS = SS_model.model_mesh3D[1].astype(float)
for key in zone_map.keys():
    Kh[Kh == key] = SS_model.parameters.param['Kh_' + zone_map[key]]['PARVAL1']
    Kv[Kv == key] = SS_model.parameters.param['Kv_' + zone_map[key]]['PARVAL1']
    Sy[Sy == key] = SS_model.parameters.param['Sy_' + zone_map[key]]['PARVAL1']
    SS[SS == key] = SS_model.parameters.param['SS_' + zone_map[key]]['PARVAL1']

SS_model.properties.assign_model_properties('Kh', Kh)
SS_model.properties.assign_model_properties('Kv', Kv)
SS_model.properties.assign_model_properties('Sy', Sy)
SS_model.properties.assign_model_properties('SS', SS)

print "************************************************************************"
print " Interpolating rainfall data to grid "

interp_rain = SS_model.interpolate_points2mesh(rain_gauges, long_term_historic_rainfall, feature_id='Name')
# Adjust rainfall to m from mm and from year to days
interp_rain = interp_rain/1000.0/365.0
# Adjust rainfall to recharge using 10% magic number

SS_model.boundaries.create_model_boundary_condition('Rainfall', 'rainfall', bc_static=True)
SS_model.boundaries.assign_boundary_array('Rainfall', interp_rain)

for i in [1,2,3,7]:
    SS_model.parameters.create_model_parameter('rch_red_'+zone_map[i], value=0.09)
    SS_model.parameters.parameter_options('rch_red_'+zone_map[i], 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=0., 
                                          PARUBND=0.9, 
                                          PARGP='rech_mult', 
                                          SCALE=1, 
                                          OFFSET=0)

    interp_rain[SS_model.model_mesh3D[1][0]==i] = interp_rain[SS_model.model_mesh3D[1][0]==i] * SS_model.parameters.param['rch_red_'+zone_map[i]]['PARVAL1']

rch = {}
rch[0] = interp_rain

print "************************************************************************"
print " Creating recharge boundary "

SS_model.boundaries.create_model_boundary_condition('Rain_reduced', 'recharge', bc_static=True)
SS_model.boundaries.assign_boundary_array('Rain_reduced', rch)

print "************************************************************************"
print " Mapping bores to grid "


SS_model.map_points_to_grid(bores_shpfile, feature_id = 'HydroCode')

bores_more_filter = []
for bores in SS_model.points_mapped["NGIS_Bores_clipped.shp"]:
    row, col = bores[0][0], bores[0][1]
    for bore in bores[1]: 
        try:
            bore_depth = bore_data_info.loc[bore, 'depth'] #[bore_data_info["HydroCode"] == HydroCode]['depth']        
        except:
            if bore in Ecology_bores:
                print 'Ecology bore not in info: ', bore
                #sys.exit('Halting model build due to bore not being found')
            if bore in Policy_bores:
                print 'Policy bore not in info: ', bore
                #sys.exit('Halting model build due to bore not being found')
            continue
        if bore_depth > SS_model.model_mesh3D[0][0][row][col]:
            #print 'Bore can't be above surface!!!   
            if bore in Ecology_bores:
                print 'Ecology bore above surf: ', bore
                #sys.exit('Halting model build due to bore not being mapped')
            if bore in Policy_bores:
                print 'Policy bore above surf: ', bore
                #sys.exit('Halting model build due to bore not being mapped')
            continue
        if bore_depth <= SS_model.model_mesh3D[0][-2][row][col]:
            #print 'Ignoring bores in bedrock!!!        
            if bore in Ecology_bores:
                print 'Ecology bore in  bedrock: ', bore
                #sys.exit('Halting model build due to bore not being mapped')
            if bore in Policy_bores:
                print 'Policy bore in bedrock: ', bore
                #sys.exit('Halting model build due to bore not being mapped')
            continue
        bores_more_filter += [bore]        

print 'Final bores within aquifers: ', len(bores_more_filter)

final_bores = final_bores[final_bores["HydroCode"].isin(bores_more_filter)]

ecology_found = [x for x in final_bores["HydroCode"] if x in Ecology_bores]                          
                          
bore_points = [[final_bores.loc[x, "Easting"], final_bores.loc[x, "Northing"]] for x in final_bores.index]

bore_points3D = final_bores[["HydroCode", "Easting", "Northing", "depth"]] # [[final_bores.loc[x, "Easting"], final_bores.loc[x, "Northing"], final_bores.loc[x, "depth"]] for x in final_bores.index]
bore_points3D = bore_points3D.set_index("HydroCode")

bores_obs_time_series = final_bores[["HydroCode", "mean level"]]

# Modify into standard format for the GWModelBuilder class
bores_obs_time_series = bores_obs_time_series.rename(columns={'HydroCode':'name', 'mean level':'value'})

SS_model.observations.set_as_observations('head', bores_obs_time_series, bore_points3D, domain='porous', obs_type='head', units='mAHD')


bores_in_layers = SS_model.map_points_to_raster_layers(bore_points, final_bores["depth"].tolist(), hu_raster_files_reproj)

# Map bores to layers to create initial head maps for different hydrogeological units
interp_heads = {}

for i in range(len(hu_raster_files_reproj)/2):
    bores_layer = np.array(bore_points)[np.array(bores_in_layers[i])]
    print 'Creating head map for: ', hu_raster_files[2*i]
    if bores_layer.shape[0] < 4: 
        #interp_heads[hu_raster_files[2*i]] = (SS_model.model_mesh3D[0][i]+SS_model.model_mesh3D[0][i+1])/2
        interp_heads[hu_raster_files[2*i]] = np.full(SS_model.model_mesh3D[1].shape[1:], np.NaN)
    else:
        bores_head_layer = np.array(final_bores["mean level"].tolist())[np.array(bores_in_layers[i])]
        unique_bores = np.unique(bores_layer) 
    
        b = np.ascontiguousarray(bores_layer).view(np.dtype((np.void, bores_layer.dtype.itemsize * bores_layer.shape[1])))
        _, idx = np.unique(b, return_index=True)
    
        unique_bores = bores_layer[idx]    
    
        interp_heads[hu_raster_files[2*i]] = SS_model.interpolate_points2mesh(bores_layer, bores_head_layer, use='griddata', method='linear')
        
initial_heads_SS = np.full(SS_model.model_mesh3D[1].shape, 0.)

for i in range(len(hu_raster_files_reproj)/2):
    initial_heads_SS[i] = SS_model.model_mesh3D[1][0] #(interp_heads[hu_raster_files[2]]) # 2*i

initial_heads_SS = np.full(SS_model.model_mesh3D[1].shape, 400.)

SS_model.initial_conditions.set_as_initial_condition("Head", initial_heads_SS)#interp_heads[hu_raster_files[0]])

# Map river polyline feature to grid including length of river in cell
print "************************************************************************"
print "Create observation wells for C14"

SS_model.map_points_to_grid(C14_points, feature_id='Bore_id')

wel = {}

# Create another column in the pandas dataframe for the C14 data for the depth
# at which the sample was taken in mAHD ... which will be calculated in the next
# for loop

df_C14['z'] = 'null'

i = 0
well_name = {}
for C14wells in SS_model.points_mapped['C14_clipped.shp']:
    row = C14wells[0][0]
    col = C14wells[0][1]
    for well in C14wells[1]: 
        try:
            well_depth = df_C14.loc[df_C14[df_C14['Bore_id'] == int(well)].index.tolist()[0], 'avg_screen(m)']
            #well_depth = df_C14.loc[df_C14['Bore_id'] == int(well), 'avg_screen(m)']
        except:
            print 'Well was excluded due to lack of information: ', int(well)            
            continue
        
        well_depth = SS_model.model_mesh3D[0][0][row][col] - well_depth

        df_C14.set_value(df_C14['Bore_id'] == int(well), 'z', well_depth)
                
        active = False
        for i in range(SS_model.model_mesh3D[1].shape[0]):
            if well_depth < SS_model.model_mesh3D[0][i][row][col] and well_depth > SS_model.model_mesh3D[0][i+1][row][col]:
                active_layer = i
                active = True
                break
        if active == False: 
            #print 'Well not placed: ', pump            
            continue

        if SS_model.model_mesh3D[1][active_layer][row][col] == -1:
            continue

        # Well sits in the mesh, so assign to well boundary condition
        well_name[i] = well
        i=i+1            
        try:
            wel[0] += [[active_layer, row, col, 0.]]
        except:
            wel[0] = [[active_layer, row, col, 0.]]

SS_model.boundaries.create_model_boundary_condition('C14_wells', 'wells', bc_static=True)
SS_model.boundaries.assign_boundary_array('C14_wells', wel)

C14_obs_time_series = df_C14.copy(deep=True) 
C14_obs_time_series = C14_obs_time_series[['Bore_id', 'a14C(pMC)']]
C14_obs_time_series.rename(columns={'Bore_id':'name', 'a14C(pMC)':'value'}, inplace=True)
C14_bore_points3D = df_C14[['Bore_id', 'zone55_easting', 'zone55_northing', 'z']]
C14_bore_points3D = C14_bore_points3D.set_index("Bore_id")
C14_bore_points3D.rename(columns={'zone55_easting':'Easting', 'zone55_northing':'Northing'}, inplace=True)

SS_model.observations.set_as_observations('C14', C14_obs_time_series, C14_bore_points3D, domain='porous', obs_type='concentration', units='pMC')

print "************************************************************************"
print " Mapping pumping wells to grid "

SS_model.map_points_to_grid(pumps_points, feature_id = 'OLD ID')

SS_model.parameters.create_model_parameter('pump_use', value=0.6)
SS_model.parameters.parameter_options('pump_use', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.2, 
                                      PARUBND=1., 
                                      PARGP='pumping', 
                                      SCALE=1, 
                                      OFFSET=0)

# Convert pumping_data to time series

#pumping_data_ts = pd.DataFrame(cols=['Works ID', 'datetime'])

# Existing data is only for 10 years from 2005 to 2015
#pump_date_index = pd.date_range(start=datetime.datetime(2005,07,01), end=datetime.datetime(2015,06,30), freq='AS-JUL')
pump_date_index = pd.date_range(start=datetime.datetime(2014,07,01), end=datetime.datetime(2015,06,30), freq='AS-JUL')

# Normalise pumping data, i.e. find individual pumps contribution to 
# total volume of water that is being pumped
total_pumping_rate = pumping_data['Use 2014/15'].mean()

wel = {}

pump_shallow = [] # Shallow (if <25m) or Deep (>= 25m)

for pump_cell in SS_model.points_mapped['pumping wells_clipped.shp']:
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
        pump_depth = SS_model.model_mesh3D[0][0][row][col] - pumping_data.loc[pump, 'Top screen depth (m)']        
        active = False
        for i in range(SS_model.model_mesh3D[0].shape[0]-1):
            if pump_depth < SS_model.model_mesh3D[0][i][row][col] and pump_depth > SS_model.model_mesh3D[0][i+1][row][col]:
                active_layer = i
                active = True
                break
        if active == False: 
            #print 'Well not placed: ', pump            
            continue
        #Get top of screen layer and calculate length of screen in layer

        # Specify if pump is shallow
        if pump_depth < 25:
            pump_shallow += [True]
        else:
            pump_shallow += [False]
        
        p14_15 = pumping_data.loc[pump, 'Use 2014/15'] / 365. * 1000.
        pump_rates = [p14_15]        
        pumping_data_ts = pd.DataFrame(pump_rates, columns=[pump], index=pump_date_index)
        pump_install = pumping_data.loc[pump, 'Construction date']
        
        if isinstance(pump_install, datetime.time):
            pump_install = datetime.date(1950,01,01)    
        pump_date_index2 = pd.date_range(start=pump_install, end=datetime.datetime(2004,06,30), freq='AS-JUL')

        #pump_allocation = pumping_data.loc[pump, 'Annual Volume'] / 365. * 1000.

        # Assume historical pumping is a percentage of lowest non-zero use for well        
        non_zero_pumping = [x for x in pump_rates if x > 0.]         
        if non_zero_pumping == []:
            pumping_rate_old = 0.
        else:
            pumping_rate_old = np.min(non_zero_pumping)

        old_pumping_ts = pd.DataFrame(index=pump_date_index2)
        old_pumping_ts[pump] = pumping_rate_old * SS_model.parameters.param['pump_use']['PARVAL1']

        # Merge the old and measured data
        pumping_data_ts = pd.concat([pumping_data_ts, old_pumping_ts])

        # Now let's resample the data on a monthly basis, and we will take the mean    
        pumping_data_ts = pumping_data_ts.resample(SS_model.model_time.t['time_step']).mean()

        # Let's also get rid of NaN data and replace with backfilling
        pumping_data_ts = pumping_data_ts.fillna(method='bfill')

        # Let's only consider times in our date range though
        date_index = pd.date_range(start=SS_model.model_time.t['start_time'], end=SS_model.model_time.t['end_time'], freq=SS_model.model_time.t['time_step'])
        pumping_data_ts = pumping_data_ts.reindex(date_index)    
        pumping_data_ts = pumping_data_ts.ix[SS_model.model_time.t['start_time']:SS_model.model_time.t['end_time']]
        pumping_data_ts = pumping_data_ts.fillna(0.0)

      
        
        # Now fill in the well dictionary with the values of pumping at relevant stress periods where Q is not 0.0
        for index, time in enumerate(pumping_data_ts.iterrows()):
            #if index >= SS_model.model_time.t['steps']: 
            #    continue
            #if time[1]['m3/day used'] != 0.0 :
            if time[1][pump] == 0:
                continue
            try:
                wel[index] += [[active_layer, row, col, -time[1][pump] / total_pumping_rate]]
            except:
                wel[index] = [[active_layer, row, col, -time[1][pump] / total_pumping_rate]]
                
print "************************************************************************"
print " Creating pumping boundary "

SS_model.boundaries.create_model_boundary_condition('licenced_wells', 'wells', bc_static=True)
SS_model.boundaries.assign_boundary_array('licenced_wells', wel)

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
   
#sw_stream_gauges = ['406214', '406219', '406201', '406224', '406218', '406202', '406265']
# Removed 406214 as it represents Axe Creek which is not included in this model
# Remvoed 406224 as it represents Mount Pleasant Creek which is not included in the model.
# Remvoed 406219 as it represents the head in the reservoir which is not modelled.
# Added in 406203 as it represents stage just downstream of the weir, NB not currently included in list of gauges online

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

#SS_model.map_points_to_grid(river_gauges, feature_id='Site_Name')
SS_model.map_points_to_grid(river_gauges, feature_id='Site_ID')

Campaspe_river_gauges = SS_model.points_mapped['processed_river_sites_stage_clipped.shp']

filter_gauges = []
for riv_gauge in Campaspe_river_gauges:
    #if riv_gauge[1][0] in use_gauges:
    if str(riv_gauge[1][0]) in Stream_gauges:
        filter_gauges += [riv_gauge]

SS_model.map_polyline_to_grid(Campaspe_river_poly)
SS_model.parameters.create_model_parameter('bed_depress', value=0.01)
SS_model.parameters.parameter_options('bed_depress', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.001, 
                                      PARUBND=0.1, 
                                      PARGP='spec_stor', 
                                      SCALE=1, 
                                      OFFSET=0)
SS_model.parameters.create_model_parameter('Kv_riv', value=5E-3)
SS_model.parameters.parameter_options('Kv_riv', 
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

# Map river from high to low
new_riv = SS_model.polyline_mapped['Campaspe_Riv_model.shp']
for index, riv_cell in enumerate(SS_model.polyline_mapped['Campaspe_Riv_model.shp']):
    row = riv_cell[0][0]
    col = riv_cell[0][1]
    new_riv[index] += [SS_model.model_mesh3D[0][0][row][col]]

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
        stages[index] = river_stage_data["Mean stage (m)"].loc[river_stage_data["Site ID"] == filter_gauges[gauge_ind[0]][1][0]]
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



for index, riv_cell in enumerate(SS_model.polyline_mapped['Campaspe_Riv_model.shp']):
    row = riv_cell[0][0]
    col = riv_cell[0][1]
    if SS_model.model_mesh3D[1][0][row][col] == -1:
        continue
    #print SS_model.model_mesh3D
    stage = stages[index] #SS_model.model_mesh3D[0][0][row][col]
    bed = beds[index] #SS_model.model_mesh3D[0][0][row][col] - SS_model.parameters.param['bed_depress']['PARVAL1']
    cond = riv_cell[1] * riv_width_avg * SS_model.parameters.param['Kv_riv']['PARVAL1'] / riv_bed_thickness
    simple_river += [[0, row, col, stage, cond, bed]]

riv = {}
riv[0] = simple_river

#print SS_model.polyline_mapped
print "************************************************************************"
print " Creating Campaspe river boundary"

SS_model.boundaries.create_model_boundary_condition('Campaspe River', 'river', bc_static=True)
SS_model.boundaries.assign_boundary_array('Campaspe River', riv)

print "************************************************************************"
print " Mapping Murray River to grid"

SS_model.map_polyline_to_grid(Murray_river_poly)

SS_model.parameters.create_model_parameter('RMstage', value=0.01)
SS_model.parameters.parameter_options('RMstage', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.001, 
                                      PARUBND=0.1, 
                                      PARGP='spec_stor', 
                                      SCALE=1, 
                                      OFFSET=0)
SS_model.parameters.create_model_parameter('Kv_RM', value=5E-3)
SS_model.parameters.parameter_options('Kv_RM', 
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
for riv_cell in SS_model.polyline_mapped['River_Murray_model.shp']:
    row = riv_cell[0][0]
    col = riv_cell[0][1]
    if SS_model.model_mesh3D[1][0][row][col] == -1:
        continue
    #print SS_model.model_mesh3D
    stage = SS_model.model_mesh3D[0][0][row][col]
    bed = SS_model.model_mesh3D[0][0][row][col] - SS_model.parameters.param['RMstage']['PARVAL1']
    cond = riv_cell[1] * riv_width_avg * SS_model.parameters.param['Kv_RM']['PARVAL1'] / riv_bed_thickness
    simple_river += [[0, row, col, stage, cond, bed]]

riv = {}
riv[0] = simple_river

#print SS_model.polyline_mapped
print "************************************************************************"
print " Creating Murray River boundary"

SS_model.boundaries.create_model_boundary_condition('Murray River', 'river', bc_static=True)
SS_model.boundaries.assign_boundary_array('Murray River', riv)

print "************************************************************************"
print " Setting up Murray River GHB boundary"

SS_model.parameters.create_model_parameter('MGHB_stage', value=0.01)
SS_model.parameters.parameter_options('MGHB_stage', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=-20.0, 
                                      PARUBND=50, 
                                      PARGP='ghb', 
                                      SCALE=1, 
                                      OFFSET=0)
SS_model.parameters.create_model_parameter('MGHBcond', value=5E-3)
SS_model.parameters.parameter_options('MGHBcond', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=1E-8, 
                                      PARUBND=50, 
                                      PARGP='ghb', 
                                      SCALE=1, 
                                      OFFSET=0)


MurrayGHB = []
for MurrayGHB_cell in SS_model.polyline_mapped['River_Murray_model.shp']:
    row = MurrayGHB_cell[0][0]
    col = MurrayGHB_cell[0][1]
    #print SS_model.model_mesh3D
    for lay in range(SS_model.model_mesh3D[1].shape[0]):    
        if SS_model.model_mesh3D[1][0][row][col] == -1:
            continue
        #MurrayGHBstage = (SS_model.model_mesh3D[0][lay+1][row][col] + SS_model.model_mesh3D[0][lay][row][col])/2. + SS_model.parameters.param['MGHB_stage']['PARVAL1']
        MurrayGHBstage = SS_model.model_mesh3D[0][0][row][col] + SS_model.parameters.param['MGHB_stage']['PARVAL1']
        if MurrayGHBstage < SS_model.model_mesh3D[0][0][row][col]:
            continue
        dx = SS_model.gridHeight
        dz = SS_model.model_mesh3D[0][lay][row][col] - SS_model.model_mesh3D[0][lay+1][row][col]
        MGHBconductance = dx * dz * SS_model.parameters.param['MGHBcond']['PARVAL1']
        MurrayGHB += [[lay, row, col, MurrayGHBstage, MGHBconductance]]

ghb = {}
ghb[0] = MurrayGHB

#print "************************************************************************"
#print " Mapping Western GW boundary to grid"
#
#SS_model.map_polyline_to_grid(WGWbound_poly)
#
#print "************************************************************************"
#print " Setting up Western GHB boundary"
#
#SS_model.parameters.create_model_parameter('WGHB_stage', value=0.01)
#SS_model.parameters.parameter_options('WGHB_stage', 
#                                      PARTRANS='log', 
#                                      PARCHGLIM='factor', 
#                                      PARLBND=-20.0, 
#                                      PARUBND=50, 
#                                      PARGP='ghb', 
#                                      SCALE=1, 
#                                      OFFSET=0)
#SS_model.parameters.create_model_parameter('WGHBcond', value=5E-3)
#SS_model.parameters.parameter_options('WGHBcond', 
#                                      PARTRANS='log', 
#                                      PARCHGLIM='factor', 
#                                      PARLBND=1E-8, 
#                                      PARUBND=50, 
#                                      PARGP='ghb', 
#                                      SCALE=1, 
#                                      OFFSET=0)
#
#
#WestGHB = []
#for WestGHB_cell in SS_model.polyline_mapped['western_head_model.shp']:
#    row = WestGHB_cell[0][0]
#    col = WestGHB_cell[0][1]
#    #print SS_model.model_mesh3D
#    for lay in range(SS_model.model_mesh3D[1].shape[0]):    
#        if SS_model.model_mesh3D[1][lay][row][col] == -1:
#            continue
#        #WestGHBstage = (SS_model.model_mesh3D[0][lay+1][row][col] + SS_model.model_mesh3D[0][lay][row][col])/2. + SS_model.parameters.param['WGHB_stage']['PARVAL1']
#        WestGHBstage = SS_model.model_mesh3D[0][0][row][col] + SS_model.parameters.param['WGHB_stage']['PARVAL1']
#        if WestGHBstage < SS_model.model_mesh3D[0][0][row][col]:
#            continue
#        dx = SS_model.gridHeight
#        dz = SS_model.model_mesh3D[0][lay][row][col] - SS_model.model_mesh3D[0][lay+1][row][col]
#        WGHBconductance = dx * dz * SS_model.parameters.param['WGHBcond']['PARVAL1']
#        WestGHB += [[lay, row, col, WestGHBstage, WGHBconductance]]
#
#ghb[0] += WestGHB
#
#
#print "************************************************************************"
#print " Mapping Eastern GW boundary to grid"
#
#SS_model.map_polyline_to_grid(EGWbound_poly)
#
#print "************************************************************************"
#print " Setting up Eastern GHB boundary"
#
#SS_model.parameters.create_model_parameter('EGHB_stage', value=0.01)
#SS_model.parameters.parameter_options('EGHB_stage', 
#                                      PARTRANS='log', 
#                                      PARCHGLIM='factor', 
#                                      PARLBND=-20, 
#                                      PARUBND=20, 
#                                      PARGP='ghb', 
#                                      SCALE=1, 
#                                      OFFSET=0)
#SS_model.parameters.create_model_parameter('EGHBcond', value=5E-3)
#SS_model.parameters.parameter_options('EGHBcond', 
#                                      PARTRANS='log', 
#                                      PARCHGLIM='factor', 
#                                      PARLBND=1E-8, 
#                                      PARUBND=50, 
#                                      PARGP='ghb', 
#                                      SCALE=1, 
#                                      OFFSET=0)
#
#
#EastGHB = []
#for EastGHB_cell in SS_model.polyline_mapped['eastern_head_model.shp']:
#    row = EastGHB_cell[0][0]
#    col = EastGHB_cell[0][1]
#    #print SS_model.model_mesh3D
#    for lay in range(SS_model.model_mesh3D[1].shape[0]):    
#        if SS_model.model_mesh3D[1][lay][row][col] == -1:
#            continue
#        #EastGHBstage = (SS_model.model_mesh3D[0][lay+1][row][col] + SS_model.model_mesh3D[0][lay][row][col])/2. + SS_model.parameters.param['EGHB_stage']['PARVAL1']
#        EastGHBstage = SS_model.model_mesh3D[0][0][row][col] + SS_model.parameters.param['EGHB_stage']['PARVAL1']
#        if EastGHBstage < SS_model.model_mesh3D[0][0][row][col]:
#            continue
#        
#        dx = SS_model.gridHeight
#        dz = SS_model.model_mesh3D[0][lay][row][col] - SS_model.model_mesh3D[0][lay+1][row][col]
#        EGHBconductance = dx * dz * SS_model.parameters.param['EGHBcond']['PARVAL1']
#        EastGHB += [[lay, row, col, EastGHBstage, EGHBconductance]]
#
#ghb[0] += EastGHB

print "************************************************************************"
print " Creating GHB boundary"

SS_model.boundaries.create_model_boundary_condition('GHB', 'general head', bc_static=True)
SS_model.boundaries.assign_boundary_array('GHB', ghb)

print "************************************************************************"
print " Collate observations"

SS_model.map_obs_loc2mesh3D(method='nearest')
#SS_model.map_obs2model_times()
#SS_model.observations.collate_observations()

print "************************************************************************"
print " Package up groundwater model builder object"

SS_model.package_model()
