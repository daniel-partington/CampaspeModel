import os
import datetime

from osgeo import osr
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt

from HydroModelBuilder.GWModelBuilder import GWModelBuilder
from HydroModelBuilder.GISInterface.GDALInterface.GDALInterface import GDALInterface

from CampaspeModel.CustomScripts import Campaspe_data
from CampaspeModel.build_common import Campaspe_mesh
from CampaspeModel.build_utils.multifrequency_resampling import resample_to_model_data_index 
from CampaspeModel.build_utils.multifrequency_resampling import resample_obs_time_series_to_model_data_index
from CampaspeModel.build_common.rainfall_recharge import prepare_transient_rainfall_data_for_model 
from CampaspeModel.build_common.groundwater_boundary import prepare_ghb_boundary_from_Murray_data
from CampaspeModel.build_common import rivers
from CampaspeModel.build_common.pumping import prepare_pumping_data_for_model
from CampaspeModel.build_common.channels import prepare_channel_data_for_model
from CampaspeModel.build_common.drains import prepare_drain_data_for_model

"""

"""
# Define basic model parameters:
Proj_CS = osr.SpatialReference()
Proj_CS.ImportFromEPSG(28355) # 28355 is the code for gda94 mga zone 55; http://spatialreference.org/ref/epsg/gda94-mga-zone-55/

Interface = GDALInterface()
Interface.projected_coordinate_system = Proj_CS 
Interface.pcs_EPSG = "EPSG:28355"

root_dir = r"C:\Workspace\part0075\MDB modelling\testbox"
tr_model = GWModelBuilder(name="02_transient_flow", 
                          data_folder=os.path.join(root_dir,r"input_data\\"),
                          model_data_folder=os.path.join(root_dir,r"00_Campaspe_Cascade\02_transient_flow\\"),
                          out_data_folder=os.path.join(root_dir,r"data_build\\"),
                          GISInterface=Interface,
                          model_type='Modflow',
                          mesh_type='structured')

# Define the units for the project for consistency and to allow converions on input data
# tr_model.length = 'm'
# tr_model.time = 'd'

Campaspe_data_folder = r"C:\Workspace\part0075\MDB modelling\Campaspe_data"
hu_raster_path = r"C:\Workspace\part0075\MDB modelling\ESRI_GRID_raw\ESRI_GRID"

custom_data = \
    Campaspe_data.process_custom_scripts_and_spatial_data(tr_model, 
                                                          Campaspe_data_folder,
                                                          verbose=True)
HGU_props = custom_data['HGU_props']
rain_gauges = custom_data['rain_gauges']
long_term_historic_weather = custom_data['long_term_historic_weather'] 
recharge_zones = custom_data['recharge_zones']
surface_raster_high_res = custom_data['surface_raster_high_res'] 
surface_raster_high_res_GSA = custom_data['surface_raster_high_res_GSA'] 
river_gauges = custom_data['river_gauges']
Campaspe_river_poly_file = custom_data['Campaspe_river_poly_file']
Murray_river_poly_file = custom_data['Murray_river_poly_file']
river_flow_data = custom_data['river_flow_data']
river_stage_data = custom_data['river_stage_data']
river_ec_data = custom_data['river_ec_data']
river_diversion_data = custom_data['river_diversions_data']
Campaspe = custom_data['Campaspe']
Campaspe_field_elevations = custom_data['Campaspe_field_elevations']
Campaspe_relevant = custom_data['Campaspe_relevant']
bores_shpfile = custom_data['bores_shpfile']
final_bores = custom_data['final_bores'] 
pumping_data = custom_data['pumping_data']
C14_points = custom_data['C14_points']
df_C14 = custom_data['df_C14']
pumps_points = custom_data['pumps_points']
FieldData = custom_data['FieldData']
bore_data_info = custom_data['bore_data_info']
bore_data_levels = custom_data['bore_data_levels']

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
start_irrigation = datetime.date(1881, 01, 01)
start_pumping = datetime.date(1966, 01, 01)
start_time_interest = datetime.date(2015, 01, 01)

end = datetime.date(2017, 05, 31)

end_post_clearance = datetime.date(1880, 12, 31)
before_pumping = datetime.date(1965, 12, 31)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@ CONSTRUCTION OF TIME PERIODS FOR MODEL @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
frequencies = ['40A', '84A', '20A', 'M']

date_index_post_clearance = pd.date_range(start=start, end=start_irrigation, freq=frequencies[0])
date_index_post_irrigation = pd.date_range(start=start_irrigation, end=start_pumping, freq=frequencies[1])
date_index_post_pumping = pd.date_range(start=start_pumping, end=start_time_interest, freq=frequencies[2])
date_index_time_of_interest = pd.date_range(start=start_time_interest, end=end, freq=frequencies[3])

date_group = [start, start_irrigation, start_pumping, \
              start_time_interest, end]

date_index_temp = date_index_post_clearance[:-1].append(date_index_post_irrigation)
date_index = date_index_temp[:-1].append(date_index_post_pumping)
date_index = date_index.append(date_index_time_of_interest)

tr_model.model_time.set_temporal_components(steady_state=False, start_time=start, end_time=end, date_index=date_index)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

HGU, hu_raster_files_reproj = Campaspe_mesh.build_mesh_and_set_properties(tr_model,
                                                  hu_raster_path,
                                                  HGU_props,
                                                  resolution=1000,
                                                  create_basement=False)

print "************************************************************************"
print " Interpolating rainfall data to grid and time steps"

interp_rain, interp_et, recharge_zone_array, rch_zone_dict = \
    prepare_transient_rainfall_data_for_model(tr_model,
                                              long_term_historic_weather,
                                              recharge_zones,
                                              long_term_historic_weather,
                                              date_index,
                                              frequencies,
                                              date_group,
                                              start,
                                              end,
                                              rain_gauges)    

    
print "************************************************************************"
print " Creating recharge boundary "

tr_model.boundaries.create_model_boundary_condition('Rain_reduced', 'recharge', bc_static=False)
tr_model.boundaries.associate_zonal_array_and_dict('Rain_reduced', recharge_zone_array, rch_zone_dict)
tr_model.boundaries.assign_boundary_array('Rain_reduced', interp_rain)

print "************************************************************************"
print " Mapping bores to grid "

tr_model.map_points_to_grid(bores_shpfile, feature_id='HydroCode')

bores_more_filter = []
bores_above_surface = []
bores_below_top_of_bedrock = []
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
            bores_above_surface += [bore]
            continue
        if bore_depth <= tr_model.model_mesh3D[0][-2][row][col]:
            #print 'Ignoring bores in bedrock!!!        
            bores_below_top_of_bedrock += [bore]
            continue
        bores_more_filter += [bore]        

print('Bores above the surface: {}'.format(len(bores_above_surface)))
print('Bores below top of bedrock: {}'.format(len(bores_below_top_of_bedrock)))
print('Final bores within aquifers: {}'.format(len(bores_more_filter)))

final_bores = final_bores[final_bores["HydroCode"].isin(bores_more_filter)]

bore_points = [[final_bores.loc[x, "Easting"], final_bores.loc[x, "Northing"]] for x in final_bores.index]

bore_points3D = final_bores[["HydroCode", "Easting", "Northing", "depth"]] # [[final_bores.loc[x, "Easting"], final_bores.loc[x, "Northing"], final_bores.loc[x, "depth"]] for x in final_bores.index]
bore_points3D = bore_points3D.set_index("HydroCode")

bores_obs_time_series = bore_data_levels[bore_data_levels["HydroCode"].isin( \
                        final_bores["HydroCode"])]

# Modify into standard format for the GWModelBuilder class
bores_obs_time_series = bores_obs_time_series.rename(columns= \
                                                     {'HydroCode':'name', \
                                                     'bore_date':'datetime', \
                                                     'result':'value'})

bores_obs_time_series['datetime'] = pd.to_datetime(bores_obs_time_series['datetime'])

# Kill all values where head observation is 0.0 indicating bad data that was missed
# by quality control codes in pre-processing function getBoreData
bores_obs_time_series = bores_obs_time_series[bores_obs_time_series['value'] > 5.]
# Remove values that are false records, i.e. bores with a reading at 30/12/1899
bores_obs_time_series = bores_obs_time_series[bores_obs_time_series['datetime'] \
                                              != datetime.datetime(1899,12,30)]
# Only consider head observations after 2010 to limit the size of the jacobian in PEST
bores_obs_time_series = bores_obs_time_series[bores_obs_time_series['datetime'] \
                                              > datetime.datetime(2010,12,30)]

# Create outlier identifier, i.e mutltiplier for standard deviations from the mean
# Using 4 here which is a rather large bound                                              
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

# Now to resample the time series to get the mean heads within a time interval
bores_obs_time_series = resample_obs_time_series_to_model_data_index(
                            bores_obs_time_series, 
                            date_index, 
                            frequencies, 
                            date_group,
                            start,
                            end,
                            fill='none')
                               
# For the weigts of observations we need to specify them as 1/sigma, where sigma is the standard deviation of measurement error
tr_model.observations.set_as_observations('head', 
                                          bores_obs_time_series, 
                                          bore_points3D, 
                                          domain='porous', 
                                          obs_type='head', 
                                          units='mAHD', 
                                          weights=1.0/0.2, 
                                          by_zone=True)

bores_in_layers = tr_model.map_points_to_raster_layers(bore_points, final_bores["depth"].tolist(), hu_raster_files_reproj)

# Initalise model with head from elevations
initial_heads_tr = np.full(tr_model.model_mesh3D[1].shape, 0.)

for i in range(len(hu_raster_files_reproj)/2):
    initial_heads_tr[i] = (tr_model.model_mesh3D[0][i] + tr_model.model_mesh3D[0][i + 1]) / 2.

tr_model.initial_conditions.set_as_initial_condition("Head", initial_heads_tr) #interp_heads[hu_raster_files[0]])

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
C14_obs_time_series = C14_obs_time_series[['Bore_id', 'a14C_corrected']]
C14_obs_time_series['datetime'] = pd.to_datetime(datetime.date(2015,12,30))
C14_obs_time_series.rename(columns={'Bore_id':'name', 'a14C_corrected':'value'}, inplace=True)
C14_bore_points3D = df_C14[['Bore_id', 'zone55_easting', 'zone55_northing', 'z']]
C14_bore_points3D = C14_bore_points3D[C14_bore_points3D['z'] != 'null']
C14_bore_points3D = C14_bore_points3D.set_index("Bore_id")
C14_bore_points3D.rename(columns={'zone55_easting':'Easting', 'zone55_northing':'Northing'}, inplace=True)

tr_model.observations.set_as_observations('C14', 
                                          C14_obs_time_series, \
                                          C14_bore_points3D, 
                                          domain='porous', \
                                          obs_type='concentration', \
                                          units='pMC', \
                                          weights=1.0 / 5.0)

print "************************************************************************"
print " Mapping pumping wells to grid "

wel = prepare_pumping_data_for_model(tr_model,
                                   pumps_points,
                                   start_pumping,
                                   start,
                                   end,
                                   date_index,
                                   pumping_data,
                                   frequencies,
                                   date_group)
                
print "************************************************************************"
print " Creating pumping boundary "

tr_model.boundaries.create_model_boundary_condition('licenced_wells', 'wells', bc_static=True)
tr_model.boundaries.assign_boundary_array('licenced_wells', wel)

print "************************************************************************"
print " Mapping Campaspe river to grid"

num_reaches = 20    
river_seg, reach_df, reach_data, known_points = \
    rivers.prepare_river_data_for_Campaspe(tr_model, 
                                    surface_raster_high_res,
                                    river_gauges,
                                    Campaspe_river_poly_file,
                                    Campaspe,
                                    Campaspe_field_elevations,
                                    num_reaches=num_reaches)

Campaspe_info = Campaspe
Campaspe_info.index = Campaspe_info['Site Id']
Campaspe_info = Campaspe_info[['Easting', 'Northing', 'Site Name', 'seg_loc']]

segment_data, seg_dict = \
    rivers.create_segment_data_transient(tr_model,
                                      river_seg,
                                      river_flow_data,
                                      Campaspe_info,
                                      river_diversion_data,
                                      FieldData,
                                      interp_et,
                                      interp_rain,
                                      date_index, 
                                      frequencies, 
                                      date_group,
                                      start,
                                      end)    
    
tr_model.save_MODFLOW_SFR_dataframes('Campaspe', reach_df, segment_data)
tr_model.river_mapping['Campaspe'] = river_seg

FieldData_info = FieldData.groupby('Name').first()[['Easting', 'Northing', 'Distance_Eppalock', 'Zone']]
field_points = [x for x in zip(FieldData_info.Easting, FieldData_info.Northing)]
river_field_seg = tr_model.get_closest_riv_segments('Campaspe', field_points)
FieldData_info.loc[:, 'seg_loc'] = river_field_seg

Camp_riv_cells = [x for x in zip(river_seg['i'], river_seg['j'])]
                      
print "************************************************************************"
print " Creating Campaspe river observations for stage and discharge at "
print " locations downstream of Lake Eppalock"


Campaspe_flow_sites_list = Campaspe['Site Id'].tolist()
 
# Vic Government stream gauges:
stage_time_series = pd.DataFrame()
for key in river_stage_data[0].keys():
    # Ignore 406207 Lake Eppalock    
    if key == 406207:
        continue
    elif key not in Campaspe_flow_sites_list:
        continue
    # end if
    site_ts = river_stage_data[0][key].copy()
    site_ts.loc[:, 'name'] = key
    site_ts['datetime'] = site_ts.index
    site_ts.index = range(site_ts.shape[0])
    site_ts.rename(columns={'Mean':'value'}, inplace=True)
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # NO FILTERING ON "Qual" AT THIS TIME SO BAD DATA ARE POSSIBLE!!!
    site_ts = site_ts[['name', 'value', 'datetime']]
    stage_time_series = pd.concat([stage_time_series, site_ts])
    
# Now to resample the time series to get the mean stage within a time interval
stage_time_series_resampled = resample_obs_time_series_to_model_data_index(
                                  stage_time_series, 
                                  date_index, 
                                  frequencies, 
                                  date_group,
                                  start, 
                                  end,
                                  fill='none')
    
tr_model.observations.set_as_observations('stage', 
                                          stage_time_series_resampled, 
                                          Campaspe_info, 
                                          domain='stream', 
                                          obs_type='stage', 
                                          units='mAHD', 
                                          weights=1.0, 
                                          real=True)

discharge_time_series = pd.DataFrame()
for key in river_flow_data.keys():
    # Ignore 406207 Lake Eppalock    
    if key == 406207:
        continue
    elif key not in Campaspe_flow_sites_list:
        continue
    # end if
    # Ignore keys not from Campaspe river
    site_ts = river_flow_data[key].copy()
    site_ts.loc[:, 'name'] = key
    site_ts['datetime'] = site_ts.index
    site_ts.index = range(site_ts.shape[0])
    site_ts.rename(columns={'Mean':'value'}, inplace=True)
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # NO FILTERING ON "Qual" AT THIS TIME SO BAD DATA ARE POSSIBLE!!!
    site_ts = site_ts[['name', 'value', 'datetime']]
    discharge_time_series = pd.concat([discharge_time_series, site_ts])

# Now to resample the time series to get the mean discharge within a time interval
discharge_time_series_resampled = resample_obs_time_series_to_model_data_index(
                                  discharge_time_series, 
                                  date_index, 
                                  frequencies, 
                                  date_group,
                                  start,
                                  end,
                                  fill='none')

tr_model.observations.set_as_observations('gflow', 
                                          discharge_time_series_resampled, 
                                          Campaspe_info, 
                                          domain='stream', 
                                          obs_type='discharge', 
                                          units='m3/d', 
                                          weights=1.0, 
                                          real=True)    


# Field data for stage and discharge
field_discharge_time_series = FieldData[['Name', 'Flows_Field']]
field_discharge_time_series.loc[:, 'datetime'] = field_discharge_time_series.index
field_discharge_time_series.index = range(field_discharge_time_series.shape[0])
field_discharge_time_series.rename(columns={'Name':'name', 'Flows_Field': 'value'}, inplace=True)
field_discharge_time_series = field_discharge_time_series.dropna()
# Ignore tributaries for now ...
field_discharge_time_series = \
    field_discharge_time_series[field_discharge_time_series['name'].str.contains('Camp')]

# Convert data units in m3/s to that of model, i.e. m3/d
field_discharge_time_series['value'] = field_discharge_time_series['value'] * 86400.

tr_model.observations.set_as_observations('fflow', 
                                          field_discharge_time_series, 
                                          FieldData_info, 
                                          domain='stream', 
                                          obs_type='discharge', 
                                          units='m3/d', 
                                          weights=1.0, 
                                          real=True)

field_stage_time_series = FieldData[['Name', 'Depth_Field']]
field_stage_time_series.loc[:, 'datetime'] = field_stage_time_series.index
field_stage_time_series.index = range(field_stage_time_series.shape[0])
field_stage_time_series.rename(columns={'Name':'name', 'Depth_Field': 'value'}, inplace=True)
field_stage_time_series = field_stage_time_series.dropna()
# Ignore tributaries for now ...
field_stage_time_series = field_stage_time_series[field_stage_time_series['name'].str.contains('Camp')]

# NEED DATAFRAME WITH COLS 'Gauge_id' and 'river_seg' 'Easting', 'Northing'

tr_model.observations.set_as_observations('fdepth', 
                                          field_stage_time_series, 
                                          FieldData_info, 
                                          domain='stream', 
                                          obs_type='depth', 
                                          units='m', 
                                          weights=1.0, 
                                          real=True)

print "************************************************************************"
print " Creating Campaspe river observations for EC at "
print " locations downstream of Lake Eppalock"
ec_time_series = pd.DataFrame()
EC_input_site = 406207
for key in river_ec_data.keys():
    # Ignore 406219 Lake Eppalock    
    if key == EC_input_site:
        continue
    elif key not in Campaspe_flow_sites_list:
        continue
    # end if
    site_ts = river_ec_data[key].copy()
    site_ts.loc[:, 'name'] = key
    site_ts['datetime'] = site_ts.index
    site_ts.index = range(site_ts.shape[0])
    site_ts.rename(columns={'Mean':'value'}, inplace=True)
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # NO FILTERING ON "Qual" AT THIS TIME SO BAD DATA ARE POSSIBLE!!!
    site_ts = site_ts[['name', 'value', 'datetime']]
    ec_time_series = pd.concat([ec_time_series, site_ts])

# Now to resample the time series to get the mean discharge within a time interval
ec_time_series_resampled = resample_obs_time_series_to_model_data_index(
                                  ec_time_series, 
                                  date_index, 
                                  frequencies, 
                                  date_group,
                                  start,
                                  end,
                                  fill='none')

tr_model.observations.set_as_observations('gstrec', 
                                          ec_time_series_resampled, 
                                          Campaspe_info, 
                                          domain='stream', 
                                          obs_type='EC', 
                                          units='uS/cm', 
                                          weights=1.0, 
                                          real=True)    


field_ec_time_series = FieldData[['Name', 'EC_Field']]
field_ec_time_series.loc[:, 'datetime'] = field_ec_time_series.index
field_ec_time_series.index = range(field_ec_time_series.shape[0])
field_ec_time_series.rename(columns={'Name':'name', 'EC_Field': 'value'}, inplace=True)
field_ec_time_series = field_ec_time_series.dropna()
# Ignore tributaries for now ...
field_ec_time_series = field_ec_time_series[field_ec_time_series['name'].str.contains('Camp')]

# NEED DATAFRAME WITH COLS 'Gauge_id' and 'river_seg' 'Easting', 'Northing'

tr_model.observations.set_as_observations('fstrec', 
                                          field_ec_time_series, 
                                          FieldData_info, 
                                          domain='stream', 
                                          obs_type='EC', 
                                          units='uS/cm', 
                                          weights=1.0, 
                                          real=True)


print "************************************************************************"
print " Creating Campaspe river observations for Radon at "
print " locations downstream of Lake Eppalock"

radon_time_series = FieldData[['Name', 'Radon']]
radon_time_series.loc[:, 'datetime'] = radon_time_series.index
radon_time_series.index = range(radon_time_series.shape[0])
radon_time_series.rename(columns={'Name':'name', 'Radon': 'value'}, inplace=True)
radon_time_series = radon_time_series.dropna()
# Ignore tributaries for now ...
radon_time_series = radon_time_series[radon_time_series['name'].str.contains('Camp')]
# Adjust units of Radon from Bq/l to mBq/l
radon_time_series['value'] = radon_time_series['value'] * 1000.
# NEED DATAFRAME WITH COLS 'Gauge_id' and 'river_seg' 'Easting', 'Northing'

tr_model.observations.set_as_observations('radon', 
                                          radon_time_series, 
                                          FieldData_info, 
                                          domain='stream', 
                                          obs_type='Radon', 
                                          units='Bq/l', 
                                          weights=1.0, 
                                          real=True)

   
   
print "************************************************************************"
print " Creating Campaspe river simulated exchange observations for data worth analysis"
   
fy_start = end - datetime.timedelta(days=365)
fy_end = end 
Years = 1

def create_obs_for_sw_gw_interaction(entries, name, start, end, freq, riv_segs, \
                                     obs_name, obs_type):
    '''
    Function to create the necessary observations for sw-gw exchange for the 
    different spatial and temporal periods
    '''
    if entries == 1:
        name_list = '{}'.format(name) 
    else:
        name_list = ['{}{}'.format(name, x) for x in range(0, entries)]
    # end if
    time_series = pd.DataFrame({'name': name_list, 
                                'value':[0.0] * entries, 
                                'datetime':pd.date_range(start=start, end=end, freq=freq)[1:]})
    # Create zero weighted observations for model predictions of interest
    tr_model.observations.set_as_observations(obs_name, 
                                              time_series, 
                                              riv_segs, 
                                              domain='stream', 
                                              obs_type=obs_type, 
                                              units='m^3/d', 
                                              weights=0.0, 
                                              real=False)

# Create river spatial groups: Whole of river, reach by gauge, reach by segment
entries = np.array([1, 4, 12]) * Years
names = ['a_swgw', 's_swgw', 'm_swgw']
swgw_exch_obs_freqs = ['A-MAY', '3M', 'M']
obs_names_whole = ['nrf_a', 'nrf_s', 'nrf_m']
obs_types = ['swgw_a', 'swgw_s', 'swgw_m']

for i in range(3):
    create_obs_for_sw_gw_interaction(entries[i], 
                                     names[i], 
                                     fy_start, 
                                     fy_end, 
                                     swgw_exch_obs_freqs[i], 
                                     river_seg['iseg'],
                                     obs_names_whole[i], 
                                     obs_types[i])

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# rrf = reach river flux
# srf = segment river flux

obs_names_g2g_reach = ['rrf_a', 'rrf_s', 'rrf_m']
# Create reaches based on locations of gauging stations:
gauge_locations = np.sort(Campaspe_info['seg_loc'].unique())
reach_no = range(len(gauge_locations))
river_seg['reach'] = np.nan    
river_seg.loc[river_seg['iseg'].isin(gauge_locations), 'reach'] = reach_no         
river_seg['reach'].fillna(method='ffill', inplace=True)
river_seg['reach'].fillna(method='bfill', inplace=True)
river_seg['reach'].astype(int, inplace=True)
river_segs_reach = [river_seg['iseg'][river_seg['reach'] == x].tolist() for x in reach_no]

obs_names_seg = ['srf_a', 'srf_s', 'srf_m']
river_segs_seg = [x for x in river_seg['iseg']]
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

for reach in reach_no:
    entries_reach = entries
    names_reach = [x + str(reach) for x in names]
    obs_names_reach = [x + str(reach) for x in obs_names_g2g_reach]
    obs_types_reach = [x + "r" for x in obs_types]
    for i in range(3):
        create_obs_for_sw_gw_interaction(entries_reach[i], 
                                         names_reach[i], 
                                         fy_start, 
                                         fy_end, 
                                         swgw_exch_obs_freqs[i], 
                                         river_segs_reach[reach],
                                         obs_names_reach[i], 
                                         obs_types_reach[i])

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$ Potential obs data $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

def create_obs_for_sim_obs(sim_obs, entries, names, start, end, freq, locs, \
                                     obs_name, obs_type, units, domain):
    '''
    Function to create the necessary observations for sw-gw exchange for the 
    different spatial and temporal periods
    '''
    # end if

    time_series = pd.DataFrame()
    for name in names:
        time_series = pd.concat((time_series, pd.DataFrame({'name': [name] * entries, 
                                    'value': [0.0] * entries, 
                                    'datetime': pd.date_range(start=start, end=end, freq=freq)[1:]})), axis=0)
    time_series.reset_index(inplace=True)
    # Create zero weighted observations for model predictions of interest
    tr_model.observations.set_as_observations(obs_name, 
                                              time_series, 
                                              locs, 
                                              domain=domain, 
                                              obs_type=obs_type, 
                                              units=units, 
                                              weights=0.0, 
                                              real=False)

sim_river_obs = ['ec_sim', 'rn_sim', 'st_sim', 'fl_sim']
units = ['uS/cm', 'Bq/l', 'mAHD', 'm^3/d']
#
for index, sim in enumerate(sim_river_obs):
    names_seg = ["{}{}".format(sim, x) for x in river_segs_seg]
    create_obs_for_sim_obs(len(river_segs_seg),
                           12,
                           names_seg,
                           fy_start,
                           fy_end,
                           'M',
                           river_segs_seg,
                           sim,
                           sim,
                           units[index],
                           'stream')

def create_cell_filter(cache):
    """
    Generates a cell filter that removes duplicate cells based on x and y position.
    Note that this filter assumes that the structure of the passed array is [z, x, y]
    
    :param cache: set object, caching object that has an `add` method.
    
    :returns: function
    """
    def filt(x):
        if (x[1], x[2]) in cache:
            return False

        cache.add((x[1], x[2]))
        return True
    # End filt()
    
    return filt
# End create_filter()

def filter_to_highest_cells(array):
    '''
    Removes cells that are below the 'highest' cells for each x and y position.
    
    :param array: ndarray, array to sort and filter in [z, x, y] format.
    '''
    cache = set()
    filt = create_cell_filter(cache)
    sorted_list = sorted(array, key=lambda x: (x[0], x[1], x[2]))

    return np.asarray(filter(filt, sorted_list))
# End filter_to_highest_cells()

def porous_potential_obs(sim_type_obs_hgu, hgus, fy_start, fy_end, freq, units, domain, skip=0):
    for index, sim_type_hgu in enumerate(sim_type_obs_hgu):
        sim_locs = np.argwhere(tr_model.model_mesh3D[1] == hgus[index])
        sim_locs = filter_to_highest_cells(sim_locs)
        sim_locs = sim_locs[0::skip + 1]
        names_seg = ["{}{}".format(sim_type_hgu, x) for x in range(len(sim_locs))]
        sim_loc_dict = {x[0]:x[1] for x in zip(names_seg, sim_locs)}
        create_obs_for_sim_obs(len(sim_locs),
                               12,
                               names_seg,
                               fy_start,
                               fy_end,
                               freq,
                               sim_loc_dict,
                               sim_type_hgu,
                               sim_type_hgu,
                               units,
                               domain)


hgus = [1, 3, 5, 6]
# Create obs for potential head observations        
sim_heads_obs_hgu = ['shcoon', 'shshep', 'shrenm', 'shcali']
porous_potential_obs(sim_heads_obs_hgu, hgus, fy_start, fy_end, 'M', 'mAHD', 'porous', skip=4)
# Create obs for potential C14 observations        
sim_c14_obs_hgu = ['c14coo', 'c14she', 'c14ren', 'c14cal']
porous_potential_obs(sim_c14_obs_hgu, hgus, fy_start, fy_end, 'M', 'pmc', 'porous', skip=4)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    
## These are the "Tableau 20" colors as RGB.    
#import matplotlib.pyplot as plt
#colors = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
#             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
#             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
#             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
#             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
#  
## Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
#for i in range(len(colors)):    
#    r, g, b = colors[i]    
#    colors[i] = (r / 255., g / 255., b / 255.)    
#    
#flatten = lambda l: [item for sublist in l for item in sublist]
#fig = plt.figure(figsize=(2.5,5))
#for index, reach in enumerate(river_segs_reach[1:]):
#    reach_river = river_seg[river_seg['iseg'].isin(reach)]
#    points = [x for x in reach_river['amalg_riv_points_collection']]
#    points = flatten(points)
#    x_points = [x[0] for x in points]
#    y_points = [y[1] for y in points]    
#    plt.plot(x_points, y_points, color=colors[index], label=str(index + 1))
#plt.legend()
#plt.tight_layout()
#plt.savefig(r"C:\Workspace\part0075\MDB modelling\testbox\PEST5000\master_Colossus3\river_reaches.png")

print "************************************************************************"
print " Creating Campaspe river boundary"

tr_model.boundaries.create_model_boundary_condition('Campaspe River', 
                                                    'river_flow', 
                                                    bc_static=False)
tr_model.boundaries.assign_boundary_array('Campaspe River', 
                                          [reach_data, seg_dict])

Eppalock_EC_ts = river_ec_data[EC_input_site] #river_ec_data[406219]
Eppalock_EC_ts = Eppalock_EC_ts.dropna()
Eppalock_EC_ts = Eppalock_EC_ts.resample('M').mean()
Eppalock_EC_ts_resampled = resample_to_model_data_index(Eppalock_EC_ts, 
                                                        date_index, 
                                                        frequencies, 
                                                        date_group,
                                                        start,
                                                        end)

tr_model.boundaries.create_model_boundary_condition('Eppalock_EC', 
                                                    'river_ec', 
                                                    bc_static=False)
tr_model.boundaries.assign_boundary_array('Eppalock_EC', 
                                          Eppalock_EC_ts_resampled)


print "************************************************************************"
print " Mapping Murray River to grid"

riv, mriver_seg_ghb = \
    rivers.prepare_river_data_for_Murray(tr_model, surface_raster_high_res_GSA,
                                         r"C:\Workspace\part0075\MDB modelling\test_model.shp", #Murray_river_poly_file,
                                         Campaspe_relevant,
                                         river_stage_data,
                                         river_seg) 

print "************************************************************************"
print " Creating Murray River boundary"

tr_model.boundaries.create_model_boundary_condition('Murray River', 'river', bc_static=True)
tr_model.boundaries.assign_boundary_array('Murray River', riv)

print "************************************************************************"
print " Setting up Murray River GHB boundary"

ghb = prepare_ghb_boundary_from_Murray_data(tr_model,
                                            mriver_seg_ghb)

print "************************************************************************"
print " Creating GHB boundary"

tr_model.boundaries.create_model_boundary_condition('GHB', 'general head', bc_static=True)
tr_model.boundaries.assign_boundary_array('GHB', ghb)

print "************************************************************************"
print " Mapping Drains to grid"

drain = prepare_drain_data_for_model(tr_model,
                                 Camp_riv_cells,
                                 start_irrigation,
                                 date_index)

print "************************************************************************"
print " Creating Drains boundary"

tr_model.boundaries.create_model_boundary_condition('Drain', 'drain', bc_static=True)
tr_model.boundaries.assign_boundary_array('Drain', drain)

print "************************************************************************"
print " Mapping Channels to grid"

#channel = prepare_channel_data_for_model(tr_model,
#                                   start_irrigation,
#                                   date_index,
#                                   Camp_riv_cells)

print "************************************************************************"
print " Creating Channel boundary"

#tr_model.boundaries.create_model_boundary_condition('Channel', 'channel', bc_static=True)
#tr_model.boundaries.assign_boundary_array('Channel', channel)

print "************************************************************************"
print " Creating parameters for transport "

# General parameters for transport
for unit in HGU:
    tr_model.parameters.create_model_parameter('por_{}'.format(unit), value=0.25)
    tr_model.parameters.parameter_options('por_{}'.format(unit), 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=0.05, 
                                          PARUBND=0.4, 
                                          PARGP='transp', 
                                          SCALE=1, 
                                          OFFSET=0)

tr_model.parameters.create_model_parameter('disp', value=0.01)
tr_model.parameters.parameter_options('disp', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=1E-5, 
                                      PARUBND=100., 
                                      PARGP='transp', 
                                      SCALE=1, 
                                      OFFSET=0)

# Parameters for the SFT model
tr_model.parameters.create_model_parameter('sfdisp', value=100.)
tr_model.parameters.parameter_options('sfdisp', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=1E-5, 
                                      PARUBND=1000., 
                                      PARGP='ec', 
                                      SCALE=1, 
                                      OFFSET=0)

# Parameters for the Radon model

# Hyporheic zone porosity
tr_model.parameters.create_model_parameter('hzporo', value=0.3)
tr_model.parameters.parameter_options('hzporo', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.1, 
                                      PARUBND=0.4, 
                                      PARGP='radon', 
                                      SCALE=1, 
                                      OFFSET=0)

# Hyporheic zone production of radon
tr_model.parameters.create_model_parameter('hzprod', value=3000.)
tr_model.parameters.parameter_options('hzprod', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=1000., 
                                      PARUBND=10000., 
                                      PARGP='radon', 
                                      SCALE=1, 
                                      OFFSET=0)

# Hyporheic zone residence time
tr_model.parameters.create_model_parameter('hz_rt', value=1.)
tr_model.parameters.parameter_options('hz_rt', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.05, 
                                      PARUBND=5., 
                                      PARGP='radon', 
                                      SCALE=1, 
                                      OFFSET=0)

# Gas transfer velocity
tr_model.parameters.create_model_parameter('gtv', value=1.)
tr_model.parameters.parameter_options('gtv', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.6, 
                                      PARUBND=1.4, 
                                      PARGP='radon', 
                                      SCALE=1, 
                                      OFFSET=0)

# Groundwater concentration of Radon
tr_model.parameters.create_model_parameter('gwconc', value=30000.)
tr_model.parameters.parameter_options('gwconc', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=10000, 
                                      PARUBND=50000, 
                                      PARGP='radon', 
                                      SCALE=1, 
                                      OFFSET=0)

# Hyporheic zone depth         
tr_model.parameters.create_model_parameter_set('hzdpth', value=0.01, \
                                               num_parameters=num_reaches)
tr_model.parameters.parameter_options_set('hzdpth', 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=0., 
                                          PARUBND=1., 
                                          PARGP='radon', 
                                          SCALE=1, 
                                          OFFSET=0)


print "************************************************************************"
print " Collate observations"
#
tr_model.map_obs_loc2mesh3D(method='nearest', ignore=[-1, 7])
tr_model.map_obs2model_times()
tr_model.observations.collate_observations()

print "************************************************************************"
print " Package up groundwater model builder object"

tr_model.package_model()
