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
river_gauges = custom_data['river_gauges']
Campaspe_river_poly_file = custom_data['Campaspe_river_poly_file']
Murray_river_poly_file = custom_data['Murray_river_poly_file']
river_flow_data = custom_data['river_flow_data']
river_stage_data = custom_data['river_stage_data']
river_ec_data = custom_data['river_ec_data']
Campaspe = custom_data['Campaspe']
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
                                                  resolution=5000,
                                                  create_basement=True)

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

# Map bores to layers to create initial head maps for different hydrogeological units
#interp_heads = {}

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
#    
#import matplotlib.pyplot as plt    
#for key in interp_heads:
#    bores_layer_df = pd.DataFrame()
#    bores_layer_df["Easting"] = [x[0] for x in bores_layer] 
#    bores_layer_df["Northing"] = [x[1] for x in bores_layer]
#    bores_layer_df["mean level"] = bores_head_layer
#    (XI, YI) = tr_model.model_mesh_centroids
#    plt.figure()
#    z_min = np.min(interp_heads[key])
#    z_max = np.max(interp_heads[key])
#    plt.pcolor(XI, YI, interp_heads[key], vmin=z_min, vmax=z_max)
#    plt.scatter([x[0] for x in bores_layer], [x[1] for x in bores_layer], 50, bores_head_layer, vmin=z_min, vmax=z_max, cmap="jet")
#    plt.colorbar()
#
#    bores_layer_df.plot(kind='scatter', x="Easting", y="Northing", c="mean level", cmap="Spectral") # , edgecolor='None'
#    plt.scatter(x=[x[0] for x in bores_layer], y=[x[1] for x in bores_layer], c=bores_head_layer)
#import sys
#sys.exit()

# Initalise model with head from elevations
initial_heads_tr = np.full(tr_model.model_mesh3D[1].shape, 0.)

for i in range(len(hu_raster_files_reproj)/2):
    initial_heads_tr[i] = (tr_model.model_mesh3D[0][i]+tr_model.model_mesh3D[0][i+1])/2

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
C14_bore_points3D = C14_bore_points3D.set_index("Bore_id")
C14_bore_points3D.rename(columns={'zone55_easting':'Easting', 'zone55_northing':'Northing'}, inplace=True)

tr_model.observations.set_as_observations('C14', 
                                          C14_obs_time_series, \
                                          C14_bore_points3D, 
                                          domain='porous', \
                                          obs_type='concentration', \
                                          units='pMC', \
                                          weights=1.0/5.0)

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
        if pumping_data.loc[pump, 'Top screen depth (m)'] == 0.: 
            #'No data to place pump at depth ... ignoring '            
            continue
        pump_depth = tr_model.model_mesh3D[0][0][row][col] - pumping_data.loc[pump, 'Top screen depth (m)']        
        active = False
        for i in range(tr_model.model_mesh3D[0].shape[0]-1):
            if pump_depth < tr_model.model_mesh3D[0][i][row][col] and \
                 pump_depth > tr_model.model_mesh3D[0][i+1][row][col]:
                active_layer = i
                active = True
                break
            # end if
        # end for
        if tr_model.model_mesh3D[1][active_layer][row][col] == -1:
            active = False
        # end if
        if active == False: 
            #print 'Well not placed: ', pump            
            continue
        # end if
        # Specify if pump is shallow
        if pump_depth < 25:
            pump_shallow += [True]
        else:
            pump_shallow += [False]
        # end if
        
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
        date_index2 = pd.date_range(start=start_pumping, end=end, freq='M')
        pumping_data_ts = pumping_data_ts.reindex(date_index2)    
        pumping_data_ts = pumping_data_ts.ix[start_pumping:end]
        pumping_data_ts = pumping_data_ts.fillna(0.0)

        resampled_pumping_data_ts = \
            resample_to_model_data_index(pumping_data_ts, date_index, 
                                         frequencies, date_group, start, end,
                                         index_report=False, fill='zero')
        
        # Now fill in the well dictionary with the values of pumping at relevant stress periods
        for index, time in enumerate(resampled_pumping_data_ts.iterrows()):
            if index >= tr_model.model_time.t['steps']: 
                print index
                continue
#            try:
#                wel[wells_start - 1 + index] += [[active_layer, row, col, -time[1][pump]]]
#            except:
#                wel[wells_start - 1 + index] = [[active_layer, row, col, -time[1][pump]]]
            try:
                wel[index] += [[active_layer, row, col, -time[1][pump]]]
            except:
                wel[index] = [[active_layer, row, col, -time[1][pump]]]
                
print "************************************************************************"
print " Creating pumping boundary "

tr_model.boundaries.create_model_boundary_condition('licenced_wells', 'wells', bc_static=True)
tr_model.boundaries.assign_boundary_array('licenced_wells', wel)

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

# Convert flows from Ml/d to m^3/d
for key in river_flow_data.keys():
    river_flow_data[key]['Mean'] = river_flow_data[key]['Mean'] * 1000.


###############
##############
########  This needs to live elsewhere but will do for now ... TODO: Fix this
##############
###############
tr_model.GISInterface.raster_reproject_by_grid(surface_raster_high_res,
                                               surface_raster_high_res[:-4] + '_reproj.tif',
                                               resample_method='min')

surface_raster_high_res = surface_raster_high_res[:-4] + '_reproj.tif'


tr_model.map_points_to_grid(river_gauges, feature_id='Site Name')

Campaspe_river_gauges = tr_model.points_mapped['processed_river_sites_stage_clipped.shp']

filter_gauges = []
for riv_gauge in Campaspe_river_gauges:
    #if riv_gauge[1][0] in use_gauges:
    if str(riv_gauge[1][0]) in use_gauges:
        filter_gauges += [riv_gauge]


tr_model.create_river_dataframe('Campaspe', Campaspe_river_poly_file, surface_raster_high_res)

# Create reach data
river_seg = tr_model.river_mapping['Campaspe']
# Parameters are ordered from upstream to downstream
num_reaches = 4

tr_model.create_pilot_points('Campaspe', linear=True)
camp_pp = tr_model.pilot_points['Campaspe']
camp_pp.set_uniform_points(river_seg['rchlen'].sum(), num_reaches)

known_points = camp_pp.points

# Define split on river for which unique values will be given to props at 
# those points which will then be interpolated along the length of the river

# Setting up river bed hydraulic conductivity values
tr_model.parameters.create_model_parameter_set('kv_riv', 
                                           value=1., 
                                           num_parameters=num_reaches)
tr_model.parameters.parameter_options_set('kv_riv', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.01, 
                                      PARUBND=10.0, 
                                      PARGP='kv_riv', 
                                      SCALE=1, 
                                      OFFSET=0)
# Setting up river bed elevation correction parameter to account for 
# uncertainty in where bed lies relative to zero gauge
tr_model.parameters.create_model_parameter_set('beddep', 
                                           value=0.01, 
                                           num_parameters=num_reaches)
tr_model.parameters.parameter_options_set('beddep', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.001, 
                                      PARUBND=1.0, 
                                      PARGP='rivbed', 
                                      SCALE=1, 
                                      OFFSET=0)
# Setting up river bed roughness values
tr_model.parameters.create_model_parameter_set('mn_riv', 
                                           value=0.01, 
                                           num_parameters=num_reaches)
tr_model.parameters.parameter_options_set('mn_riv', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.001, 
                                      PARUBND=0.1, 
                                      PARGP='rough', 
                                      SCALE=1, 
                                      OFFSET=0)
# Setting up river width values
tr_model.parameters.create_model_parameter_set('rivwdth', 
                                           value=10.0, 
                                           num_parameters=num_reaches)
tr_model.parameters.parameter_options_set('rivwdth', 
                                      PARTRANS='fixed', 
                                      PARCHGLIM='factor', 
                                      PARLBND=4., 
                                      PARUBND=40., 
                                      PARGP='rivwdt', 
                                      SCALE=1, 
                                      OFFSET=0)
# Setting up riverbed thickness values
tr_model.parameters.create_model_parameter_set('bedthck', 
                                           value=0.10, 
                                           num_parameters=num_reaches)
tr_model.parameters.parameter_options_set('bedthck', 
                                      PARTRANS='fixed', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.01, 
                                      PARUBND=1., 
                                      PARGP='bedthk', 
                                      SCALE=1, 
                                      OFFSET=0)


strcond_val = [tr_model.parameters.param['kv_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
river_seg['strhc1'] = np.interp(river_seg['Cumulative Length'].tolist(), 
                                known_points, strcond_val)
strthick_val = [tr_model.parameters.param['bedthck{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
river_seg['strthick'] = np.interp(river_seg['Cumulative Length'].tolist(), known_points, strthick_val)

amalg_riv_points = []
for row in river_seg[['i', 'j']].iterrows():
    amalg_riv_points += [[row[1]['j'], row[1]['i']]]

# The depths in the column at row j and col i can be obtained using:
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

min_length = 500.
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
            if curr['rchlen'] < min_length:
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
tr_model.river_mapping['Campaspe'] = river_seg
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

river_seg.loc[already_defined, 'strhc1'] = 0.0

new_k = []

for row in river_seg.iterrows():
    j_mesh = row[1]['i'] 
    i_mesh = row[1]['j']
    strtop = row[1]['strtop']
    strbot = row[1]['strtop'] - row[1]['strthick'] 
    new_k += [find_layer(strbot, tr_model.model_mesh3D[0][:, j_mesh, i_mesh])]

river_seg['k'] = new_k
       
# Remove any stream segments for which the elevation could not be mapped to a layer
#river_seg.dropna(inplace=True)

river_seg['ireach'] = 1
river_seg['iseg'] = [x + 1 for x in range(river_seg.shape[0])]

FieldData_info = FieldData.groupby('Name').first()[['Easting', 'Northing', 'Distance_Eppalock', 'Zone']]
                      
# Set up bed elevations based on the gauge zero levels:
gauge_points = [x for x in zip(Campaspe.Easting, Campaspe.Northing)]
field_points = [x for x in zip(FieldData_info.Easting, FieldData_info.Northing)]

river_gauge_seg = tr_model.get_closest_riv_segments('Campaspe', gauge_points)
river_field_seg = tr_model.get_closest_riv_segments('Campaspe', field_points)
river_seg.loc[:, 'bed_from_gauge'] = np.nan

Campaspe['new_gauge'] = Campaspe[['Gauge Zero (Ahd)', 'Cease to flow level', 'Min value']].max(axis=1) 
Campaspe['seg_loc'] = river_gauge_seg 
FieldData_info['seg_loc'] = river_field_seg        
Campaspe_gauge_zero = Campaspe[Campaspe['new_gauge'] > 10.]
# There are two values at the Campaspe weir, while it would be ideal to split the
# reach here it will cause problems for the segment
Campaspe_gauge_zero2 = Campaspe_gauge_zero[Campaspe_gauge_zero['Site Id'] != 406203]

river_seg.loc[river_seg['iseg'].isin(Campaspe_gauge_zero2['seg_loc'].tolist()), 'bed_from_gauge'] = sorted(Campaspe_gauge_zero2['new_gauge'].tolist(), reverse=True)
river_seg['bed_from_gauge'] = river_seg.set_index(river_seg['Cumulative Length'])['bed_from_gauge'].interpolate(method='values', limit_direction='both').tolist()
river_seg['bed_from_gauge'] = river_seg['bed_from_gauge'].bfill()

new_k = []
surface_layers = {}
bottom_layer = []
for row in river_seg.iterrows():
    j_mesh = row[1]['i'] 
    i_mesh = row[1]['j']
    strbot = row[1]['bed_from_gauge'] - row[1]['strthick']
    new_k += [find_layer(strbot, tr_model.model_mesh3D[0][:, j_mesh, i_mesh])]
    k = find_layer(strbot, tr_model.model_mesh3D[0][:, j_mesh, i_mesh])
    bottom_layer += [tr_model.model_mesh3D[0][k + 1, j_mesh, i_mesh]] 
    for layer in range(7):
        try:
            surface_layers[layer] += [tr_model.model_mesh3D[0][layer, j_mesh, i_mesh]]
        except:
            surface_layers[layer] = [tr_model.model_mesh3D[0][layer, j_mesh, i_mesh]]

for layer in range(7):
    river_seg["surf{}".format(layer)] = surface_layers[layer]

river_seg['bottom_layer'] = bottom_layer

river_seg.plot(x='Cumulative Length', y=['bed_from_gauge'] + ["surf{}".format(x) for x in range(7)])
river_seg.plot(x='Cumulative Length', y=['bed_from_gauge', 'bottom_layer'])


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

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@ Create temporal segment data @@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

nseg = river_seg['iseg'].tolist()
icalc = [1] * len(nseg)
outseg = river_seg['iseg'] + 1
outseg = outseg.tolist()
outseg[-1] = 0
iupseg = [0] * len(nseg)
iprior = [0] * len(nseg)
nstrpts = [0] * len(nseg)
# Inflows to the system ... ignoring tribs and only considering Eppalock
Eppalock_inflow = river_flow_data[406207]
Eppalock_inflow_resampled = resample_to_model_data_index(Eppalock_inflow, 
                                                         date_index, 
                                                         frequencies, 
                                                         date_group,
                                                         start,
                                                         end,
                                                         df_freq='D',
                                                         fill='stats',
                                                         stat='25%')
flow = [0] * len(nseg)
runoff = [0] * len(nseg)
etsw = [0] * len(nseg)
pptsw = [0] * len(nseg)

# Set the roughness for the channel
roughch_val = [tr_model.parameters.param['mn_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
roughch = np.interp(river_seg['Cumulative Length'].tolist(), 
                                known_points, roughch_val)
# Set the roughness for the banks
roughbk_val = [tr_model.parameters.param['mn_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
roughbk = np.interp(river_seg['Cumulative Length'].tolist(), 
                                known_points, roughbk_val)
river_seg['roughch'] = roughch
river_seg['roughbk'] = roughbk

cdpth = [0] * len(nseg)
fdpth = [0] * len(nseg)
awdth = [0] * len(nseg)
bwdth = [0] * len(nseg)

width1_val = [tr_model.parameters.param['rivwdth{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
width1 = width2 = np.interp(river_seg['Cumulative Length'].tolist(), 
                                known_points, width1_val)
river_seg['width2'] = river_seg['width1'] = width1

segment_data = {}
segment_data1 = {}

# Prepare the inflow data from Lake Eppalock by resampling it to fit the model
# stress periods and backfilling missing data with the 25th percentile of flow

seg_dict = {}
for per, date in enumerate(date_index[:-1]):

    flow[0] = Eppalock_inflow_resampled.loc[date]['Mean']
    etsw = [interp_et[per][x[1]['i'], x[1]['j']] for x in river_seg.iterrows()] 
    pptsw = [interp_rain[per][x[1]['i'], x[1]['j']] for x in river_seg.iterrows()] 
        
    segment_data[per] = pd.DataFrame({'nseg':nseg, 'icalc':icalc, 'outseg':outseg, 'iupseg':iupseg, 'iprior':iprior, 'nstrpts':nstrpts, \
                                 'flow':flow, 'runoff':runoff, 'etsw':etsw, 'pptsw':pptsw, 'roughch':roughch, 'roughbk':roughbk, \
                                 'cdpth':cdpth, 'fdpth':fdpth, 'awdth':awdth, 'bwdth':bwdth, 'width1':width1, 'width2':width2})
    cols_ordered = ['nseg', 'icalc', 'outseg', 'iupseg', 'iprior', 'nstrpts', \
                    'flow', 'runoff', 'etsw', 'pptsw', 'roughch', 'roughbk', \
                    'cdpth', 'fdpth', 'awdth', 'bwdth', 'width1', 'width2']
    segment_data[per] = segment_data[per][cols_ordered]
    segment_data1[per] = segment_data[per].to_records(index=False)

tr_model.save_MODFLOW_SFR_dataframes('Campaspe', reach_df, segment_data)

tr_model.river_mapping['Campaspe'] = river_seg
                      
seg_dict = segment_data1

Camp_riv_cells = [x for x in zip(river_seg['i'], river_seg['j'])]
                      
print "************************************************************************"
print " Creating Campaspe river observations for stage and discharge at "
print " locations downstream of Lake Eppalock"

Campaspe_info = Campaspe
Campaspe_info.index = Campaspe_info['Site Id']
Campaspe_info = Campaspe_info[['Easting', 'Northing', 'Site Name', 'seg_loc']]

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
for key in river_ec_data.keys():
    # Ignore 406219 Lake Eppalock    
    if key == 406219:
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
entries = [1, 4, 12]
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
# TODO: Implement spatial scales too ...
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
river_segs_seg = [[x] for x in river_seg['iseg']]
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
    
#swgw_exch_obs_spatial = 3 # whole of river, long reach, cell reach
#
#for swgw_exch_obs_freq in swgw_exch_obs_freqs:
#    for swgw_exch_obs_space in swgw_exch_obs_spatial:
#        create_obs_for_sw_gw_interaction

# These are the "Tableau 20" colors as RGB.    
import matplotlib.pyplot as plt
colors = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(colors)):    
    r, g, b = colors[i]    
    colors[i] = (r / 255., g / 255., b / 255.)    
    
fig = plt.figure(figsize=(2.5,5))
for index, reach in enumerate(river_segs_reach[1:]):
    reach_river = river_seg[river_seg['iseg'].isin(reach)]
    points = [x for x in reach_river['amalg_riv_points_collection']]
    points = flatten(points)
    x_points = [x[0] for x in points]
    y_points = [y[1] for y in points]    
    plt.plot(x_points, y_points, color=colors[index], label=str(index + 1))
plt.legend()
plt.tight_layout()
plt.savefig(r"C:\Workspace\part0075\MDB modelling\testbox\PEST5000\master_Colossus3\river_reaches.png")

print "************************************************************************"
print " Creating Campaspe river boundary"

tr_model.boundaries.create_model_boundary_condition('Campaspe River', 
                                                    'river_flow', 
                                                    bc_static=False)
tr_model.boundaries.assign_boundary_array('Campaspe River', 
                                          [reach_data, seg_dict])

Eppalock_EC_ts = river_ec_data[406219]
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

#tr_model.map_polyline_to_grid(Murray_river_poly)

# Parameter to modify the stage, thus accounting for errors in values specified for stage
tr_model.parameters.create_model_parameter('rmstage', value=0.01)
tr_model.parameters.parameter_options('rmstage', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.001, 
                                      PARUBND=0.1, 
                                      PARGP='murr_riv', 
                                      SCALE=1, 
                                      OFFSET=0)
# Parameter to all shifting the location of the bed which is only estimated based on assumed depth below zero gauge
tr_model.parameters.create_model_parameter('rmbed', value=0.01)
tr_model.parameters.parameter_options('rmbed', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.001, 
                                      PARUBND=0.1, 
                                      PARGP='murr_riv', 
                                      SCALE=1, 
                                      OFFSET=0)
# Parameter for River Murray bed thickness
tr_model.parameters.create_model_parameter('rmbedthk', value=0.01)
tr_model.parameters.parameter_options('rmbedthk', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.001, 
                                      PARUBND=0.5, 
                                      PARGP='murr_riv', 
                                      SCALE=1, 
                                      OFFSET=0)
# Parameter for the vertical hydraulic conductivity of the River Murray
tr_model.parameters.create_model_parameter('kv_rm', value=5E-3)
tr_model.parameters.parameter_options('kv_rm', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=1E-8, 
                                      PARUBND=20, 
                                      PARGP='murr_riv', 
                                      SCALE=1, 
                                      OFFSET=0)
# Parameter for the width of the River Murray
tr_model.parameters.create_model_parameter('rmwdth', value=30)
tr_model.parameters.parameter_options('rmwdth', 
                                      PARTRANS='fixed', 
                                      PARCHGLIM='factor', 
                                      PARLBND=20, 
                                      PARUBND=50, 
                                      PARGP='murr_riv', 
                                      SCALE=1, 
                                      OFFSET=0)

tr_model.create_river_dataframe('Murray', Murray_river_poly_file, surface_raster_high_res)

mriver_seg = tr_model.river_mapping['Murray']
mriver_seg['strthick'] = tr_model.parameters.param['rmbedthk']['PARVAL1']

# Set up bed elevations based on the gauge zero levels:
Murray = Campaspe_relevant[Campaspe_relevant['Site Name'].str.contains("MURRAY RIVER")]
gauge_points = [x for x in zip(Murray.Easting, Murray.Northing)]
mriver_gauge_seg = tr_model.get_closest_riv_segments('Murray', gauge_points)
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
    k = find_layer(strtop, tr_model.model_mesh3D[0][:, j_mesh, i_mesh])
    new_k += [k]
    bottom_layer += [tr_model.model_mesh3D[0][k+1, j_mesh, i_mesh]] 
    active += [tr_model.model_mesh3D[1][k, j_mesh, i_mesh]]
    for layer in range(7):
        try:
            surface_layers[layer] += [tr_model.model_mesh3D[0][layer, j_mesh, i_mesh]]
        except:
            surface_layers[layer] = [tr_model.model_mesh3D[0][layer, j_mesh, i_mesh]]

for layer in range(7):
    mriver_seg["surf{}".format(layer)] = surface_layers[layer]

mriver_seg['bottom_layer'] = bottom_layer

mriver_seg.plot(x='Cumulative Length', y=['bed_from_gauge'] + ["surf{}".format(x) for x in range(7)])
mriver_seg.plot(x='Cumulative Length', y=['bed_from_gauge', 'bottom_layer'])

mriver_seg['k'] = new_k
mriver_seg['active'] = active
      
# Remove any stream segments for which the elevation could not be mapped to a layer
mriver_seg[mriver_seg['active'] == -1] = np.nan
mriver_seg.dropna(inplace=True)
tr_model.river_mapping['Murray'] = mriver_seg

mriver_seg['strtop'] = mriver_seg['stage_from_gauge']                      
                      
mriver_seg['strhc1'] = tr_model.parameters.param['kv_rm']['PARVAL1']                      

mriver_seg['width1'] = tr_model.parameters.param['rmwdth']['PARVAL1']

mriver_seg['stage'] = mriver_seg['strtop'] + tr_model.parameters.param['rmstage']['PARVAL1']

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

mriver_seg['cell_loc_tuple'] = [(x[1]['k'], x[1]['i'], x[1]['j']) for x in mriver_seg.iterrows()]
mriver_seg = mriver_seg.groupby(by='cell_loc_tuple').mean()
mriver_seg.index = range(mriver_seg.shape[0])

tr_model.river_mapping['Murray'] = mriver_seg


simple_river = []
for row in mriver_seg.iterrows():
    row = row[1]
    simple_river += [[row['k'], row['i'], row['j'], row['stage'], \
                      row['strhc1'] * row['rchlen'] * row['width1'], \
                      row['strtop']]]

riv = {}
riv[0] = simple_river
   
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
                                      PARLBND=-20.0, 
                                      PARUBND=50, 
                                      PARGP='ghb', 
                                      SCALE=1, 
                                      OFFSET=0)
tr_model.parameters.create_model_parameter('mghbk', value=10)
tr_model.parameters.parameter_options('mghbk', 
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
    for lay in range(tr_model.model_mesh3D[1].shape[0]):    
        if [lay, row, col] in checked:
            continue
        checked += [lay, row, col]
        if tr_model.model_mesh3D[1][0][row][col] == -1:
            continue
        MurrayGHBstage = mrow['stage'] + tr_model.parameters.param['mghb_stage']['PARVAL1']
        #if MurrayGHBstage < tr_model.model_mesh3D[0][lay][row][col]:
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
zone = tr_model.model_mesh3D[1]
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

#for MurrayGHB_cell in tr_model.polyline_mapped['River_Murray_model.shp']:
for index, MurrayGHB_cell in enumerate(Final_MurrayGHB_cells):

    lay = MurrayGHB_cell[0]
    row = MurrayGHB_cell[1]
    col = MurrayGHB_cell[2]
        
    MurrayGHBstage = mriver_seg['stage'].loc[Murray_df_ind2[index]] + tr_model.parameters.param['mghb_stage']['PARVAL1']
    #if MurrayGHBstage < tr_model.model_mesh3D[0][0][row][col]:
    #    continue
    dx = tr_model.gridHeight
    dz = tr_model.model_mesh3D[0][lay][row][col] - tr_model.model_mesh3D[0][lay + 1][row][col]
    MGHBconductance = dx * dz * tr_model.parameters.param['mghbk']['PARVAL1']
    MurrayGHB += [[lay, row, col, MurrayGHBstage, MGHBconductance]]

ghb = {}
ghb[0] = MurrayGHB

print "************************************************************************"
print " Creating GHB boundary"

tr_model.boundaries.create_model_boundary_condition('GHB', 'general head', bc_static=True)
tr_model.boundaries.assign_boundary_array('GHB', ghb)

print "************************************************************************"
print " Mapping Drains to grid"

drain_poly = tr_model.read_poly("Drain_Clip.shp", path=r"C:\Workspace\part0075\MDB modelling\testbox\input_data\SW\\") 
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
    if (row, col) in Camp_riv_cells:
        continue
    #print tr_model.model_mesh3D
    drain_bed = tr_model.model_mesh3D[0][0][row][col] - tr_model.parameters.param['drain_drop']['PARVAL1']
    drain_cond = drain_cell[1] * drain_width_avg * tr_model.parameters.param['kv_drain']['PARVAL1'] / drain_bed_thickness
    simple_drain += [[0, row, col, drain_bed, drain_cond]]

drain_start = findInterval(start_irrigation, date_index)
drain = {}
drain[drain_start + 1] = simple_drain

print "************************************************************************"
print " Creating Drains boundary"

tr_model.boundaries.create_model_boundary_condition('Drain', 'drain', bc_static=True)
tr_model.boundaries.assign_boundary_array('Drain', drain)


print "************************************************************************"
print " Mapping Channels to grid"

channel_poly = tr_model.read_poly("Channel_Clip.shp", path=r"C:\Workspace\part0075\MDB modelling\testbox\input_data\SW\\") 
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
channel_width_avg = 5.0 #m
channel_bed_thickness = 0.10 #m
for channel_cell in tr_model.polyline_mapped['Channel_Clip_model.shp']:
    row = channel_cell[0][0]
    col = channel_cell[0][1]
    if tr_model.model_mesh3D[1][0][row][col] == -1:
        continue
    if (row, col) in Camp_riv_cells:
        continue
    channel_stage = tr_model.model_mesh3D[0][0][row][col]
    channel_bed = tr_model.model_mesh3D[0][0][row][col] - tr_model.parameters.param['chan_drop']['PARVAL1']
    channel_cond = channel_cell[1] * channel_width_avg * tr_model.parameters.param['kv_chan']['PARVAL1'] / channel_bed_thickness
    simple_channel += [[0, row, col, channel_stage, channel_cond, channel_bed]]

channel_start = findInterval(start_irrigation, date_index)

channel = {}
channel[channel_start + 1] = simple_channel

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

# Parameters for the SFT model
tr_model.parameters.create_model_parameter('sfdisp', value=100)
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
tr_model.parameters.create_model_parameter('hz_poro', value=0.3)
tr_model.parameters.parameter_options('hz_poro', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=0.1, 
                                      PARUBND=0.4, 
                                      PARGP='radon', 
                                      SCALE=1, 
                                      OFFSET=0)

# Hyporheic zone production of radon
tr_model.parameters.create_model_parameter('hz_prod', value=3000.)
tr_model.parameters.parameter_options('hz_prod', 
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
tr_model.parameters.create_model_parameter('gw_conc', value=30000.)
tr_model.parameters.parameter_options('gw_conc', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=10000, 
                                      PARUBND=50000, 
                                      PARGP='radon', 
                                      SCALE=1, 
                                      OFFSET=0)

# Hyporheic zone depth         
tr_model.parameters.create_model_parameter_set('hz_dpth', value=0.01, \
                                               num_parameters=num_reaches)
tr_model.parameters.parameter_options_set('hz_dpth', 
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
