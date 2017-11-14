from osgeo import osr
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt


from HydroModelBuilder.GWModelBuilder import GWModelBuilder
from HydroModelBuilder.GISInterface.GDALInterface.GDALInterface import GDALInterface

from CampaspeModel.CustomScripts import Campaspe_data
from CampaspeModel.build_common import Campaspe_mesh

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

Campaspe_data_folder = r"C:\Workspace\part0075\MDB modelling\Campaspe_data"
hu_raster_path = r"C:\Workspace\part0075\MDB modelling\ESRI_GRID_raw\ESRI_GRID"

custom_data = \
    Campaspe_data.process_custom_scripts_and_spatial_data(SS_model, 
                                                          Campaspe_data_folder,
                                                          verbose=True)
HGU_props = custom_data['HGU_props']
rain_gauges = custom_data['rain_gauges']
long_term_historic_rainfall = custom_data['long_term_historic_rainfall'] 
recharge_zones = custom_data['recharge_zones']
surface_raster_high_res = custom_data['surface_raster_high_res'] 
river_gauges = custom_data['river_gauges']
Campaspe_river_poly_file = custom_data['Campaspe_river_poly_file']
Murray_river_poly_file = custom_data['Murray_river_poly_file']
river_flow_data = custom_data['river_flow_data']
river_stage_data = custom_data['river_stage_data']
Campaspe = custom_data['Campaspe']
Campaspe_relevant = custom_data['Campaspe_relevant']

#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************

print '########################################################################'
print '########################################################################'
print '## Mesh specific model building '
print '########################################################################'
print '########################################################################'

HGU, hu_raster_files_reproj = Campaspe_mesh.build_mesh_and_set_properties(SS_model,
                                                  hu_raster_path,
                                                  HGU_props,
                                                  resolution=5000,
                                                  create_basement=True)

print "************************************************************************"
print " Interpolating rainfall data to grid "

interp_rain = SS_model.interpolate_points2mesh(rain_gauges, long_term_historic_rainfall, feature_id='Name', method='linear')
# Adjust rainfall to m from mm and from year to day
interp_rain = interp_rain / 1000.0 / 365.0
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
num_reaches = 20
base_guesses_for_k_bed_x = [0., 0.33, 0.66, 1.0]
base_guesses_for_k_bed_y = [1., 1., 0.01, 0.001]
interp_guesses_for_k_bed_x = np.linspace(0., 1., num_reaches)
interp_guesses_for_k_bed_y = np.interp(interp_guesses_for_k_bed_x, 
                                       base_guesses_for_k_bed_x,
                                       base_guesses_for_k_bed_y) 

SS_model.create_pilot_points('Campaspe', linear=True)
camp_pp = SS_model.pilot_points['Campaspe']
camp_pp.set_uniform_points(river_seg['rchlen'].sum(), num_reaches)

known_points = camp_pp.points

# Define split on river for which unique values will be given to props at 
# those points which will then be interpolated along the length of the river
#for reach in range(num_reaches):
# Setting up river bed hydraulic conductivity values
SS_model.parameters.create_model_parameter_set('kv_riv', 
                                           value=list(interp_guesses_for_k_bed_y), 
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
river_seg['strthick'] = np.interp(river_seg['Cumulative Length'].tolist(), known_points, strthick_val)

# The depths in the column at row j and col i can be obtained using:
# SS_model.model_mesh3D[0][:,0,1]
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
def find_layer(elev, col_vals):
    for index, val in enumerate(col_vals):
        if elev > val:
            if index == 0:
                return index
            else:
                return index - 1
        #end if

        
from operator import itemgetter
from itertools import groupby

river_seg['amalg_riv_points_tuple'] = river_seg['amalg_riv_points'].apply(lambda x: (x[0], x[1]))    

def merge_collocated_stream_reaches(river_segment):
    # Sort out collocated stream reaches to avoid short circuiting:
    
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
    river_seg2 = river_segment.copy()
    
    max_length = 3000.
    merge_row = []
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
            #def loc_tup(row):
            #    return (row['i'], row['j'])
            if prev['amalg_riv_points_tuple'] == nexx['amalg_riv_points_tuple']:
                if curr['rchlen'] < max_length:
                    merge_row += [ind]
                
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
        merge_group = merge_group + [index_list[index_dict[merge_group[-1]] + 1]] 
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
        
        river_seg2.drop(merge_group[1:], inplace=True)
    
    river_seg2.index = range(river_seg2.shape[0])

    ###########################################################################
    # BREAK THIS NEXT BIT INTO NEW FUNCTION
    
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

    #river_seg = river_seg2
    return river_seg2
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

river_seg = merge_collocated_stream_reaches(river_seg)
SS_model.river_mapping['Campaspe'] = river_seg

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
#Campaspe_gauge_zero2 = Campaspe_gauge_zero[Campaspe_gauge_zero['Site Id'] != 406203]
if Campaspe_gauge_zero[Campaspe_gauge_zero['Site Id'] == 406203]['seg_loc'].tolist() == \
   Campaspe_gauge_zero[Campaspe_gauge_zero['Site Id'] == 406218]['seg_loc'].tolist():
    Campaspe_gauge_zero.at[Campaspe_gauge_zero['Site Id'] == 406203, 'seg_loc'] = \
        Campaspe_gauge_zero[Campaspe_gauge_zero['Site Id'] == 406203]['seg_loc'].tolist()[0] + 1
Campaspe_gauge_zero2 = Campaspe_gauge_zero

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
    k = find_layer(strbot, SS_model.model_mesh3D[0][:, j_mesh, i_mesh])
    new_k += [k]
    bottom_layer += [SS_model.model_mesh3D[0][k + 1, j_mesh, i_mesh]] 
    for layer in range(7):
        try:
            surface_layers[layer] += [SS_model.model_mesh3D[0][layer, j_mesh, i_mesh]]
        except:
            surface_layers[layer] = [SS_model.model_mesh3D[0][layer, j_mesh, i_mesh]]

for layer in range(7):
    river_seg["surf{}".format(layer)] = surface_layers[layer]

river_seg['k'] = new_k
river_seg['strtop'] = river_seg['bed_from_gauge']
river_seg['bottom_layer'] = bottom_layer

river_seg.plot(x='Cumulative Length', y=['bed_from_gauge'] + ["surf{}".format(x) for x in range(7)] + ['bottom_layer'])
river_seg.plot(x='Cumulative Length', y=['bed_from_gauge'] + ['bottom_layer'])

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

#SS_model.map_polyline_to_grid(Murray_river_poly)

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
#mriver_seg.loc[:, 'bed_from_gauge'] = np.nan
mriver_seg.loc[:, 'bed_from_gauge'] = mriver_seg['strtop']
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

#values_to_edit = ['bed_from_gauge', 'stage_from_gauge']
values_to_edit = ['stage_from_gauge']
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

mriver_seg['k'] = new_k
mriver_seg['active'] = active
      
# Remove any stream segments for which the elevation could not be mapped to a layer
mriver_seg[mriver_seg['active'] == -1] = np.nan
mriver_seg.dropna(inplace=True)
SS_model.river_mapping['Murray'] = mriver_seg

mriver_seg['strtop'] = mriver_seg['bed_from_gauge']                      
                      
mriver_seg['strhc1'] = SS_model.parameters.param['kv_rm']['PARVAL1']                      

mriver_seg['width1'] = SS_model.parameters.param['rmwdth']['PARVAL1']

mriver_seg['stage'] = mriver_seg['stage_from_gauge'] #+ SS_model.parameters.param['rmstage']['PARVAL1']

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
mriver_seg_ghb = mriver_seg.copy()

#mriver_seg[mriver_seg['cell_used'] == 0] = np.nan
#mriver_seg.dropna(inplace=True)

#mriver_seg['cell_loc_tuple'] = [(x[1]['k'], x[1]['i'], x[1]['j']) for x in mriver_seg.iterrows()]
#mriver_seg['cell_loc_tuple'] = [(x[1]['i'], x[1]['j']) for x in mriver_seg.iterrows()]
#mriver_seg = mriver_seg.groupby(by='cell_loc_tuple').mean()
#mriver_seg.index = range(mriver_seg.shape[0])

mriver_seg['amalg_riv_points_tuple'] = mriver_seg['amalg_riv_points'].apply(lambda x: (x[0], x[1]))    
mriver_seg3 = mriver_seg.iloc[::-1]
mriver_seg3.index = range(mriver_seg3.shape[0])
mriver_seg = merge_collocated_stream_reaches(mriver_seg3)
SS_model.river_mapping['Murray'] = mriver_seg


mriver_seg.plot(x='Cumulative Length', y=['bed_from_gauge'] + ["surf{}".format(x) for x in range(7)])
mriver_seg.plot(x='Cumulative Length', y=['bed_from_gauge', 'bottom_layer', 'strtop', 'surf0'])
mriver_seg.plot(x='Cumulative Length', y=['strtop', 'stage'])

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
for mrow in mriver_seg_ghb.iterrows():
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
        
    MurrayGHBstage = mriver_seg_ghb['stage'].loc[Murray_df_ind2[index]] + SS_model.parameters.param['mghb_stage']['PARVAL1']
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
