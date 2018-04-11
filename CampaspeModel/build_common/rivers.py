import numpy as np
import pandas as pd

from CampaspeModel.build_utils import river_df_tools
from CampaspeModel.build_utils.multifrequency_resampling import resample_to_model_data_index

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
def find_layer(elev, col_vals):
    for index, val in enumerate(col_vals):
        if elev > val:
            if index == 0:
                return index
            else:
                return index - 1
            # end if
        # end if
    # end for

def prepare_river_data_for_campaspe(model_builder_object, 
                                    surface_raster_high_res,
                                    river_gauges,
                                    campaspe_river_poly_file,
                                    campaspe,
                                    campaspe_field_elevations,
                                    num_reaches=4,
                                    plot=False,
                                    riv_boundary_only=False):
    
    mbo = model_builder_object
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
    


    mbo.GISInterface.raster_reproject_by_grid(surface_raster_high_res,
                                              surface_raster_high_res[:-4] + '_reproj.tif',
                                              resample_method='min')
    
    surface_raster_high_res = surface_raster_high_res[:-4] + '_reproj.tif'
    
    
    mbo.map_points_to_grid(river_gauges, feature_id='Site Name')
    #SS_model.map_points_to_grid(river_gauges, feature_id='Site_ID')
    
    Campaspe_river_gauges = mbo.points_mapped['processed_river_sites_stage_clipped.shp']
    
    filter_gauges = []
    for riv_gauge in Campaspe_river_gauges:
        #if riv_gauge[1][0] in use_gauges:
        if str(riv_gauge[1][0]) in use_gauges:
            filter_gauges += [riv_gauge]
    
    
    mbo.create_river_dataframe('Campaspe', campaspe_river_poly_file, surface_raster_high_res)
    
    # Create reach data
    river_seg = mbo.river_mapping['Campaspe']
    # Parameters are ordered from upstream to downstream
    num_reaches = num_reaches
    base_guesses_for_k_bed_x = [0.0, 0.33, 0.66, 1.0]
    base_guesses_for_k_bed_y = [1., 1., 0.01, 0.001]
    interp_guesses_for_k_bed_x = np.linspace(0., 1., num_reaches)
    interp_guesses_for_k_bed_y = np.interp(interp_guesses_for_k_bed_x, 
                                           base_guesses_for_k_bed_x,
                                           base_guesses_for_k_bed_y) 
    
    mbo.create_pilot_points('Campaspe', linear=True)
    camp_pp = mbo.pilot_points['Campaspe']
    camp_pp.set_uniform_points(river_seg['rchlen'].sum(), num_reaches)
    
    known_points = camp_pp.points
    
    # Define split on river for which unique values will be given to props at 
    # those points which will then be interpolated along the length of the river
    #for reach in range(num_reaches):
    # Setting up river bed hydraulic conductivity values
    mbo.parameters.create_model_parameter_set('kv_riv', 
                                               value=list(interp_guesses_for_k_bed_y), 
                                               num_parameters=num_reaches)
    
    mbo.parameters.parameter_options_set('kv_riv', 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=0.0001, 
                                          PARUBND=10.0, 
                                          PARGP='kv_riv', 
                                          SCALE=1, 
                                          OFFSET=0)

    # Setting up river bed elevation correction parameter to account for 
    # uncertainty in where bed lies relative to zero gauge
    mbo.parameters.create_model_parameter_set('beddep', 
                                               value=0.01, 
                                               num_parameters=num_reaches)
    mbo.parameters.parameter_options_set('beddep', 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=0.001, 
                                          PARUBND=1.0, 
                                          PARGP='rivbed', 
                                          SCALE=1, 
                                          OFFSET=0)
    if not riv_boundary_only:    
        # Setting up river bed roughness values
        mbo.parameters.create_model_parameter_set('mn_riv', 
                                                   value=0.001, 
                                                   num_parameters=num_reaches)
        mbo.parameters.parameter_options_set('mn_riv', 
                                              PARTRANS='log', 
                                              PARCHGLIM='factor', 
                                              PARLBND=0.001, 
                                              PARUBND=0.1, 
                                              PARGP='rough', 
                                              SCALE=1, 
                                              OFFSET=0)
    # Setting up river width values
    mbo.parameters.create_model_parameter_set('rivwdth', 
                                               value=20.0, 
                                               num_parameters=num_reaches)
    
    mbo.parameters.parameter_options_set('rivwdth', 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=4., 
                                          PARUBND=40., 
                                          PARGP='rivwdt', 
                                          SCALE=1, 
                                          OFFSET=0)
    # Setting up riverbed thickness values
    mbo.parameters.create_model_parameter_set('bedthck', 
                                               value=0.10, 
                                               num_parameters=num_reaches)
    
    mbo.parameters.parameter_options_set('bedthck', 
                                          PARTRANS='fixed', 
                                          PARCHGLIM='factor', 
                                          PARLBND=0.01, 
                                          PARUBND=1., 
                                          PARGP='bedthk', 
                                          SCALE=1, 
                                          OFFSET=0)
    
    
    strcond_val = [mbo.parameters.param['kv_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    river_seg['strhc1'] = np.interp(river_seg['Cumulative Length'].tolist(), 
                                    known_points, strcond_val)
    strthick_val = [mbo.parameters.param['bedthck{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    river_seg['strthick'] = np.interp(river_seg['Cumulative Length'].tolist(), known_points, strthick_val)
    
    river_seg['amalg_riv_points_tuple'] = river_seg['amalg_riv_points'].apply(lambda x: (x[0], x[1]))    
    
    river_seg = river_df_tools.merge_collocated_stream_reaches(river_seg, max_length=500.)
    river_seg = river_df_tools.merge_very_short_stream_reaches(river_seg, min_length=289.6) #200.)
    
    mbo.river_mapping['Campaspe'] = river_seg
    
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
        new_k += [find_layer(strbot, mbo.model_mesh3D[0][:, j_mesh, i_mesh])]
    
    river_seg['k'] = new_k
           
    # Remove any stream segments for which the elevation could not be mapped to a layer
    #river_seg.dropna(inplace=True)
    
    river_seg['ireach'] = 1
    river_seg['iseg'] = [x + 1 for x in range(river_seg.shape[0])]
    
    campaspe_field_elevations = campaspe_field_elevations[campaspe_field_elevations.index != 'TribB']
    # Set up bed elevations based on the gauge zero levels:
    gauge_points = [x for x in zip(campaspe.Easting, campaspe.Northing)]
    field_points = [x for x in zip(campaspe_field_elevations.Easting, campaspe_field_elevations.Northing)]
    river_gauge_seg = mbo.get_closest_riv_segments('Campaspe', gauge_points)
    river_field_seg = mbo.get_closest_riv_segments('Campaspe', field_points)
    river_seg.loc[:, 'bed_from_gauge'] = np.nan
    
    campaspe['new_gauge'] = campaspe[['Gauge Zero (Ahd)', 'Cease to flow level', 'Min value']].max(axis=1) 
    campaspe['seg_loc'] = river_gauge_seg         
    campaspe_field_elevations['seg_loc'] = river_field_seg         
    campaspe_gauge_zero = campaspe[campaspe['new_gauge'] > 10.]
    
    # There are two values at the Campaspe weir, while it would be ideal to split the
    # reach here it will cause problems for the segment
    #Campaspe_gauge_zero2 = Campaspe_gauge_zero[Campaspe_gauge_zero['Site Id'] != 406203]
    if campaspe_gauge_zero[campaspe_gauge_zero['Site Id'] == 406203]['seg_loc'].tolist() == \
       campaspe_gauge_zero[campaspe_gauge_zero['Site Id'] == 406218]['seg_loc'].tolist():
        campaspe_gauge_zero.at[campaspe_gauge_zero['Site Id'] == 406203, 'seg_loc'] = \
            campaspe_gauge_zero[campaspe_gauge_zero['Site Id'] == 406203]['seg_loc'].tolist()[0] + 1
    campaspe_gauge_zero2 = campaspe_gauge_zero
    
    if len(campaspe_field_elevations['seg_loc'].unique()) < len(campaspe_field_elevations['seg_loc']):
        campaspe_field_elevations = campaspe_field_elevations.drop_duplicates(subset='seg_loc')
    
    river_seg.loc[river_seg['iseg'].isin(campaspe_gauge_zero2['seg_loc'].tolist()), 'bed_from_gauge'] = sorted(campaspe_gauge_zero2['new_gauge'].tolist(), reverse=True)
    river_seg.loc[river_seg['iseg'].isin(campaspe_field_elevations['seg_loc'].tolist()), 'bed_from_gauge'] = sorted(campaspe_field_elevations['Elevation'].tolist(), reverse=True)

    river_seg['bed_from_gauge'] = river_seg.set_index(river_seg['Cumulative Length'])['bed_from_gauge'].interpolate(method='values', limit_direction='both').tolist()
    river_seg['bed_from_gauge'] = river_seg['bed_from_gauge'].bfill() 
    
    new_k = []
    surface_layers = {}
    bottom_layer = []
    for row in river_seg.iterrows():
        j_mesh = row[1]['i'] 
        i_mesh = row[1]['j']
        strbot = row[1]['bed_from_gauge'] - row[1]['strthick']
        k = find_layer(strbot, mbo.model_mesh3D[0][:, j_mesh, i_mesh])
        new_k += [k]
        bottom_layer += [mbo.model_mesh3D[0][k + 1, j_mesh, i_mesh]] 
        for layer in range(7):
            try:
                surface_layers[layer] += [mbo.model_mesh3D[0][layer, j_mesh, i_mesh]]
            except:
                surface_layers[layer] = [mbo.model_mesh3D[0][layer, j_mesh, i_mesh]]
    
    for layer in range(7):
        river_seg["surf{}".format(layer)] = surface_layers[layer]

    river_seg.loc[:, 'k'] = new_k

    if not riv_boundary_only:
        river_seg.loc[:, 'strtop'] = river_seg['bed_from_gauge']
        river_seg.loc[:, 'bottom_layer'] = bottom_layer
    else:
        river_seg.loc[:, 'strhc1'] = mbo.parameters.param['kv_riv']['PARVAL1']                      
        river_seg.loc[:, 'width1'] = mbo.parameters.param['rivwdth']['PARVAL1']
        river_seg.loc[:, 'stage'] = river_seg['stage_from_gauge']

    
    if plot:
        river_seg.plot(x='Cumulative Length', y=['bed_from_gauge'] + ["surf{}".format(x) for x in range(7)] + ['bottom_layer'])
        river_seg.plot(x='Cumulative Length', y=['bed_from_gauge'] + ['bottom_layer'])
    
    # For stream reaches that didn't map properly to the mesh for z elevation we 
    # can still include by setting to layer 0 with a bed hydraulic conductivity of 0
    #inds = np.where(river_seg['k'].isnull())[0]
    #river_seg['strhc1'].loc[inds] = 0.0
    #river_seg['k'].loc[inds] = 0
    river_seg.dropna(inplace=True)
              
    river_seg['iseg'] = [x + 1 for x in range(river_seg.shape[0])]
             
    def slope_corrector(x):
        if  x  < 0.0001:
            return 0.0001
        elif x > 1E5:
            return 0.0001
        else:
            return x
        # end if
        
    river_seg['slope'] = river_seg['slope'].apply(lambda x: slope_corrector(x))
    
    reach_df = river_seg[['k','i','j','iseg','ireach','rchlen','strtop','slope','strthick','strhc1']]
    reach_data = reach_df.to_records(index=False)

    return river_seg, reach_df, reach_data, known_points
    
def create_segment_data(model_builder_object, river_seg, river_flow_data):
    mbo = model_builder_object
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
    
    num_reaches = mbo.pilot_points['Campaspe'].num_points
    known_points = mbo.pilot_points['Campaspe'].points
    # Set the roughness for the channel
    roughch_val = [mbo.parameters.param['mn_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    roughch = np.interp(river_seg['Cumulative Length'].tolist(), 
                                    known_points, roughch_val)
    # Set the roughness for the banks
    roughbk_val = [mbo.parameters.param['mn_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    roughbk = np.interp(river_seg['Cumulative Length'].tolist(), 
                                    known_points, roughbk_val)
    river_seg.loc[:, 'roughch'] = roughch
    river_seg.loc[:, 'roughbk'] = roughbk
    
    cdpth = [0] * len(nseg)
    fdpth = [0] * len(nseg)
    awdth = [0] * len(nseg)
    bwdth = [0] * len(nseg)
    
    width1_val = [mbo.parameters.param['rivwdth{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    width1 = width2 = np.interp(river_seg['Cumulative Length'].tolist(), 
                                    known_points, width1_val)
    river_seg.loc[:, 'width2'] = river_seg.loc[:, 'width1'] = width1
    
    segment_data = pd.DataFrame({'nseg':nseg, 'icalc':icalc, 'outseg':outseg, 'iupseg':iupseg, 'iprior':iprior, 'nstrpts':nstrpts, \
                                 'flow':flow, 'runoff':runoff, 'etsw':etsw, 'pptsw':pptsw, 'roughch':roughch, 'roughbk':roughbk, \
                                 'cdpth':cdpth, 'fdpth':fdpth, 'awdth':awdth, 'bwdth':bwdth, 'width1':width1, 'width2':width2})
    cols_ordered = ['nseg', 'icalc', 'outseg', 'iupseg', 'iprior', 'nstrpts', \
                    'flow', 'runoff', 'etsw', 'pptsw', 'roughch', 'roughbk', \
                    'cdpth', 'fdpth', 'awdth', 'bwdth', 'width1', 'width2']
    segment_data = segment_data[cols_ordered]
    
    segment_data1 = segment_data.to_records(index=False)
    seg_dict = {0: segment_data1}

    return segment_data, seg_dict

def create_segment_data_transient(model_builder_object,
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
                                  end,
                                  include_et=True,
                                  include_extractions=True):
    
    mbo = model_builder_object
    num_reaches = mbo.pilot_points['Campaspe'].num_points
    known_points = mbo.pilot_points['Campaspe'].points
    
    # Convert flows from Ml/d to m^3/d
    for key in river_flow_data.keys():
        river_flow_data[key]['Mean'] = river_flow_data[key]['Mean'] * 1000.

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

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #@@ Diversions from the Campaspe @@
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    # There are three types of diversion, input to the Campaspe for which there
    # is data.
    #
    # 1. Private diversions which can be broken into three reaches:
    #       a. Eppalock to Campaspe Weir (406207 to 406203)
    #       b. Campaspe Weir to Waranga Western Channel (406203 to 406202)
    #       c. Waranga Western Channel to Murray River (406202 to end)
    #
    # 2. Diversion to channel 1 and 2 just upstream of the Campaspe Weir @ 406218
    #
    # 3. Pumping into and out of the Waranga Western Channel @ 406202
    #
    
    rdd = river_diversion_data
    rdd_resampled = {}
    rdd_resampled['CID Diversions'] = resample_to_model_data_index(rdd['CID Diversions'].sum(axis=1),
                                                      date_index, 
                                                      frequencies, 
                                                      date_group,
                                                      start,
                                                      end,
                                                      df_freq='M',
                                                      fill='zero')

    rdd_resampled['Campaspe River Use'] = resample_to_model_data_index(rdd['Campaspe River Use'],
                                                      date_index, 
                                                      frequencies, 
                                                      date_group,
                                                      start,
                                                      end,
                                                      df_freq='M',
                                                      fill='none',
                                                      retain_na=True)

    rdd_resampled['Campaspe River Use'].ffill(inplace=True)
    rdd_resampled['Campaspe River Use'].fillna(0., inplace=True)
    
    rdd_resampled['WWC'] = resample_to_model_data_index(rdd['WWC'],
                                                      date_index, 
                                                      frequencies, 
                                                      date_group,
                                                      start,
                                                      end,
                                                      df_freq='M',
                                                      fill='none',
                                                      retain_na=True)

    rdd_resampled['WWC'].ffill(inplace=True)
    rdd_resampled['WWC'].fillna(0., inplace=True)

    
    #!!!! 
    
    # Create reaches based on locations of gauging stations:
    div_gauge_locations = np.sort(Campaspe_info[Campaspe_info.index.isin([406207, 406203, 406202])]['seg_loc'].unique())
    div_reach_no = range(len(div_gauge_locations))
    river_seg['div_reach'] = np.nan    
    river_seg.loc[river_seg['iseg'].isin(div_gauge_locations), 'div_reach'] = div_reach_no         
    river_seg['div_reach'].fillna(method='ffill', inplace=True)
    river_seg['div_reach'].fillna(method='bfill', inplace=True)
    river_seg['div_reach'].astype(int, inplace=True)
    # Remove diversions from 406207 and upstream if it didn't map to first segment
    river_seg.loc[river_seg['iseg'].isin(range(1, div_gauge_locations[0] + 1)), 'div_reach'] = -99
    river_segs_reach = [river_seg['iseg'][river_seg['div_reach'] == x].tolist() for x in div_reach_no]
    river_segs_reach_len = [float(len(river_segs_reach[ind])) for ind, x in enumerate(river_segs_reach)]
                        
    CID_gauge_loc = np.sort(Campaspe_info[Campaspe_info.index == 406218]['seg_loc'])[0]
    WWC_gauge_loc = np.sort(Campaspe_info[Campaspe_info.index == 406202]['seg_loc'])[0]
                        
    flow = [0] * len(nseg)
    runoff = [0] * len(nseg)
    etsw = [0] * len(nseg)
    pptsw = [0] * len(nseg)
    
    # Set the roughness for the channel
    roughch_val = [mbo.parameters.param['mn_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    roughch = np.interp(river_seg['Cumulative Length'].tolist(), 
                                    known_points, roughch_val)
    # Set the roughness for the banks
    roughbk_val = [mbo.parameters.param['mn_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    roughbk = np.interp(river_seg['Cumulative Length'].tolist(), 
                                    known_points, roughbk_val)
    river_seg['roughch'] = roughch
    river_seg['roughbk'] = roughbk
    
    cdpth = [0] * len(nseg)
    fdpth = [0] * len(nseg)
    awdth = [0] * len(nseg)
    bwdth = [0] * len(nseg)
    
    width1_val = [mbo.parameters.param['rivwdth{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
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

        # Private diversions:
        cols = ['Eppalock to C/Weir', 'Campaspe Weir to WWC', 'WWC to Murray']
        for ind, div_reach in enumerate(river_segs_reach):
            flow_temp = - rdd_resampled['Campaspe River Use'][cols[ind]].loc[date] / river_segs_reach_len[ind]
            for div_seg in div_reach:
                flow[div_seg - 1] = flow_temp

        if include_extractions:
            # CID diversions:
            flow[CID_gauge_loc] = flow[CID_gauge_loc] - rdd_resampled['CID Diversions'].loc[date]
            # WWC pumping:            
            flow[WWC_gauge_loc] = flow[WWC_gauge_loc] \
                                 + rdd_resampled['WWC']['WWC TO CAMPASPE (RO317)'].loc[date] \
                                 + rdd_resampled['WWC']['WWC TO CAMPASPE (RO317).1'].loc[date] \
                                 - rdd_resampled['WWC']['CAMPASPE PUMP TO WWC -REG'].loc[date] \
                                 - rdd_resampled['WWC']['CAMPASPE PUMP TO WWC -UNREG'].loc[date] \

        if include_et:
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
    
    
    seg_dict = segment_data1

    return segment_data, seg_dict
    
def prepare_river_data_for_murray(model_builder_object, 
                                    surface_raster_high_res,
                                    Murray_river_poly_file,
                                    Campaspe_relevant,
                                    river_stage_data,
                                    river_seg,
                                    pilot_points_YX=False,
                                    plot=False):
    
    mbo = model_builder_object
    # Parameter to modify the stage, thus accounting for errors in values specified for stage
    
    if not pilot_points_YX:
        mbo.parameters.create_model_parameter('rmstag', value=0.01)
        mbo.parameters.parameter_options('rmstag', 
                                              PARTRANS='log', 
                                              PARCHGLIM='factor', 
                                              PARLBND=0.001, 
                                              PARUBND=0.1, 
                                              PARGP='murriv', 
                                              SCALE=1, 
                                              OFFSET=0)
        # Parameter to shift the location of the bed which is only estimated based on assumed depth below zero gauge
        mbo.parameters.create_model_parameter('rmbed', value=0.01)
        mbo.parameters.parameter_options('rmbed', 
                                              PARTRANS='log', 
                                              PARCHGLIM='factor', 
                                              PARLBND=0.001, 
                                              PARUBND=0.1, 
                                              PARGP='murriv', 
                                              SCALE=1, 
                                              OFFSET=0)
        # Parameter for River Murray bed thickness
        mbo.parameters.create_model_parameter('rmbdtk', value=0.01)
        mbo.parameters.parameter_options('rmbdtk', 
                                              PARTRANS='log', 
                                              PARCHGLIM='factor', 
                                              PARLBND=0.001, 
                                              PARUBND=0.5, 
                                              PARGP='murriv', 
                                              SCALE=1, 
                                              OFFSET=0)
        # Parameter for the vertical hydraulic conductivity of the River Murray
        mbo.parameters.create_model_parameter('kv_rm', value=5E-3)
        mbo.parameters.parameter_options('kv_rm', 
                                              PARTRANS='log', 
                                              PARCHGLIM='factor', 
                                              PARLBND=1E-8, 
                                              PARUBND=20, 
                                              PARGP='murriv', 
                                              SCALE=1, 
                                              OFFSET=0)
        # Parameter for the width of the River Murray
        mbo.parameters.create_model_parameter('rmwdth', value=30)
        mbo.parameters.parameter_options('rmwdth', 
                                              PARTRANS='fixed', 
                                              PARCHGLIM='factor', 
                                              PARLBND=20, 
                                              PARUBND=50, 
                                              PARGP='murriv', 
                                              SCALE=1, 
                                              OFFSET=0)
    elif pilot_points_YX:
        mbo.parameters.create_model_parameter('rmstag', value=0.01)
        mbo.parameters.parameter_options('rmstag', 
                                              PARTRANS='fixed', 
                                              PARCHGLIM='factor', 
                                              PARLBND=0.001, 
                                              PARUBND=0.1, 
                                              PARGP='murriv', 
                                              SCALE=1, 
                                              OFFSET=0)
        # Parameter to allow shifting the location of the bed which is only estimated based on assumed depth below zero gauge
        mbo.parameters.create_model_parameter('rmbed', value=0.01)
        mbo.parameters.parameter_options('rmbed', 
                                              PARTRANS='fixed', 
                                              PARCHGLIM='factor', 
                                              PARLBND=0.001, 
                                              PARUBND=0.1, 
                                              PARGP='murriv', 
                                              SCALE=1, 
                                              OFFSET=0)
        # Parameter for River Murray bed thickness
        mbo.parameters.create_model_parameter('rmbdtk', value=0.01)
        mbo.parameters.parameter_options('rmbdtk', 
                                              PARTRANS='fixed', 
                                              PARCHGLIM='factor', 
                                              PARLBND=0.001, 
                                              PARUBND=0.5, 
                                              PARGP='murriv', 
                                              SCALE=1, 
                                              OFFSET=0)
        # Parameter for the vertical hydraulic conductivity of the River Murray
        mbo.parameters.create_model_parameter('kv_rm', value=5E-3)
        mbo.parameters.parameter_options('kv_rm', 
                                              PARTRANS='log', 
                                              PARCHGLIM='factor', 
                                              PARLBND=1E-8, 
                                              PARUBND=20, 
                                              PARGP='murriv', 
                                              SCALE=1, 
                                              OFFSET=0)
        # Parameter for the width of the River Murray
        mbo.parameters.create_model_parameter('rmwdth', value=30)
        mbo.parameters.parameter_options('rmwdth', 
                                              PARTRANS='fixed', 
                                              PARCHGLIM='factor', 
                                              PARLBND=20, 
                                              PARUBND=50, 
                                              PARGP='murriv', 
                                              SCALE=1, 
                                              OFFSET=0)
        
    print(" ** Creating river dataframe")
    mbo.create_river_dataframe('Murray', Murray_river_poly_file, surface_raster_high_res, verbose=True)
    
    mriver_seg = mbo.river_mapping['Murray']

    mriver_seg.loc[:, 'strthick'] = mbo.parameters.param['rmbdtk']['PARVAL1']
    
    # Set up bed elevations based on the gauge zero levels:
    Murray = Campaspe_relevant[Campaspe_relevant['Site Name'].str.contains("MURRAY RIVER")]
    gauge_points = [x for x in zip(Murray.Easting, Murray.Northing)]
    mriver_gauge_seg = mbo.get_closest_riv_segments('Murray', gauge_points)
    
    mriver_seg.loc[:, 'bed_from_gauge'] = np.nan
    #mriver_seg.loc[:, 'bed_from_gauge'] = mriver_seg['strtop']
    mriver_seg.loc[:, 'stage_from_gauge'] = np.nan
    mriver_seg.loc[:, 'strtop_offset'] = np.nan

    Murray.loc[:, 'new_gauge'] = Murray[['Gauge Zero (Ahd)', 'Cease to flow level', 'Min value']].max(axis=1) 
    Murray.loc[:, 'seg_loc'] = mriver_gauge_seg         
    mriver_seg.loc[:, 'iseg'] = [mriver_seg.shape[0] - x for x in range(mriver_seg.shape[0])]
    #Murray_gauge_zero['Cumulative Length'] = mriver_seg.loc[Murray_gauge_zero['seg_loc'].tolist(), 'Cumulative Length'].tolist()
    
    Murray = pd.merge(Murray, river_stage_data[1], on='Site Name', how='inner', suffixes=('','_r'))
    Murray = Murray[[x for x in Murray.columns if '_r' not in x]]
    Murray_gauge_zero = Murray[Murray['new_gauge'] > 10.]
  
    def values_from_gauge(column, to, interp=True):
        mriver_seg.loc[mriver_seg['iseg'].isin( \
            Murray_gauge_zero['seg_loc'].tolist()), column] = sorted( \
            Murray_gauge_zero[to].tolist(), reverse=True)
        
        #strtop_stage_from_gauge_diff = mriver_seg[mriver_seg['stage_from_gauge']]
        if interp:
            mriver_seg[column] = mriver_seg.set_index( \
                mriver_seg['Cumulative Length'])[column].interpolate( \
                method='values', limit_direction='both').tolist()
            mriver_seg[column].fillna(method='bfill', inplace=True)
    
    values_from_gauge('stage_from_gauge', 'Mean stage (m)', interp=False)
    
    # Setup the Murray river bed level from the zero gauge data
    values_from_gauge('bed_from_gauge', 'new_gauge', interp=False)        

    # Find the offset between the strtop from elevation data to help inform bed levels 
    # where gauge data isn't present
    inds = mriver_seg['bed_from_gauge'].notnull()
    inds2 = mriver_seg['bed_from_gauge'].isnull()
    mriver_seg.loc[inds, 'strtop_offset'] = mriver_seg[inds]['bed_from_gauge'] - \
                      mriver_seg[inds]['strtop']

    # Interpolate between gauges for the offset and then fill the ends of the
    # df with the closest value
    mriver_seg['strtop_offset'] = mriver_seg.set_index( \
        mriver_seg['Cumulative Length'])['strtop_offset'].interpolate( \
        method='values').where(mriver_seg.bfill().notnull()).tolist()
    mriver_seg.loc[:, 'strtop_offset'] = mriver_seg['strtop_offset'].bfill()
    mriver_seg.loc[:, 'strtop_offset'] = mriver_seg['strtop_offset'].ffill()
                    
    mriver_seg.loc[inds2, 'bed_from_gauge'] = mriver_seg.loc[inds2, 'strtop'] + \
                                              mriver_seg.loc[inds2, 'strtop_offset']

    # DO the same but for the river stage
    mriver_seg.loc[inds, 'stage_offset'] = \
        mriver_seg[inds]['stage_from_gauge'] - \
        mriver_seg[inds]['bed_from_gauge']

    # Interpolate between gauges for the offset and then fill the ends of the
    # df with the closest value
    mriver_seg['stage_offset'] = mriver_seg.set_index( \
        mriver_seg['Cumulative Length'])['stage_offset'].interpolate( \
        method='values').where(mriver_seg.bfill().notnull()).tolist()
    mriver_seg.loc[:, 'stage_offset'] = mriver_seg['stage_offset'].bfill()
    mriver_seg.loc[:, 'stage_offset'] = mriver_seg['stage_offset'].ffill()
                    
    mriver_seg.loc[inds2, 'stage_from_gauge'] = \
        mriver_seg.loc[inds2, 'bed_from_gauge'] + \
        mriver_seg.loc[inds2, 'stage_offset']

                                 
    new_k = []
    active = []
    surface_layers = {}
    bottom_layer = []
    
    for row in mriver_seg.iterrows():
        j_mesh = int(row[1]['i'])
        i_mesh = int(row[1]['j'])
        strtop = row[1]['bed_from_gauge']
        k = find_layer(strtop, mbo.model_mesh3D[0][:, j_mesh, i_mesh])
        new_k += [k]
        bottom_layer += [mbo.model_mesh3D[0][k + 1, j_mesh, i_mesh]] 
        active += [mbo.model_mesh3D[1][k, j_mesh, i_mesh]]
        for layer in range(7):
            try:
                surface_layers[layer] += [mbo.model_mesh3D[0][layer, j_mesh, i_mesh]]
            except:
                surface_layers[layer] = [mbo.model_mesh3D[0][layer, j_mesh, i_mesh]]
            # end try
        # end for
    # end for
        
    for layer in range(7):
        mriver_seg["surf{}".format(layer)] = surface_layers[layer]
    
    mriver_seg.loc[:, 'bottom_layer'] = bottom_layer
    
    mriver_seg.loc[:, 'k'] = new_k
    mriver_seg.loc[:, 'active'] = active
          
    # Remove any stream segments for which the elevation could not be mapped to a layer
    mriver_seg[mriver_seg['active'] == -1] = np.nan
    mriver_seg.dropna(inplace=True)
    mbo.river_mapping['Murray'] = mriver_seg
    
    mriver_seg.loc[:, 'strtop'] = mriver_seg['bed_from_gauge']                      
                          
    mriver_seg.loc[:, 'strhc1'] = mbo.parameters.param['kv_rm']['PARVAL1']                      
    
    mriver_seg.loc[:, 'width1'] = mbo.parameters.param['rmwdth']['PARVAL1']
    
    mriver_seg.loc[:, 'stage'] = mriver_seg['stage_from_gauge'] #+ SS_model.parameters.param['rmstage']['PARVAL1']
    
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
    
    
    mriver_seg.loc[:, 'amalg_riv_points_tuple'] = mriver_seg['amalg_riv_points'].apply(lambda x: (x[0], x[1]))    
    mriver_seg = mriver_seg.iloc[::-1]
    mriver_seg.index = range(mriver_seg.shape[0])
    
    print(" ** Merging collocated stream reaches")
    mriver_seg = river_df_tools.merge_collocated_stream_reaches(mriver_seg, max_length=2000.)
    print(" ** Merging short stream reaches")
    mriver_seg = river_df_tools.merge_very_short_stream_reaches(mriver_seg, min_length=300.)

    print(" ** Rechecking for collocated Campaspe and Murray cells")
    #mriver_seg['cell_loc_tuple'] = [(x[1]['k'], x[1]['i'], x[1]['j']) for x in mriver_seg.iterrows()]
    #mriver_seg['cell_loc_tuple'] = [(x[1]['i'], x[1]['j']) for x in mriver_seg.iterrows()]
    #mriver_seg = mriver_seg.groupby(by='cell_loc_tuple').mean()
    #mriver_seg.index = range(mriver_seg.shape[0])
    
    cells_overlapping = is_in_other_river(mriver_seg, river_seg)
    mriver_seg.loc[:, 'cell_used'] = cells_overlapping
    mriver_seg_ghb = mriver_seg.copy()

    mriver_seg[mriver_seg['cell_used'] == 0] = np.nan
    mriver_seg.dropna(inplace=True)
    
    
    mbo.river_mapping['Murray'] = mriver_seg
    
    if plot:
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

    return riv, mriver_seg_ghb

def create_riv_data_transient(model_builder_object,
                                  river_seg,
                                  river_stage_data,
                                  Campaspe_info,
                                  date_index, 
                                  frequencies, 
                                  date_group,
                                  start,
                                  end):
    
    mbo = model_builder_object
    num_reaches = mbo.pilot_points['Campaspe'].num_points
    known_points = mbo.pilot_points['Campaspe'].points
    
    ###########################################################################
    # Convert stages to depths for assistance in the interpolation allowing use of bed elevations
    #
    # Estimate gauge zeros where there is no data
    #
    #
    #
    #
    Campaspe_stage = river_stage_data[1]
    Campaspe_stage = Campaspe_stage[~Campaspe_stage['Site Name'].str.contains('MURRAY')]
    filtered = Campaspe_stage['Site ID'].tolist()
    gauge_points = [x for x in zip(Campaspe_stage.Easting, Campaspe_stage.Northing)]
    river_gauge_seg = mbo.get_closest_riv_segments('Campaspe', gauge_points)
    Campaspe_stage.loc[:, 'new_gauge'] = Campaspe_stage[['Gauge Zero (Ahd)', 'Low stage (m)']].max(axis=1) 
    Campaspe_stage.loc[:, 'seg_loc'] = river_gauge_seg         

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #@@ Create temporal segment data @@
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
    stage_resampled = {}
    depth_resampled = {}
    for gauge in filtered:
        stage = river_stage_data[0][gauge]
        #Tidy up stage data:
        stage = stage[stage['Mean'] > 0.]
        stage_resampled[gauge] = resample_to_model_data_index(stage, 
                                                              date_index, 
                                                              frequencies, 
                                                              date_group,
                                                              start,
                                                              end,
                                                              df_freq='D',
                                                              fill='mean'
                                                              )                                                             
        depth_resampled[gauge] = stage_resampled[gauge]['Mean'] - Campaspe_stage[Campaspe_stage['Site ID'] == gauge]['new_gauge'].tolist()[0]


    width1_val = [mbo.parameters.param['rivwdth{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    width1 = np.interp(river_seg['Cumulative Length'].tolist(), 
                                    known_points, width1_val)
    river_seg['width1'] = width1
    
    # Prepare the stage data by resampling it to fit the model
    # stress periods and backfilling missing data with the 25th percentile of flow
    
    riv = {}
    for per, date in enumerate(date_index[:-1]):
        Campaspe_stage.loc[:, 'depth'] = np.nan
        for gauge in depth_resampled:
            Campaspe_stage.loc[Campaspe_stage['Site ID'] == gauge, 'depth'] = depth_resampled[gauge].loc[date] 

        river_seg.loc[:, 'depth'] = np.nan
        river_seg.loc[river_seg['iseg'].isin(Campaspe_stage['seg_loc'].tolist()), 'depth'] = sorted(Campaspe_stage['depth'].tolist(), reverse=True)
        river_seg.loc[:, 'depth'] = river_seg.set_index(river_seg['Cumulative Length'])['depth'].interpolate(method='values', limit_direction='both').tolist()
        river_seg.loc[:, 'depth'] = river_seg['depth'].bfill() 
        river_seg.loc[:, 'stage'] = river_seg['strtop'] + river_seg['depth']
        
        simple_river = []
        for row in river_seg.iterrows():
            row = row[1]
            simple_river += [[row['k'], row['i'], row['j'], row['stage'], \
                              row['strhc1'] * row['rchlen'] * row['width1'], \
                              row['strtop']]]
    
        riv[per] = simple_river

    gauge_points = [x for x in zip(Campaspe_stage.Easting, Campaspe_stage.Northing)]
    river_gauge_seg = mbo.get_closest_riv_segments('Campaspe', gauge_points)
    river_seg.loc[:, 'gauge_id'] = 'none'
    river_seg.loc[river_seg['iseg'].isin(Campaspe_stage['seg_loc'].tolist()),
                  'gauge_id'] = \
    Campaspe_stage.sort_values('Mean stage (m)', ascending=False)['Site ID'].astype(str).tolist()
   

    return riv