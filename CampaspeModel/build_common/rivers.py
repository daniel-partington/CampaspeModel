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
        #end if
    # end for
    

def prepare_river_data_for_Campaspe(ModelBuilderObject, 
                                    surface_raster_high_res,
                                    river_gauges,
                                    Campaspe_river_poly_file,
                                    Campaspe,
                                    num_reaches=4,
                                    plot=False):
    
    MBO = ModelBuilderObject
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
    


    MBO.GISInterface.raster_reproject_by_grid(surface_raster_high_res,
                                                   surface_raster_high_res[:-4] + '_reproj.tif',
                                                   resample_method='min')
    
    surface_raster_high_res = surface_raster_high_res[:-4] + '_reproj.tif'
    
    
    MBO.map_points_to_grid(river_gauges, feature_id='Site Name')
    #SS_model.map_points_to_grid(river_gauges, feature_id='Site_ID')
    
    Campaspe_river_gauges = MBO.points_mapped['processed_river_sites_stage_clipped.shp']
    
    filter_gauges = []
    for riv_gauge in Campaspe_river_gauges:
        #if riv_gauge[1][0] in use_gauges:
        if str(riv_gauge[1][0]) in use_gauges:
            filter_gauges += [riv_gauge]
    
    
    MBO.create_river_dataframe('Campaspe', Campaspe_river_poly_file, surface_raster_high_res)
    
    # Create reach data
    river_seg = MBO.river_mapping['Campaspe']
    # Parameters are ordered from upstream to downstream
    num_reaches = num_reaches
    base_guesses_for_k_bed_x = [0., 0.33, 0.66, 1.0]
    base_guesses_for_k_bed_y = [1., 1., 0.01, 0.001]
    interp_guesses_for_k_bed_x = np.linspace(0., 1., num_reaches)
    interp_guesses_for_k_bed_y = np.interp(interp_guesses_for_k_bed_x, 
                                           base_guesses_for_k_bed_x,
                                           base_guesses_for_k_bed_y) 
    
    MBO.create_pilot_points('Campaspe', linear=True)
    camp_pp = MBO.pilot_points['Campaspe']
    camp_pp.set_uniform_points(river_seg['rchlen'].sum(), num_reaches)
    
    known_points = camp_pp.points
    
    # Define split on river for which unique values will be given to props at 
    # those points which will then be interpolated along the length of the river
    #for reach in range(num_reaches):
    # Setting up river bed hydraulic conductivity values
    MBO.parameters.create_model_parameter_set('kv_riv', 
                                               value=list(interp_guesses_for_k_bed_y), 
                                               num_parameters=num_reaches)
    
    MBO.parameters.parameter_options_set('kv_riv', 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=0.01, 
                                          PARUBND=10.0, 
                                          PARGP='kv_riv', 
                                          SCALE=1, 
                                          OFFSET=0)

    # Setting up river bed elevation correction parameter to account for 
    # uncertainty in where bed lies relative to zero gauge
    MBO.parameters.create_model_parameter_set('beddep', 
                                               value=0.01, 
                                               num_parameters=num_reaches)
    MBO.parameters.parameter_options_set('beddep', 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=0.001, 
                                          PARUBND=1.0, 
                                          PARGP='rivbed', 
                                          SCALE=1, 
                                          OFFSET=0)
    # Setting up river bed roughness values
    MBO.parameters.create_model_parameter_set('mn_riv', 
                                               value=0.001, 
                                               num_parameters=num_reaches)
    MBO.parameters.parameter_options_set('mn_riv', 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=0.001, 
                                          PARUBND=0.1, 
                                          PARGP='rough', 
                                          SCALE=1, 
                                          OFFSET=0)
    # Setting up river width values
    MBO.parameters.create_model_parameter_set('rivwdth', 
                                               value=20.0, 
                                               num_parameters=num_reaches)
    
    MBO.parameters.parameter_options_set('rivwdth', 
                                          PARTRANS='fixed', 
                                          PARCHGLIM='factor', 
                                          PARLBND=4., 
                                          PARUBND=40., 
                                          PARGP='rivwdt', 
                                          SCALE=1, 
                                          OFFSET=0)
    # Setting up riverbed thickness values
    MBO.parameters.create_model_parameter_set('bedthck', 
                                               value=0.10, 
                                               num_parameters=num_reaches)
    
    MBO.parameters.parameter_options_set('bedthck', 
                                          PARTRANS='fixed', 
                                          PARCHGLIM='factor', 
                                          PARLBND=0.01, 
                                          PARUBND=1., 
                                          PARGP='bedthk', 
                                          SCALE=1, 
                                          OFFSET=0)
    
    
    strcond_val = [MBO.parameters.param['kv_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    river_seg['strhc1'] = np.interp(river_seg['Cumulative Length'].tolist(), 
                                    known_points, strcond_val)
    strthick_val = [MBO.parameters.param['bedthck{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    river_seg['strthick'] = np.interp(river_seg['Cumulative Length'].tolist(), known_points, strthick_val)
    
    river_seg['amalg_riv_points_tuple'] = river_seg['amalg_riv_points'].apply(lambda x: (x[0], x[1]))    
    
    river_seg = river_df_tools.merge_collocated_stream_reaches(river_seg, max_length=3000.)
    river_seg = river_df_tools.merge_very_short_stream_reaches(river_seg, min_length=500.)
    
    MBO.river_mapping['Campaspe'] = river_seg
    
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
        new_k += [find_layer(strbot, MBO.model_mesh3D[0][:, j_mesh, i_mesh])]
    
    river_seg['k'] = new_k
           
    # Remove any stream segments for which the elevation could not be mapped to a layer
    #river_seg.dropna(inplace=True)
    
    river_seg['ireach'] = 1
    river_seg['iseg'] = [x + 1 for x in range(river_seg.shape[0])]
    
                          
    # Set up bed elevations based on the gauge zero levels:
    gauge_points = [x for x in zip(Campaspe.Easting, Campaspe.Northing)]
    river_gauge_seg = MBO.get_closest_riv_segments('Campaspe', gauge_points)
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
        k = find_layer(strbot, MBO.model_mesh3D[0][:, j_mesh, i_mesh])
        new_k += [k]
        bottom_layer += [MBO.model_mesh3D[0][k + 1, j_mesh, i_mesh]] 
        for layer in range(7):
            try:
                surface_layers[layer] += [MBO.model_mesh3D[0][layer, j_mesh, i_mesh]]
            except:
                surface_layers[layer] = [MBO.model_mesh3D[0][layer, j_mesh, i_mesh]]
    
    for layer in range(7):
        river_seg["surf{}".format(layer)] = surface_layers[layer]
    
    river_seg['k'] = new_k
    river_seg['strtop'] = river_seg['bed_from_gauge']
    river_seg['bottom_layer'] = bottom_layer
    
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
        else:
            return x
        # end if
        
    river_seg['slope'] = river_seg['slope'].apply(lambda x: slope_corrector(x))
    
    reach_df = river_seg[['k','i','j','iseg','ireach','rchlen','strtop','slope','strthick','strhc1']]
    reach_data = reach_df.to_records(index=False)
    
    return river_seg, reach_df, reach_data, known_points
    
def create_segment_data(ModelBuilderObject, river_seg, river_flow_data):
    MBO = ModelBuilderObject
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
    
    num_reaches = MBO.pilot_points['Campaspe'].num_points
    known_points = MBO.pilot_points['Campaspe'].points
    # Set the roughness for the channel
    roughch_val = [MBO.parameters.param['mn_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    roughch = np.interp(river_seg['Cumulative Length'].tolist(), 
                                    known_points, roughch_val)
    # Set the roughness for the banks
    roughbk_val = [MBO.parameters.param['mn_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    roughbk = np.interp(river_seg['Cumulative Length'].tolist(), 
                                    known_points, roughbk_val)
    river_seg.loc[:, 'roughch'] = roughch
    river_seg.loc[:, 'roughbk'] = roughbk
    
    cdpth = [0] * len(nseg)
    fdpth = [0] * len(nseg)
    awdth = [0] * len(nseg)
    bwdth = [0] * len(nseg)
    
    width1_val = [MBO.parameters.param['rivwdth{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
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

def create_segment_data_transient(ModelBuilderObject,
                                  river_seg,
                                  river_flow_data,
                                  FieldData,
                                  interp_et,
                                  interp_rain,
                                  date_index, 
                                  frequencies, 
                                  date_group,
                                  start,
                                  end):
    
    MBO = ModelBuilderObject
    num_reaches = MBO.pilot_points['Campaspe'].num_points
    known_points = MBO.pilot_points['Campaspe'].points
    
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
    flow = [0] * len(nseg)
    runoff = [0] * len(nseg)
    etsw = [0] * len(nseg)
    pptsw = [0] * len(nseg)
    
    # Set the roughness for the channel
    roughch_val = [MBO.parameters.param['mn_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    roughch = np.interp(river_seg['Cumulative Length'].tolist(), 
                                    known_points, roughch_val)
    # Set the roughness for the banks
    roughbk_val = [MBO.parameters.param['mn_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    roughbk = np.interp(river_seg['Cumulative Length'].tolist(), 
                                    known_points, roughbk_val)
    river_seg['roughch'] = roughch
    river_seg['roughbk'] = roughbk
    
    cdpth = [0] * len(nseg)
    fdpth = [0] * len(nseg)
    awdth = [0] * len(nseg)
    bwdth = [0] * len(nseg)
    
    width1_val = [MBO.parameters.param['rivwdth{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
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
    
    
    seg_dict = segment_data1

    return segment_data, seg_dict
    
def prepare_river_data_for_Murray(ModelBuilderObject, 
                                    surface_raster_high_res,
                                    Murray_river_poly_file,
                                    Campaspe_relevant,
                                    river_stage_data,
                                    river_seg,
                                    plot=False):
    
    MBO = ModelBuilderObject
    # Parameter to modify the stage, thus accounting for errors in values specified for stage
    MBO.parameters.create_model_parameter('rmstage', value=0.01)
    MBO.parameters.parameter_options('rmstage', 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=0.001, 
                                          PARUBND=0.1, 
                                          PARGP='murr_riv', 
                                          SCALE=1, 
                                          OFFSET=0)
    # Parameter to all shifting the location of the bed which is only estimated based on assumed depth below zero gauge
    MBO.parameters.create_model_parameter('rmbed', value=0.01)
    MBO.parameters.parameter_options('rmbed', 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=0.001, 
                                          PARUBND=0.1, 
                                          PARGP='murr_riv', 
                                          SCALE=1, 
                                          OFFSET=0)
    # Parameter for River Murray bed thickness
    MBO.parameters.create_model_parameter('rmbedthk', value=0.01)
    MBO.parameters.parameter_options('rmbedthk', 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=0.001, 
                                          PARUBND=0.5, 
                                          PARGP='murr_riv', 
                                          SCALE=1, 
                                          OFFSET=0)
    # Parameter for the vertical hydraulic conductivity of the River Murray
    MBO.parameters.create_model_parameter('kv_rm', value=5E-3)
    MBO.parameters.parameter_options('kv_rm', 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=1E-8, 
                                          PARUBND=20, 
                                          PARGP='murr_riv', 
                                          SCALE=1, 
                                          OFFSET=0)
    # Parameter for the width of the River Murray
    MBO.parameters.create_model_parameter('rmwdth', value=30)
    MBO.parameters.parameter_options('rmwdth', 
                                          PARTRANS='fixed', 
                                          PARCHGLIM='factor', 
                                          PARLBND=20, 
                                          PARUBND=50, 
                                          PARGP='murr_riv', 
                                          SCALE=1, 
                                          OFFSET=0)
    
    MBO.create_river_dataframe('Murray', Murray_river_poly_file, surface_raster_high_res)
    
    mriver_seg = MBO.river_mapping['Murray']
    mriver_seg.loc[:, 'strthick'] = MBO.parameters.param['rmbedthk']['PARVAL1']
    
    # Set up bed elevations based on the gauge zero levels:
    Murray = Campaspe_relevant[Campaspe_relevant['Site Name'].str.contains("MURRAY RIVER")]
    gauge_points = [x for x in zip(Murray.Easting, Murray.Northing)]
    mriver_gauge_seg = MBO.get_closest_riv_segments('Murray', gauge_points)
    #mriver_seg.loc[:, 'bed_from_gauge'] = np.nan
    mriver_seg.loc[:, 'bed_from_gauge'] = mriver_seg['strtop']
    mriver_seg.loc[:, 'stage_from_gauge'] = np.nan
    
    Murray.loc[:, 'new_gauge'] = Murray[['Gauge Zero (Ahd)', 'Cease to flow level', 'Min value']].max(axis=1) 
    Murray.loc[:, 'seg_loc'] = mriver_gauge_seg         
    Murray_gauge_zero = Murray[Murray['new_gauge'] > 10.]
    mriver_seg.loc[:, 'iseg'] = [x + 1 for x in range(mriver_seg.shape[0])]
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
        k = find_layer(strtop, MBO.model_mesh3D[0][:, j_mesh, i_mesh])
        new_k += [k]
        bottom_layer += [MBO.model_mesh3D[0][k+1, j_mesh, i_mesh]] 
        active += [MBO.model_mesh3D[1][k, j_mesh, i_mesh]]
        for layer in range(7):
            try:
                surface_layers[layer] += [MBO.model_mesh3D[0][layer, j_mesh, i_mesh]]
            except:
                surface_layers[layer] = [MBO.model_mesh3D[0][layer, j_mesh, i_mesh]]
    
    for layer in range(7):
        mriver_seg["surf{}".format(layer)] = surface_layers[layer]
    
    mriver_seg.loc[:, 'bottom_layer'] = bottom_layer
    
    mriver_seg.loc[:, 'k'] = new_k
    mriver_seg.loc[:, 'active'] = active
          
    # Remove any stream segments for which the elevation could not be mapped to a layer
    mriver_seg[mriver_seg['active'] == -1] = np.nan
    mriver_seg.dropna(inplace=True)
    MBO.river_mapping['Murray'] = mriver_seg
    
    mriver_seg.loc[:, 'strtop'] = mriver_seg['bed_from_gauge']                      
                          
    mriver_seg.loc[:, 'strhc1'] = MBO.parameters.param['kv_rm']['PARVAL1']                      
    
    mriver_seg.loc[:, 'width1'] = MBO.parameters.param['rmwdth']['PARVAL1']
    
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
    
    cells_overlapping = is_in_other_river(mriver_seg, river_seg)
    mriver_seg.loc[:, 'cell_used'] = cells_overlapping
    mriver_seg_ghb = mriver_seg.copy()
    
    #mriver_seg[mriver_seg['cell_used'] == 0] = np.nan
    #mriver_seg.dropna(inplace=True)
    
    #mriver_seg['cell_loc_tuple'] = [(x[1]['k'], x[1]['i'], x[1]['j']) for x in mriver_seg.iterrows()]
    #mriver_seg['cell_loc_tuple'] = [(x[1]['i'], x[1]['j']) for x in mriver_seg.iterrows()]
    #mriver_seg = mriver_seg.groupby(by='cell_loc_tuple').mean()
    #mriver_seg.index = range(mriver_seg.shape[0])
    
    mriver_seg.loc[:, 'amalg_riv_points_tuple'] = mriver_seg['amalg_riv_points'].apply(lambda x: (x[0], x[1]))    
    mriver_seg3 = mriver_seg.iloc[::-1]
    mriver_seg3.index = range(mriver_seg3.shape[0])
    
    mriver_seg3 = river_df_tools.merge_collocated_stream_reaches(mriver_seg3, max_length=3000.)
    mriver_seg = river_df_tools.merge_very_short_stream_reaches(mriver_seg3, min_length=500.)
    
    MBO.river_mapping['Murray'] = mriver_seg
    
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