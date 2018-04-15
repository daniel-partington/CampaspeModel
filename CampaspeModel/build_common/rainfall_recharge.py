import numpy as np

from CampaspeModel.build_utils.multifrequency_resampling import resample_to_model_data_index 

def weather2model_grid(MBO, df, rain_gauges):
    weather_dict = {}
    for step, month in enumerate(df.iterrows()):
        weather_dict[step] = MBO.interpolate_points2mesh(rain_gauges, month[1], feature_id='Name')
        # Adjust rainfall to m from mm 
        weather_dict[step] = weather_dict[step] / 1000.0
    return weather_dict

def prepare_transient_rainfall_data_for_model(ModelBuilderObject, 
                                    recharge_zones,
                                    recharge_info,
                                    long_term_historic_weather,
                                    date_index,
                                    frequencies,
                                    date_group,
                                    start,
                                    end,
                                    rain_gauges,
                                    pilot_points_YX=False):
    MBO = ModelBuilderObject
    resampled_weather = resample_to_model_data_index(long_term_historic_weather, \
                                                     date_index, frequencies, \
                                                     date_group, start, end)
    
    resampled_rain = resampled_weather[[x for x in resampled_weather.columns if '_ET' not in x]]
    resampled_et = resampled_weather[[x for x in resampled_weather.columns if '_ET' in x]]
    resampled_et.columns = [x.replace('_ET','') for x in resampled_et.columns]
    
    interp_rain = weather2model_grid(MBO, resampled_rain, rain_gauges)
    interp_et = weather2model_grid(MBO, resampled_et, rain_gauges)
    
    MBO.boundaries.create_model_boundary_condition('Rainfall', 'rainfall', bc_static=True)
    MBO.boundaries.assign_boundary_array('Rainfall', interp_rain)
    
    MBO.boundaries.create_model_boundary_condition('ET', 'pet', bc_static=True)
    MBO.boundaries.assign_boundary_array('ET', interp_et)
    
    # Need to make copies of all rainfall arrays
    interp_rain2 = {}
    for key in interp_rain.keys():
        interp_rain2[key] = np.copy(interp_rain[key])
    interp_rain = interp_rain2
    
    recharge_zone_array = MBO.map_raster_to_regular_grid_return_array(recharge_zones)
    rch_zone_dict = {i:x for i, x in enumerate(np.unique(recharge_zone_array))}
    rch_zones = len(rch_zone_dict.keys())
    
    # Adjust rainfall to recharge using rainfall reduction
    MBO.parameters.create_model_parameter_set('ssrch', 
                                              value=0.01,
                                              num_parameters=rch_zones)
    MBO.parameters.parameter_options_set('ssrch', 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=1.0E-3, 
                                          PARUBND=0.5, 
                                          PARGP='ssrch', 
                                          SCALE=1, 
                                          OFFSET=0)
    
    
    MBO.parameters.create_model_parameter_set('rchred',
                                              value=0.05,
                                              num_parameters=rch_zones)
    if pilot_points_YX:
        PARTRANS = 'fixed'
    else:
        PARTRANS = 'log'
    # end if
    
    MBO.parameters.parameter_options_set('rchred', 
                                          PARTRANS=PARTRANS, 
                                          PARCHGLIM='factor', 
                                          PARLBND=1E-3, 
                                          PARUBND=0.9, 
                                          PARGP='rchred', 
                                          SCALE=1, 
                                          OFFSET=0)

    for i in range(rch_zones - 1):
        df_rch_row = recharge_info[recharge_info['Zone code'] == rch_zone_dict[i + 1]]
        if df_rch_row.empty:
            continue
        # End if
        MBO.parameters.param['rchred{}'.format(i)]['PARVAL1'] = (df_rch_row['Mean'] / df_rch_row['P']).tolist()[0]
        MBO.parameters.param['rchred{}'.format(i)]['PARLBND'] = max(0.001, (df_rch_row['Min'] / df_rch_row['P']).tolist()[0])
        MBO.parameters.param['rchred{}'.format(i)]['PARUBND'] = (df_rch_row['Max'] / df_rch_row['P']).tolist()[0]
        if MBO.parameters.param['rchred{}'.format(i)]['PARLBND'] < 0.:
             MBO.parameters.param['rchred{}'.format(i)]['OFFSET'] = MBO.parameters.param['rchred{}'.format(i)]['PARLBND'] - 0.1

    
    for key in interp_rain.keys():
        for i in range(rch_zones - 1):
            interp_rain[key][recharge_zone_array == rch_zone_dict[i + 1]] = \
                interp_rain[key][recharge_zone_array == rch_zone_dict[i + 1]] * \
                MBO.parameters.param['rchred{}'.format(i)]['PARVAL1']
    
        interp_rain[key][recharge_zone_array==rch_zone_dict[0]] = \
            interp_rain[key][recharge_zone_array == rch_zone_dict[0]] * 0.0
        interp_rain[key][MBO.model_mesh3D[1][0] == -1] = 0.
    
    return interp_rain, interp_et, recharge_zone_array, rch_zone_dict
    
def prepare_static_rainfall_data_for_model(model_builder_object,
                                           recharge_zones,
                                           recharge_info,
                                           long_term_historic_weather,
                                           rain_gauges):
    '''
    UNTESTED and not documented
    '''
    
    long_term_historic_rainfall = long_term_historic_weather
    mbo = model_builder_object
    interp_rain = mbo.interpolate_points2mesh(rain_gauges, long_term_historic_rainfall, feature_id='Name', method='linear')
    # Adjust rainfall to m from mm and from year to day
    interp_rain = interp_rain / 1000.0 / 365.0
    
    mbo.boundaries.create_model_boundary_condition('Rainfall', 'rainfall', bc_static=True)
    mbo.boundaries.assign_boundary_array('Rainfall', interp_rain)
    
    # Replace interp_rain with a copy to prevent alteration of assigned boundary array
    interp_rain = np.copy(interp_rain)
    
    recharge_zone_array = mbo.map_raster_to_regular_grid_return_array(recharge_zones)
    
    rch_zone_dict = {i:x for i, x in enumerate(np.unique(recharge_zone_array))}
    rch_zones = len(rch_zone_dict.keys())
    
    mbo.parameters.create_model_parameter_set('ssrch', 
                                                   value=0.01,
                                                   num_parameters=rch_zones - 1)
    
    mbo.parameters.parameter_options_set('ssrch', 
                                              PARTRANS='log', 
                                              PARCHGLIM='factor', 
                                              PARLBND=1.0E-3, 
                                              PARUBND=0.2, 
                                              PARGP='ssrch', 
                                              SCALE=1, 
                                              OFFSET=0)
    
    for i in range(rch_zones - 1):
        interp_rain[recharge_zone_array == rch_zone_dict[i + 1]] = \
            interp_rain[recharge_zone_array == rch_zone_dict[i + 1]] * \
            mbo.parameters.param['ssrch{}'.format(i)]['PARVAL1']
    
    # Ensure model recharge is 0 over areas where the domain is inactive or where zonal array is a NaN value such as over lakes                                  
    interp_rain[recharge_zone_array==rch_zone_dict[0]] = interp_rain[recharge_zone_array == rch_zone_dict[0]] * 0.0
    interp_rain[mbo.model_mesh3D[1][0] == -1] = 0.
        
    rch = {}
    rch[0] = interp_rain    