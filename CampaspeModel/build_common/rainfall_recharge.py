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
                                    resampled_weather, 
                                    recharge_zones,
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
    
    for key in interp_rain.keys():
        for i in range(rch_zones - 1):
            interp_rain[key][recharge_zone_array == rch_zone_dict[i + 1]] = \
                interp_rain[key][recharge_zone_array == rch_zone_dict[i + 1]] * \
                MBO.parameters.param['rchred{}'.format(i)]['PARVAL1']
    
        interp_rain[key][recharge_zone_array==rch_zone_dict[0]] = \
            interp_rain[key][recharge_zone_array == rch_zone_dict[0]] * 0.0
        interp_rain[key][MBO.model_mesh3D[1][0] == -1] = 0.
    
    return interp_rain, interp_et, recharge_zone_array, rch_zone_dict