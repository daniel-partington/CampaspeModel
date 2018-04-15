import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from CampaspeModel.build_utils.multifrequency_resampling import resample_to_model_data_index

def prepare_pumping_data_for_model(model_builder_object,
                                   pumps_points,
                                   start_pumping,
                                   start,
                                   end,
                                   date_index,
                                   pumping_data,
                                   frequencies,
                                   date_group,
                                   plot=False):

    '''
    NOTE: This function assumes monthly time_steps for the resampling of pumping
    data, which is resampled again based on date_index for the model
    '''
    
    mbo = model_builder_object
    
    mbo.map_points_to_grid(pumps_points, feature_id='OLD ID')
    
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
    
    pump_shallow = [] # Shallow (if <25m) or Deep (>= 25m)

    for pump_cell in mbo.points_mapped['pumping wells_clipped.shp']:
        row = pump_cell[0][0]
        col = pump_cell[0][1]
        for pump in pump_cell[1]: 
            if pumping_data.loc[pump, 'Top screen depth (m)'] == 0.: 
                #'No data to place pump at depth ... ignoring '            
                continue
            pump_depth = mbo.model_mesh3D[0][0][row][col] - pumping_data.loc[pump, 'Top screen depth (m)']        
            active = False
            for i in range(mbo.model_mesh3D[0].shape[0]-1):
                if pump_depth < mbo.model_mesh3D[0][i][row][col] and \
                     pump_depth > mbo.model_mesh3D[0][i+1][row][col]:
                    active_layer = i
                    active = True
                    break
                # end if
            # end for
            if mbo.model_mesh3D[1][active_layer][row][col] == -1:
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
            if pumping_data_ts.iloc[-1][pump] > 0.:
                pumping_data_ts = pumping_data_ts.fillna(0.0)
                pumping_data_ts = pumping_data_ts.reindex(date_index2)    
                pumping_data_ts = pumping_data_ts.ix[start_pumping:end]
                pumping_data_ts = pumping_data_ts.ffill()
                pumping_data_ts = pumping_data_ts.fillna(0.0)
            else:
                pumping_data_ts = pumping_data_ts.fillna(0.0)
                pumping_data_ts = pumping_data_ts.reindex(date_index2)    
                pumping_data_ts = pumping_data_ts.ix[start_pumping:end]
                pumping_data_ts = pumping_data_ts.fillna(0.0)

            resampled_pumping_data_ts = \
                resample_to_model_data_index(pumping_data_ts, date_index, 
                                             frequencies, date_group, start, end,
                                             index_report=False, fill='zero')
            # Now fill in the well dictionary with the values of pumping at relevant stress periods
            for index, time in enumerate(resampled_pumping_data_ts.iterrows()):
                if index >= mbo.model_time.t['steps']: 
                    print index
                    continue
                try:
                    wel[index] += [[active_layer, row, col, -time[1][pump]]]
                except:
                    wel[index] = [[active_layer, row, col, -time[1][pump]]]
        
    return wel