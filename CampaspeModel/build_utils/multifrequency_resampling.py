import pandas as pd
import numpy as np

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@ RESAMPLING TEMPORALLY WITH MULTIPLE FREQUENCIES @@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

def _fill_in_time_series_nan(df, fill='mean', stat='50%'):
    if fill == 'mean':
        df = df.fillna(df.mean())
    elif fill == 'stats':
        df = df.fillna(df.describe().loc[stat])
    elif fill == 'zero':
        df = df.fillna(0.)
    elif fill == 'none':
        pass
    #end if
    return df

def resample_to_model_data_index(df, date_index, frequencies, date_group, \
                                 start, end, \
                                 fill='mean', stat='50%', df_freq=None, 
                                 index_report=True, label='left', debug=False,
                                 retain_na=False):

    pd_dt = pd.to_datetime

    if len(frequencies) != len(date_group) - 1:
        print("Frequencies list must have one less item than the date_group list")
        return

    if df_freq != None:
        df = df.resample(df_freq).mean()
    # end if
    df = df.loc[start:end]

    #However if the time period for the model is longer we need to reindex the dataseries
    if df_freq == None:
        df_freq = pd.infer_freq(df.index)

    # Create temporary date_index 
    date_index_temp = pd.date_range(start=date_index[0], end=date_index[-1], \
                                    freq=df_freq)
    df = df.reindex(date_index_temp)
    # Then we have to fill in the missing values with mean or some other descriptor
    df = _fill_in_time_series_nan(df, fill=fill, stat=stat)
    # Create empty list for placing the resampled parts of the dataframe
    df_resamples = []

    len_frequencies = len(frequencies)
    for index, frequency in enumerate(frequencies):
        #print(frequency)
        p_start, p_end = date_group[index], date_group[index + 1]
        #resample = df[df.index.isin(pd.date_range(p_start, p_end))] \
        if index < len_frequencies - 1:
            resample = df[(df.index >= pd_dt(p_start)) & (df.index < pd_dt(p_end))] \
                          .resample(frequency, label=label).mean()
        elif len_frequencies == 1:
            resample = df[(df.index >= pd_dt(p_start))] \
                          .resample(frequency, label='right').mean()
        else:
            resample = df[(df.index >= pd_dt(p_start))] \
                          .resample(frequency, label=label).mean()
            
        if debug: print resample.index
        if index < len_frequencies - 1:
            if label == 'left':
                df_resamples += [resample.iloc[1:]]
            elif label == 'right':
                df_resamples += [resample.iloc[:-1]]
            # end if
        elif len_frequencies == 1:
            df_resamples += [resample]
        else:
            df_resamples += [resample.iloc[1:]]
        # end if
    # end for
    df_concat = pd.concat(df_resamples)
    if index_report:
        # TODO: Report if any of the df_concat indices are not in date_index
        if len_frequencies > 1:
            if np.all(np.in1d(df_concat.index, date_index[:-1])): #np.array_equal(df_concat.index, date_index):
                print("Successful match of date indices for model and resampled df")
            else:
                print("*** Failed match of some date indices for model and resampled df \n {0} \n {1}".format(df_concat.index, date_index))
                import sys        
                sys.exit("*** Failed match of some date indices for model and resampled df \n {0} \n {1}".format(df_concat.index, date_index))        
            # end if
        else:
            if np.all(np.in1d(df_concat.index, date_index)): #np.array_equal(df_concat.index, date_index):
                print("Successful match of date indices for model and resampled df")
            else:
                print("*** Failed match of some date indices for model and resampled df \n {0} \n {1}".format(df_concat.index, date_index))
                import sys        
                sys.exit("*** Failed match of some date indices for model and resampled df \n {0} \n {1}".format(df_concat.index, date_index))        
            # end if
        # end if
    # end if
    
    # Remove the dead rows from the dataframe if there was no filling
    if fill == 'none' and not retain_na:
        df_concat = df_concat.dropna()
    # end if
    
    return df_concat

def resample_obs_time_series_to_model_data_index(df_obs, date_index, \
                                                 frequencies, date_group, \
                                                 start, end, \
                                                 fill='mean', stat='50%',
                                                 df_freq='M'):
    '''
    Function to resample obs time series to the model time periods
    '''    
    obs_names = df_obs['name'].unique()

    count = 0
    obs_by_name_temp = []
    for obs_name in obs_names:
        df_obs_name = df_obs[df_obs['name'] == obs_name]
        df_obs_name.index = df_obs_name['datetime'] 
        df_obs_name = resample_to_model_data_index(df_obs_name, date_index, 
                                              frequencies, date_group, 
                                              start, end,
                                              fill='none', df_freq=df_freq,
                                              index_report=False)
        # Restore the datetime column
        if df_obs_name.empty: 
            continue
        
        df_obs_name.loc[:, 'datetime'] = df_obs_name.index
        # Restore the name column
        df_obs_name.loc[:,'name'] = obs_name                    
        df_obs_name.index = range(count, count + df_obs_name.shape[0])
        count += df_obs_name.shape[0] + 1
        obs_by_name_temp += [df_obs_name]

    df_obs = pd.concat(obs_by_name_temp)

    # TEST: Report if any of the df_obs datetime are not in date_index
    test_index = df_obs['datetime']
    if len(frequencies) > 1:
        if np.all(np.in1d(test_index, date_index[:-1])): #np.array_equal(df_concat.index, date_index):
            print("Successful match of date indices for model and resampled df")
        else:
            print("*** Failed match of some date indices for model and resampled df \n {0} \n {1}".format(test_index, date_index))
            import sys        
            sys.exit("*** Failed match of some date indices for model and resampled df \n {0} \n {1}".format(test_index, date_index))        
        # end if
    else:
        if np.all(np.in1d(test_index, date_index)):
            print("Successful match of date indices for model and resampled df")
        else:
            print("*** Failed match of some date indices for model and resampled df \n {0} \n {1}".format(test_index, date_index))
            import sys        
            sys.exit("*** Failed match of some date indices for model and resampled df \n {0} \n {1}".format(test_index, date_index))        
        # end if
        
    return df_obs