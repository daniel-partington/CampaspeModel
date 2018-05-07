import os
import numpy as np


def write_model_predictions(modflow_model, m):
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    #^^^ MODEL PREDICTIONS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # a = annual
    # s = seasonal
    # m = monthly
    #
    temporal = ['a', 's', 'm']
    freqs = ['A', 'Q', 'M']
    swgw_obs_groups = ['nrf_{}'.format(x) for x in temporal]
    reach_swgw_obs_groups_base = ['rrf_{}'.format(x) for x in temporal]
    reaches = range(11) # This needs to come automatically out of build!!
    reach_freqs = freqs * len(reaches)
    reach_swgw_obs_groups = []
    for reach in reaches:
        reach_swgw_obs_groups += [x + str(reach) for x in reach_swgw_obs_groups_base]
    
    all_swgw_obs_groups = swgw_obs_groups + reach_swgw_obs_groups
    all_freqs = freqs + reach_freqs
    sfr_df = modflow_model.importSfrOut()
      
    for index, obs_group in enumerate(all_swgw_obs_groups):
        swgw_obs = m.observations.obs_group[obs_group]
        swgw_obs_ts = swgw_obs['time_series']
        obs = m.observations.obs_group
        # Some locations from observations object are defined as a dataframe 
        # and others as list, so need two ways to handle ...
        try:
            sfr_location = obs[obs_group]['locations']['seg_loc']
        except:
            sfr_location = obs[obs_group]['locations']
        # end try            
        with open(os.path.join(modflow_model.data_folder, 
                               'observations_{}.txt'.format(obs_group))
                               , 'w') as f:
    
            sfr = sfr_df.copy()
            col_of_interest = 'Qaquifer'
            #dt = observation['datetime']
            dateindex = m.model_time.t['dateindex'][1:]
            sfr = sfr[sfr['segment'].isin(sfr_location)][[col_of_interest, 'time']]
            sfr = sfr.groupby('time').sum()
            sfr.index = dateindex
            dt = swgw_obs_ts['datetime'].tolist()[-1]
            month = dt.strftime('%b').upper()
            freq = all_freqs[index]
            if freq in ['A', 'Q']:
                freq = "{}-{}".format(freq, month)
            # end if
        
            sfr = sfr.loc[sfr.index > "2014-1-1"]
            sfr = sfr.resample(freq).mean()
            for observation in swgw_obs_ts.iterrows():
                sim_obs = sfr.loc[observation[1]['datetime']]
                f.write('%f\n' % sim_obs)                


def write_potential_observations(modflow_model):

    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    #^^^ POTENTIAL OBS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # Set model output arrays to None to initialise
    head = None
    sfr_df = None
    stream_options = ['st_sim', 'fl_sim']
    head_options = ['shshal', 'shdeep']
    # Write observation to file
    for obs_set in modflow_model.model_data.observations.obs_group.keys():
        obs_type = modflow_model.model_data.observations.obs_group[obs_set]['obs_type']
        # Import the required model outputs for processing
        if obs_type in head_options:
            # Check if model outputs have already been imported and if not import
            if head == None:
                headobj = modflow_model.importHeads()
                head = headobj.get_alldata()
        elif obs_type in stream_options:
            try:
                sfr_df = modflow_model.sfr_df
            except:
                sfr_df = modflow_model.importSfrOut()
            # End except
        else:
            continue
        # End if

        obs_df = modflow_model.model_data.observations.obs_group[obs_set]['time_series']
        obs_df = obs_df[obs_df['active'] == True]
        sim_map_dict = modflow_model.model_data.observations.obs_group[obs_set]['mapped_observations']

        if obs_type in stream_options:
            sfr_location = modflow_model.model_data.observations.obs_group[obs_set]['locations']#['seg_loc']
            for zone in obs_df['zone'].unique():
                if len(obs_df['zone'].unique()) == 1:
                    zone_txt = obs_set
                else:
                    zone_txt = obs_set + zone
                # End if
                with open(os.path.join(modflow_model.data_folder, 'observations_' + zone_txt + '.txt'), 'w') as f:
                    obs_df_zone = obs_df[obs_df['zone'] == zone]
                    for observation in obs_df_zone.index:
                        interval = int(obs_df_zone['interval'].loc[observation])
                        name = obs_df_zone['name'].loc[observation]
                        seg = int(sfr_location.loc[name])
                        sfr = sfr_df
                        col_of_interest = obs_type
                        if obs_type == 'fl_sim':
                            col_of_interest = 'Qout'
                        elif obs_type == 'st_sim':
                            col_of_interest = 'stage'
                        # End if
                        sim_obs = sfr[(sfr['segment'] == seg) &
                                      (sfr['time'] == interval)][col_of_interest]
                        f.write('%f\n' % sim_obs)
                    # End for
                # End with
            # End for
        # End if

        if obs_type in head_options:
            zone_txt = obs_set
            with open(os.path.join(modflow_model.data_folder, 'observations_' + zone_txt + '.txt'), 'w') as f:
                for observation in obs_df.index:
                    interval = int(obs_df['interval'].loc[observation])
                    name = obs_df['name'].loc[observation]
                    (x_cell, y_cell) = modflow_model.model_data.mesh2centroid2Dindex[
                        (sim_map_dict[name][1], sim_map_dict[name][2])]
                    (lay, row, col) = [sim_map_dict[name][0],
                                       sim_map_dict[name][1], sim_map_dict[name][2]]

                    sim_heads = [head[interval][lay][row][col]]

                    sim_head = np.mean(sim_heads)
                    f.write('%f\n' % sim_head)
                # End for
            # End with
        # End if
    # End for
# End write_potential_observations()