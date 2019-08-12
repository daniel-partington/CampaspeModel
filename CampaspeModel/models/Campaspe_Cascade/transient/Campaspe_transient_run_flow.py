import sys
import os
import time
import datetime
import numpy as np
import flopy.utils.binaryfile as bf

from HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface

# Configuration Loader
from HydroModelBuilder.Utilities.Config.ConfigLoader import ConfigLoader


def run(model_folder, data_folder, mf_exe, param_file="", verbose=False):
    
    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder, "02_transient_flow_packaged.pkl"))
    name = list(MM.GW_build.keys())[0]
    m = MM.GW_build[name]
    # Load in the new parameters based on parameters.txt or dictionary of new parameters
 
    if param_file != "":
        m.updateModelParameters(os.path.join(data_folder, 'parameters.txt'), verbose=verbose)
    
    if verbose:
        print("************************************************************************")
        print(" Updating HGU parameters ")
    
        # This needs to be automatically generated with the map_raster2mesh routine ...
    zone_map = {1: 'qa', 2: 'utb', 3: 'utqa', 4: 'utam', 5: 'utaf', 6: 'lta', 7: 'bse'}

    default_array = m.model_mesh3D[1].astype(float)
    Zone = np.copy(default_array)
    Kh = np.copy(default_array)
    Kv = np.copy(default_array)
    Sy = np.copy(default_array)
    SS = np.copy(default_array)
    
    def create_pp_points_dict(zone_map, Zone, prop_array, prop_name, m):
        points_values_dict = {}
        for index, key in enumerate(zone_map.keys()):
            for index2, param in enumerate(m.parameters.param_set[prop_name + zone_map[key]]):
                if index2 == 0:
                    points_values_dict[index] = [m.parameters.param[param]['PARVAL1']]
                else: 
                    points_values_dict[index] += [m.parameters.param[param]['PARVAL1']]
        return points_values_dict    
        
    def update_pilot_points(zone_map, Zone, prop_array, par_name, prop_name, prop_folder, m, prop_array_fname):
        points_values_dict = create_pp_points_dict(zone_map, Zone, prop_array, prop_name, m)
        p = m.pilot_points[par_name]
        zones = len(list(zone_map.keys()))
        p.output_directory = os.path.join(data_folder, prop_folder)
        p.update_pilot_points_files_by_zones(zones, points_values_dict)
        time.sleep(3)
        p.run_pyfac2real_by_zones(zones) 
        #p.save_mesh3D_array(filename=os.path.join(data_folder, prop_array_fname))
        return p.val_array

    Kh = update_pilot_points(zone_map, Zone, Kh, 'hk', 'kh', 'hk_pilot_points',
                             m, 'hk_val_array')  
    m.save_array(os.path.join(data_folder, 'Kh'), Kh)

    print(("Erroneous K pilot cells: {}".format(len(Kh[Kh > 10000.]))))
    Kh[Kh > 10000.] = 25.
    Kv = Kh * 0.1
    m.save_array(os.path.join(data_folder, 'Kv'), Kv)

    Sy = update_pilot_points(zone_map, Zone, Sy, 'sy', 'sy', 'sy_pilot_points',
                             m, 'sy_val_array')
    print(("Erroneous Sy pilot cells: {}".format(len(Sy[Sy > 0.5]))))
    Sy[Sy > 0.5] = 0.25
    m.save_array(os.path.join(data_folder, 'Sy'), Sy)
    
    SS = update_pilot_points(zone_map, Zone, SS, 'ss', 'ss', 'ss_pilot_points',
                             m, 'ss_val_array')
    print(("Erroneous Ss pilot cells: {}".format(len(SS[SS > 0.01]))))
    SS[SS > 0.01] = 1E-5
    m.save_array(os.path.join(data_folder, 'SS'), SS)
    
    m.properties.update_model_properties('Kh', Kh)
    m.properties.update_model_properties('Kv', Kv)
    m.properties.update_model_properties('Sy', Sy)
    m.properties.update_model_properties('SS', SS)
    
    if verbose:
        print("************************************************************************")
        print(" Updating river parameters ")
    
    reach_df = m.mf_sfr_df['Campaspe']['reach_df'] 
    segment_data = m.mf_sfr_df['Campaspe']['seg_df']

    num_reaches = m.pilot_points['Campaspe'].num_points #4
    known_points = m.pilot_points['Campaspe'].points
    # Create reach data
    river_seg = m.river_mapping['Campaspe']
    
    strcond_val = [m.parameters.param['kv_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    river_seg.loc[river_seg['strhc1'] != 0.0, 'strhc1'] = np.interp(
            river_seg[river_seg['strhc1'] != 0.0]['Cumulative Length'].tolist(), 
            known_points, strcond_val)
    
    strthick_val = [m.parameters.param['bedthck{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    river_seg['strthick'] = np.interp(river_seg['Cumulative Length'].tolist(), known_points, strthick_val)
    
    reach_df = river_seg[['k','i','j','iseg','ireach','rchlen','strtop','slope','strthick','strhc1']]
    reach_data = reach_df.to_records(index=False)
    
    # Set the roughness for the channel
    roughch_val = [m.parameters.param['mn_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    roughch = np.interp(river_seg['Cumulative Length'].tolist(), 
                                    known_points, roughch_val)
    # Set the roughness for the banks
    roughbk_val = [m.parameters.param['mn_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    roughbk = np.interp(river_seg['Cumulative Length'].tolist(), 
                                    known_points, roughbk_val)
    river_seg['roughch'] = roughch
    river_seg['roughbk'] = roughbk
    
    width1_val = [m.parameters.param['rivwdth{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    width1 = np.interp(river_seg['Cumulative Length'].tolist(), 
                                    known_points, width1_val)

    segment_data1 = {}
    for key in list(segment_data.keys()):
        segment_data[key]['width2'] = segment_data[key]['width1'] = width1
        segment_data1[key] = segment_data[key].to_records(index=False)
    
    seg_dict = segment_data1    

    if verbose:
        print("************************************************************************")
        print(" Updating Campaspe river boundary")
    
    m.boundaries.update_boundary_array('Campaspe River', [reach_data, seg_dict])
    
    if verbose:
        print("************************************************************************")
        print(" Updating Murray River boundary")
    
    mriver_seg = m.river_mapping['Murray']
    mriver_seg['strhc1'] = m.parameters.param['kv_rm']['PARVAL1']                      
    
    simple_river = []
    for row in mriver_seg.iterrows():
        row = row[1]
        simple_river += [[row['k'], row['i'], row['j'], row['stage'], \
                          row['strhc1'] * row['rchlen'] * row['width1'], row['strtop']]]
    
    riv = {}
    riv[0] = simple_river
    m.boundaries.update_boundary_array('Murray River', riv)
    
    if verbose:
        print("************************************************************************")
        print(" Updating recharge boundary ")

    # Adjust rainfall to recharge using zoned rainfall reduction parameters
    # Need to make copies of all rainfall arrays
    interp_rain = m.boundaries.bc['Rainfall']['bc_array']
    for key in list(interp_rain.keys()):
        interp_rain[key] = np.copy(interp_rain[key])

    recharge_zone_array = m.boundaries.bc['Rain_reduced']['zonal_array']
    rch_zone_dict = m.boundaries.bc['Rain_reduced']['zonal_dict']

    rch_zones = len(list(rch_zone_dict.keys()))

    par_rech_vals = [m.parameters.param['rchred{}'.format(i)]['PARVAL1'] \
                     for i in range(rch_zones - 1)]

    def update_recharge(vals):
        for key in list(interp_rain.keys()):
            for i in range(rch_zones - 1):
                interp_rain[key][recharge_zone_array == rch_zone_dict[i + 1]] = \
                    interp_rain[key][recharge_zone_array == rch_zone_dict[i + 1]] * \
                    vals[i]

            interp_rain[key][recharge_zone_array == rch_zone_dict[0]] = \
                interp_rain[key][recharge_zone_array == rch_zone_dict[0]] * 0.0
            interp_rain[key][m.model_mesh3D[1][0] == -1] = 0.
        return interp_rain

    interp_rain = update_recharge(par_rech_vals)
    rch = interp_rain

    
    m.boundaries.update_boundary_array('Rain_reduced', rch)

    if verbose:
        print("************************************************************************")
        print(" Updating Murray River GHB boundary")
    
    MurrayGHB = []

    MurrayGHB_cells = [[x[0], x[1], x[2], x[3]] for x in m.boundaries.bc['GHB']['bc_array'][0]]
    for MurrayGHB_cell in MurrayGHB_cells:
        lay, row, col = MurrayGHB_cell[:3]
        MurrayGHBstage = MurrayGHB_cell[3] # m.parameters.param['mghbst']['PARVAL1']
        # Don't use GHB cell is stage is below cell bottom.
        if MurrayGHBstage < m.model_mesh3D[0][lay + 1][row][col]:
            continue
       
        dx = m.gridHeight
        dz = m.model_mesh3D[0][lay][row][col] - \
            m.model_mesh3D[0][lay + 1][row][col]
        MGHBconductance = dx * dz * m.parameters.param['mghbk']['PARVAL1']
        MurrayGHB += [[lay, row, col, MurrayGHBstage, MGHBconductance]]

    ghb = {}
    ghb[0] = MurrayGHB

    if verbose:
        print("************************************************************************")
        print(" Updating GHB boundary")

    m.boundaries.update_boundary_array('GHB', ghb)

    if verbose:
        print("************************************************************************")
        print(" Updating Drains boundary")
    
    mapped_drains = m.polyline_mapped['Drain_Clip_model.shp']
    
    simple_drain = []
    drain_width_avg = 3.0 #m
    drain_bed_thickness = 0.10 #m
    for drain_cell in mapped_drains:
        row = drain_cell[0][0]
        col = drain_cell[0][1]
        if m.model_mesh3D[1][0][row][col] == -1:
            continue
        #print m.model_mesh3D
        drain_bed = m.model_mesh3D[0][0][row][col] - m.parameters.param['drndrp']['PARVAL1']
        drain_cond = drain_cell[1] * drain_width_avg * m.parameters.param['kv_drn']['PARVAL1'] / drain_bed_thickness
        simple_drain += [[0, row, col, drain_bed, drain_cond]]
    
    drain = {}
    drain[0] = simple_drain
    
    m.boundaries.assign_boundary_array('Drain', drain)
    
    if verbose:
        print("************************************************************************")
        print(" Updating Channels boundary")
    
#    simple_channel = []
#    channel_width_avg = 10.0 #m
#    channel_bed_thickness = 0.10 #m
#    for channel_cell in m.polyline_mapped['Channel_Clip_model.shp']:
#        row = channel_cell[0][0]
#        col = channel_cell[0][1]
#        if m.model_mesh3D[1][0][row][col] == -1:
#            continue
#        #print m.model_mesh3D
#        channel_stage = m.model_mesh3D[0][0][row][col]
#        channel_bed = m.model_mesh3D[0][0][row][col] - m.parameters.param['chan_drop']['PARVAL1']
#        channel_cond = channel_cell[1] * channel_width_avg * m.parameters.param['kv_chan']['PARVAL1'] / channel_bed_thickness
#        simple_channel += [[0, row, col, channel_stage, channel_cond, channel_bed]]
#    
#    channel = {}
#    channel[0] = simple_channel
#    
#    m.boundaries.assign_boundary_array('Channel', channel)

    if verbose:
        print("************************************************************************")
        print(" Updating pumping boundary")

    #pumpy = m.boundaries.bc['licenced_wells']['bc_array']
    #wel = {key: [[b[0], b[1], b[2], b[3] * pumping] for b in a] for key, a in pumpy.iteritems()}

    #m.boundaries.assign_boundary_array('licenced_wells', wel)

    if verbose:
        print("************************************************************************")
        print(" Check for boundary condition updating")
        m.generate_update_report()
    
    if verbose:
        print("************************************************************************")
        print(" Set initial head ")

    path=os.path.join(data_folder,"model_01_steady_state")
    fname="01_steady_state"
    head = flopyInterface.get_previous_conditions(os.path.join(path, fname))
    m.initial_conditions.set_as_initial_condition("Head", head)
    
    if verbose:
        print("************************************************************************")
        print(" Build and run MODFLOW model ")
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ## Currently using flopyInterface directly rather than running from the ModelManager ...

    modflow_model = flopyInterface.ModflowModel(m, data_folder=os.path.join(data_folder, "model_" + m.name))

    modflow_model.executable = mf_exe_folder

    modflow_model.headtol = 1E-5
    modflow_model.fluxtol = 1.

    modflow_model.buildMODFLOW(transport=True, write=True, verbose=False, check=False)

    converge = modflow_model.runMODFLOW(silent=True)

    if converge:
        if verbose:
            print('model converged')
            #break
    
        modflow_model.writeObservations()

        #modflow_model.waterBalanceTS()
    
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
        reaches = list(range(11)) # This needs to come automatically out of build!!
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
        
        #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        #^^^ POTENTIAL OBS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        def writePotentialObservations():
            
            # Set model output arrays to None to initialise
            head = None
            sfr_df = None
            stream_options = ['st_sim', 'fl_sim']
            head_options = ['shshal', 'shdeep']
            # Write observation to file
            for obs_set in list(modflow_model.model_data.observations.obs_group.keys()):
                obs_type = modflow_model.model_data.observations.obs_group[obs_set]['obs_type']
                # Import the required model outputs for processing
                if obs_type in head_options:
                    # Check if model outputs have already been imported and if not import
                    if head == None:
                        headobj = modflow_model.headobj
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
        # End writeObservations()
        
        writePotentialObservations() 

        #modflow_model.compareAllObs()

    return modflow_model

if __name__ == "__main__":

    verbose = False
                    
    args = sys.argv
    if len(args) > 1:
        model_folder = sys.argv[1]
        data_folder = sys.argv[2]
        mf_exe_folder = sys.argv[3]
        if len(args) > 4:
            param_file = sys.argv[4]
        else:
            param_file = ""
    else:
        # Get general model config information
        CONFIG = ConfigLoader('../../../config/model_config.json')\
                        .set_environment("02_transient_flow")
        model_config = CONFIG.model_config
        model_folder = model_config['model_folder'] + model_config['grid_resolution'] + os.path.sep
        data_folder = model_config['data_folder']
        mf_exe_folder = model_config['mf_exe_folder']
        param_file = model_config['param_file']
        

    if param_file:
        run = run(model_folder, data_folder, mf_exe_folder, param_file=param_file, verbose=verbose)
    else:
        run = run(model_folder, data_folder, mf_exe_folder, verbose=verbose)

       
