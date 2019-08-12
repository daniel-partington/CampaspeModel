"""Update parameters and run transient flow model for Campaspe."""

import os
import sys

import numpy as np

#sys.path.append('C:\Workspace\part0075\GIT_REPOS')
from HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface
from HydroModelBuilder.Utilities.model_assessment import plot_obs_vs_sim

## Configuration Loader
#from HydroModelBuilder.Utilities.Config.ConfigLoader import CONFIG

# Configuration Loader
from HydroModelBuilder.Utilities.Config.ConfigLoader import ConfigLoader

def run(model_folder, data_folder, mt_exe_folder, param_file=None, verbose=True, plots=False):
    """Model Runner."""

    # MM is short for model manager
    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder, "02_transient_flow_packaged.pkl"))

    name = list(MM.GW_build.keys())[0]
    m = MM.GW_build[name]

    # Load in the new parameters based on parameters.txt or dictionary of new parameters

    if param_file:
        MM.GW_build[name].updateModelParameters(os.path.join(data_folder, 'parameters.txt'), verbose=verbose)

    if verbose:
        print("************************************************************************")
        print(" Instantiate MODFLOW model ")

    ###########################################################################
    ###########################################################################
    ###########################################################################
    # Currently using flopyInterface directly rather than running from the ModelManager ...
    modflow_model = flopyInterface.ModflowModel(MM.GW_build[name], 
                                                data_folder=os.path.join(data_folder, "model_" + m.name))
    modflow_model.buildMODFLOW(transport=True, write=False, verbose=False)


    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
    #@@@ SIMPLE RADON MODEL @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 

    num_reaches = m.pilot_points['Campaspe'].num_points #4
    known_points = m.pilot_points['Campaspe'].points

    sfr_df = modflow_model.importSfrOut()
    sfr_info = m.river_mapping['Campaspe']
    cum_len = sfr_info['Cumulative Length'].tolist()
    sfr_df.loc[:, 'Cumulative Length'] = cum_len * (sfr_df['time'].max() + 1)
    
    # Hyporheic zone depth         
    hz_depth_vals = [m.parameters.param['hzdpth{}'.format(x)]['PARVAL1'] for \
                                        x in range(num_reaches)] 
    R_depth_HZ = np.interp(sfr_info['Cumulative Length'].tolist(), 
                           known_points, 
                           hz_depth_vals)

    df_size = sfr_info.shape[0]
    t_steps = sfr_df['time'].max() + 1
    # Hyporheic zone porosity
    sfr_df.loc[:, 'HZ_poro'] = [m.parameters.param['hzporo']['PARVAL1']] \
                               * df_size * t_steps
    # Hyporheic zone production of radon
    sfr_df.loc[:, 'HZ_Prod_Rate'] = [m.parameters.param['hzprod']['PARVAL1']] \
                                    * df_size * t_steps
    # Hyporheic zone residence time
    sfr_df.loc[:, 'HZ_RTime'] = [m.parameters.param['hz_rt']['PARVAL1']] \
                                * df_size * t_steps
    # Depth of the hyporheic zone
    sfr_df.loc[:, 'R_depth_HZ'] = R_depth_HZ.tolist() * t_steps              
    # Gas transfer velocity
    sfr_df.loc[:, 'GTV'] = [m.parameters.param['gtv']['PARVAL1']] \
                           * df_size * t_steps
    # Groundwater radon concentration
    sfr_df.loc[:, 'GW_Rn_conc'] = [m.parameters.param['gwconc']['PARVAL1']] \
                                  * df_size * t_steps
    # Groundwater EC
    sfr_df.loc[:, 'GW_EC'] = [m.parameters.param['gwecconc']['PARVAL1']] \
                                  * df_size * t_steps

    # EC of the inflowing tributary water if present
    sfr_df.loc[:, 'Trib_EC'] = [10.] * df_size * t_steps # ARBITRARY!!!! 
    # Radon concentration of inflowing tributary water if present
    sfr_df.loc[:, 'Trib_Rn'] = [10.] * df_size * t_steps # ARBITRARY!!!!
    # Reach lengths
    sfr_df.loc[:, 'dx'] = sfr_info['rchlen'].tolist() * t_steps 

    Radon_obs = m.observations.obs_group['radon']
    Radon_obs2 = m.observations.obs_group['rn_sim']
    Radon_obs_ts = Radon_obs['time_series']
    Radon_obs_ts2 = Radon_obs2['time_series']

    intervals_of_interest = np.unique(Radon_obs_ts['interval'].unique().tolist() + \
                                      Radon_obs_ts2['interval'].unique().tolist())

    ec_groups = ['fstrec', 'gstrec', 'ec_sim']
    for ec_group in ec_groups:
        ec_obs = m.observations.obs_group[ec_group]
        ec_obs_ts = ec_obs['time_series']
    
        intervals_of_interest = np.unique(intervals_of_interest.tolist() + \
                                      ec_obs_ts['interval'].unique().tolist())

    radon_df_dict = {}
    for i in intervals_of_interest:
        df = sfr_df[sfr_df['time'] == i]
        Ini_cond = (df.iloc[0]['Qin'], 0., 300.)
        df.loc[:, 'Flow'], df.loc[:, 'Rn'], df.loc[:, 'EC'] = \
              modflow_model.Calculate_Rn_from_SFR_with_simple_model(df, 
                                                                    Ini_cond)
        radon_df_dict[i] = df

    def write_obs(obs_set, col_interest):
        obs_df = m.observations.obs_group[obs_set]['time_series']
        obs_df = obs_df[obs_df['active'] == True]
        sfr_location = m.observations.obs_group[obs_set]['locations']['seg_loc']
        with open(os.path.join(modflow_model.data_folder, 'observations_{}.txt'.format(obs_set))
                  , 'w') as f:
            for observation in obs_df.index:
                interval = int(obs_df['interval'].loc[observation])
                name = obs_df['name'].loc[observation]
                seg = sfr_location.loc[name]
                col_of_interest = col_interest
                df_radon = radon_df_dict[interval] 
                sim_obs = df_radon[df_radon['segment'] == seg] \
                                  [col_of_interest]
                f.write('%f\n' % sim_obs)                
                
    radon_obs_sets = ['radon', 'rn_sim']
    for obs_set in radon_obs_sets:
        write_obs(obs_set, 'Rn')

    for obs_set in ec_groups:
        write_obs(obs_set, 'EC')
        
    post = flopyInterface.MT3DPostProcess(modflow_model, 
                                              mt_name='Radon')
        
    if plots:
        def compareObs(post, obs_set, col_interest):
            obs_df = m.observations.obs_group[obs_set]['time_series']
            obs_df = obs_df[obs_df['active'] == True]
            sfr_location = m.observations.obs_group[obs_set]['locations']['seg_loc']
            obs_sim_zone_all = []
            for observation in obs_df.index:
                interval = int(obs_df['interval'].loc[observation])
                name = obs_df['name'].loc[observation]
                obs = obs_df['value'].loc[observation]
                seg = sfr_location.loc[name]
                col_of_interest = col_interest
                df_radon = radon_df_dict[interval] 
                                
                sim_obs = df_radon[df_radon['segment'] == seg] \
                                  [col_of_interest].tolist()[0]
                obs_sim_zone_all += [[obs, sim_obs, seg]]                
                    
            plot_obs_vs_sim(obs_set, obs_sim_zone_all, unc=2)
        # End compareObsRn
        compareObs(post, 'radon', 'Rn')    

        for ec_sim_obs_plot in ec_groups[:-1]:
            compareObs(post, ec_sim_obs_plot, 'EC')  

if __name__ == "__main__":

    verbose = False

    args = sys.argv
    if len(args) > 1:
        model_folder = sys.argv[1]
        data_folder = sys.argv[2]
        mt_exe_folder = sys.argv[3]
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
        mt_exe_folder = model_config['mt_exe_folder']
        param_file = model_config['param_file']

    if param_file:
        run = run(model_folder, data_folder, mt_exe_folder, param_file=param_file, verbose=verbose, plots=verbose)
    else:
        run = run(model_folder, data_folder, mt_exe_folder, verbose=verbose)