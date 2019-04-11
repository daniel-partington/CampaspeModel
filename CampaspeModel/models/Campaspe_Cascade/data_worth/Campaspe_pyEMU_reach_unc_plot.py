import os
import datetime

import numpy as np
import pandas as pd
import pickle
from itertools import chain, combinations
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import colors as mpl_colors
from matplotlib import gridspec   
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
     
from PIL import Image

import flopy
from flopy.utils.sfroutputfile import SfrFile
from flopy.utils import sfroutputfile
import pyemu

from HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface

# Configuration Loader
from HydroModelBuilder.Utilities.Config.ConfigLoader import ConfigLoader

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set up some nice colours for figures
# These are the "Tableau 20" colors as RGB.    

colors = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(colors)):    
    r, g, b = colors[i]    
    colors[i] = (r / 255., g / 255., b / 255.)  

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def unix2dos(path):
    text = open(path, "U").read() 
    text = text.replace("\n", "\r\n") 
    open(path, "wb").write(text)
    
def load_obj(filename):
    """
    :param filename:
    """

    if filename.endswith('.pkl'):
        with open(filename, 'rb') as f:
            print "Loading: ", f, filename
            p = pickle.load(f)
            return p
        # End with
    else:
        raise TypeError('File type not recognised as "pkl": {}'.format(filename))
    # End if
# End load_obj()
    
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
def update_forecasts_in_pest_file(pst, forecasts=None):
    # Forecasts here are those with name containing nrf = net river swgw flux, and rrf = reach river swgw flux
    # The other parts of these obs names refers to the temporal aspect, i.e. a = annual, s = seasonal (0, 1, 2, 3), and m = monthly range(12)
    
    with open('pest_emu.pst', 'r') as f:
        txt = f.read()
        if "++forecasts" in txt:
            print("forecasts written already")
        else:
            df = pst.observation_data
            forecasts = df[df['obgnme'].str.contains('rrf') | df['obgnme'].str.contains('nrf')]['obsnme'].tolist()
            with open('pest_emu.pst', 'a') as f:
                f.write("++forecasts({})".format(",".join(forecasts)))
    if forecasts:
        return forecasts

def forecast_names(pst_df):
    '''Set up conversion of obsnme to useable prediction name'''
    forecasts_df = pst_df[pst_df['obgnme'].str.contains('rrf') | pst_df['obgnme'].str.contains('nrf')]
    forecasts_df_count_obgnme = forecasts_df.groupby('obgnme').count()
    forecasts_df_count_obgnme = forecasts_df_count_obgnme[forecasts_df_count_obgnme['obsnme'] > 1]
    for obgp in forecasts_df_count_obgnme.index:
        new_names = ["{}_{}".format(obgp, x) for x in range(forecasts_df_count_obgnme.loc[obgp,'obsnme'])]
        forecasts_df.loc[forecasts_df['obgnme'] == obgp, 'unique_name'] = new_names
    unique_null = pd.isnull(forecasts_df['unique_name'])
    forecasts_df.loc[unique_null, 'unique_name'] = forecasts_df.loc[unique_null, 'obgnme']        
    return forecasts_df
        
def adjust_potential_by_existing(pst, pot_existing_map={'shshal':['head1', 'head3'], 'shdeep':['head5', 'head6'],
     'st_sim':['stage'], 'fl_sim':['gflow'],
     'c14shal':['c14'], 'c14deep':['c14'], 'rn_sim':['radon'],
     'ec_sim':['gstrec', 'fstrec']}):    
    
    for pot in pot_existing_map:
        pst.observation_data.loc[pst.observation_data['obgnme'] == pot, 'weight'] = \
            np.array([pst.observation_data[pst.observation_data['obgnme'] == exist]['weight'].mean() for exist in pot_existing_map[pot]]).mean()
    print pst.observation_data
    return pst

def load_pest_file_and_jacobian(model_folder, res_file=None):
    jco = os.path.join(model_folder, "pest_emu.jco")
    pst_file = jco.replace(".jco", ".pst")
    #unc_file = jco.replace(".jco", ".unc")
    unix2dos(pst_file)
    # SETUP PST OBJECT
    pst = pyemu.Pst(pst_file)             

    if res_file:
        pst.observation_data.loc[pst.observation_data['obgnme'] == 'radon', 'weight'] = 20.
        pst.adjust_weights_resfile(resfile=os.path.join(model_folder, res_file))
        pst = adjust_potential_by_existing(pst)
        
    return jco, pst
        
def load_pest_file_and_jacobian_into_Schur(model_folder, forecasts=None, res_file=None, pst_obj=None):         
    
    if not pst_obj:
        jco, pst = load_pest_file_and_jacobian(model_folder, res_file=res_file)
        
    if forecasts:
        # force weights of forecasts to 0
        pst.observation_data.loc[pst.observation_data.index.isin(forecasts), 'weight'] = 0.0
        la = pyemu.Schur(jco=jco, pst=pst, forecasts=forecasts)       
    else:
        forecasts = update_forecasts_in_pest_file(pst)                        
        la = pyemu.Schur(jco=jco, pst=pst, forecasts=forecasts)       
    
    return la    
    
def _zone_array2layers(zone_array, plots=False):
    '''
    Function to generate 2D masked layer arrays for each zone
    '''
    zones = [x + 1 for x in range(len(np.unique(zone_array)) - 1)]
    layers = zone_array.shape[0]
    zone_mask2D = {}
    zone_top2D = {}

    for index, zone in enumerate(zones):
        zone_mask = np.ma.masked_array(zone_array, 
                                       zone_array == zone).mask
    
        zone_mask2D[index] = np.full_like(zone_mask[0], False, dtype=bool)
        zone_top2D[index] = np.full_like(zone_mask[0], 0, dtype=int)
        for layer in range(layers):
            zone_mask2D[index] |= zone_mask[layer] 
    
        for layer in range(layers, 0, -1):
            zone_top2D[index][zone_mask[layer - 1]] = layer - 1

    if plots:
        import matplotlib.pyplot as plt
        for index, zone in enumerate(zones):
            plt.figure()
            plt.imshow(np.ma.masked_array(zone_top2D[index], ~zone_mask2D[index]), interpolation='none')

    return zone_mask2D, zone_top2D    

def get_model(config_path=None):
    if config_path: 
        CONFIG = ConfigLoader(os.path.join(config_path,'model_config.json'))\
                        .set_environment("02_transient_flow")
    else:
        CONFIG = ConfigLoader('../../../config/model_config.json')\
                        .set_environment("02_transient_flow")
    model_config = CONFIG.model_config
    model_folder = model_config['model_folder'] + model_config['grid_resolution'] + os.path.sep
    data_folder = model_config['data_folder']
    param_file = model_config['param_file']
    
    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder, "02_transient_flow_packaged.pkl"))
    name = MM.GW_build.keys()[0]
    m = MM.GW_build[name]
    return m, data_folder
    
def get_model_build(m, data_folder):
    modflow_model = flopyInterface.ModflowModel(m, data_folder=data_folder)
    modflow_model.buildMODFLOW(transport=True, write=False)
    return modflow_model

def get_heads_for_sfr(sfr_df, modflow_model, dry_nan=False):
    headobj = modflow_model.import_heads_from_file(path=os.path.join(data_folder,"model_02_transient_flow"), name="02_transient_flow")
    heads = headobj.get_alldata()
    sfr_df.loc[:, 'head'] = sfr_df.apply(lambda x: heads[x['time']][x['layer'] - 1][x['row'] - 1][x['column'] - 1], axis=1)
    if dry_nan:
        sfr_df['head'][sfr_df['head'] < -999.] = np.nan 
    return sfr_df
    
def get_stream_info_and_sim(m, modflow_model, data_folder):    
    sfr = sfroutputfile.SfrFile(os.path.join(os.path.join(data_folder,"model_02_transient_flow"), modflow_model.sfr.file_name[0] + '.out'))
    sfr_df = sfr.get_dataframe()
    #fn = lambda x: x.iloc[0]
    #sfr_df_group = sfr_df.groupby(by='time').agg({'Qaquifer':np.sum,'Qin':fn, 'Qovr':np.sum, 'Qprecip':np.sum, 'Qet':np.sum})
    #sfr_df_group.plot(x=date_index[1:], y=['Qaquifer', 'Qin', 'Qovr', 'Qprecip', 'Qet'])
    sfr_info = m.river_mapping['Campaspe']
    cum_len = sfr_info['Cumulative Length'].tolist()
    sfr_df.loc[:, 'Cumulative Length'] = cum_len * (sfr_df['time'].max() + 1)
    return sfr_info, sfr_df
    
if __name__ == '__main__':

    # Setup all of the requisite data for Data Worth Analysis and pretty plotting
    p_j = os.path.join

    model_folder = r"C:\Workspace\part0075\MDB modelling\testbox\PEST1000\HPC_Feb2018"
    save_folder = p_j(model_folder, r"adjust")
    observations_folder = r"C:\Workspace\part0075\MDB modelling\testbox\PEST1000\master\model_02_transient_flow"
    pst_file = "pest_emu.pst"
    os.chdir(model_folder)    

    # To recalculate for the spatio analysis
    force_recalc = False

    if not os.path.exists(save_folder):
        os.mkdir(save_folder)
    # end if

    obs_colour_map={'head':'#1b9e77',
                    'stage':'#d95f02',
                    'flow':'#7570b3',
                    'radon':'#e7298a',
                    'c14':'#66a61e',
                    'EC':'#e6ab02'}

    obs_marker_map={'head':'^',
                    'stage':'o',
                    'flow':'o',
                    'radon':'o',
                    'c14':'*',
                    'EC':'P'}

    # Some useful bits from the simulation outputs

    field_sampling = [datetime.datetime(2016,03,31),
                  datetime.datetime(2016,12,31),
                  datetime.datetime(2017,04,30)]

    m, data_folder = get_model(config_path=r"C:/Workspace/part0075/GIT_REPOS/CampaspeModel/CampaspeModel/config")
    date_index = m.model_time.t['dateindex']
    modflow_model = get_model_build(m, data_folder)
    mesh3D = load_obj("mesh_stuff.pkl")
    zone2D_info = _zone_array2layers(mesh3D[1])
    sfr_info, sfr_df = get_stream_info_and_sim(m, modflow_model, data_folder)
    get_heads_for_sfr(sfr_df, modflow_model, dry_nan=True)
    sfr_df.loc[:,'reach_dw'] = sfr_info['reach'].tolist() * len(sfr_df['time'].unique())
    sfr_df.loc[:,'rchlen'] = sfr_info['rchlen'].tolist() * len(sfr_df['time'].unique())
    sfr_df.loc[:,'rchlen_perc'] = sfr_df.loc[:,'rchlen'] / sfr_info['rchlen'].sum() 
    sfr_df.loc[:,'area'] = sfr_df['width'] * sfr_df['rchlen']
    # Sort by reaches:
    sfr_df_reaches = pd.DataFrame()
    for time in sfr_df['time'].unique():
        sfr_df_gp = sfr_df[sfr_df['time'] == time].groupby('reach_dw').agg({'rchlen':'sum', 'rchlen_perc':'sum', 'Qaquifer':'sum', 'area':'sum', 'head':'mean', 'stage':'mean', 'Cumulative Length':'max'})
        sfr_df_gp.loc[:, 'time'] = time
        sfr_df_gp.loc[:, 'reach_dw'] = sfr_df_gp.index
        #print sfr_df_gp
        if sfr_df_reaches.empty:
            sfr_df_reaches = sfr_df_gp
        sfr_df_reaches = pd.concat([sfr_df_reaches, sfr_df_gp], ignore_index=True)     
        
    sfr_df_reaches_relevant = sfr_df_reaches[sfr_df_reaches['time'].isin(sfr_df_reaches['time'].unique()[-12:])]    
    sfr_df_reaches_relevant.loc[:, 'Qaquifer_lin'] = sfr_df_reaches_relevant['Qaquifer'] / sfr_df_reaches_relevant['rchlen']
    sfr_df_reaches_relevant.loc[:, 'Qaquifer_lin2'] = sfr_df_reaches_relevant['Qaquifer'] / sfr_df_reaches_relevant['area']

    # Lake Eppalock Inflow:
    Eppalock_inflow = pd.read_csv("Eppalock_inflow.csv", skiprows=2, index_col='Date', parse_dates=True, dayfirst=True)
    qual_codes_to_ignore = [8, 9, 21, 100, 101, 120, 149, 150, 151, 152,
                                153, 154, 155, 156, 160, 161, 165, 180, 190,
                                200, 201, 237, 250, 254, 255]

    Eppalock_inflow = Eppalock_inflow[~Eppalock_inflow['Qual'].isin(qual_codes_to_ignore)]

    flow_high_quantile = 0.8
    flow_low_quantile = 0.35

    Eppalock_inflow_monthly = Eppalock_inflow.resample('M').mean()
    q_lower = Eppalock_inflow['Mean'].quantile(flow_low_quantile)
    q_upper = Eppalock_inflow['Mean'].quantile(flow_high_quantile)
    Eppalock_inflow_monthly.loc[:, 'flow_group'] = np.nan    
    Eppalock_inflow_monthly.loc[(Eppalock_inflow_monthly['Mean'] < q_lower).tolist(), 'flow_group'] = 'low'
    Eppalock_inflow_monthly.loc[(Eppalock_inflow_monthly['Mean'] > q_upper).tolist(), 'flow_group'] = 'high'
    Eppalock_inflow_monthly.loc[pd.isnull(Eppalock_inflow_monthly['flow_group']), 'flow_group'] = 'regular'

    flow_type_colours = {'low':'orangered', 'high':'dodgerblue', 'regular':'mediumseagreen'}                                
    #flow_type_colours = {'low':'#a1dab4', 'high':'#2c7fb8', 'regular':'#41b6c4'}                                
    Eppalock_inflow_monthly.loc[:, 'colour'] = Eppalock_inflow_monthly.apply(lambda x: flow_type_colours[x['flow_group']], axis=1)
    Eppalock_inflow_monthly = Eppalock_inflow_monthly.append(pd.DataFrame(data={'Mean': [Eppalock_inflow_monthly['Mean'][-1]], 
                                                            'Qual':[0], 
                                                            'colour':[Eppalock_inflow_monthly['colour'][-1]],
                                                            'flow_group':[Eppalock_inflow_monthly['flow_group'][-1]]}, 
                                                      index=[pd.datetime(2017,5,31)]))
                                      
    date_range = pd.date_range(start=datetime.datetime(2016,6,1), end=datetime.datetime(2017,5,30), freq='D')
    date_range_month = pd.date_range(start=datetime.datetime(2016,6,1), end=datetime.datetime(2017,5,31), freq='M')
    
    '''
    pyEMU time ...
    '''    
    
    la = load_pest_file_and_jacobian_into_Schur(model_folder, res_file='pest.rei.11')
    #la2 = load_pest_file_and_jacobian_into_Schur(model_folder, res_file='pest.rei.11')
    forecasts_df = forecast_names(la.pst.observation_data)
    forecasts_df.loc[:, 'var'] = [la.posterior_forecast[row] for row in forecasts_df.index]
    forecasts_df.loc[:, 'unc'] = forecasts_df['var'].apply(np.sqrt)

    # Exclude seasonal predictions:
    pred_categories = ['nrf_a', 'nrf_m', 'rrf_a', 'rrf_m']
    marker_dict_predtype = {'nrf_a': 'D', 'nrf_m': 'D', 'rrf_a': 'o' , 'rrf_m':'o'}
    markerfill_dict_predtype = {'nrf_a': True, 'nrf_m': False, 'rrf_a': True , 'rrf_m':False}
    tick_offsets = {'nrf_a': -0.26, 'nrf_m': -0.16, 'rrf_a': 0.06, 'rrf_m':0.16}
    pred_flow_types = Eppalock_inflow_monthly.loc[date_range_month, 'flow_group'].tolist()    
    pred_colours = Eppalock_inflow_monthly.loc[date_range_month, 'colour'].tolist()    

    rows_name2 = ['Annual', 'Jun, 16', 'Jul, 16', 'Aug, 16',
                 'Sep, 16', 'Oct, 16', 'Nov, 16', 'Dec, 16', 'Jan, 17', 'Feb, 17', 'Mar, 17', 'Apr, 17', 'May, 17']

    months = rows_name2[1:]
    reaches = 11

    pred_exchange = pd.DataFrame(index= ["rrf_m{}_{}".format(x, y) for x, y in zip(range(reaches), range(12))] +\
                                          ["rrf_a{}".format(x) for x in range(reaches)] +\
                                          ['nrf_m_{}'.format(x) for x in range(12)] +\
                                          ['nrf_a'] , columns=['Exchange'])
    pred_exchange.loc['nrf_a', 'Exchange'] = sfr_df_reaches_relevant.groupby(by='time').mean()['Qaquifer_lin2'].mean()
    for i in range(reaches):
        pred_exchange.loc["rrf_a{}".format(i), 'Exchange'] = \
             sfr_df_reaches_relevant[sfr_df_reaches_relevant['reach_dw'] == i]['Qaquifer_lin2'].mean()
    # end for 
    for i in range(12):
        pred_exchange.loc['nrf_m_{}'.format(i), 'Exchange'] = \
             sfr_df_reaches_relevant[sfr_df_reaches_relevant['time'] == 21 + i]['Qaquifer_lin2'].mean()
    # end for
    for i in range(reaches):
        for j in range(12):
            sfr_df_reaches_relevant_reach = \
                 sfr_df_reaches_relevant[sfr_df_reaches_relevant['reach_dw'] == i]
            pred_exchange.loc["rrf_m{}_{}".format(i, j), 'Exchange'] = \
                 sfr_df_reaches_relevant_reach[sfr_df_reaches_relevant_reach['time'] == 21 + j]['Qaquifer_lin2'].mean()
        # end for
    # end for    
    #exchange_type = {'gaining': , 'neutral':, 'losing':}
    
    max_exchange = pred_exchange.max()
    min_exchange = pred_exchange.max()
    mean_exchange = pred_exchange.mean()
    neutral_criterion = -0.004 #0.2 * mean_exchange.tolist()[0]
    losing = pred_exchange[pred_exchange['Exchange'] > abs(neutral_criterion)]
    gaining = pred_exchange[pred_exchange['Exchange'] < neutral_criterion]
    neutral = pred_exchange[(pred_exchange['Exchange'] < abs(neutral_criterion)) & (pred_exchange['Exchange'] > neutral_criterion)]

    #_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
    # Uncertatinty for all the forecast observations 

    #forecasts_df[forecasts_df['obgnme'].str.contains('rrf_m')].head()
    sfr_df_reaches_relevant.loc[:, 'unc'] = sfr_df_reaches_relevant.apply(
        lambda x: forecasts_df[forecasts_df['unique_name'] == "rrf_m{}_{}".format(
        int(x['reach_dw']), int(x['time'] - 21))]['unc'].tolist()[0], axis=1)
    
    sfr_df_reaches_relevant.loc[:, 'unc_lin'] = sfr_df_reaches_relevant['unc'] / sfr_df_reaches_relevant['rchlen']
    sfr_df_reaches_relevant.loc[:, 'Cumulative Length (km)'] = sfr_df_reaches_relevant['Cumulative Length'] / 1000.
    # Plot the boxplots up
        
    import brewer2mpl
    bmap = brewer2mpl.get_map('Set1', 'Qualitative', 5)
    palette = bmap.mpl_colors
    
    #plt.close('all')         
    fig2, ax = plt.subplots(2,1, figsize=(12, 6), facecolor='w', edgecolor='k')
    fig2.tight_layout(rect=[0.02,0.05,1,0.99])
    plt.subplots_adjust(hspace = 0.1)
    counter = 0
    
    def gen_synthetic_data(mu, sigma, size):
        return np.random.normal(mu, sigma, size)
    
    Camp_weir_seg = m.observations.obs_group['gflow']['locations'].loc[406203, 'seg_loc']
    sfr_df.loc[:, 'up_down'] = sfr_df.apply(lambda x: 0 if x['segment'] < Camp_weir_seg else 1, axis=1)
    sfr_df_reaches_campweir = pd.DataFrame()
    for time in sfr_df['time'].unique():
        sfr_df_gp = sfr_df[sfr_df['time'] == time].groupby('up_down').agg({'rchlen':'sum', 'rchlen_perc':'sum', 'Qaquifer':'sum', 'area':'sum', 'head':'mean', 'stage':'mean', 'Cumulative Length':'max'})
        sfr_df_gp.loc[:, 'time'] = time
        sfr_df_gp.loc[:, 'reach_ud'] = sfr_df_gp.index
        #print sfr_df_gp
        if sfr_df_reaches_campweir.empty:
            sfr_df_reaches_campweir = sfr_df_gp
        sfr_df_reaches_campweir = pd.concat([sfr_df_reaches_campweir, sfr_df_gp], ignore_index=True)     
        
    sfr_df_reaches_campweir_relevant = sfr_df_reaches_campweir[sfr_df_reaches_campweir['time'].isin(sfr_df_reaches_campweir['time'].unique()[-12:])]    
    sfr_df_reaches_campweir_relevant.loc[:, 'Qaquifer_lin'] = sfr_df_reaches_campweir_relevant['Qaquifer'] / sfr_df_reaches_campweir_relevant['rchlen']
    sfr_df_reaches_campweir_relevant.loc[:, 'Qaquifer_lin2'] = sfr_df_reaches_campweir_relevant['Qaquifer'] / sfr_df_reaches_campweir_relevant['area']
    
    with open(p_j(save_folder, "Up and downstream of Campapspe weir avg exchange.txt"), 'w') as f:    
        f.write("Lake Eppalock to Campaspe Weir: {} m^3/m/dat".format(
                sfr_df_reaches_campweir_relevant[sfr_df_reaches_campweir_relevant['reach_ud'] == 0]['Qaquifer_lin'].mean()))
        f.write("Campaspe Weir to Murray River: {} m^3/m/dat".format(
                sfr_df_reaches_campweir_relevant[sfr_df_reaches_campweir_relevant['reach_ud'] == 1]['Qaquifer_lin'].mean()))
        
    for i in range(0,33): #sfr_df['time'].max()):
        if date_index[i] in field_sampling[1:]:    
            ax[counter].axvline(x=73.4, color='black', lw=2) 
            ttes_data = []
            mean_GW_Inpt = []
            compt_width = 0
            data = sfr_df_reaches_relevant[sfr_df_reaches_relevant['time'] == i]
            cum_length = [0.] + data['Cumulative Length (km)'].tolist()
            positions = [(cum_length[j + 1]
                           + cum_length[j]) / 2. 
                          for j in range(0, len(cum_length) - 1)]
            widths = [(cum_length[j + 1]
                      - cum_length[j]) * 0.96
                          for j in range(0, len(cum_length) - 1)]
            for k, row in data.iterrows():
                ttes_data += [gen_synthetic_data(-row['Qaquifer_lin'], row['unc_lin'], 1000)]
                #mean_GW_Inpt += [np.array(row.dropna(), dtype=float).mean() * widths[compt_width]]
                compt_width+=1
            box_param  = ax[counter].boxplot(ttes_data, positions=positions, patch_artist=True, widths=widths, showfliers=False)
            for box in box_param['boxes']:
                box.set(facecolor=palette[counter + 1], alpha=0.3)
            for median in box_param['medians']:
                median.set(color=palette[counter + 1], linewidth=2)
            ax[counter].set_ylabel('GW Input (m$^3$/m/day)')
            #ax[counter].set_ylim([-5, 5])
            ax[counter].set_xticks(np.arange(0, 160, 20))
            ax[counter].set_xticklabels([])
            ax[counter].set_xticklabels(np.arange(0, 160, 20))
            counter += 1
    
    ax[-1].set_xlabel('Distance from Lake Eppalock (km)')
    handles = []
    labels = ['03-04 / 2016', '12 / 2016', '05 / 2017', '01 / 2018']
    for Camp in labels[1:3]:
        handles.append(mpatches.Patch(color=palette[labels.index(Camp)], alpha=0.5))
    legend = fig2.legend(handles = handles, labels = labels,  ncol = 4, 
                         prop = {'size': 12}, frameon = True, numpoints=1,
                         loc = 'center', bbox_to_anchor=(0.52, 0.02))
    plt.savefig(p_j(save_folder, 'SWGWexch_and_unc.png'), dpi = 300)