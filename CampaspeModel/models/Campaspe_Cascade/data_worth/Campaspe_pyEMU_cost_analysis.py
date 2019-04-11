import os
import numpy as np
import pandas as pd
import pickle
from itertools import chain, combinations
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

import pyemu

from HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface

# Configuration Loader
from HydroModelBuilder.Utilities.Config.ConfigLoader import ConfigLoader

from helper_functions import *

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

    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    new_weights_dict = {
        # Existing monitoring network
        'head6':  1./1.0, #0.2, # This is representative only of measurement error and then some
        'head5':  1./1.0, #0.2, 
        'head3':  1./1.0, #0.2, 
        'head1':  1./1.0, #0.2, 
        'stage':  1./1.0, #0.1, 
        'gstrec': 1./100., 
        # Field work
        'fstrec': 1./100., 
        'fdepth': 1./0.5, 
        'radon':  1./200.,
        'c14':    1./10., #5.,
        # Potential Observations
        'c14shal': 1./10., #5.,
        'c14deep': 1./10., #5.,
        'shshal': 1./1.0, #0.2,
        'shdeep': 1./1.0, #0.2,
        'st_sim': 1./1.0, #0.1,
        'ec_sim': 1./100.,
        'rn_sim': 1./200.
        }

    new_weights_perc_dict = {
        # Existing monitoring network
        'gflow':  0.4,
        # Field work
        'fflow':  0.4, 
        # Potential Observations
        'fl_sim': 0.4,
        }

    new_par_bounds = {
        # Campaspe River parameter groups for SFR package
        'rivwdt': [], # width
        'kv_riv': [], # bed hydraulic conductivity
        'bedthk': [], # bed thickness
        'rivbed': [], # River bed
        'rough': [],# River bed roughness
        # Murray River parameter groups for RIV package
        'murriv': [],
        # Transport group parameters for C14
        'transp': [],
        # Transport group parameters for EC
        'ec': [],
        # Transport group parameters for Radon
        'radon': [],
        # Recharge parameters
        'rchred': [], # Transient sim recharge params
        'ssrch': [], # Steady state sim recharge params
        # Specific yield
        'syqa': [],
        'syutqa': [],
        'syutaf': [],
        'syutam': [],
        'syutb': [],
        'sylta': [],
        'sybse': [],
        # Drain parameters group
        'drain': [],
        # GHB parameters group
        'ghb': [],
        # Specific storage parameters groups
        'ssqa': [],
        'ssutqa': [],
        'ssutaf': [],
        'ssutam': [],
        'ssutb': [],
        'sslta': [],
        'ssbse': [],
        # Hydraulic conductivity horizontal (vertical was set at 10% of horizontal)
        'kqa': [],
        'kutqa': [],
        'kutaf': [],
        'kutam': [],
        'kutb': [],
        'klta': [],
        'kbse': [],
                  }           
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    #--D-A-T-A-W-O-R-T-H---T-I-M-E----------------------------------------------
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    
    '''
    pyEMU time ...
    '''    
    
    la = load_pest_file_and_jacobian_into_Schur(model_folder, res_file='pest.rei.11')
    forecasts_df = forecast_names(la.pst.observation_data)
    forecasts_df.loc[:, 'var'] = [la.posterior_forecast[row] for row in forecasts_df.index]
    forecasts_df.loc[:, 'unc'] = forecasts_df['var'].apply(np.sqrt)

    # Some usefule definitions of what is old, new and potential data    
    old_data = ['head1', 'head3', 'head5', 'head6', 'stage', 'gflow', 'gstrec']
    new_field_data = ['fstrec', 'radon', 'c14'] #ignore ['fflow, 'fdepth'] 
    potential_data = ['st_sim', 'fl_sim', 'rn_sim', 'ec_sim', 'shshal', 'shdeep', 'c14shal', 'c14deep']
    potential_data_stream = ['st_sim', 'fl_sim', 'rn_sim', 'ec_sim']
    potential_data_ground = ['shshal', 'shdeep', 'c14shal', 'c14deep']

    #_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
    # ANALYSE SIMPLE COMBOS IN TRADITIONAL VS NOVEL DATA TYPES
    
    obs_groups_name = [
                       "h", 
                       "h_s",
                       "h_s_f", 
                       "h_s_f_c14", 
                       "h_s_f_radon", 
                       "h_s_f_stream_ec",
                       "all" 
                       ]
    
    obs_groups = [
                  ['head1', 'head3', 'head5', 'head6'], 
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'fdepth'],
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'fdepth', 'gflow', 'fflow'], 
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'fdepth', 'gflow', 'fflow', 'c14'], 
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'fdepth', 'gflow', 'fflow', 'radon'], 
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'fdepth', 'gflow', 'fflow', 'gstrec', 'fstrec'], 
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'fdepth', 'gflow', 'fflow', 'c14', 'radon', 'gstrec', 'fstrec']]

    combo_obs_dict = create_combos(la, obs_groups_name, obs_groups)

    df_worth_add_combo, df_worth_add_combo_unc, df_unc_perc_add_combo = process_combos(la, combo_obs_dict, obs_groups_name)
    ob_interest = 'nrf_a'

    #_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
    # ANALYSE ALL COMBOS IN THE CONTEXT OF COST
    
    # Let's look at ALL combos for DATA additions

    obs_groups_types = ['head', 'stage', 
                        'flow', 'c14', 'radon', 
                        'strec'] 
    
    obs_groups_types_map = {'head': ['head1', 'head3', 'head5', 'head6'], #, 'shshal', 'shdeep'], 
                            'stage': ['stage'], #, 'st_sim'], 
                            'flow': ['gflow'],#, 'fl_sim'], 
                            'c14': ['c14'], #, 'c14shal', 'c14deep'], 
                            'radon': ['radon'],#, 'rn_sim'], 
                            'strec': ['gstrec']}#, 'ec_sim']}
    
    obs_groups_types_abbrev_map = {'head': 'h', 
                                    'stage': 's', 
                                    'flow': 'f', 
                                    'c14': 'c', 
                                    'radon': 'r', 
                                    'strec': 'e'}
    
    obs_combos_df = all_combos(la, obs_groups_types, obs_groups_types_map, obs_groups_types_abbrev_map) 
    
    shared_costs = {'h':['c'],
                    'sh':['c'],
                    'dh':['c'],
                    'f':['s', 'e', 'r'],
                    's':['f', 'e', 'r'],
                    'e':['f', 's', 'r'],
                    'r':['f', 's', 'e'],
                    'c':['h']}

    def calculate_costs(obs_to_cost):
        '''
        1. Sort observations by type:
        2. Loop over types and calculate total cost
        3. 
        '''
        standard_costs_per_day = {'personnel':1500,
                                  'vehicle':250,
                                  'accom_meals':250}
        measurements_per_day = {'SW':10,
                                'GW':3}      
        maintenance_stream_gauge_per_day = 2000 / (365/2)                                

    obs_costs = {# Measurement costs
                 'base':0,
                 'h':100,
                 'sh':100,
                 'dh':100,
                 'f':100,
                 's':40,
                 'e':10,
                 'r':170,
                 'c':577,
                 # Infrastructure costs
                 'h_inf':20000,
                 'sh_inf':5000,
                 'dh_inf':50000,
                 'f_inf':100,
                 's_inf':100,
                 'e_inf':3500,
                 'r_inf':200, # If using external lab to analyse
                 'c_inf':100  # If using external lab to analyse
                 }
    obs_costs_unc = {'base':0,
                     'h':60,
                     'sh':60,
                     'dh':60,
                     'f':80,
                     's':30,
                     'e':5,
                     'r':50,
                     'c':200
                     }
    
    obs_map =   {'h':'Head',
                 'f':'Flow',
                 's':'Stage',
                 'e':'EC',
                 'r':'$^{222}$Rn',
                 'c':'$^{14}$C'
                 }
                 
    combo_df = obs_combos_df
    combo_df_percred = perc_red_unc(combo_df, 'add')
    combo_df = combo_df.apply(np.sqrt)
    combo_df = combo_df[combo_df.index != 'base']
    combo_df.loc[:, 'cost'] = combo_df.apply(get_cost, axis=1)
    combo_df.loc[:, 'cost_unc'] = combo_df.apply(get_cost_unc, axis=1)

    combo_df_percred.loc[:, 'cost'] = combo_df.apply(get_cost, axis=1)
    combo_df_percred.loc[:, 'cost_unc'] = combo_df.apply(get_cost_unc, axis=1)
    
    combo_df_lite = combo_df
    combo_cost_unc_plot(combo_df_lite, 'nrf_a')

    combo_df_lite_percred = combo_df_percred
    combo_cost_unc_plot(combo_df_lite_percred, 'nrf_a')
    
    # Let's look at clustering the combos by uncertainty and cost using
    # unsupervised machine learning via K-means clustering
    clusters = 12
    combo_df_pred, combo_df_grouped = cluster_combos_kmeans(combo_df, obs_map, clusters=clusters, normalise=False)
    scatter_unc_cost(combo_df_pred, combo_df_grouped, clusters=clusters, method='sklearn', 
                     title="", xax_title="Annual whole of river SW-GW exchange uncertainty [m$^3$/d]")
#    clusters = 20
#    combo_df_pred_percred, combo_df_grouped_percred = cluster_combos_kmeans(combo_df_percred, obs_map, clusters=clusters, normalise=True)
#    scatter_unc_cost(combo_df_pred_percred, combo_df_grouped_percred, 
#                     clusters=clusters, method='sklearn',
#                     title="", append_text='perc', xlim=(90.0, 65.0), 
#                     ylim=(0,500),
#                     xax_title="SW-GW exchange uncertainty reduction [%]")
    
    #combo_df_pred2, combo_df_grouped2 = cluster_combos_dbscan(combo_df, eps=5000, min_samples=2)
    #scatter_unc_cost(combo_df_pred2, combo_df_grouped2, method='DBSCAN')

    #_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
    # ANALYSE ALL COMBOS IN THE CONTEXT OF COST
    
    # Let's look at ALL combos for DATA additions

    obs_groups_types_p = ['head', 'stage', 
                        'flow', 'c14', 'radon', 
                        'strec'] 
    
    obs_groups_types_map_p = {'head': ['head1', 'head3', 'head5', 'head6', 'shshal', 'shdeep'], 
                            'stage': ['stage', 'st_sim'], 
                            'flow': ['gflow', 'fl_sim'], 
                            'c14': ['c14','c14shal', 'c14deep'], 
                            'radon': ['radon', 'rn_sim'], 
                            'strec': ['gstrec', 'ec_sim']}
    
    obs_groups_types_abbrev_map_p = {'head': 'h', 
                                    'stage': 's', 
                                    'flow': 'f', 
                                    'c14': 'c', 
                                    'radon': 'r', 
                                    'strec': 'e'}
    
    obs_combos_df_p = all_combos(la, obs_groups_types_p, obs_groups_types_map_p, obs_groups_types_abbrev_map_p) 
        
    obs_map_p = {'h':'Head',
                 'f':'Flow',
                 's':'Stage',
                 'e':'EC',
                 'r':'$^{222}$Rn',
                 'c':'$^{14}$C'
                 }
                 
    combo_df_p = obs_combos_df_p
    combo_df_p = combo_df_p[combo_df_p.index != 'base']
    combo_df_p = combo_df_p.apply(np.sqrt)
    combo_df_p.loc[:, 'cost'] = combo_df_p.apply(get_cost, axis=1)
    combo_df_p.loc[:, 'cost_unc'] = combo_df_p.apply(get_cost_unc, axis=1)
    
    combo_df_lite_p = combo_df_p
    combo_cost_unc_plot(combo_df_lite_p, 'nrf_a',
                        fname='Bar_unc_vs_cost_{}_potential',
                        ax_text='h=head \ns=stage \nf=flow \ne=ec \nr=radon \nc=$^{14}$C')
    
    # Let's look at clustering the combos by uncertainty and cost using
    # unsupervised machine learning via K-means clustering
    clusters = 10
    combo_df_pred_p, combo_df_grouped_p = cluster_combos_kmeans(combo_df_p, obs_map_p, clusters=clusters, normalise=False)
    scatter_unc_cost(combo_df_pred_p, combo_df_grouped_p, clusters=clusters, method='sklearn',
                     title="Data worth (potential data) for SW-GW exchange in the Campaspe River",
                     append_text="_potential", xlim=(2500, 24000))

    #combo_df_pred2, combo_df_grouped2 = cluster_combos_dbscan(combo_df, eps=5000, min_samples=2)
    #scatter_unc_cost(combo_df_pred2, combo_df_grouped2, method='DBSCAN')    

    #_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
    # ANALYSE ALL COMBOS IN THE CONTEXT OF COST
    
    # Let's look at ALL combos for DATA additions

    obs_groups_types_h = ['shallow head', 'deep head', 'stage', 
                        'flow', 'c14', 'radon', 
                        'strec'] 
    
    obs_groups_types_map_h = {'shallow head': ['head1', 'head3'], 
                            'deep head':  ['head5', 'head6'], 
                            'stage': ['stage'], 
                            'flow': ['gflow'], 
                            'c14': ['c14'], 
                            'radon': ['radon'], 
                            'strec': ['gstrec']}
    
    obs_groups_types_abbrev_map_h = {'shallow head': 'sh', 
                                    'deep head': 'dh',
                                    'stage': 's', 
                                    'flow': 'f', 
                                    'c14': 'c', 
                                    'radon': 'r', 
                                    'strec': 'e'}
    
    obs_combos_df_h = all_combos(la, obs_groups_types_h, obs_groups_types_map_h, obs_groups_types_abbrev_map_h) 
    
    obs_map_h = {'sh':'Shallow Head',
                 'dh':'Deep Head',
                 'f':'Flow',
                 's':'Stage',
                 'e':'EC',
                 'r':'$^{222}$Rn',
                 'c':'$^{14}$C'
                 }
                 
    combo_df_h = obs_combos_df_h
    combo_df_h = combo_df_h[combo_df_h.index != 'base']
    combo_df_h = combo_df_h.apply(np.sqrt)
    combo_df_h.loc[:, 'cost'] = combo_df_h.apply(get_cost, axis=1)
    combo_df_h.loc[:, 'cost_unc'] = combo_df_h.apply(get_cost_unc, axis=1)
    
    combo_df_lite_h = combo_df_h
    combo_cost_unc_plot(combo_df_lite_h, 'nrf_a',
                        fname='Bar_unc_vs_cost_{}_head_split',
                        ax_text='sh=shallow head \ndh=deep head \ns=stage \nf=flow \ne=ec \nr=radon \nc=$^{14}$C')
    
    # Let's look at clustering the combos by uncertainty and cost using
    # unsupervised machine learning via K-means clustering
    clusters = 14
    combo_df_pred_h, combo_df_grouped_h = cluster_combos_kmeans(combo_df_h, obs_map_h, clusters=clusters, normalise=True)
    scatter_unc_cost(combo_df_pred_h, combo_df_grouped_h, clusters=clusters, method='sklearn',
                     title="Data worth (existing data) for SW-GW exchange in the Campaspe River",
                     append_text="_head_split", xlim=(5000, 25000))



    # next most important obs

    #existing_obs_groups = ['head1', 'head3', 'head5', 'head6', 'stage', 'fdepth', 'gflow', 'fflow', 'c14', 'radon', 'gstrec', 'fstrec']
    #next_most = la.next_most_important_added_obs(
    #    forecast=pst_df[pst_df['obgnme'] == 'nrf_a']['obsnme'].tolist()[0], 
    #    base_obslist=pst_df[pst_df['obgnme'].isin(existing_obs_groups)]['obsnme'].tolist())
