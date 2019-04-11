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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
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
    zone2D_info = zone_array2layers(mesh3D[1])
    sfr_info, sfr_df = get_stream_info_and_sim(m, modflow_model, data_folder)
    get_heads_for_sfr(sfr_df, modflow_model, data_folder, dry_nan=True)
    sfr_df, sfr_df_reaches_relevant = process_sfr_df_and_sfr_info_into_reaches(sfr_df, sfr_info)

    # Lake Eppalock Inflow:
    Eppalock_inflow = pd.read_csv("Eppalock_inflow.csv", skiprows=2, index_col='Date', parse_dates=True, dayfirst=True)
    qual_codes_to_ignore = [8, 9, 21, 100, 101, 120, 149, 150, 151, 152,
                                153, 154, 155, 156, 160, 161, 165, 180, 190,
                                200, 201, 237, 250, 254, 255]

    Eppalock_inflow = Eppalock_inflow[~Eppalock_inflow['Qual'].isin(qual_codes_to_ignore)]
    flow_type_colours = {'low':'orangered', 'high':'dodgerblue', 'regular':'mediumseagreen'}                                
    flow_high_quantile = 0.8
    flow_low_quantile = 0.35

    Eppalock_inflow_monthly = process_Eppalock_inflow(Eppalock_inflow,
                                                      qual_codes_to_ignore=[8, 9, 21, 100, 101, 120, 149, 150, 151, 152,
                                                                            153, 154, 155, 156, 160, 161, 165, 180, 190,
                                                                            200, 201, 237, 250, 254, 255],
                                                      flow_high_quantile=flow_high_quantile,
                                                      flow_low_quantile=flow_low_quantile)                                    

    date_range = pd.date_range(start=datetime.datetime(2016,6,1), end=datetime.datetime(2017,5,30), freq='D')
    date_range_month = pd.date_range(start=datetime.datetime(2016,6,1), end=datetime.datetime(2017,5,31), freq='M')
    
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
        'rn_sim': 1./50.
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
    #la2 = load_pest_file_and_jacobian_into_Schur(model_folder, res_file='pest.rei.11')
    forecasts_df = forecast_names(la.pst.observation_data)
    forecasts_df.loc[:, 'var'] = [la.posterior_forecast[row] for row in forecasts_df.index]
    forecasts_df.loc[:, 'unc'] = forecasts_df['var'].apply(np.sqrt)

    # Some usefule definitions of what is old, new and potential data    
    old_data = ['head1', 'head3', 'head5', 'head6', 'stage', 'gflow', 'gstrec']
    new_field_data = ['fstrec', 'radon', 'c14'] #ignore ['fflow, 'fdepth'] 
    potential_data = ['st_sim', 'fl_sim', 'rn_sim', 'ec_sim', 'shshal', 'shdeep', 'c14shal', 'c14deep']
    potential_data_stream = ['st_sim', 'fl_sim', 'rn_sim', 'ec_sim']
    potential_data_ground = ['shshal', 'shdeep', 'c14shal', 'c14deep']

    if os.path.exists(os.path.join(save_folder, 'multi_df_add_subtract.pkl')) and not force_recalc:
        df_worth, df_perc, df_worth_added, df_perc_add = load_obj(os.path.join(save_folder, 'multi_df_add_subtract.pkl'))
    else:
        df_worth, df_perc, df_worth_added, df_perc_add = compare_adding_and_subtraction_of_obs_groups(la)
        save_obj([df_worth, df_perc, df_worth_added, df_perc_add], os.path.join(save_folder, 'multi_df_add_subtract.pkl'))    
    #plot_add_sub_comp(df_worth, df_perc, df_worth_added, df_perc_add)

    #_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
    # ANALYSE OBS TYPE ALONE IN TRADITIONAL VS NOVEL DATA TYPES
    
    obs_groups_name_al = [
                       "h", 
                       "s",
                       "f", 
                       "c14", 
                       "radon", 
                       "ec",
                       "all", 
                       "less_h",
                       "less_s",
                       "less_f",
                       "less_c14",
                       "less_radon",
                       "less_ec"
                       ]
    
    obs_groups_al = [
                  ['head1', 'head3', 'head5', 'head6'], 
                  ['stage'], #, 'fdepth'],
                  ['gflow'], #, 'fflow'], 
                  ['c14'], 
                  ['radon'], 
                  ['gstrec'], #, 'fstrec'], 
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'gflow', 'c14', 'radon', 'gstrec'],
                  ['stage', 'gflow', 'c14', 'radon', 'gstrec'],
                  ['head1', 'head3', 'head5', 'head6', 'gflow', 'c14', 'radon', 'gstrec'],
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'c14', 'radon', 'gstrec'],
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'gflow', 'radon', 'gstrec'],
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'gflow', 'c14', 'gstrec'],
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'gflow', 'c14', 'radon']] # 'fdepth', 'fflow', , 'fstrec'

    combo_obs_dict_al = create_combos(la, obs_groups_name_al, obs_groups_al)

    df_worth_add_combo_al, df_worth_add_combo_unc_al, df_unc_perc_add_combo_al = process_combos(la, combo_obs_dict_al, obs_groups_name_al)

    la.pst.observation_data[la.pst.observation_data['weight'] == 0]['obgnme'].unique()
    
    df_unc_perc_add_combo_al = df_unc_perc_add_combo_al.reindex([
                       "h", 
                       "less_h",
                       "s",
                       "less_s",
                       "f", 
                       "less_f",
                       "c14", 
                       "less_c14",
                       "radon", 
                       "less_radon",
                       "ec",
                       "less_ec",
                       "all", 
                       ])

    use_cols = [x for x in df_unc_perc_add_combo_al.columns if "s" not in x]
    
    df_unc_perc_add_combo_al = df_unc_perc_add_combo_al[use_cols]
    
    fig = plt.figure(figsize=(5.2, 6))
    ax = fig.add_subplot(211)

    box_vals = []
    for i in ["h", "s", "f", "c14", "radon", "ec"]: 
        box_vals += [df_unc_perc_add_combo_al.loc[i, :].tolist()]
    arxe = ax.boxplot(box_vals, positions=list(np.array(range(1, 7, 1))), widths=0.4, patch_artist=True)
    custom_boxplot(arxe, boxcol='silver', boxfillcol='silver', lw=1, 
                   whiskcol='silver', capcol='silver', mediancol='black')

#    box_vals = []
#    for i in ["less_h", "less_s", "less_f", "less_c14", "less_radon", "less_ec"]:
#        box_vals += [df_unc_perc_add_combo_al.loc[i, :].tolist()]
#    arxe = ax.boxplot(box_vals, positions=list(np.array(range(1, 7, 1)) + 0.2), widths=0.4, patch_artist=True)
#    custom_boxplot(arxe, boxcol='teal', boxfillcol='teal', lw=1, 
#                   whiskcol='teal', capcol='teal', mediancol='black')

    ax.text(0.5, 95, "a) Addition of observation group", fontsize=10, color='black', verticalalignment='center')
    ax.text(-1.25, -7, "With only:", fontsize=10, color='black', fontweight='bold', verticalalignment='center')

    ax.text(6.5, 89, "All observations", fontsize=10, color='r', verticalalignment='center', fontweight='bold')

    # Plot markers for annual exchange reduncs
    ax.plot(range(1, 7, 1), df_unc_perc_add_combo_al.loc[["h", "s", "f", "c14", "radon", "ec"]]['nrf_a'].tolist(), 
            linewidth=0, marker='o', color='teal')
    
    all_lst = np.array(df_unc_perc_add_combo_al.loc['all', :].tolist())
    ax.axhline(all_lst.max(), color='r', linestyle='--', alpha=0.6)
    ax.text(6.5, all_lst.max(), "Maximum", fontsize=9, color='r', verticalalignment='center')
    ax.axhline(np.median(all_lst), color='r', linestyle='-', alpha=0.6)
    ax.text(6.5, np.median(all_lst), "Median", fontsize=9, color='r', verticalalignment='center')
    IQR = np.percentile(all_lst, 75) - np.percentile(all_lst, 25)
    ax.axhline(np.percentile(all_lst, 75), color='r', linestyle=':', alpha=0.6)
    ax.text(6.5, np.percentile(all_lst, 75), "Upper quartile", fontsize=9, color='r', verticalalignment='center')
    ax.axhline(np.percentile(all_lst, 25), color='r', linestyle=':', alpha=0.6)
    ax.text(6.5, np.percentile(all_lst, 25), "Lower quartile", fontsize=9, color='r', verticalalignment='center')
    ax.axhline(np.min(all_lst), color='r', linestyle='--', alpha=0.6)
    ax.text(6.5, np.min(all_lst), "Minimum", fontsize=9, color='r', verticalalignment='center')

    ax.tick_params(direction='out', top=False, right=False)
    ax.set_ylim(0, 100) 
    ax.set_xlim(0.4, 6.4) 
    ax.set_xticks(range(1,7,1))
    ax.set_yticklabels([0, 20, 40, 60, 80, 100], fontsize=10)
    ax.set_xticklabels(["Head", "Stage", "Flow", "$^{14}$C", "$^{222}$Rn", "EC"], fontsize=10)
    ax.set_ylabel("% change in predictive uncertainty\n for SW-GW exchange", fontsize=10)
    fig.subplots_adjust(bottom=0.08, left=0.14, right=0.78, top=0.92, hspace=0.45)
    
    ax = fig.add_subplot(212)

    box_vals = []
    for i in ["less_h", "less_s", "less_f", "less_c14", "less_radon", "less_ec"]: 
        box_vals += [df_unc_perc_add_combo_al.loc["all", :] - df_unc_perc_add_combo_al.loc[i, :]]
    arxe = ax.boxplot(box_vals, positions=list(np.array(range(1, 7, 1))), widths=0.4, patch_artist=True)
    custom_boxplot(arxe, boxcol='silver', boxfillcol='silver', lw=1, 
                   whiskcol='silver', capcol='silver', mediancol='black')

    ax.plot(range(1, 7, 1), [df_unc_perc_add_combo_al.loc["all", 'nrf_a'] - x for x in df_unc_perc_add_combo_al.loc[["less_h", "less_s", "less_f", "less_c14", "less_radon", "less_ec"]]['nrf_a'].tolist()], 
            linewidth=0, marker='o', color='teal')
    
    
#    box_vals = []
#    for i in ["less_h", "less_s", "less_f", "less_c14", "less_radon", "less_ec"]:
#        box_vals += [df_unc_perc_add_combo_al.loc[i, :].tolist()]
#    arxe = ax.boxplot(box_vals, positions=list(np.array(range(1, 7, 1)) + 0.2), widths=0.4, patch_artist=True)
#    custom_boxplot(arxe, boxcol='teal', boxfillcol='teal', lw=1, 
#                   whiskcol='teal', capcol='teal', mediancol='black')

#    ax.text(8.5, 50, "All observations", fontsize=10, color='r', verticalalignment='center', rotation=-90)

#    all_lst = np.array(df_unc_perc_add_combo_al.loc['all', :].tolist())
#    ax.axhline(all_lst.max(), color='r', linestyle='--', alpha=0.6)
#    ax.text(6.5, all_lst.max(), "Maximum", fontsize=9, color='r', verticalalignment='center')
#    ax.axhline(np.median(all_lst), color='r', linestyle='-', alpha=0.6)
#    ax.text(6.5, np.median(all_lst), "Median", fontsize=9, color='r', verticalalignment='center')
#    IQR = np.percentile(all_lst, 75) - np.percentile(all_lst, 25)
#    ax.axhline(np.percentile(all_lst, 75), color='r', linestyle=':', alpha=0.6)
#    ax.text(6.5, np.percentile(all_lst, 75), "Upper quartile", fontsize=9, color='r', verticalalignment='center')
#    ax.axhline(np.percentile(all_lst, 25), color='r', linestyle=':', alpha=0.6)
#    ax.text(6.5, np.percentile(all_lst, 25), "Lower quartile", fontsize=9, color='r', verticalalignment='center')
#    ax.axhline(np.min(all_lst), color='r', linestyle='--', alpha=0.6)
#    ax.text(6.5, np.min(all_lst), "Minimum", fontsize=9, color='r', verticalalignment='center')

    ax.text(0.5, 18, "b) Removal of each observation group\n from all groups combined", fontsize=10, color='black', verticalalignment='center')
    ax.text(-1.10, -1.4, "Without:", fontsize=10, color='black', fontweight='bold', verticalalignment='center')

    ax.tick_params(direction='out', top=False, right=False)
    ax.set_ylim(0, 20) 
    ax.set_xlim(0.4, 6.4) 
    ax.set_xticks(range(1,7,1))
    ax.yaxis.set_tick_params(labelsize=10)
    #ax.set_yticklabels([0, 20, 40, 60, 80, 100], fontsize=10)
    ax.set_xticklabels(["Head", "Stage", "Flow", "$^{14}$C", "$^{222}$Rn", "EC"], fontsize=10)
    ax.set_ylabel("% change in predictive uncertainty\n for SW-GW exchange", fontsize=10)

    nrf_a_patch = Line2D([0], [0], marker='o', linestyle='none',  color='teal', label='Whole of river\nannual exchange')
    handles = [nrf_a_patch]
    plt.legend(handles=handles, ncol=1, fontsize=9, numpoints=1,
               frameon=False, bbox_to_anchor=(0.95, 1.7), loc=2, 
               handletextpad=0.1, title='')

    
    fig.subplots_adjust(bottom=0.08, left=0.16, right=0.74, top=0.97, hspace=0.18)    
    
    plt.savefig(os.path.join(save_folder, "With_or_without_data_type.png"), dpi=300)
    
    ob_interest_al = 'nrf_a'

    #_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
    # ANALYSE SIMPLE COMBOS IN TRADITIONAL VS NOVEL DATA TYPES
    
    obs_groups_name = [
                       "hydraulic", 
                       "h_s_f_c14", 
                       "h_s_f_radon", 
                       "h_s_f_stream_ec",
                       "all" 
                       ]
    
    obs_groups = [
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'fdepth', 'gflow', 'fflow'], 
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'fdepth', 'gflow', 'fflow', 'c14'], 
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'fdepth', 'gflow', 'fflow', 'radon'], 
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'fdepth', 'gflow', 'fflow', 'gstrec', 'fstrec'], 
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'fdepth', 'gflow', 'fflow', 'c14', 'radon', 'gstrec', 'fstrec']]

    combo_obs_dict = create_combos(la, obs_groups_name, obs_groups)

    df_worth_add_combo, df_worth_add_combo_unc, df_unc_perc_add_combo = process_combos(la, combo_obs_dict, obs_groups_name)

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    # Combined plot of individual and cumulative 
    def plot_boxplot_points(df, ax, annual_color, pred_categories, 
                            pred_flow_types, pred_colours, pred_exchange,
                            marker_dict_predtype, markerfill_dict_predtype, alpha,
                            gain_lose_indicate=False,
                            flow_type_indicate=True):
        max_exchange = pred_exchange.max()
        min_exchange = pred_exchange.max()
        mean_exchange = pred_exchange.mean()
        neutral_criterion = 0.1 * mean_exchange
        pred_exchange[pred_exchange['Exchange'] > 0]

        for cat in pred_categories:
            for num, ind in enumerate(df.index):
                relevant_columns = [col for col in df.columns if cat in col]
                y_cat_ind = df[relevant_columns].loc[ind, :]
                y_cat_ind_ = pred_exchange.loc[relevant_columns, :]
                y_cat_ind2 = y_cat_ind_[y_cat_ind_['Exchange'] > 0.004]
                y_cat_ind3 = y_cat_ind_[y_cat_ind_['Exchange'] < -0.004]
                y_cat_ind4 = y_cat_ind_[(y_cat_ind_['Exchange'] <= 0.004) & (y_cat_ind_['Exchange'] >= -0.004)]
                relevant_columns2 = y_cat_ind2.index
                relevant_columns3 = y_cat_ind3.index
                relevant_columns4 = y_cat_ind4.index
                y_cat_ind2 = df[relevant_columns2].loc[ind, :]
                y_cat_ind3 = df[relevant_columns3].loc[ind, :]
                y_cat_ind4 = df[relevant_columns4].loc[ind, :]
                x_cat_ind = [num + 1 + tick_offsets[cat]] * len(y_cat_ind)
                x_cat_ind2 = [num + 1 + tick_offsets[cat]] * len(y_cat_ind2)
                x_cat_ind3 = [num + 1 + tick_offsets[cat]] * len(y_cat_ind3)
                x_cat_ind4 = [num + 1 + tick_offsets[cat]] * len(y_cat_ind4)
                if 'a' in cat:
                    markerfacecolor = 'none' #annual_color
                    edgecolor = annual_color
                else:
                    markerfacecolor =  [pred_colours[int(val.split('_')[-1])] for val in relevant_columns]
                    edgecolor = 'none' #[pred_colours[int(val.split('_')[-1])] for val in relevant_columns]
                    extra_adjust = [pred_flow_types[int(val.split('_')[-1])] for val in relevant_columns]
                    extra_adjust2 = []
                    extra_adjust3 = [pred_flow_types[int(val.split('_')[-1])] for val in relevant_columns2]
                    extra_adjust4 = []
                    extra_adjust5 = [pred_flow_types[int(val.split('_')[-1])] for val in relevant_columns3]
                    extra_adjust6 = []
                    extra_adjust7 = [pred_flow_types[int(val.split('_')[-1])] for val in relevant_columns4]
                    extra_adjust8 = []
                    for val in extra_adjust:
                        if val == 'low':
                            extra_adjust2 += [0]
                        elif val == 'regular':
                            extra_adjust2 += [0.05]
                        elif val == 'high':
                            extra_adjust2 += [0.1]
                    for val in extra_adjust3:
                        if val == 'low':
                            extra_adjust4 += [0]
                        elif val == 'regular':
                            extra_adjust4 += [0.05]
                        elif val == 'high':
                            extra_adjust4 += [0.1]
                    for val in extra_adjust5:
                        if val == 'low':
                            extra_adjust6 += [0]
                        elif val == 'regular':
                            extra_adjust6 += [0.05]
                        elif val == 'high':
                            extra_adjust6 += [0.1]
                    for val in extra_adjust7:
                        if val == 'low':
                            extra_adjust8 += [0]
                        elif val == 'regular':
                            extra_adjust8 += [0.05]
                        elif val == 'high':
                            extra_adjust8 += [0.1]
                    x_cat_ind = np.array(x_cat_ind) + np.array(extra_adjust2)
                    x_cat_ind2 = np.array(x_cat_ind2) + np.array(extra_adjust4)
                    x_cat_ind3 = np.array(x_cat_ind3) + np.array(extra_adjust6)
                    x_cat_ind4 = np.array(x_cat_ind4) + np.array(extra_adjust8)
                # end if
                if flow_type_indicate:
                    ax.scatter(x=x_cat_ind, y=y_cat_ind, marker=marker_dict_predtype[cat], 
                            s=12, linewidths=1, c=markerfacecolor, alpha=alpha, 
                            edgecolors=edgecolor, zorder=2)   
                if gain_lose_indicate:
                    ax.scatter(x=x_cat_ind2, y=y_cat_ind2, marker='+', 
                            s=10, linewidths=0.7, c='orange', alpha=alpha, 
                            zorder=2)   
                    ax.scatter(x=x_cat_ind3, y=y_cat_ind3, marker='x', 
                            s=10, linewidths=0.7, c='green', alpha=alpha, 
                            zorder=2)   
                    ax.scatter(x=x_cat_ind4, y=y_cat_ind4, marker='*', 
                            s=10, linewidths=0.7, c='blue', alpha=alpha, 
                            zorder=2)   

    def combined_alone_and_combo_boxplots_for_select_preds(df_select_al, df_select, identifier, 
                                                       pred_categories,
                                                       marker_dict_predtype,
                                                       markerfill_dict_predtype,
                                                       tick_offsets,
                                                       pred_flow_types,
                                                       pred_colours,
                                                       colour_dict,
                                                       pred_exchange,
                                                       alpha=0.7,
                                                       titles=['a', 'b'],
                                                       annual_color='black',
                                                       flow_type_indicate=True,
                                                       gain_lose_indicate=False):

        #fig = plt.figure(figsize=(4.92,5))
        fig = plt.figure(figsize=(7,7))
        ax = fig.add_subplot(2, 1, 1)
        ticklabels_al = ("Head", "Stage", "Flow", "$^{14}$C", "$^{222}$Rn", 
                         "EC", "All data")
        #arxe = df_select_al.transpose().plot(kind='box', ax=ax, patch_artist=True)#bar_chart(df_unc_perc_add_combo_al[ob_interest], ticklabels_al, ax=ax, save=False)
        arxe = ax.boxplot(df_select_al.transpose().as_matrix(), patch_artist=True)

        ax.text(0.2, 1.05, 'Hydraulic', horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes, fontsize=12)
        ax.text(0.65, 1.05, 'Chemical', horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes, fontsize=12)
        #ax.text(0.93, 1.05, 'Combined', horizontalalignment='center',
        #        verticalalignment='center', transform=ax.transAxes, fontsize=12)
        ax.fill_between([3.5, 6.5], [0, 0], [100, 100], color='lightgrey', alpha=0.15) #, interpolate=True)

        width = 0.5
        ax.plot([6.5, 6.5], [0, 100], color='grey', linestyle='--')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.tick_params(direction='out', top=False, right=False)
        ax.set_ylabel('% reduction in predictive\n uncertainty for net annual exchange flux')
        ax.yaxis.set_label_coords(-0.13,-0.1)
        ax.set_xticklabels(ticklabels_al, rotation=45)    
        ax.set_ylim(0, 100)
        ax.text(0.05, 0.9, 'a.', horizontalalignment='center',
                     verticalalignment='center', transform=ax.transAxes, fontsize=12)    

        custom_boxplot(arxe, boxcol='lightgrey', boxfillcol='lightgrey', lw=1, 
                   whiskcol='lightgrey', capcol='lightgrey', mediancol='black')

        plot_boxplot_points(df_select_al, ax, annual_color, pred_categories, 
                            pred_flow_types, pred_colours, pred_exchange,
                            marker_dict_predtype, markerfill_dict_predtype, alpha,
                            flow_type_indicate=flow_type_indicate,
                            gain_lose_indicate=gain_lose_indicate)

        if flow_type_indicate:
            nrf_a_patch = Line2D([0], [0], marker=marker_dict_predtype['nrf_a'], linestyle='none', markerfacecolor='none', markeredgewidth=1, color='black', label='Whole annual')
            nrf_m_patch = Line2D([0], [0], marker=marker_dict_predtype['nrf_m'], linestyle='none',  color='black', label='Whole monthly')
            rrf_a_patch = Line2D([0], [0], marker=marker_dict_predtype['rrf_a'], linestyle='none', markerfacecolor='none', markeredgewidth=1, color='black', label='Reach annual')
            rrf_m_patch = Line2D([0], [0], marker=marker_dict_predtype['rrf_m'], linestyle='none',  color='black', label='Reach monthly')
            handles = [nrf_a_patch, nrf_m_patch, rrf_a_patch, rrf_m_patch]
            plt.legend(handles=handles, ncol=1, fontsize=10, numpoints=1,
                       frameon=False, bbox_to_anchor=(0.97, 1.05), loc=2, 
                       handletextpad=0.1, title='Prediction type')
        if gain_lose_indicate:
            gain_patch = Line2D([0], [0], marker='x', linestyle='none', color='black', label='Gaining')
            lose_patch = Line2D([0], [0], marker='+', linestyle='none', color='black', label='Neutral')
            neutral_patch = Line2D([0], [0], marker='*', linestyle='none',  color='black', label='Losing')
            handles = [gain_patch, neutral_patch, lose_patch]
            plt.legend(handles=handles, ncol=1, fontsize=10, numpoints=1,
                       frameon=False, bbox_to_anchor=(0.97, 1.05), loc=2, 
                       handletextpad=0.1, title='Exchange type')
        
        ax = fig.add_subplot(2, 1, 2)
        ticklabels_c = ("Head", "+ Stage", "+ Flow", "Hydraulic \n+ $^{14}$C",
                        "Hydraulic \n+ $^{222}$Rn", "Hydraulic \n+ EC", "All data")
        arxe2 = ax.boxplot(df_select.transpose().as_matrix(), patch_artist=True)
        ax.fill_between([3.5, 6.5], [0, 0], [100, 100], color='lightgrey', alpha=0.15) #, interpolate=True)

        ax.plot([6.5, 6.5], [0, 100], color='grey', linestyle='--')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.tick_params(direction='out', top=False, right=False)
        #ax.ticks['right'].set_visible(False)
        #ax.ticks['top'].set_visible(False)
        ax.set_xticklabels(ticklabels_c, rotation=45)    
        ax.set_ylim(0, 100)
        ax.text(0.05, 0.9, 'b.', horizontalalignment='center',
                     verticalalignment='center', transform=ax.transAxes, fontsize=12)    

        custom_boxplot(arxe2, boxcol='lightgrey', boxfillcol='lightgrey', lw=1, 
                   whiskcol='lightgrey', capcol='lightgrey', mediancol='black')

        plot_boxplot_points(df_select, ax, annual_color, pred_categories, 
                            pred_flow_types, pred_colours, pred_exchange,
                            marker_dict_predtype, markerfill_dict_predtype, alpha,
                            flow_type_indicate=flow_type_indicate,
                            gain_lose_indicate=gain_lose_indicate)
        
        if flow_type_indicate:
            silver_patch = mpatches.Rectangle((0,0), 0.1, 0.1, color=annual_color, label='mean', alpha=alpha)
            red_patch = mpatches.Patch(color=flow_type_colours['low'], label='low', alpha=alpha)
            orange_patch = mpatches.Patch(color=flow_type_colours['regular'], label='regular', alpha=alpha)
            blue_patch = mpatches.Patch(color=flow_type_colours['high'], label='high', alpha=alpha)
            handles = [silver_patch, red_patch, orange_patch, blue_patch]
            plt.legend(handles=handles, ncol=1, fontsize=10, 
                       frameon=False, bbox_to_anchor=(0.99, 1.9), loc=2, title='Flow type')                

        fig.subplots_adjust(bottom=0.13, left=0.15, right=0.81, top=0.96, hspace=0.29)
        plt.savefig(p_j(save_folder, 
            'Predictive uncert individual vs cumulative w hydraulic vs chem_{}.png'.format(identifier)), dpi=300)

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
        
                            
    df_noseason_al = df_unc_perc_add_combo_al[[col for col in df_unc_perc_add_combo_al.columns if 's' not in col]]
    df_noseason = df_unc_perc_add_combo[[col for col in df_unc_perc_add_combo.columns if 's' not in col]]
    combined_alone_and_combo_boxplots_for_select_preds(df_noseason_al, df_noseason, 'ALL',
                                                       pred_categories,
                                                       marker_dict_predtype,
                                                       markerfill_dict_predtype,
                                                       tick_offsets,
                                                       pred_flow_types,
                                                       pred_colours,
                                                       flow_type_colours,
                                                       pred_exchange,
                                                       alpha=0.7,
                                                       gain_lose_indicate=False,
                                                       flow_type_indicate=True)

    #_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
    # Uncertatinty for all the forecast observations 
    df = la.obscov.to_dataframe()
    obscov = pd.Series(np.diag(df), index=[df.index, df.columns])

    sfr_df_reaches_relevant.loc[:, 'unc'] = sfr_df_reaches_relevant.apply(
        lambda x: forecasts_df[forecasts_df['unique_name'] == "rrf_m{}_{}".format(
        int(x['reach_dw']), int(x['time'] - 21))]['unc'].tolist()[0], axis=1)
    
    sfr_df_reaches_relevant.loc[:, 'unc_lin'] = sfr_df_reaches_relevant['unc'] / sfr_df_reaches_relevant['rchlen']
    sfr_df_reaches_relevant.loc[:, 'Cumulative Length (km)'] = sfr_df_reaches_relevant['Cumulative Length'] / 1000.

    #_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
    # ANALYSE SAME COMBOS ACROSS VARIOUS SPATIOTEMPORAL RANGES
    columns_name = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 'nrf']

    rows_name_orig = ['Annual', 'June', 'July', 'August',
                 'September', 'October', 'November', 'December', 'January', 'February', 'March', 'April', 'May', 
                     'Quarter1', 'Quarter2', 'Quarter3', 'Quarter4', ]

    rows_name = ['Annual', 'Quarter1', 'Quarter2', 'Quarter3', 'Quarter4', 'June', 'July', 'August',
                 'September', 'October', 'November', 'December', 'January', 'February', 'March', 'April', 'May']

    rows_name2 = ['Annual', 'Jun, 16', 'Jul, 16', 'Aug, 16',
                 'Sep, 16', 'Oct, 16', 'Nov, 16', 'Dec, 16', 'Jan, 17', 'Feb, 17', 'Mar, 17', 'Apr, 17', 'May, 17']

    columns = {}

    novel_data_sets_old = ['h', 'h_s'    , 'h_s_f'      , 'h_s_f_c14', 'h_s_f_radon', 'h_s_f_stream_ec', 'all']
    novel_data_sets_alias_old = ['Head', '+ Stage', '+ Flow', 'Hydraulic + $^{14}$C' , 'Hydraulic + $^{222}$Rn', 'Hydraulic + EC',  'All data']
    novel_data_sets = ['hydraulic', 'h_s_f_c14', 'h_s_f_radon', 'h_s_f_stream_ec', 'all']
    novel_data_sets_alias = ['Hydraulic', 'Hydraulic + $^{14}$C' , 'Hydraulic + $^{222}$Rn', 'Hydraulic + EC',  'All data']

    novel_data_sets_al = ['all', 'h', 's', 'f', 'c14', 'radon', 'ec']
    novel_data_sets_alias_al = ['All', 'Head', 'Stage' , 'Flow', '$^{14}$C',  '$^{222}$Rn', 'EC']

    novel_data_set_to_alias = {key:val for key, val in zip(novel_data_sets, novel_data_sets_alias)}
    novel_data_set_to_alias_al = {key:val for key, val in zip(novel_data_sets_al, novel_data_sets_alias_al)}

    month_and_annual = [x for x in df_worth_add_combo_unc.columns if 's' not in x]

    df_worth_add_combo_unc_filter = df_worth_add_combo_unc[month_and_annual] 
    df_worth_add_combo_unc_filter = df_worth_add_combo_unc_filter[df_worth_add_combo_unc_filter.index != 'base']
    df_worth_add_combo_unc_filter.index = [novel_data_set_to_alias[x] for x in df_worth_add_combo_unc_filter.index]
    df_perc_add_combo_unc_filter = df_unc_perc_add_combo[month_and_annual] 
    df_perc_add_combo_unc_filter.index = [novel_data_set_to_alias[x] for x in df_perc_add_combo_unc_filter.index]

    df_worth_add_combo_unc_filter_al = df_worth_add_combo_unc_al[month_and_annual] 
    df_worth_add_combo_unc_filter_al = df_worth_add_combo_unc_filter_al[df_worth_add_combo_unc_filter_al.index != 'base']
    df_worth_add_combo_unc_filter_al = df_worth_add_combo_unc_filter_al.loc[novel_data_sets_al, :]
    df_worth_add_combo_unc_filter_al.index = [novel_data_set_to_alias_al[x] for x in df_worth_add_combo_unc_filter_al.index if x in novel_data_set_to_alias_al.keys()]
    df_perc_add_combo_unc_filter_al = df_unc_perc_add_combo_al[month_and_annual] 
    df_perc_add_combo_unc_filter_al = df_unc_perc_add_combo_al.loc[novel_data_sets_al, :] 
    df_perc_add_combo_unc_filter_al.index = [novel_data_set_to_alias_al[x] for x in df_perc_add_combo_unc_filter_al.index]


    months = rows_name2[1:]
    reaches = 11

    ###########################################################################
    #
    # 
    #
    ###########################################################################
    
    #plot_spatiotemporal_individual(novel_data_sets_alias, df_worth_add_combo_unc_filter, 
    #                              adjust=365. * 1000. / 1000000000. , cbar_text='[Gl/yr]', unc_type='',
    #                              vlim=(0., 10.))
        
    #plot_spatiotemporal_individual(novel_data_sets_alias, df_perc_add_combo_unc_filter, 
    #                              adjust=1.0, cbar_text='[%]', unc_type='perc_red',
    #                              vlim=(0., 100.))

    plot_spatiotemporal_unc_top_plus6(novel_data_sets_alias, 
                                      df_perc_add_combo_unc_filter, 
                                      months, reaches, 
                                      m,
                                      zone2D_info,
                                      save_folder,
                                      Eppalock_inflow, 
                                      date_range, 
                                      rows_name2,
                                      unc_type='',
                                      flow_low_quantile=flow_low_quantile,
                                      flow_high_quantile=flow_high_quantile,
                                      flow_type_colours=flow_type_colours
                                      )

    plot_spatiotemporal_unc_top_plus6(novel_data_sets_alias_al, 
                                      df_perc_add_combo_unc_filter_al, 
                                      months, reaches, 
                                      m,
                                      zone2D_info,
                                      save_folder,
                                      Eppalock_inflow, 
                                      date_range, 
                                      rows_name2,
                                      unc_type='individual',
                                      flow_low_quantile=flow_low_quantile,
                                      flow_high_quantile=flow_high_quantile,
                                      flow_type_colours=flow_type_colours
                                      )

    #_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
    #_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
    # Potential Observations

    obs_groups_name2 = []
    for ob_gp in obs_groups_name_al:
        obs_groups_name2 += [ob_gp, ob_gp + '_new']

    obs_groups2 = [
                  ['head1', 'head3', 'head5', 'head6'], 
                  ['head1', 'head3', 'head5', 'head6', 'shshal', 'shdeep'], 
                  ['stage', 'fdepth'],
                  ['stage', 'fdepth', 'st_sim'],
                  ['gflow', 'fflow'], 
                  ['gflow', 'fflow', 'fl_sim'], 
                  ['c14'], 
                  ['c14', 'c14shal', 'c14deep'], 
                  ['radon'], 
                  ['radon', 'rn_sim'], 
                  ['gstrec', 'fstrec'], 
                  ['gstrec', 'fstrec', 'ec_sim'], 
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'fdepth', 'gflow', 'fflow', 'c14', 'radon', 'gstrec', 'fstrec'],
                  ['head1', 'head3', 'head5', 'head6', 'shshal', 'shdeep', 'stage', 'fdepth', 'st_sim', 'gflow', 'fflow', 'fl_sim', 'c14', 'c14shal', 'c14deep', 'radon', 'rn_sim', 'gstrec', 'fstrec', 'ec_sim']]

    combo_obs_dict2 = create_combos(la, obs_groups_name2, obs_groups2)
    df_worth_add_combo2, df_worth_add_combo_unc2, df_unc_perc_add_combo2 = process_combos(la, combo_obs_dict2, obs_groups_name2)
    
    pst_df = la.pst.observation_data
 
    df_ob = df_unc_perc_add_combo2['nrf_a']        
    existing = df_ob[~df_ob.index.str.contains('new')].tolist() 
    new = df_ob[df_ob.index.str.contains('new')].tolist()

    #
    ###
    ######
    #########
    ######
    ###
    #

    #fig = plt.figure(figsize=(4.92,3))
    fig = plt.figure(figsize=(4.92*0.8,3*0.8))

    df_ob = df_unc_perc_add_combo2#['nrf_a']        
    existing = df_ob[~df_ob.index.str.contains('new')]#.tolist() 
    new = df_ob[df_ob.index.str.contains('new')]#.tolist()
    
    ind = np.arange(len(existing))  # the x locations for the groups
    width = 0.35  # the width of the bars

    ax = fig.add_subplot(111)
    ax.text(0.2, 1.05, 'Hydraulic', horizontalalignment='center',
         verticalalignment='center', transform=ax.transAxes, fontsize=10)
    ax.text(0.64, 1.05, 'Chemical', horizontalalignment='center',
         verticalalignment='center', transform=ax.transAxes, fontsize=10)
    #ax.text(0.93, 1.05, 'Combined', horizontalalignment='center',
    #     verticalalignment='center', transform=ax.transAxes, fontsize=10)
    ax.fill_between([2.5, 5.45], [0, 0], [100, 100], color='lightgrey', alpha=0.6) #, interpolate=True)
    ax.plot([5.45, 5.45], [0, 100], color='grey', linestyle='--')
    #arxe3 = new.transpose().plot(kind='box', ax=ax, positions=ind + width/2, widths=width, patch_artist=True)
    #arxe2 = existing.transpose().plot(kind='box', ax=ax, positions=ind - width/2, widths=width, patch_artist=True)

    arxe2 = ax.boxplot(existing.transpose().as_matrix(), positions=ind - width/2, widths=width, patch_artist=True)
    arxe3 = ax.boxplot(new.transpose().as_matrix(), positions=ind + width/2, widths=width, patch_artist=True)

    custom_boxplot(arxe2)
    custom_boxplot(arxe3, boxfillcol='teal')
            
    #rects1 = ax.bar(ind - width/2, existing, width, 
    #                color='IndianRed', label='Existing', edgecolor='none')
    #rects2 = ax.bar(ind + width/2, new, width, 
    #                color='SkyBlue', label='Potential', edgecolor='none')
    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('% reduction in predictive\n uncertainty', fontsize=10)
    ax.set_ylim(0, 100)        
    ax.set_xlim(0 - width*1.3, len(existing) - width)        
    #ax.set_title('Observation types reduction in uncertainty')
    #ax.set_xticks([x + width/2 for x in ind])
    comp_ticklabels = ("Head", "Stage", "Flow", "$^{14}$C",
                        "$^{222}$Rn", "EC", "All data")
    ax.set_xticks(range(0, 7))
    ax.tick_params(axis='y', labelsize=10)
    ax.set_xticklabels(comp_ticklabels, rotation=45, size=10)
    ax.tick_params(direction='out', top=False, right=False)

    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    #ax.legend(loc='top center', fontsize=12)  
    #ax.text(0.03, 0.37, 'Existing', transform=ax.transAxes, fontsize=10, rotation=90)#, color='white')
    #ax.text(0.084, 0.40, 'Potential', transform=ax.transAxes, fontsize=10, rotation=90)

    fig.subplots_adjust(bottom=0.28, left=0.18, right=0.95, top=0.92, hspace=0.45)
    plt.savefig(p_j(save_folder, 'Comp_existing_potential_combos_all.png'), dpi=300)    

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    

    novel_data_set_to_alias2 = {'all_new': 'All',
                                 'h_new': 'Head',
                                 's_new': 'Stage',
                                 'f_new': 'Flow',
                                 'c14_new': '$^{14}$C',
                                 'radon_new': '$^{222}$Rn',
                                 'ec_new': 'EC'}

    novel_data_set_to_alias3 = {'all_diff': 'All',
                                 'h_diff': 'Head',
                                 's_diff': 'Stage',
                                 'f_diff': 'Flow',
                                 'c14_diff': '$^{14}$C',
                                 'radon_diff': '$^{222}$Rn',
                                 'ec_diff': 'EC'}

    #df_worth_add_combo_unc_filter2 = df_unc_perc_add_combo2[month_and_annual] 
    #df_worth_add_combo_unc_filter2 = df_worth_add_combo_unc_filter2[df_worth_add_combo_unc_filter2.index != 'base']
    #df_worth_add_combo_unc_filter2 = df_worth_add_combo_unc_filter2.loc[[col for col in df_worth_add_combo_unc_filter2.index if 'new' in col], :]    
    #df_worth_add_combo_unc_filter2.index = [novel_data_set_to_alias2[x] for x in df_worth_add_combo_unc_filter2.index]
    df_perc_add_combo_unc_filter2 = df_unc_perc_add_combo2[month_and_annual] 
    df_perc_add_combo_unc_filter2 = df_perc_add_combo_unc_filter2.loc[[col for col in df_perc_add_combo_unc_filter2.index if 'new' in col], :]    
    df_perc_add_combo_unc_filter2 = df_perc_add_combo_unc_filter2[~df_perc_add_combo_unc_filter2.index.str.contains('less')]    
    df_perc_add_combo_unc_filter2.index = [novel_data_set_to_alias2[x] for x in df_perc_add_combo_unc_filter2[~df_perc_add_combo_unc_filter2.index.str.contains('less')].index]

    plot_spatiotemporal_unc_top_plus6(novel_data_sets_alias_al, 
                                      df_perc_add_combo_unc_filter2, 
                                      months, reaches, 
                                      m,
                                      zone2D_info,
                                      save_folder,
                                      Eppalock_inflow, 
                                      date_range, 
                                      rows_name2,
                                      unc_type='pot',
                                      flow_low_quantile=flow_low_quantile,
                                      flow_high_quantile=flow_high_quantile,
                                      flow_type_colours=flow_type_colours
                                      )

    df_perc_add_combo_unc_filter3 = df_unc_perc_add_combo2[month_and_annual] 
    for i in range(0, len(df_perc_add_combo_unc_filter3.index), 2):
        base = df_perc_add_combo_unc_filter3.index[i]
        df_perc_add_combo_unc_filter3.loc['{}_diff'.format(base), :] = df_perc_add_combo_unc_filter3.loc['{}_new'.format(base), :] - df_perc_add_combo_unc_filter3.loc['{}'.format(base), :] 
    df_perc_add_combo_unc_filter3 = df_perc_add_combo_unc_filter3.loc[[col for col in df_perc_add_combo_unc_filter3.index if 'diff' in col], :]    
    df_perc_add_combo_unc_filter3 = df_perc_add_combo_unc_filter3[~df_perc_add_combo_unc_filter3.index.str.contains('less')]    
    df_perc_add_combo_unc_filter3.index = [novel_data_set_to_alias3[x] for x in df_perc_add_combo_unc_filter3.index]

    plot_spatiotemporal_unc_top_plus6(novel_data_sets_alias_al, 
                                      df_perc_add_combo_unc_filter3, 
                                      months, reaches, 
                                      m,
                                      zone2D_info,
                                      save_folder,
                                      Eppalock_inflow, 
                                      date_range, 
                                      rows_name2,
                                      unc_type='pot_diff',
                                      flow_low_quantile=flow_low_quantile,
                                      flow_high_quantile=flow_high_quantile,
                                      flow_type_colours=flow_type_colours,
                                      vlim=(0, 40),
                                      title_prefix='Improvement in reduction from existing to potential: ',
                                      cbar_label='[% difference]'
                                      )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    

    novel_data_set_to_alias2_al = {'all_new': 'All data',
                                 'h_new': 'Head',
                                 's_new': 'Stage',
                                 'f_new': 'Flow',
                                 'c14_new': '$^{14}$C',
                                 'radon_new': '$^{222}$Rn',
                                 'ec_new': 'EC'}

    novel_data_set_to_alias3_al = {'all_diff': 'All data',
                                 'h_diff': 'Head',
                                 'h_s_diff': '+ Stage',
                                 'h_s_f_diff': '+ Flow',
                                 'h_s_f_c14_diff': 'Hydraulic + $^{14}$C',
                                 'h_s_f_radon_diff': 'Hydraulic + $^{222}$Rn',
                                 'h_s_f_stream_ec_diff': 'Hydraulic + EC'}

    #df_worth_add_combo_unc_filter2 = df_unc_perc_add_combo2[month_and_annual] 
    #df_worth_add_combo_unc_filter2 = df_worth_add_combo_unc_filter2[df_worth_add_combo_unc_filter2.index != 'base']
    #df_worth_add_combo_unc_filter2 = df_worth_add_combo_unc_filter2.loc[[col for col in df_worth_add_combo_unc_filter2.index if 'new' in col], :]    
    #df_worth_add_combo_unc_filter2.index = [novel_data_set_to_alias2[x] for x in df_worth_add_combo_unc_filter2.index]
    df_perc_add_combo_unc_filter2_al = df_unc_perc_add_combo2_al[month_and_annual] 
    df_perc_add_combo_unc_filter2_al = df_perc_add_combo_unc_filter2_al.loc[[col for col in df_perc_add_combo_unc_filter2_al.index if 'new' in col], :]    
    df_perc_add_combo_unc_filter2_al.index = [novel_data_set_to_alias2_al[x] for x in df_perc_add_combo_unc_filter2_al.index]

    plot_spatiotemporal_unc_top_plus6(novel_data_sets_alias_al, 
                                      df_perc_add_combo_unc_filter2_al, 
                                      months, reaches, 
                                      m,
                                      zone2D_info,
                                      save_folder,
                                      Eppalock_inflow, 
                                      date_range, 
                                      rows_name2,
                                      unc_type='pot_al',
                                      flow_low_quantile=flow_low_quantile,
                                      flow_high_quantile=flow_high_quantile,
                                      flow_type_colours=flow_type_colours
                                      )

    df_perc_add_combo_unc_filter3 = df_unc_perc_add_combo2[month_and_annual] 
    for i in range(0, len(df_perc_add_combo_unc_filter3.index), 2):
        base = df_perc_add_combo_unc_filter3.index[i]
        df_perc_add_combo_unc_filter3.loc['{}_diff'.format(base), :] = df_perc_add_combo_unc_filter3.loc['{}_new'.format(base), :] - df_perc_add_combo_unc_filter3.loc['{}'.format(base), :] 
    df_perc_add_combo_unc_filter3 = df_perc_add_combo_unc_filter3.loc[[col for col in df_perc_add_combo_unc_filter3.index if 'diff' in col], :]    
    df_perc_add_combo_unc_filter3.index = [novel_data_set_to_alias3[x] for x in df_perc_add_combo_unc_filter3.index]

    plot_spatiotemporal_unc_top_plus6(novel_data_sets_alias, 
                                      df_perc_add_combo_unc_filter3, 
                                      months, reaches, 
                                      m,
                                      zone2D_info,
                                      save_folder,
                                      Eppalock_inflow, 
                                      date_range, 
                                      rows_name2,
                                      unc_type='pot_diff',
                                      flow_low_quantile=flow_low_quantile,
                                      flow_high_quantile=flow_high_quantile,
                                      flow_type_colours=flow_type_colours,
                                      vlim=(0, 26),
                                      title_prefix='Improvement in reduction from existing to potential: ',
                                      cbar_label='[% difference]'
                                      )
    
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Look at things spatiotemporally    

    observations = load_obj('observations.pkl')

    model_mesh3D = m.model_mesh3D
    model_boundary = m.model_boundary
    nrow = model_mesh3D[0].shape[1]
    ncol = model_mesh3D[0].shape[2]
    delr = m.gridHeight
    delc = m.gridWidth
    xul = model_boundary[0]
    yul = model_boundary[3]

    print("Estimate for all stream potential obs: {:.2f} hrs".format(
        pst_df[pst_df['obgnme'].isin(
            potential_data_stream)].shape[0] * 3.5 / 3600.))

    if os.path.exists(os.path.join(save_folder, 'All_potential_data_stream.csv')) and not force_recalc:
        stream_pot_obs = pd.read_csv(os.path.join(save_folder, 'All_potential_data_stream.csv'), index_col=0)
    else:
        stream_pot_obs = la.get_added_obs_importance(
            {ob:[ob] for ob in pst_df[pst_df['obgnme'].isin(potential_data_stream)].index.tolist()})
        stream_pot_obs.to_csv(os.path.join(save_folder, 'All_potential_data_stream.csv'))

    existing_groups = old_data + new_field_data    
    existing_obs = pst_df[pst_df['obgnme'].isin(existing_groups)].index.tolist()    
    
    if os.path.exists(os.path.join(save_folder, 'All_potential_data_stream_add.csv')) and not force_recalc:
        stream_pot_obs_add = pd.read_csv(os.path.join(save_folder, 'All_potential_data_stream_add.csv'), index_col=0)
    else:
        stream_pot_obs_add = la.get_added_obs_importance(obslist_dict= \
            {ob:[ob] for ob in pst_df[pst_df['obgnme'].isin(potential_data_stream)].index.tolist()}, 
            base_obslist=existing_obs)
        stream_pot_obs_add.to_csv(os.path.join(save_folder, 'All_potential_data_stream_add.csv'))
        
        
    stream_pot_obs_unc = stream_pot_obs[~stream_pot_obs.index.isin(['base'])].apply(
        lambda x:(1 - x / stream_pot_obs.loc["base", :]) * 100., axis=1)

    stream_pot_obs_unc_add = stream_pot_obs_add[~stream_pot_obs_add.index.isin(['base'])].apply(
        lambda x:(1 - x / stream_pot_obs_add.loc["base", :]) * 100., axis=1)
    
    river_seg = pd.read_csv('river_top.csv')
    riv_pp_df = pd.read_csv('river_pp.csv')

    #for stream_ob_type in potential_data_stream:
    #    mega_plot_river_predunc_red(stream_pot_obs_unc, pst_df, stream_ob_type, 'ob3454', v_bounds=(0, 25))

    mf = flopy.modflow.Modflow.load(os.path.join(observations_folder, '02_transient_flow.nam'), 
                                    version='mfnwt')        

    sfr_seg_rec = mf.get_package('sfr').segment_data
    sfr_seg = pd.DataFrame.from_records(sfr_seg_rec[0])
    sfr_orig = pd.DataFrame.from_records(mf.get_package('sfr').reach_data)
    pd.DataFrame.cumsum(sfr_orig)['rchlen'].tolist()
    sfr_orig.loc[:, 'Cumulative Length'] = pd.DataFrame.cumsum(sfr_orig)['rchlen']
    sfr_orig.loc[:, 'Cumulative Length km'] = pd.DataFrame.cumsum(sfr_orig)['rchlen'] / 1000.0

    sfr_df.loc[:, 'x'] = [xul + col * delc for col in sfr_df['column']]
    sfr_df.loc[:, 'y'] = [yul - row * delr for row in sfr_df['row']]
    sfr_df.loc[:, 'Qaquifer_norm'] = \
    sfr_df[sfr_df['time'] == 32].reset_index(
        range(sfr_df.shape[0]))['Qaquifer'] / (sfr_orig['rchlen'] * sfr_seg['width2'])

    stream_ob_types = ['stage', 'gflow', 'radon', 'gstrec']
    stream_ob_dic = {}
    for stream_ob_type in stream_ob_types:    
        existing_flow_observations = observations.obs_group[stream_ob_type]
        eoo = existing_flow_observations
        eoo_ts = eoo['time_series']
        eoo_ts = eoo_ts[eoo_ts['active'] == True]
        eoo_ts['name'].unique()
        eoo_loc = eoo['locations']
        try:
            eoo_loc_used = eoo_loc[eoo_loc.index.isin([int(x) for x in eoo_ts['name'].unique()])]    
        except:
            eoo_loc_used = eoo_loc[eoo_loc.index.isin(eoo_ts['name'].unique())]    
            eoo_loc_used = eoo_loc_used.sort_values('seg_loc')
        eoo_loc_used.loc[:, 'Cumulative Length km'] = \
            sfr_orig[sfr_orig['iseg'].isin(eoo_loc_used['seg_loc'])]['Cumulative Length km'].tolist()
        stream_ob_dic[stream_ob_type] = eoo_loc_used

    obs_colour_map={'head':'#1b9e77',
                    'stage':'#d95f02',
                    'gflow':'#7570b3',
                    'radon':'#e7298a',
                    'c14':'#e6ab02',
                    'fstrec':'#66a61e',
                    'gstrec':'#66a61e'}

    obs_facecolour_map={'head':'#1b9e77',
                    'stage':'none',
                    'gflow':'none',
                    'radon':'#e7298a',
                    'c14':'#e6ab02',
                    'fstrec':'#66a61e',
                    'gstrec':'#66a61e'}

    obs_marker_map={'head':'^',
                    'stage':'o',
                    'gflow':'o',
                    'radon':'o',
                    'c14':'*',
                    'fstrec':'+',
                    'gstrec':'+'}

    obs_size_map={'head':8,
                  'stage':40,
                  'gflow':140,
                  'radon':25,
                  'c14':8,
                  'fstrec':85,
                  'gstrec':85}

    obs_lw_map={'head':2,
                'stage':2,
                'gflow':1.5,
                'radon':2,
                'c14':2,
                'fstrec':2,
                'gstrec':2}

    fig = plt.figure(figsize=(7, 1.2))
    ax = fig.add_subplot(111) 
    for index, key in enumerate(stream_ob_types):
        ax.scatter(x=stream_ob_dic[key]['Cumulative Length km'], y=[index] * stream_ob_dic[key].shape[0],
                   marker=obs_marker_map[key], facecolors=obs_facecolour_map[key], edgecolors=obs_colour_map[key],
                   s=obs_size_map[key] * 0.5, lw=obs_lw_map[key])  
    ax.set_ylim(-0.5, 3.5)
    ax.set_xlim(0., 141)
    ax.set_ylabel('Existing data')
    ax.set_yticks([k for k in range(4)])
    ax.set_yticklabels(['Stage', 'Flow', '$^{222}$Rn', 'EC'], minor=False, rotation=0)
    ax.set_xticklabels('')
    fig.subplots_adjust(right=0.8)

    plt.savefig(p_j(save_folder, 'Existing_stream_data_locs_linear.png'), dpi=300)

#    for month in [1]:#range(12):
#        mega_plot_river_predunc_red_month(month, stream_pot_obs_unc, pst_df, 
#                                          potential_data_stream, 'ob3454', 
#                                          v_bounds=(0,25))

    #select_months = [4, 5, 11] 
    select_months = [1, 4, 5] 
    stream_pot_obs_unc_add.columns = [forecasts_df.loc[x, 'unique_name'] for x in stream_pot_obs_unc_add.columns]    
    mega_plot_river_predunc_red_month_axes(select_months, stream_pot_obs_unc_add, pst_df, 
                                           potential_data_stream, ['nrf_m_1', 'nrf_m_4', 'nrf_m_5'], sfr_df, sfr_orig,
                                           riv_pp_df, save_folder, v_bounds=(0, 1),
                                           color_bar_exchange=[0.85, 0.663, 0.03, 0.163])    

    if os.path.exists(os.path.join(save_folder, 'All_potential_data_ground.csv')) and not force_recalc:
        ground_pot_obs = pd.read_csv(os.path.join(save_folder, 'All_potential_data_ground.csv'), index_col=0)
    else:
        ground_pot_obs = la.get_added_obs_importance(
            {ob:[ob] for ob in pst_df[pst_df['obgnme'].isin(potential_data_ground)].index.tolist()})
        ground_pot_obs.to_csv(os.path.join(save_folder, 'All_potential_data_ground.csv'))
    # end if

    if os.path.exists(os.path.join(save_folder, 'All_potential_data_ground_add.csv')) and not force_recalc:
        ground_pot_obs_add = pd.read_csv(os.path.join(save_folder, 'All_potential_data_ground_add.csv'), index_col=0)
    else:
        ground_pot_obs_add = la.get_added_obs_importance(obslist_dict= \
            {ob:[ob] for ob in pst_df[pst_df['obgnme'].isin(potential_data_ground)].index.tolist()}, 
            base_obslist=existing_obs)
        ground_pot_obs_add.to_csv(os.path.join(save_folder, 'All_potential_data_ground_add.csv'))

    ground_pot_obs_unc = ground_pot_obs[~ground_pot_obs.index.isin(['base'])].apply(
        lambda x:(1 - x / ground_pot_obs.loc["base", :]) * 100., axis=1)

    ground_pot_obs_unc_add = ground_pot_obs_add[~ground_pot_obs_add.index.isin(['base'])].apply(
        lambda x:(1 - x / ground_pot_obs_add.loc["base", :]) * 100., axis=1)
    
    existing_head_observations = observations.obs_group['head']
    eho = existing_head_observations
    eho_ts = eho['time_series']
    eho_loc = eho['locations']
    eho_ts_shal_bores = eho_ts[eho_ts['zone'].isin(['head1', 'head3'])]['name'].unique()
    eho_ts_deep_bores = eho_ts[eho_ts['zone'].isin(['head5', 'head6'])]['name'].unique()
    eho_shal_loc = eho_loc[eho_loc.index.isin(eho_ts_shal_bores)]
    eho_deep_loc = eho_loc[eho_loc.index.isin(eho_ts_deep_bores)]

    observations.obs_group.keys()

    existing_c14_observations = observations.obs_group['C14']
    eco = existing_c14_observations
    eco_ts = eco['time_series']
    eco_loc = eco['locations']
    eco_ts_shal_bores = eco_ts[eco_ts['zone'].isin(['C141', 'C143'])]['name'].unique()
    eco_ts_deep_bores = eco_ts[eco_ts['zone'].isin(['C145', 'C146'])]['name'].unique()
    eco_shal_loc = eco_loc[eco_loc.index.isin(eco_ts_shal_bores)]
    eco_deep_loc = eco_loc[eco_loc.index.isin(eco_ts_deep_bores)]


    shallow_xyz_dict = observations.obs_group['shshal']['locations']
    deep_xyz_dict = observations.obs_group['shdeep']['locations']

    shallow_obs_ts = observations.obs_group['shshal']['time_series']
    deep_obs_ts = observations.obs_group['shdeep']['time_series']

    shallow_obs_ts['x'] = [shallow_xyz_dict[ob][2] for ob in shallow_obs_ts['name']]
    shallow_obs_ts['y'] = [shallow_xyz_dict[ob][1] for ob in shallow_obs_ts['name']]
    shallow_obs_ts['z'] = [shallow_xyz_dict[ob][0] for ob in shallow_obs_ts['name']]
    deep_obs_ts['x'] = [deep_xyz_dict[ob][2] for ob in deep_obs_ts['name']]
    deep_obs_ts['y'] = [deep_xyz_dict[ob][1] for ob in deep_obs_ts['name']]
    deep_obs_ts['z'] = [deep_xyz_dict[ob][0] for ob in deep_obs_ts['name']]

    shallow_obs_ts['x_m'] = [xul + shallow_xyz_dict[ob][2] * delc for ob in shallow_obs_ts['name']]
    shallow_obs_ts['y_m'] = [yul - shallow_xyz_dict[ob][1] * delr for ob in shallow_obs_ts['name']]
    deep_obs_ts['x_m'] = [xul + deep_xyz_dict[ob][2] * delc for ob in deep_obs_ts['name']]
    deep_obs_ts['y_m'] = [yul - deep_xyz_dict[ob][1] * delr for ob in deep_obs_ts['name']]

    df = ground_pot_obs_unc['ob3454']
    shallow_obs_ts['perc_red_unc'] = [df.loc[ob] for ob in shallow_obs_ts['obs_map']]
    deep_obs_ts['perc_red_unc'] = [df.loc[ob] for ob in deep_obs_ts['obs_map']]

    ground_pot_obs_unc.describe().transpose().sort_values(by='max', ascending=False).head(10)['max']
    ground_pot_obs_unc.describe().transpose().sort_values(by='min').head(100)['min']
    ground_pot_obs_unc.describe().transpose()['mean']

    hgus = [[0, 2], [4, 5]]
    hgus_new = [[1, 3], [5, 6]]
    def plot_ground_pot_obs_predunc_axes(fig, obtype, dfs, hgus, vlim=(0, 60)):
        plot_titles = ["{} {}".format(x, obtype) for x in ['Shallow', 'Deep']]
        for index in range(2):
            hgu = hgus[index]
            ax = fig.add_subplot(1, 3, index + 2, aspect='equal')
            ax.set_title("Potential\n {}".format(plot_titles[index]), fontsize=10)
            #ax.axis('off')
            ax.set_ylabel("")
            plot_stream_reaches_basic(m, ax, zone2D_info, zones=hgu)            
            ax2 = ax.scatter(dfs[index]['x_m'], dfs[index]['y_m'], c=dfs[index]['perc_red_unc'], 
                       cmap='viridis', linewidth=0.0, vmin=vlim[0], vmax=vlim[1])
            bounds = m.model_boundary
            ax.set_xlim(bounds[0], bounds[1])
            ax.set_ylim(bounds[2], bounds[3])    
            ax.yaxis.set_ticklabels("")
            ax.set_ylabel("")
            plt.xticks(fontsize=10)
            plt.yticks(fontsize=10)
            x_offset_text = ax.xaxis.get_offset_text()
            x_offset_text.set_size(10)

            if index == 1:
                fig.subplots_adjust(right=0.8)
                cbar_ax = fig.add_axes([0.85, 0.30, 0.03, 0.6])
                fig.colorbar(ax2, cax=cbar_ax)    
                cbar_ax.get_yaxis().labelpad = 18
                cbar_ax.set_ylabel('Uncertainty reduction (%)', fontsize=10, rotation=270)
                cbar_ax.tick_params(labelsize=10)
        return

    for date_of_interest in shallow_obs_ts['datetime'].unique()[0:1]:#[0]
        dfs = [shallow_obs_ts[shallow_obs_ts['datetime'] == date_of_interest],
               deep_obs_ts[deep_obs_ts['datetime'] == date_of_interest]]

        fig = plt.figure(figsize=(6,3.6))
        ax = fig.add_subplot(1, 3, 1, aspect='equal')
        plot_stream_reaches_basic(m, ax, zone2D_info, bounds=m.model_boundary)
        x_offset_text = ax.xaxis.get_offset_text()
        x_offset_text.set_size(10)
        y_offset_text = ax.yaxis.get_offset_text()
        y_offset_text.set_size(10)
        ax.scatter(eho_shal_loc['Easting'], eho_shal_loc['Northing'], marker='^', edgecolor='', c=obs_colour_map['head'], label='shallow') 
        ax.scatter(eho_deep_loc['Easting'], eho_deep_loc['Northing'], marker='^', edgecolor='', c=obs_colour_map['stage'], label='deep') 
        ax.legend(scatterpoints=1, fontsize=8, borderpad=0.2)
        ax.set_title("Existing\n Head", fontsize=10)

        plot_ground_pot_obs_predunc_axes(fig, 'Head', dfs, hgus, vlim=[0, 60])
        fig.subplots_adjust(bottom=0.01, left=0.1, right=0.83, top=1.00, hspace=0.40, wspace=0.1)
        plt.savefig(p_j(save_folder, 'Potential_head_shallow_and_deep.png'), dpi=300)

#    for date_of_interest in shallow_obs_ts['datetime'].unique():#[0]
#        dfs = [shallow_obs_ts[shallow_obs_ts['datetime'] == date_of_interest],
#               deep_obs_ts[deep_obs_ts['datetime'] == date_of_interest]]
#        plot_ground_pot_obs_predunc(zone2D_info, 'Head', dfs, hgus, vlim=[0, 50])

    #
    # C14 potential observations:
    #

    hgus = [[0, 2], [4, 5]]
    #date_of_interest = '2017-05-31'

    c14_shallow_xyz_dict = observations.obs_group['c14shal']['locations']
    c14_deep_xyz_dict = observations.obs_group['c14deep']['locations']

    c14_shallow_obs_ts = observations.obs_group['c14shal']['time_series']
    c14_deep_obs_ts = observations.obs_group['c14deep']['time_series']        

    c14_shallow_obs_ts['perc_red_unc'] = [df.loc[ob] for ob in c14_shallow_obs_ts['obs_map']]
    c14_deep_obs_ts['perc_red_unc'] = [df.loc[ob] for ob in c14_deep_obs_ts['obs_map']]

    c14_shallow_obs_ts['x'] = [c14_shallow_xyz_dict[ob][2] for ob in c14_shallow_obs_ts['name']]
    c14_shallow_obs_ts['y'] = [c14_shallow_xyz_dict[ob][1] for ob in c14_shallow_obs_ts['name']]
    c14_deep_obs_ts['x'] = [c14_deep_xyz_dict[ob][2] for ob in c14_deep_obs_ts['name']]
    c14_deep_obs_ts['y'] = [c14_deep_xyz_dict[ob][1] for ob in c14_deep_obs_ts['name']]  
    c14_shallow_obs_ts['x_m'] = [xul + c14_shallow_xyz_dict[ob][2] * delc for ob in c14_shallow_obs_ts['name']]
    c14_shallow_obs_ts['y_m'] = [yul - c14_shallow_xyz_dict[ob][1] * delr for ob in c14_shallow_obs_ts['name']]
    c14_deep_obs_ts['x_m'] = [xul + c14_deep_xyz_dict[ob][2] * delc for ob in c14_deep_obs_ts['name']]
    c14_deep_obs_ts['y_m'] = [yul - c14_deep_xyz_dict[ob][1] * delr for ob in c14_deep_obs_ts['name']]  

    hgus = [[0, 2], [4, 5]]
    #date_of_interest = '2017-05-31'
    for date_of_interest in c14_shallow_obs_ts['datetime'].unique():
        dfs = [c14_shallow_obs_ts[c14_shallow_obs_ts['datetime'] == date_of_interest],
               c14_deep_obs_ts[c14_deep_obs_ts['datetime'] == date_of_interest]]

    for date_of_interest in shallow_obs_ts['datetime'].unique()[0:1]:#[0]
        #dfs = [c14_shallow_obs_ts.groupby(by='name').mean(),
        #       c14_deep_obs_ts.groupby(by='name').mean()]

        dfs = [c14_shallow_obs_ts[c14_shallow_obs_ts['datetime'] == date_of_interest],
               c14_deep_obs_ts[c14_deep_obs_ts['datetime'] == date_of_interest]]

        fig = plt.figure(figsize=(6,3.6))
        ax = fig.add_subplot(1, 3, 1, aspect='equal')
        plot_stream_reaches_basic(m, ax, zone2D_info, bounds=m.model_boundary)
        x_offset_text = ax.xaxis.get_offset_text()
        x_offset_text.set_size(10)
        y_offset_text = ax.yaxis.get_offset_text()
        y_offset_text.set_size(10)
        ax.scatter(eco_shal_loc['Easting'], eco_shal_loc['Northing'], marker='*', edgecolor='', c=obs_colour_map['c14'], s=50, label='shallow') 
        ax.scatter(eco_deep_loc['Easting'], eco_deep_loc['Northing'], marker='*', edgecolor='', c=obs_colour_map['radon'], s=50, label='deep') 
        ax.legend(scatterpoints=1, fontsize=8, borderpad=0.2)
        ax.set_title("Existing\n $^{14}$C", fontsize=10)

        plot_ground_pot_obs_predunc_axes(fig, '$^{14}$C', dfs, hgus, vlim=[0, 60])
        fig.subplots_adjust(bottom=0.01, left=0.1, right=0.83, top=1.00, hspace=0.40, wspace=0.1)
        plt.savefig(p_j(save_folder, 'Potential_c14_shallow_and_deep.png'), dpi=300)

    # Create shapefile for shallow obs
    import geopandas
    from shapely.geometry import Point
    pot_flow_observations = observations.obs_group['fl_sim']['locations']
    sfr_df_pot = sfr_df[sfr_df['time'] == sfr_df['time'].unique()[0]]   
    sfr_df_pot = sfr_df_pot[sfr_df_pot['segment'].isin(pot_flow_observations['seg_loc'])]
    sfr_df_pot.loc[:, 'geometry'] = sfr_df_pot.apply(lambda x: Point(x.x, x.y), axis=1)
    sfr_df_pot_geo = geopandas.GeoDataFrame(sfr_df_pot, geometry='geometry')                        
    sfr_df_pot_geo.crs = {'init': 'epsg:28355'}                        
    sfr_df_pot_geo.to_file('stream_pot_locs.shp', driver='ESRI Shapefile')


    shallow_obs_pot = shallow_obs_ts.groupby(by='name').mean()        
    shallow_obs_pot.loc[:, 'geometry'] = shallow_obs_pot.apply(lambda x: Point(x.x_m, x.y_m), axis=1)
    shallow_obs_pot_geo = geopandas.GeoDataFrame(shallow_obs_pot, geometry='geometry')                        
    shallow_obs_pot_geo.crs = {'init': 'epsg:28355'}                        
    shallow_obs_pot_geo.drop('active', axis=1, inplace=True)
    shallow_obs_pot_geo.to_file('gw_pot_locs.shp', driver='ESRI Shapefile')

    deep_obs_pot = deep_obs_ts.groupby(by='name').mean()        
    deep_obs_pot.loc[:, 'geometry'] = deep_obs_pot.apply(lambda x: Point(x.x_m, x.y_m), axis=1)
    deep_obs_pot_geo = geopandas.GeoDataFrame(deep_obs_pot, geometry='geometry')                        
    deep_obs_pot_geo.crs = {'init': 'epsg:28355'}                        
    deep_obs_pot_geo.drop('active', axis=1, inplace=True)
    deep_obs_pot_geo.to_file('gw_pot_locs_deep.shp', driver='ESRI Shapefile')


    # next most important obs

    obs_groups_name_pot = [
                       "new_head_shallow",
                       "new_head_deep",
                       "new_stage",
                       "new_flow", 
                       "new_c14", 
                       "new_radon", 
                       "new_ec",
                       ]
                       
    obs_groups_pot = [
                  ['shshal'],
                  ['shdeep'], 
                  ['st_sim'],
                  ['fl_sim'],
                  ['c14shal', 'c14deep'], 
                  ['rn_sim'], 
                  ['ec_sim']]

    combo_obs_dict_pot = create_combos(la, obs_groups_name_pot, obs_groups_pot)
    
    existing_obs_groups = ['head1', 'head3', 'head5', 'head6', 'stage', 'fdepth', 'gflow', 'fflow', 'c14', 'radon', 'gstrec', 'fstrec']
    
    if os.path.exists(os.path.join(save_folder, 'next_most_pot.csv')) and not force_recalc:
        next_most = pd.read_csv(os.path.join(save_folder, 'next_most_pot.csv'), index_col=0)
    else:
        print("Calculating next most important obs to add")
        next_most = la.next_most_important_added_obs(
            forecast=pst_df[pst_df['obgnme'] == 'nrf_a']['obsnme'].tolist()[0], 
            base_obslist=pst_df[pst_df['obgnme'].isin(existing_obs_groups)]['obsnme'].tolist(),
            obslist_dict=combo_obs_dict_pot, niter=7)
        next_most.to_csv(os.path.join(save_folder, 'next_most_pot.csv'))
    # end if

    potential_data_mod = [i for i in potential_data if i not in ['st_sim', 'c14shal', 'c14deep']]
    pot_all_list = all_obs(potential_data_mod, la.pst.observation_data)
    len(pot_all_list)
    pot_all_dict = {x:x for x in pot_all_list} 
    
    if os.path.exists(os.path.join(save_folder, 'next_most_pot_all.csv')) and not force_recalc:
        next_most2 = pd.read_csv(os.path.join(save_folder, 'next_most_pot_all.csv'), index_col=0)
    else:
        print("Calculating next most important obs to add")
        next_most2 = la.next_most_important_added_obs(
            forecast=pst_df[pst_df['obgnme'] == 'nrf_a']['obsnme'].tolist()[0], 
            base_obslist=pst_df[pst_df['obgnme'].isin(existing_obs_groups)]['obsnme'].tolist(),
            obslist_dict=pot_all_dict, niter=5)
        next_most2.to_csv(os.path.join(save_folder, 'next_most_pot_all.csv'))
    # end if

    sfr_orig.loc[:, 'x_m'] = [xul + ob * delc for ob in sfr_orig['j']]
    sfr_orig.loc[:, 'y_m'] = [yul - ob * delr for ob in sfr_orig['i']]
    
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    sfr_orig[['x_m', 'y_m']].plot(kind='scatter', x='x_m', y='y_m', ax=ax, )

    for ob in next_most2.index:
        data_type = la.pst.observation_data.loc[ob, 'obgnme']
        seg = int(observations.obs_group['rn_sim']['time_series'][observations.obs_group['rn_sim']['time_series']['obs_map'] == ob]['name'].tolist()[0].replace(data_type, ''))
        time = observations.obs_group['rn_sim']['time_series'][observations.obs_group['rn_sim']['time_series']['obs_map'] == ob]['interval'].tolist()[0]
        ax.scatter(sfr_orig[sfr_orig['iseg'] == seg]['x_m'].tolist()[0], sfr_orig[sfr_orig['iseg'] == seg]['y_m'].tolist()[0], c='red')         
