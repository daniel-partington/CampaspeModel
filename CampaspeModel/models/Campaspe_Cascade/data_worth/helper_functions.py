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

from flopy.utils import sfroutputfile
import pyemu

from HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface

# Configuration Loader
from HydroModelBuilder.Utilities.Config.ConfigLoader import ConfigLoader

p_j = os.path.join

def unix2dos(path):
    text = open(path, "U").read() 
    text = text.replace("\n", "\r\n") 
    open(path, "wb").write(text)

def save_obj(obj, filename):
        filename.replace('.pkl', '')
        with open(filename + '.pkl', 'wb') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)    
    
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
        pst.observation_data.loc[pst.observation_data['obgnme'] == 'radon', 'weight'] = 1./50.
        temp = pst.observation_data.copy()
        pst.adjust_weights_resfile(resfile=os.path.join(model_folder, res_file))
        pst = adjust_potential_by_existing(pst)
        pst.observation_data.loc[:, 'weight_diff'] = pst.observation_data['weight'] - temp['weight']
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

def perc_red(df, mode):
    if mode == 'add':
        df_perc = df[~df.index.isin(['base', 'nrf_m'])].apply(
           lambda x:(1 - x / df.loc["base", :]) * 100., axis=1)             
    elif mode == 'subtract':
        df_perc = df[~df.index.isin(['base', 'nrf_m'])].apply(
           lambda x:(x - df.loc["base", :]) / df.loc["base", :] * 100., axis=1)             
    else:
        print "Mode not recognised, try 'add' or 'subtract'"

    return df_perc

def perc_red_unc(df, mode):
    df_unc = df.apply(np.sqrt)
    df_perc = perc_red(df_unc, mode)             
    return df_perc

def plot_add_sub_comp(df_worth, df_perc, df_worth_added, df_perc_add, forecast_of_interest=None):
    if not forecast_of_interest:
        forecast_of_interest = la.pst.observation_data[la.pst.observation_data['obgnme'] == 'nrf_a']['obsnme']
    # end if
    df_unc_perc = perc_red_unc(df_worth, 'subtract')             
    df_unc_perc_add = perc_red_unc(df_worth_added, 'add')
    for index, df in enumerate([(df_perc, df_perc_add), (df_unc_perc, df_unc_perc_add)]):
        fig = plt.figure()
        ax = fig.add_subplot(311)
        df[0][forecast_of_interest].plot(kind='bar', color='red', alpha=0.1, ax=ax, label='subtract')
        ax = fig.add_subplot(312)
        df[1][forecast_of_interest].plot(kind='bar', alpha=0.1, ax=ax, label='add')
        ax.set_title('Comparison of addition and subtraction of data')
        if index == 0:
            label = 'posterior'
        elif index == 1:
            label = 'uncertainty'
        # end if
        ax.set_ylabel('Percent reduction in {}'.format(label))
        df_add_sub_diff = df[0][forecast_of_interest] - df[1][forecast_of_interest]
        ax = fig.add_subplot(313)
        df_add_sub_diff.plot(kind='bar', color='red', label='difference', ax=ax)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(p_j(save_folder, 'Comparison of addition and subtraction of data on {}.png'.format(label)), dpi=300)

def bar_chart(df_ob, ticklabels, colors, title=True, ylabel=True, ax=None, 
              legend=True, save=True, spines=True):
    
    if not ax:
        ax = fig.add_subplot(1, 1, 1)
    
    if title:
        ax.text(0.2, 1.05, 'Hydraulic', horizontalalignment='center',
             verticalalignment='center', transform=ax.transAxes, fontsize=12)
        ax.text(0.65, 1.05, 'Chemical', horizontalalignment='center',
             verticalalignment='center', transform=ax.transAxes, fontsize=12)
        #ax.text(0.93, 1.05, 'Combined', horizontalalignment='center',
        #     verticalalignment='center', transform=ax.transAxes, fontsize=12)
    ax.fill_between([2.25, 5.25], [0, 0], [100, 100], color='lightgrey', alpha=0.6) #, interpolate=True)

    width = 0.5
    ax.plot([5.25, 5.25], [0, 100], color='grey', linestyle='--')
    ax.bar(np.array(range(len(df_ob.index))) - width, df_ob, width, color=colors[0], edgecolor='none')
    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    if ylabel:
        ax.set_ylabel('% reduction in predictive\n uncertainty for net annual exchange flux')
        #ax.yaxis.set_label_coords(-0.1,1.02)
    ax.set_ylim(0, 100)        
    ax.set_xlim(0 - 2 * width + 0.25, len(df_ob.index) - width - 0.25)        
    #ax.set_title('Observation types reduction in uncertainty')
    ax.set_xticks([x - 0.5 * width for x in np.array(range(len(df_ob.index)))])
    ax.set_xticklabels(ticklabels)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    if spines:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    if legend:
        ax.legend(loc='top center', fontsize=12)
    if save:
        plt.savefig(p_j(save_folder, 'Predictive uncert cumulative hydraulic vs chem_alone.png'), dpi=300)
    return ax        
        
def compare_adding_and_subtraction_of_obs_groups(la, forecast_of_interest=None, 
                                                 plot_difference=False):
    df_worth = la.get_removed_obs_group_importance()
    df_worth_added = la.get_added_obs_group_importance()
    df_perc = perc_red(df_worth, 'subtract')             
    df_perc_add = perc_red(df_worth_added, 'add')
    if plot_difference:
        plot_add_sub_comp(df_worth, df_perc, df_worth_added, df_perc_add, forecast_of_interest)
    # end if
    return df_worth, df_perc, df_worth_added, df_perc_add
    
def all_obs(obsgnme, obs_exist):
    obs_list = obs_exist[obs_exist['obgnme'].isin(obsgnme)].index.tolist()
    return obs_list

def random_obs(obsgnme, num, obs_exist):
    obs_list = obs_exist[obs_exist['obgnme'].isin(obsgnme)].index.tolist()
    obs_list = [obs_list[i] for i in np.random.choice(len(obs_list), num)]
    return obs_list

def create_combos(la, obs_groups_name, obs_groups):
    # Let's look at some combos
    combo_obs_dict = {}
    obs_exist = la.pst.observation_data
    
    for ogn, og in zip(obs_groups_name, obs_groups):
        combo_obs_dict[ogn] = all_obs(og, obs_exist)
    
    return combo_obs_dict 

def box_plot_for_the_punters(mean_exchange, df, unit_convert=1.0):
    obs_groups = df.index.tolist()
    print obs_groups
    mu = mean_exchange
    sigma = {}

    #colors = ['white', 'gray', 'b', 'r', 'y', 'g', ]

    width = 8
    height = 7
    multiplier = 1.
    fig = plt.figure(figsize=(width * multiplier, height * multiplier))

    data = []
    for index, item in enumerate(obs_groups):
        sigma[item] = df.loc[item] * unit_convert
        data += [np.random.normal(mu, sigma[item], 1E6)]

    ax = fig.add_subplot(1, 1, 1)

    bp = plt.boxplot(data, 0, '', showmeans=False, patch_artist=True) #, whis='range') #, colors=colors)

    xtickNames = plt.setp(ax, xticklabels= ['Head', '+ Stage', '+ Discharge', 'also C14', 
                                            'also Radon', 'also Stream EC',  'all obs']) 

    y_lims = ax.get_ylim()
    y_lims = [y_lims[0], y_lims[1]]
    #y_lims = [y_lims[0] * 0.9, y_lims[1] * 1.1]
    plt.plot([3.5, 3.5], y_lims, color='gray', linestyle='--', linewidth=1.)

    plt.setp(xtickNames, rotation=45, fontsize=11)

    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['medians'], color='black')
    plt.setp(bp['means'], color='black')

    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)

    ax.spines["top"].set_visible(False)    
    #ax.spines["bottom"].set_visible(False)    
    ax.spines["right"].set_visible(False)    
    #ax.spines["left"].set_visible(False)         
    ax.spines["bottom"].set_color('gray')        
    ax.spines["left"].set_color('gray')        

    #ax.set_title('Worth of various data for net annual SW-GW exchange')
    ax.set_ylabel('Net Annual SW-GW exchange [Gl/yr]')

    ax = plt.gca()
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    
    plt.tight_layout()
    plt.savefig(p_j(save_folder, 'predunc_annual_swgw_{}.png'.format(df.name)), dpi=300)  

def create_spatiotemporal_unc_array(df, months, reaches):
    heat = pd.DataFrame(columns=["r{}".format(x) for x in range(1, reaches + 1)] + ['nrf'], index=['Annual'] + months)
    heat.loc['Annual', 'nrf'] = df['nrf_a']
    for i in range(reaches):
        heat.loc['Annual', "r{}".format(i + 1)] = df['rrf_a{}'.format(i)]
    for i, m in enumerate(months):
        heat.loc[m, 'nrf'] = df['nrf_m_{}'.format(i)]
    for i in range(reaches):
        for j, m in enumerate(months):
            heat.loc[m, "r{}".format(i + 1)] = df['rrf_m{}_{}'.format(i, j)]
    return heat
    
def process_combos(la, combo_obs_dict, obs_groups_name, reset_col_names_by_obgnme=True):    
    df_worth_add_combo = la.get_added_obs_importance(obslist_dict=combo_obs_dict)
    if reset_col_names_by_obgnme:
        forecasts_df = forecast_names(la.pst.observation_data)
        df_worth_add_combo.columns = [forecasts_df.loc[x, 'unique_name'] for x in df_worth_add_combo.columns]
    # end if
    df_worth_add_combo = df_worth_add_combo.reindex([
                   "base"] + obs_groups_name)
    df_worth_add_combo_unc = df_worth_add_combo.apply(np.sqrt)
    df_unc_perc_add_combo = df_worth_add_combo_unc[~df_worth_add_combo_unc.index.isin(['nrf_m'])] \
       .apply(lambda x:(1 - x / df_worth_add_combo_unc.loc["base", :]) * 100., axis=1)
    df_unc_perc_add_combo = df_unc_perc_add_combo[df_unc_perc_add_combo.index != 'base']
    return df_worth_add_combo, df_worth_add_combo_unc, df_unc_perc_add_combo

def process_combos_remove(la, combo_obs_dict, obs_groups_name, reset_col_names_by_obgnme=True):    
    df_worth_remove_combo = la.get_removed_obs_importance(obslist_dict=combo_obs_dict)
    if reset_col_names_by_obgnme:
        forecasts_df = forecast_names(la.pst.observation_data)
        df_worth_remove_combo.columns = [forecasts_df.loc[x, 'unique_name'] for x in df_worth_remove_combo.columns]
    # end if
    df_worth_remove_combo = df_worth_remove_combo.reindex([
                   "base"] + obs_groups_name)
    df_worth_remove_combo_unc = df_worth_remove_combo.apply(np.sqrt)
    df_unc_perc_remove_combo = df_worth_remove_combo_unc[~df_worth_remove_combo_unc.index.isin(['nrf_m'])] \
       .apply(lambda x:((x - df_worth_remove_combo_unc.loc["base", :]) / df_worth_remove_combo_unc.loc["base", :]) * 100., axis=1)
    df_unc_perc_remove_combo = df_unc_perc_remove_combo[df_unc_perc_remove_combo.index != 'base']
    return df_worth_remove_combo, df_worth_remove_combo_unc, df_unc_perc_remove_combo
    
def plot_combos(obs_interest):
    plt.figure()
    ax = df_unc_perc_add_combo[ob_interest].plot(kind='bar', colors=colors)
    ax.set_title('Worth of combinations of data types for "{}"'.format(ob_interest))
    ax.set_ylabel('Percent reduction in uncertainty')
    plt.savefig(p_j(save_folder, 'Worth of combinations of data types for {}.png'.format(ob_interest)), dpi=300)


def plot_spatiotemporal_individual(novel_data_sets_alias, df_core, adjust=1.0, 
                                  cbar_text='', unc_type='default', vlim=None):
    
    for dset in novel_data_sets_alias:
        df = df_core.loc[dset]
        heat = create_spatiotemporal_unc_array(df, months, reaches) 
        heat = heat.astype(np.float64)
        heat = heat * adjust
        heat = heat.transpose()
        heat = heat.reindex(heat.index.tolist()[::-1])
        plt.figure(figsize=(9,5))
        if not vlim:
            vlim = (heat.min(), heat.max())
        ax = sns.heatmap(heat, vmin=vlim[0], vmax=vlim[1], cmap='viridis')
        ax.set_yticklabels(heat.index, rotation=0)    
        plt.title("{}".format(dset))
        plt.text(14.7, 9.2, cbar_text)
        plt.tight_layout()
        plt.savefig(p_j(save_folder, 'Spatiotemporal_{}_analysis_for_{}.png'.format(unc_type, dset)), dpi=300)         
    
def plot_spatiotemporal_unc_top_plus6(novel_data_sets_alias, 
                                      df_active, 
                                      months, reaches,
                                      m,
                                      zone2D_info,
                                      save_folder,
                                      Eppalock_inflow, 
                                      date_range, 
                                      rows_name2,
                                      adjust=1.0, unc_type='', vlim=(0, 100),
                                      flow_low_quantile=0.25, flow_high_quantile=0.75,
                                      flow_type_colours=None,
                                      title_prefix="Reduction in uncertainty in SW-GW exchange for ",
                                      cbar_label='[% reduction]'):
    
    fig = plt.figure(figsize=(2*9/2., 3*5/2.))
    
    labels = {0:'b.', 1:'c.', 2:'d.', 3:'e.', 4:'f.', 5:'g.', 6:'h.'}
    
    num_groups = len(novel_data_sets_alias[1:])
    for index, dset in enumerate(novel_data_sets_alias[1:]):
        df = df_active.loc[dset]
        heat = create_spatiotemporal_unc_array(df, months, reaches) 
        heat = heat.astype(np.float64)
        heat = heat * adjust
        heat = heat.transpose()
        heat = heat.reindex(heat.index.tolist()[::-1])
        ax = fig.add_subplot(3, 2, index + 1)
        sns.heatmap(heat, vmin=vlim[0], vmax=vlim[1], cmap='viridis', ax=ax, 
                    cbar=False, yticklabels=1)
        if index < num_groups - 2:
            ax.set_xticklabels("")
        if index in [1, 3, 5, 7, 9]:
            ax.set_yticklabels("")    
        else:            
            ax.set_yticklabels(heat.index, rotation=0)   
        ax.text(0, -0.6, labels[index], fontsize=14)
        plt.title("{}".format(dset))

    plt.savefig(p_j(save_folder, 'Spatiotemporal_{}_analysis_matrix_for_rest_6.png'.format(unc_type)), dpi=300)         
    #
    # Top of figure
    #
    
    red_fac = 0.7
    fig = plt.figure(figsize=(9., 8 * red_fac))  
    gs = gridspec.GridSpec(2, 3,
                       width_ratios=[0.85, 3, 0.1],
                       height_ratios=[4, 1.5]
                       )
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    plot_stream_reaches(m, ax1, zone2D_info, x_offset=8000)

    df = df_active.loc[novel_data_sets_alias[0]]
    heat = create_spatiotemporal_unc_array(df, months, reaches) 
    heat = heat.astype(np.float64)
    heat = heat * adjust
    heat = heat.transpose()
    heat = heat.reindex(heat.index.tolist()[::-1])
    sns.heatmap(heat, vmin=vlim[0], vmax=vlim[1], cmap='viridis', ax=ax2, cbar=False, yticklabels=1)
    ax2.set_yticklabels(heat.index, rotation=0)    
    plt.title("a. {}{}".format(title_prefix, novel_data_sets_alias[0]))
    ax3 = plt.subplot(gs[2])
    norm = mpl.colors.Normalize(vmin=vlim[0], vmax=vlim[1])
    cb1 = mpl.colorbar.ColorbarBase(ax3, cmap='viridis',
                                    norm=norm,
                                    orientation='vertical')
    cb1.set_label(cbar_label)   
    cb1.outline.set_visible(False)
    ax4 = plt.subplot(gs[4])
    plot_inflow_w_average(ax4, Eppalock_inflow, date_range, '30-06-2016', rows_name2,
                          lower=flow_low_quantile, upper=flow_high_quantile,
                          colour_dict=flow_type_colours)
    fig.subplots_adjust(bottom=0.02, left=0.08, right=0.90, top=0.94, hspace=0.40, wspace=0.21)
    plt.savefig(p_j(save_folder, 'Spatiotemporal_head_riverlay_flows_{}.png'.format(unc_type)), dpi=300)

    top_image = 'Spatiotemporal_head_riverlay_flows_{}.png'.format(unc_type)
    bottom_image = 'Spatiotemporal_{}_analysis_matrix_for_rest_6.png'.format(unc_type)
    img_top = Image.open(p_j(save_folder, top_image))
    img_bot = Image.open(p_j(save_folder, bottom_image))
    width = max(img_top.width, img_bot.width)
    height = img_top.height + img_bot.height
    new_img = Image.new('RGBA', (width, height), 'white')
    new_img.paste(img_top, (0, 0, width, img_top.height))    
    new_img.paste(img_bot, (0, img_top.height, width, img_top.height + img_bot.height))    
    #new_img.show()
    new_img.save(p_j(save_folder, 'Spatiotemporal_heat_maps_all_merged_{}.png'.format(unc_type)))    
#______________________________________________________________________________
# All data combos    
def powerset(iterable):
    s = list(set(iterable))
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def flatten_list(lst):
    flat = []
    for item in lst:
        flat += item
    return flat

def list_to_list_with_dict_expansion(lst, dicto):
    expanded = [dicto[x] for x in lst]
    expansion = flatten_list(expanded)
    return expansion

def get_cost(row, obs_costs):
    ind = row.name
    ind = ind.split('_')
    cost = 0.
    for letter in ind:
        cost += obs_costs[letter]
    return cost

def get_cost_unc(row, obs_costs_unc):
    ind = row.name
    ind = ind.split('_')
    cost_unc = 0.
    for letter in ind:
        cost_unc += obs_costs_unc[letter]
    return cost_unc

def all_combos(la, obs_groups_types, obs_groups_types_map, obs_groups_types_abbrev_map, save_folder):
    obs_exist = la.pst.observation_data
    combo_all_dict = {}
    for i, combo in enumerate(powerset(obs_groups_types), 1):
        if i == 0:
            # The first set will be empty so skip over this one
            continue
        
        #combo_all_sub_dict = {}
        combo_expansion = list_to_list_with_dict_expansion(combo, obs_groups_types_map)
        
        print('combo #{}: {}'.format(i, combo))
        print('combo #{}: {}'.format(i, combo_expansion))
    
        OBS_GROUPS_NAME = ["base"]
        temp = ""
        for x in combo:
            if temp == "":
                temp = obs_groups_types_abbrev_map[x]
            else:
                temp += "_{}".format(obs_groups_types_abbrev_map[x])
            # end if
            OBS_GROUPS_NAME += [temp]
    
        OBS_GROUPS = []
        for i2 in range(len(combo) + 1):
            OBS_GROUPS += [list_to_list_with_dict_expansion(combo[0:i2], obs_groups_types_map)]
    
        for ogn, og in zip(OBS_GROUPS_NAME, OBS_GROUPS):
            #combo_all_sub_dict[ogn] = all_obs(og, obs_exist)
            combo_all_dict[ogn] = all_obs(og, obs_exist)
        
        #combo_all_dict[i] = combo_all_sub_dict
    combo_all_dict.pop('base')
    
    obs_combos_df = la.get_added_obs_importance(obslist_dict=combo_all_dict)
    obs_combos_df.columns = [la.pst.observation_data.loc[x, 'obgnme'] for x in obs_combos_df.columns]
    
    return obs_combos_df    

def combo_cost_unc_plot(df, col_interest, flow_units='m3/d',
                        fname='Bar_unc_vs_cost_{}',
                        ax_text='h=head \ns=stage \nf=flow \ne=ec \nr=radon \nc=$^{14}$C'):
    fname = fname.format(col_interest)
    combo_df_lite = df[[col_interest, 'cost']]
    combo_df_lite[col_interest] = combo_df_lite[col_interest].apply(np.sqrt)
    if flow_units == 'Ml/yr':
        combo_df_lite['cost'] = combo_df_lite['cost'] 
    combo_df_lite.columns = ['uncertainty', 'cost']
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(1,1,1)
    combo_df_lite.sort_values('uncertainty').plot(kind='bar', y='uncertainty', ax=ax)    
    ax2 = combo_df_lite.sort_values('uncertainty').plot(kind='bar', y='cost', color='r', alpha = 0.3, secondary_y=True, ax=ax)    
    ax.set_ylabel("Uncertainty [m$^3$/d]")   
    ax2.set_ylabel("Cost [$]")
    ax.set_title('Comparison of uncertainty and cost for combinations of data types: {}'.format(col_interest))
    ax.text(-4, -10000, ax_text, style='italic',
            bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
    plt.tight_layout()
    plt.savefig(p_j(save_folder, fname), dpi=300)

def normalise_series(df):
    df = (df - df.min()) / (df.max() - df.min())
    return df

def cluster_combos_kmeans(combo_df, obs_map, col_interest=None, clusters=12, normalise=False):
    from sklearn.cluster import KMeans
    # Convert DataFrame to matrix
    if not col_interest:
        col_interest = 'nrf_a'

    combo_df_pred = combo_df[[col_interest, 'cost']]#, 'cost_unc']]
    combo_df_pred.columns = ['uncertainty', 'cost']#, 'cost_unc']
    if normalise:
        combo_df_pred.loc[:, 'cost_n'] = normalise_series(combo_df_pred['cost'])
        combo_df_pred.loc[:, 'uncertainty_n'] = normalise_series(combo_df_pred['uncertainty'])
        mat = combo_df_pred[['uncertainty_n', 'cost_n']].as_matrix()
    else:
        mat = combo_df_pred[['uncertainty', 'cost']].as_matrix()
    # end if
    # Using sklearn
    km = KMeans(n_clusters=clusters)
    km.fit(mat)
    # Get cluster assignment labels
    labels = km.labels_
    # Format results as a DataFrame
    results = pd.DataFrame(data=labels, columns=['cluster'])
    
    combo_df_pred.loc[:,'cluster'] = results['cluster'].tolist()
    
    cluster_name = []
    undefined = 0
    for i in range(clusters):
        par_comb_cluster = combo_df_pred[combo_df_pred['cluster'] == i].index.tolist()
        p = [x.split('_') for x in par_comb_cluster]
        result = set(p[0])
        for s in p[1:]:
            result.intersection_update(s)
        if len(result) == 0:
            cluster_name += ["Mixed: {}".format(['_'.join(t) for t in p])]
            undefined += 1
        else:
            #print par_comb_cluster
            new_result = [obs_map[x] for x in list(result)]
            cluster_name += [', '.join(list(new_result))]
        # end if
    
    combo_df_grouped = combo_df_pred.groupby('cluster').mean()
    combo_df_grouped.index = cluster_name
    return combo_df_pred, combo_df_grouped

def cluster_combos_dbscan(combo_df, eps=300, min_samples=12):
    from sklearn.cluster import DBSCAN
    # Convert DataFrame to matrix
    col_interest = 'nrf_a'
    combo_df_pred = combo_df[[col_interest, 'cost']]#, 'cost_unc']]
    combo_df_pred.columns = ['uncertainty', 'cost']#, 'cost_unc']
    mat = combo_df_pred[['uncertainty', 'cost']].as_matrix()
    # Using sklearn
    db = DBSCAN(eps=eps, min_samples=min_samples)
    db.fit(mat)
    # Get cluster assignment labels
    labels = db.labels_
    # Format results as a DataFrame
    results = pd.DataFrame(data=labels, columns=['cluster'])
    
    combo_df_pred.loc[:,'cluster'] = results['cluster'].tolist()
    print combo_df_pred
    cluster_name = []
    undefined = 0
    for i in range(len(results['cluster'].unique())):
        par_comb_cluster = combo_df_pred[combo_df_pred['cluster'] == i].index.tolist()
        p = [x.split('_') for x in par_comb_cluster]
        result = set(p[0])
        for s in p[1:]:
            result.intersection_update(s)
        if len(result) == 0:
            cluster_name += ["Mixed: {}".format(['_'.join(t) for t in p])]
            undefined += 1
        else:
            #print par_comb_cluster
            new_result = [obs_map[x] for x in list(result)]
            cluster_name += [', '.join(list(new_result))]
        # end if
    
    combo_df_grouped = combo_df_pred.groupby('cluster').mean()
    combo_df_grouped.index = cluster_name
    return combo_df_pred, combo_df_grouped

def scatter_unc_cost(combo_df_pred, combo_df_grouped, clusters=12, method='ml',
                     title="Data worth (existing data) for SW-GW exchange in the Campaspe River",
                     append_text='', xlim=(5000.0, 25000.0), ylim=(0,900),
                     xax_title="SW-GW exchange uncertainty [m$^3$/d]"):
    fig = plt.figure(figsize=(7,5))
    ax = fig.add_subplot(1,1,1)
    for cluster in range(clusters):
        combo_df_pred[combo_df_pred['cluster'] == cluster].plot(kind='scatter', x='uncertainty', 
                                                      y='cost', 
                                                      color=colors[cluster], 
                                                      alpha=0.1, ax=ax)
    
    for ind, cluster in enumerate(combo_df_grouped.index.tolist()):
        s = 50#combo_df_grouped[combo_df_grouped.index == cluster]['cost_unc']
        combo_df_grouped[combo_df_grouped.index == cluster].plot(kind='scatter', x='uncertainty', 
                                            y='cost', s=s, ax=ax,
                                            color=colors[ind])
    
    ax.set_xlabel(xax_title)   
    ax.set_ylabel("Cost [$]")
    ax.set_title(title)
    
    for ind, kv in enumerate(combo_df_grouped.iterrows()):
        k, v = kv
        if 'Mixed' in  k:
            k = 'Mixed'
        ax.annotate(k, xy=(v['uncertainty'], v['cost']), 
                    xytext=(5,2), textcoords='offset points',
                    size=10, color='black', alpha=0.9)
    
    #ax.text(3E8, 1000, 'Cluster annotations signifiy that \nobservation type is present in \nall points within cluster', style='italic',
    #        bbox={'facecolor':'gray', 'alpha':0.5, 'pad':10})
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    plt.tight_layout()
    plt.savefig(p_j(save_folder, "Cost_vs_predunc_scatter_{}{}".format(method, append_text)), dpi=300)    
    
#______________________________________________________________________________
#______________________________________________________________________________
#______________________________________________________________________________
# Potential obs
    
def mega_plot_river_predunc_red(df, pst_df, name, ob_interest, v_bounds=(0, 100)):
    months = ['June', 'July', 'August',
                 'September', 'October', 'November', 'December', 'January', 'February', 'March', 'April', 'May']

    fig = plt.figure(figsize=(8,3))
    ax = fig.add_subplot(111)
    nice_name = {'rn_sim': '${^{222}}$Rn',
                 'st_sim': 'Stage',
                 'ec_sim': 'EC',
                 'fl_sim': 'Flow'}
    for i in range(12):
        #ax = fig.add_subplot(12, 1, i + 1)
        ilocs = [12 * x + i for x in range(79)]
        c = (df[df.index.isin(pst_df[pst_df['obgnme'].isin([name])].index.tolist())][ob_interest].sort_index()).iloc[ilocs].tolist()
        #riv_pp_df.plot(kind='scatter', x='Cumulative Length', y='strtop', ax=ax, c=c, s=50, cmap='viridis', linewidth=0, vmin=0, vmax=max(c))
        riv_pp_df.loc[:, 'y_static'] = [i] * len(riv_pp_df)
        #riv_pp_df.plot(kind='scatter', x='Cumulative Length', y='y_static', ax=ax, c=c, s=50, cmap='viridis', linewidth=0, vmin=0, vmax=100.0) #max(c))
        ax2 = ax.scatter(x=riv_pp_df['Cumulative Length'], y=riv_pp_df['y_static'], c=c, s=50, cmap='viridis', linewidth=0, vmin=v_bounds[0], vmax=v_bounds[1]) #max(c))
        #fig.delaxes(fig.axes[1]) 
        
    ax.set_title("Reduction in predictive uncertainty for net annual exchange through {}".format(nice_name[name]))
    ax.set_xlabel('Chainage (m)')
    fig.colorbar(ax2)
    ax.set_ylim((-0.2, 11.2))
    ax.set_xlim((0, 143000))
    ax.set_ylabel('')
    ax.set_yticks([k for k in range(12)])
    ax.set_yticklabels(months, minor=False, rotation=0)
    #fig = plt.figure(figsize=(12,5))
        #ax = fig.add_subplot(111)
        #(Rad_test_unc[Rad_test_unc.index.isin(pst_df[pst_df['obgnme'].isin(['rn_sim'])].index.tolist())]['ob3454'].sort_index())[79*i:79*(i+1)].plot(kind='bar', ax=ax)
    plt.subplots_adjust(wspace=0, hspace=0)
    #plt.savefig('{}_pot_obs_nrf_a_all.png'.format(name))
    
def mega_plot_river_predunc_red_month(month, df, pst_df, names, ob_interest, v_bounds=(0, 100)):
    months = ['June', 'July', 'August',
                 'September', 'October', 'November', 'December', 'January', 'February', 'March', 'April', 'May']

    riv_pp_df.loc[:, 'Cumulative Length km'] = riv_pp_df['Cumulative Length'] / 1000.0
    fig = plt.figure(figsize=(7,3.5))
    ax = fig.add_subplot(2,1,1)
    nice_name = {'rn_sim': '${^{222}}$Rn',
                 'st_sim': 'Stage',
                 'ec_sim': 'EC',
                 'fl_sim': 'Flow'}
    for index, name in enumerate(names):
        #ax = fig.add_subplot(12, 1, i + 1)
        ilocs = [12 * x + month for x in range(79)]
        c = (df[df.index.isin(pst_df[pst_df['obgnme'].isin([name])].index.tolist())][ob_interest].sort_index()).iloc[ilocs].tolist()
        #riv_pp_df.plot(kind='scatter', x='Cumulative Length', y='strtop', ax=ax, c=c, s=50, cmap='viridis', linewidth=0, vmin=0, vmax=max(c))
        riv_pp_df.loc[:, 'y_static'] = [index] * len(riv_pp_df)
        #riv_pp_df.plot(kind='scatter', x='Cumulative Length', y='y_static', ax=ax, c=c, s=50, cmap='viridis', linewidth=0, vmin=0, vmax=100.0) #max(c))
        ax2 = ax.scatter(x=riv_pp_df['Cumulative Length km'], y=riv_pp_df['y_static'], c=c, s=50, cmap='viridis', linewidth=0, vmin=v_bounds[0], vmax=v_bounds[1]) #max(c))
        #fig.delaxes(fig.axes[1]) 
        
    ax.set_title("{}".format(months[month]))
    #fig.colorbar(ax2)
    ax.set_ylim((-0.2, float(len(names) - 1) + 0.2))
    ax.set_xlim((0, 143))
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_yticks([k for k in range(len(names))])
    ax.set_yticklabels([nice_name[name] for name in names], minor=False, rotation=0)
    ax.set_xticklabels('')
    ax.invert_yaxis()    
    
    # Adding the exchange flux
    ax = fig.add_subplot(2,1,2)
    month_num = 32 - 11 + month
    exchange = sfr_df[sfr_df['time'] == month_num].reset_index(
                range(sfr_df.shape[0]))['Qaquifer'] / (sfr_orig['rchlen'] * sfr_seg['width2']).tolist()
    c = np.array(sfr_df[sfr_df['time'] == month_num].reset_index(
                range(sfr_df.shape[0]))['Qaquifer'].tolist()) #[0 if i > 0 else 1 for i in exchange]
    max_val = np.max(c)
    min_val = np.min(c)
    if max_val > abs(min_val):
        min_val = -max_val
    else:
        max_val = -min_val
    # End if

    #ax.fill_between((0, 143), 0, -max(abs(min(exchange)), max(exchange)), facecolor='grey', alpha=0.4)
    ax.axhline(0.0, linestyle='--', color='grey', alpha=0.8, linewidth=1)
    ax3 = ax.scatter(sfr_orig['Cumulative Length km'].tolist(), 
               exchange, c=c, cmap='seismic', linewidth=0.5, alpha=1.0, vmin=min_val, vmax=max_val)
    ax.axvline(73.4, linestyle='-', color='black', alpha=0.8, linewidth=1.5)
    ax.set_xlabel('Chainage (km)')
    ax.set_xlim((0, 143))
    ax.set_ylabel('Exchange (m/d)')
    ax.text(0.8, 0.8, 'Losing', color='red', #horizontalalignment='center',
         transform=ax.transAxes)
    ax.text(0.8, 0.1, 'Gaining', color='blue', #horizontalalignment='center',
         transform=ax.transAxes)
    yticks = ax.get_yticks().tolist() # get list of ticks
    yticks[-1] = ''                          # set last tick to empty string
    #yticks[-2] = ''                          # set last tick to empty string
    ax.set_yticklabels(yticks)               # set the labels

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.50, 0.03, 0.4])
    fig.colorbar(ax2, cax=cbar_ax)    
    cbar_ax.get_yaxis().labelpad = 18
    cbar_ax.set_ylabel('Uncertainty reduction (%)', rotation=270)

    cbar_ax2 = fig.add_axes([0.85, 0.12, 0.03, 0.35])
    fig.colorbar(ax3, cax=cbar_ax2)    
    cbar_ax2.get_yaxis().labelpad = 16
    cbar_ax2.set_ylabel('Reach Exchange (m${^3}$/d)', rotation=270)
    
    #fig = plt.figure(figsize=(12,5))
        #ax = fig.add_subplot(111)
        #(Rad_test_unc[Rad_test_unc.index.isin(pst_df[pst_df['obgnme'].isin(['rn_sim'])].index.tolist())]['ob3454'].sort_index())[79*i:79*(i+1)].plot(kind='bar', ax=ax)
    plt.subplots_adjust(wspace=0, hspace=0, bottom=0.1)
    plt.savefig(p_j(save_folder, 'Instream_pot_obs_nrf_a_by_{}.png'.format(months[month])), dpi=300)    
    
def mega_plot_river_predunc_red_month_axes(select_months, df, pst_df, names, 
                                           ob_interest, sfr_df, sfr_orig,
                                           riv_pp_df, save_folder, v_bounds=(0, 100),
                                           color_bar_exchange=[0.85, 0.663, 0.03, 0.286]):
    months = ['June 2016', 'July 2016', 'August 2016',
                 'September 2016', 'October 2016', 'November 2016', 'December 2016', 
                 'January 2017', 'February 2017', 'March 2017', 'April 2017', 'May 2017']

    fig = plt.figure(figsize=(7, len(select_months) * 3.5))

    for index, select_month in enumerate(select_months):
        month_num = 32 - 11 + select_month
        exchange = sfr_df[sfr_df['time'] == month_num].reset_index(
                    range(sfr_df.shape[0]))['Qaquifer'] / (sfr_orig['rchlen']).tolist() # * sfr_seg['width2']).tolist()
        c = np.array(sfr_df[sfr_df['time'] == month_num].reset_index(
                    range(sfr_df.shape[0]))['Qaquifer'].tolist()) #[0 if i > 0 else 1 for i in exchange]
        if index == 0:
            max_val = np.max(c)
            min_val = np.min(c)
        else:
            max_val = max(max_val, np.max(c))
            min_val = min(min_val, np.min(c))
    # End for
    if max_val > abs(min_val):
        min_val = -max_val
    else:
        max_val = -min_val
    # End if

    colors_ = ['orangered', 'dodgerblue', 'mediumseagreen']
    extra_text = [' (low)', ' (high)', ' (regular)']    

    for index, select_month in enumerate(select_months):
        axes = [fig.add_subplot(len(select_months) * 2, 1, 1 + 2 * index), 
                fig.add_subplot(len(select_months) * 2, 1, 2  + 2 * index)]
        if index == 0:
            legend = True
        else:
            legend = False
        if index == len(select_months) - 1:
            skip_xlabels = False
        else:
            skip_xlabels = True
                 
        riv_pp_df.loc[:, 'Cumulative Length km'] = riv_pp_df['Cumulative Length'] / 1000.0
        ax = axes[0]
        nice_name = {'rn_sim': '${^{222}}$Rn',
                     'st_sim': 'Stage',
                     'ec_sim': 'EC',
                     'fl_sim': 'Flow'}
        markers = ['*', 'x', '+', 'D']    
        size = [5, 5, 5, 2]         
        #colours = ['red', 'blue', 'black', 'green']# ['#f0f9e8', '#bae4bc', '#7bccc4', '#2b8cbe']#plt.get_cmap('viridis', 4)
        colours=['#d95f02', '#7570b3', '#e7298a', '#66a61e']

        for index2, name in enumerate(names):
            #ax = fig.add_subplot(12, 1, i + 1)
            ilocs = [12 * x + select_month for x in range(79)]
            if type(ob_interest) == list:
                c = (df[df.index.isin(pst_df[pst_df['obgnme'].isin([name])].index.tolist())][ob_interest[index]].sort_index()).iloc[ilocs].tolist()
            else:
                c = (df[df.index.isin(pst_df[pst_df['obgnme'].isin([name])].index.tolist())][ob_interest].sort_index()).iloc[ilocs].tolist()
            if len(riv_pp_df['Cumulative Length km'].tolist()) > len(c):
                tmp = riv_pp_df['Cumulative Length km'].tolist()[:-1]
            else:
                tmp = riv_pp_df['Cumulative Length km'].tolist()
            ax.plot(tmp, c, color=colours[index2], 
                          marker=markers[index2], ms=size[index2], mew=2.0, 
                          mec=colours[index2], label=nice_name[name]) #max(c))
            #fig.delaxes(fig.axes[1]) 
            
        ax.axvline(73.4, linestyle='-', color='black', alpha=0.3, linewidth=1.0)
        font = {#'family': 'serif',
                'color':  colors_[index],
                'weight': 'heavy',
                'size': 12, 
            }
        ax.text(-0.175, 0.05, "{}{}".format(months[select_month], extra_text[index]), 
                rotation=90, fontdict=font, horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes)
        #fig.colorbar(ax2)
        #if v_bounds:
        #    ax.set_ylim((v_bounds[0], v_bounds[1]))
        #else:
        #    ax.set_ylim((-2, 60))
        ax.set_xlim((0, 143))
        ax.set_xlabel('')
        ax.set_ylabel('Uncertainty\nreduction (%)')
        #ax.set_yticks([k for k in range(len(names))])
        #ax.set_yticklabels([nice_name[name] for name in names], minor=False, rotation=0)
        ax.set_xticklabels('')
        #ax.invert_yaxis()    
        if legend:
            ax.legend(numpoints=1, fontsize=10, 
                      bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, 
                      mode="expand", borderaxespad=0.1)
        
        # Adding the exchange flux
        ax = axes[1]
        month_num = 32 - 11 + select_month
        exchange = sfr_df[sfr_df['time'] == month_num].reset_index(
                    range(sfr_df.shape[0]))['Qaquifer'] / (sfr_orig['rchlen']).tolist() # * sfr_seg['width2']
        c = np.array(sfr_df[sfr_df['time'] == month_num].reset_index(
                    range(sfr_df.shape[0]))['Qaquifer'].tolist()) #[0 if i > 0 else 1 for i in exchange]
#        max_val = np.max(c)
#        min_val = np.min(c)
#        if max_val > abs(min_val):
#            min_val = -max_val
#        else:
#            max_val = -min_val
#        # End if
    
        #ax.fill_between((0, 143), 0, -max(abs(min(exchange)), max(exchange)), facecolor='grey', alpha=0.4)
        ax.axhline(0.0, linestyle='--', color='grey', alpha=0.8, linewidth=1)

        sfr_df.loc[:, 'Cumulative Length (km)'] = sfr_df['Cumulative Length'] / 1000.
        sfr_df['Qaquifer_adj'] = sfr_df['Qaquifer'] / (sfr_df['rchlen']) # * sfr_df['width'])
        sfr_df_mean = sfr_df.groupby('Cumulative Length').mean()
        #for t in sfr_df['time'].unique()[10:]:
        #    df_t = sfr_df[sfr_df['time'] == t]
        #    ax.plot(df_t['Cumulative Length (km)'], df_t['Qaquifer_adj'], color='gray', alpha=0.5)    
        ax.plot(sfr_df_mean['Cumulative Length (km)'], sfr_df_mean['Qaquifer_adj'], color='black', linestyle='-', alpha=0.5)    

        ax3 = ax.scatter(sfr_orig['Cumulative Length km'].tolist(), 
                   exchange, c=c, cmap='seismic', linewidth=0.5, alpha=1.0, vmin=min_val, vmax=max_val)
        ax.axvline(73.4, linestyle='-', color='black', alpha=0.3, linewidth=1.0)
        ax.set_xlabel('Distance from lake Eppalock (km)')
        ax.set_xlim((0, 143))
        ax.set_ylabel('Exchange (m$^2$/d)')
        #ax.text(0.8, 0.8, 'Losing', color='red', #horizontalalignment='center',
        #     transform=ax.transAxes)
        #ax.text(0.8, 0.1, 'Gaining', color='blue', #horizontalalignment='center',
        #     transform=ax.transAxes)
        yticks = ax.get_yticks().tolist() # get list of ticks
        yticks[-1] = ''                          # set last tick to empty string
        #yticks[-2] = ''                          # set last tick to empty string
        ax.set_yticklabels(yticks)               # set the labels
        if skip_xlabels:
            ax.set_xticklabels('')
            ax.set_xlabel('')

    fig.subplots_adjust(right=0.8)
    cbar_ax2 = fig.add_axes(color_bar_exchange)
    fig.colorbar(ax3, cax=cbar_ax2)    
    cbar_ax2.get_yaxis().labelpad = 16
    cbar_ax2.set_ylabel('Reach Exchange (m${^3}$/d)', rotation=270)
    cbar_ax2.text(0.1, 0.85, 'Losing', color='white', #horizontalalignment='center',
                  transform=cbar_ax2.transAxes, rotation=270)
    cbar_ax2.text(0.1, 0.3, 'Gaining', color='white', #horizontalalignment='center',
                  transform=cbar_ax2.transAxes, rotation=270)
    
    plt.subplots_adjust(wspace=0, hspace=0.1, bottom=0.06, top=0.95, left=0.14, right=0.82)
    plt.savefig(p_j(save_folder, 'Instream_pot_obs_nrf_a_stacked.png'), dpi=300)        

def plot_SWGW_exchange(sfr_df, save_folder, ax=None, show_gauges=False, fontsize=8,
                       inflow_data=None, colour_dict=None, date_index=None,
                       plot_only=None, adjust_left=0.14,  adjust_right=0.97, linewidth=0.8, bbox_to_anchor=(1.05, 1),
                       alpha=0.5):
    
    sfr_df.loc[:, 'Cumulative Length (km)'] = sfr_df['Cumulative Length'] / 1000.
    sfr_df['Qaquifer_adj'] = sfr_df['Qaquifer'] / (sfr_df['rchlen'])# * sfr_df['width'])
    sfr_df_mean = sfr_df.groupby('Cumulative Length').mean()
    if not ax:
        if colour_dict is None:
            fig = plt.figure(figsize=(5, 2))
        else:
            fig = plt.figure(figsize=(6, 2))
            
        ax = fig.add_subplot(111)
    # end if
    #ax2 = fig.add_subplot(212)

    for t in sfr_df['time'].unique()[6:]:
        df_t = sfr_df[sfr_df['time'] == t - 1]
        if colour_dict is None:
            ax.plot(df_t['Cumulative Length (km)'], df_t['Qaquifer_adj'], color='gray', alpha=0.5)    
        else:
            flow_type = inflow_data.loc[date_index[int(t)], 'flow_group']
            if plot_only is not None:
                if flow_type != plot_only:
                    continue
                print(date_index[int(t)])
            ax.plot(df_t['Cumulative Length (km)'], df_t['Qaquifer_adj'], color=colour_dict[flow_type], alpha=alpha)    
            
        #ax2.plot(df_t['Cumulative Length (km)'], df_t['stage'], color='blue', alpha=0.5)    
        #ax2.plot(df_t['Cumulative Length (km)'], df_t['head'], color='green', alpha=0.5)    
    #ax.axvline(74.4, color='black', linestyle='--')
    ax.axhline(0.0, color='black', linestyle=':', alpha=9.7)
    #ax2.axvline(74.4, color='black')
    ax.plot(sfr_df_mean['Cumulative Length (km)'], sfr_df_mean['Qaquifer_adj'], color='black', linestyle='-', alpha=1.0)
    ax.set_xlim(0, 141)    
    ax.set_ylabel('SW-GW exchange [m$^2$/d]', fontsize=fontsize)
    ax.set_xlabel('Distance from Lake Eppalock (km)', fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)#ax2.set_xlim(0, 140)    
    if show_gauges:
        x_stream_gauges = sfr_df[sfr_df['time'] == 0].groupby(by='reach_dw').sum()['rchlen']
        x_stream_gauges = x_stream_gauges / 1000.
        x_stream_gauges =  x_stream_gauges.cumsum().tolist()
        for ind, gauge in enumerate(x_stream_gauges):
            ax.axvline(gauge, color='gray', alpha=0.8, linestyle=':')
            if ind == 0:
                text_loc = gauge / 2.
            elif ind == 9:
                text_loc = (gauge + x_stream_gauges[ind - 1]) / 2. - 3
            else:
                text_loc = (gauge + x_stream_gauges[ind - 1]) / 2.
            # end if
            ax.text(text_loc - 2, 
                    ax.get_ylim()[1] + 0.01 * (ax.get_ylim()[1] - ax.get_ylim()[0]), 
                    'r{}'.format(ind + 1), fontsize=8)
    #ax2.set_ylabel('')
    if colour_dict:
        silver_patch = Line2D([0], [0], color='black', label='mean', linewidth=linewidth)
        red_patch = Line2D([0], [0], color=colour_dict['low'], label='low', alpha=alpha, linewidth=linewidth)
        orange_patch = Line2D([0], [0], color=colour_dict['regular'], label='regular', alpha=alpha, linewidth=linewidth)
        blue_patch = Line2D([0], [0], color=colour_dict['high'], label='high', alpha=alpha, linewidth=linewidth)
        handles = [silver_patch, red_patch, orange_patch, blue_patch]
        leg = plt.legend(handles=handles, ncol=1, fontsize=8, 
                   frameon=False, bbox_to_anchor=bbox_to_anchor, loc=2, title='Inflow type')
        leg.set_title('Flow conditions', prop={'size':8})
        
    plt.subplots_adjust(hspace=0.1, top=0.93, bottom=0.22, left=adjust_left, right=adjust_right)

    plt.savefig(p_j(save_folder, 'SWGW_exchange_all.png'), dpi=300)

def zone_array2layers(zone_array, plots=False):
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

def plot_ground_pot_obs_predunc(zone2D_info, obtype, dfs, hgus, vlim=(0, 60)):
    cmap_grey_white = mpl_colors.ListedColormap(['white', 'lightgrey'])
    fig = plt.figure()
    plot_titles = ["{} {}".format(x, obtype) for x in ['Shallow', 'Deep']]
    for index in range(2):
        hgu = hgus[index]
        sim_locs_bool = zone2D_info[0][hgu[0]] | zone2D_info[0][hgu[1]]
        ax = fig.add_subplot(1, 2, index + 1, aspect='equal')
        ax.set_title(plot_titles[index])
        ax.axis('off')
        plt.imshow(sim_locs_bool, interpolation='none', cmap=cmap_grey_white)
        #x_dw = [i[2] for i in sim_locs]
        #y_dw = [i[1] for i in sim_locs]
        #c_dw = [i[0] for i in sim_locs]
        #plt.scatter(x_dw, y_dw, c=c, cmap='viridis', linewidth=0.0) #, c=c_dw, cmap='viridis')
        #dfs[index].plot(kind='scatter', x='x', y='y', c=dfs[index]['perc_red_unc'], cmap='viridis', linewidth=0.0, ax=ax)
        ax.scatter(dfs[index]['x'], dfs[index]['y'], c=dfs[index]['perc_red_unc'], 
                   cmap='viridis', linewidth=0.0, vmin=vlim[0], vmax=vlim[1])
    plt.savefig(p_j(save_folder, 'Shallow and deep perc red for {}'.format(obtype)), dpi=300)

def plot_ground_pot_obs_predunc_stream(obtype, dfs, hgus, m, axes, zone2D_info, vlim=None):
    plot_titles = ["{} {}".format(x, obtype) for x in ['Shallow', 'Deep']]
    for index in range(2):
        hgu = hgus[index]
        sim_locs_bool = zone2D_info[0][hgu[0]] | zone2D_info[0][hgu[1]]
        ax = axes[index]
        plot_stream_reaches_basic(m, ax, zone2D_info)
        ax.set_title(plot_titles[index])
        #ax.axis('off')
        plt.imshow(sim_locs_bool, interpolation='none', cmap=cmap_grey_white)
        ax.scatter(dfs[index]['x'], dfs[index]['y'], c=dfs[index]['perc_red_unc'], 
                   cmap='viridis', linewidth=0.0, vmin=vlim[0], vmax=vlim[1])
    plt.savefig(p_j(save_folder, 'Shallow and deep perc red w existing locs for {}'.format(obtype)), dpi=300)
    
    
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

def plot_inflow(ax, inflow_df, date_range, start_date, rows_name2):
    Eppalock_inflow = inflow_df
    #Eppalock_inflow = Eppalock_inflow.reindex(date_range)
    Eppalock_inflow['Mean'].fillna(Eppalock_inflow['Mean'].mean(), inplace=True)
    Eppalock_inflow[Eppalock_inflow.index > start_date].resample('M').mean().plot(
        kind='bar', y='Mean', label='Monthly mean', ax=ax, color='IndianRed', legend=False, linewidth=0) #, marker='o'
    #Eppalock_inflow[Eppalock_inflow.index > start_date].resample('D').mean().plot(y='Mean', ax=ax, label='Daily mean')
    ax.set_ylabel("Eppalock\n Inflow [Gl/yr]")
    ax.set_yscale("log", nonposy='clip')
    #ax.legend(loc='center left', bbox_to_anchor=(-0.5, 0.5), prop={'size': 9})
    ax.set_xlabel('')

    #minorLocator = mpl.ticker.MultipleLocator(300)
    #ax.xaxis.set_minor_locator(minorLocator)
    ax.set_xticklabels([""] * len(rows_name2), minor=False, rotation=90)
    ax.grid(True)
    #ax.set_yticks([k for k in range(12)])
    #ax.set_yticklabels(months, minor=False, rotation=0)

    return ax

def plot_inflow_w_average(ax, inflow_df, date_range, start_date, rows_name2, 
                          upper=0.75, lower=0.25, colour_dict=None, plot_quantiles=False):
    Eppalock_inflow = inflow_df
    Eppalock_inflow['Mean'].fillna(Eppalock_inflow['Mean'].mean(), inplace=True)

    Eppalock_inflow_monthly = Eppalock_inflow.resample('M').mean()
    q_lower = Eppalock_inflow['Mean'].quantile(lower)
    q_upper = Eppalock_inflow['Mean'].quantile(upper)
    
    Eppalock_inflow_monthly.loc[:, 'flow_group'] = np.nan    
    Eppalock_inflow_monthly.loc[(Eppalock_inflow_monthly['Mean'] < q_lower).tolist(), 'flow_group'] = 'low'
    Eppalock_inflow_monthly.loc[(Eppalock_inflow_monthly['Mean'] > q_upper).tolist(), 'flow_group'] = 'high'
    Eppalock_inflow_monthly.loc[pd.isnull(Eppalock_inflow_monthly['flow_group']), 'flow_group'] = 'regular'

    if colour_dict is None:
        flow_type_colours = {'low':'red', 'high':'blue', 'regular':'orange'}                                
    else:
        flow_type_colours = colour_dict
    # end if 
    
    Eppalock_inflow_monthly.loc[:, 'colour'] = Eppalock_inflow_monthly.apply(lambda x: flow_type_colours[x['flow_group']], axis=1)

    #Mean_Eppalock = Eppalock_inflow[Eppalock_inflow.index > start_date].resample('M').mean()
    Mean_Eppalock = Eppalock_inflow_monthly[Eppalock_inflow_monthly.index >= start_date]
    # Hack to add inflow applied in model which was just the last value padded forward
    Mean_Eppalock = Mean_Eppalock.append(pd.DataFrame(data={'Mean': [Mean_Eppalock['Mean'][-1]], 
                                                            'Qual':[0], 
                                                            'colour':[Mean_Eppalock['colour'][-1]]}, 
                                                      index=[pd.datetime(2017,5,31)]))
    # Insert the average at the start of series to assist in the plotting of the bar chart
    AverageEppalock = pd.DataFrame(data={'Mean': [Mean_Eppalock['Mean'].mean()], 'Qual':[0]}, index=[pd.datetime(1900,1,1)])
    Mean_Eppalock = pd.concat([AverageEppalock, Mean_Eppalock])

    colours_list = ['Silver'] + Mean_Eppalock['colour'].tolist()[1:]
    
    #Mean_Eppalock.plot(kind='bar', y='Mean', label='Monthly mean', ax=ax, 
    #                   color=['Silver'] + ['Teal'] * 12, legend=False, linewidth=0) #, marker='o'
    Mean_Eppalock.plot(kind='bar', y='Mean', label='Monthly mean', ax=ax, 
                       color=colours_list, legend=False, linewidth=0) #, marker='o'
    if plot_quantiles:
        ax.axhline(q_lower, linestyle='--', color='grey', alpha=0.5)
        ax.axhline(q_upper, linestyle='--', color='grey', alpha=0.5)
    # end if
    ax.set_ylabel("Eppalock\n Inflow [Gl/yr]")
    ax.set_yscale("log", nonposy='clip')
    ax.set_xlabel('')
    
    silver_patch = mpatches.Patch(color='silver', label='mean')
    red_patch = mpatches.Patch(color=flow_type_colours['low'], label='low')
    orange_patch = mpatches.Patch(color=flow_type_colours['regular'], label='regular')
    blue_patch = mpatches.Patch(color=flow_type_colours['high'], label='high')
    handles = [silver_patch, red_patch, orange_patch, blue_patch]
    plt.legend(handles=handles, ncol=1, fontsize=10, 
               frameon=False, bbox_to_anchor=(1.05, 1), loc=2, title='Flow type')
    
    ax.set_xticklabels([""] * len(rows_name2), minor=False, rotation=90)
    ax.grid(True)
    return ax
    
def plot_stream_reaches(m, ax, zone2D_info, new_fig=False, x_offset=5000., c2=['red', 'blue']):
    model_mesh3D = m.model_mesh3D
    model_boundary = m.model_boundary
    river_seg = m.river_mapping['Campaspe']
    
    river_segs_reach = [river_seg['iseg'][river_seg['reach'] == x].tolist() for x in river_seg['reach'].unique()]
    nrow = model_mesh3D[0].shape[1]
    ncol = model_mesh3D[0].shape[2]
    delr = m.gridHeight
    delc = m.gridWidth
    #top = model_mesh3D[0][0]
    #botm = self.model_data.model_mesh3D[0][1:]
    xul = model_boundary[0]
    yul = model_boundary[3]
    
    x = np.linspace(xul, xul + ncol * delc, ncol)
    y = np.linspace(yul - nrow * delr, yul, nrow)
    X, Y = np.meshgrid(x, y)

    #ax.pcolormesh(X, Y, np.flipud(zone2D_info[0][6]))
    flatten = lambda l: [item for sublist in l for item in sublist]
    cmap_grey_white = mpl_colors.ListedColormap(['white', 'lightgrey'])
    
    if new_fig:
        fig = plt.figure(figsize=(3.5,7))
        ax = fig.add_subplot(1,1,1, aspect='equal')

    ax.pcolormesh(X, Y, np.flipud(zone2D_info[0][6]), cmap=cmap_grey_white)
    
    for index, reach in enumerate(river_segs_reach[:]):
        reach_river = river_seg[river_seg['iseg'].isin(reach)]
        points = [i for i in reach_river['amalg_riv_points_collection']]
        points = flatten(points)
        reach_dist = 0
        for ind, point in enumerate(points[:-1]):
            reach_dist += np.linalg.norm(np.array(point) - np.array(points[ind + 1]))
        x_points = [i[0] for i in points]
        y_points = [j[1] for j in points]    
        ax.plot(x_points, y_points, c=c2[index % 2])
        midx = x_points[len(x_points)/2] - x_offset
        midy = y_points[len(y_points)/2]
        #ax.scatter(midx, midy)
        ax.text(midx, midy, 'r{}'.format(index + 1), color=c2[index % 2],  fontsize=9)
        ax.text(midx + 1.2 * x_offset, midy, '{0:.1f} km'.format(reach_dist / 1000.), color=c2[index % 2],  fontsize=9)
    #Campaspe_info.plot(kind='scatter', x='Easting', y='Northing', ax=ax)#, label='Site Id')    
    start_ax, end_ax = ax.get_xlim()
    start_ax = start_ax // 1000 * 1000 + 1000
    end_ax = end_ax // 1000 * 1000 - 1000
    ax.xaxis.set_ticks(np.arange(start_ax, end_ax, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
    #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Reach Number')
    ax.set_xlim(2.65E5, xul + ncol * delc) #(xul, xul + ncol * delc)
    ax.set_ylim(5915000, 6005000)#(yul - nrow * delr, 6010000)#yul)    
    ax.set_xlabel('Easting')
    ax.set_ylabel('Northing')
    return ax

def plot_stream_reaches_basic(m, ax, zone2D_info, new_fig=False, bounds=None, zones=None):
    model_mesh3D = m.model_mesh3D
    model_boundary = m.model_boundary
    river_seg = m.river_mapping['Campaspe']
    
    river_segs_reach = [river_seg['iseg'][river_seg['reach'] == x].tolist() for x in river_seg['reach'].unique()]
    nrow = model_mesh3D[0].shape[1]
    ncol = model_mesh3D[0].shape[2]
    delr = m.gridHeight
    delc = m.gridWidth
    xul = model_boundary[0]
    yul = model_boundary[3]
    
    x = np.linspace(xul, xul + ncol * delc, ncol)
    y = np.linspace(yul - nrow * delr, yul, nrow)
    X, Y = np.meshgrid(x, y)

    flatten = lambda l: [item for sublist in l for item in sublist]
    cmap_grey_white = mpl_colors.ListedColormap(['white', 'lightgrey'])
    
    if new_fig:
        fig = plt.figure(figsize=(3.5,7))
        ax = fig.add_subplot(1,1,1, aspect='equal')

    if not zones:
        ax.pcolormesh(X, Y, np.flipud(zone2D_info[0][6]), cmap=cmap_grey_white)
    else:
        for i, zone in enumerate(zones):
            if i == 0:
                zone_temp = zone2D_info[0][zone]
            else:
                zone_temp = zone_temp | zone2D_info[0][zone]
        ax.pcolormesh(X, Y, np.flipud(zone_temp), cmap=cmap_grey_white)
        
    for index, reach in enumerate(river_segs_reach[:]):
        reach_river = river_seg[river_seg['iseg'].isin(reach)]
        points = [i for i in reach_river['amalg_riv_points_collection']]
        points = flatten(points)
        reach_dist = 0
        for ind, point in enumerate(points[:-1]):
            reach_dist += np.linalg.norm(np.array(point) - np.array(points[ind + 1]))
        x_points = [i[0] for i in points]
        y_points = [j[1] for j in points]    
        ax.plot(x_points, y_points, c='blue', alpha=0.5)
    #Campaspe_info.plot(kind='scatter', x='Easting', y='Northing', ax=ax)#, label='Site Id')    
    start_ax, end_ax = ax.get_xlim()
    start_ax = start_ax // 1000 * 1000 + 1000
    end_ax = end_ax // 1000 * 1000 - 1000
    ax.xaxis.set_ticks(np.arange(start_ax, end_ax, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
    if bounds:
        ax.set_xlim(bounds[0], bounds[1])
        ax.set_ylim(bounds[2], bounds[3])    
    ax.set_xlabel('Easting', fontsize=10)
    ax.set_ylabel('Northing', fontsize=10)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    return ax
    
def get_heads_for_sfr(sfr_df, modflow_model, data_folder, dry_nan=False):
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

def process_sfr_df_and_sfr_info_into_reaches(sfr_df, sfr_info):
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
    sfr_df_reaches_relevant.loc[:, 'Cumulative Length (km)'] = sfr_df_reaches_relevant['Cumulative Length'] / 1000.

    return sfr_df, sfr_df_reaches_relevant

def process_Eppalock_inflow(Eppalock_inflow,
                          qual_codes_to_ignore=[8, 9, 21, 100, 101, 120, 149, 150, 151, 152,
                                                153, 154, 155, 156, 160, 161, 165, 180, 190,
                                                200, 201, 237, 250, 254, 255],
                          flow_high_quantile = 0.8,
                          flow_low_quantile = 0.35):

    Eppalock_inflow = Eppalock_inflow[~Eppalock_inflow['Qual'].isin(qual_codes_to_ignore)]
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
    return Eppalock_inflow_monthly
    
def get_simple_transport_setup(m, modflow_model, sfr_df, sfr_info):
    num_reaches = m.pilot_points['Campaspe'].num_points
    known_points = m.pilot_points['Campaspe'].points
    
    # Hyporheic zone depth         
    hz_depth_vals = [m.parameters.param['hzdpth{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    R_depth_HZ = np.interp(sfr_info['Cumulative Length'].tolist(), known_points, hz_depth_vals)

    df_size = sfr_info.shape[0]
    t_steps = sfr_df['time'].max() + 1
    # Hyporheic zone porosity
    sfr_df.loc[:, 'HZ_poro'] = [m.parameters.param['hzporo']['PARVAL1']] * df_size * t_steps
    # Hyporheic zone production of radon
    sfr_df.loc[:, 'HZ_Prod_Rate'] = [m.parameters.param['hzprod']['PARVAL1']] * df_size * t_steps
    # Hyporheic zone residence time
    sfr_df.loc[:, 'HZ_RTime'] = [m.parameters.param['hz_rt']['PARVAL1']] * df_size * t_steps
    # Depth of the hyporheic zone
    sfr_df.loc[:, 'R_depth_HZ'] = R_depth_HZ.tolist() * t_steps              
    # Gas transfer velocity
    sfr_df.loc[:, 'GTV'] = [m.parameters.param['gtv']['PARVAL1']] * df_size * t_steps
    # Groundwater radon concentration
    sfr_df.loc[:, 'GW_Rn_conc'] = [m.parameters.param['gwconc']['PARVAL1']] * df_size * t_steps
    # Groundwater EC
    sfr_df.loc[:, 'GW_EC'] = [10.] * df_size * t_steps # ARBITRARY!!!!
    # EC of the inflowing tributary water if present
    sfr_df.loc[:, 'Trib_EC'] = [10.] * df_size * t_steps # ARBITRARY!!!! 
    # Radon concentration of inflowing tributary water if present
    sfr_df.loc[:, 'Trib_Rn'] = [10.] * df_size * t_steps # ARBITRARY!!!!
    # Reach lengths
    sfr_df.loc[:, 'dx'] = sfr_info['rchlen'].tolist() * t_steps 

    return sfr_df

def create_sfr_sim_plots(m, modflow_model):
    
    date_index = m.model_time.t['dateindex']
    df_list = []
    # March 2016, Dec 2016, May 2017
    import matplotlib.pyplot as plt
    
    field_sampling = [datetime.datetime(2016,03,31),
                      datetime.datetime(2016,12,31),
                      datetime.datetime(2017,04,30)]

    import brewer2mpl
    bmap = brewer2mpl.get_map('Set1', 'Qualitative', 5)
    palette = bmap.mpl_colors
    handles = []
    labels = ['03-04 / 2016', '12 / 2016', '05 / 2017', '01 / 2018']

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(0,33): #sfr_df['time'].max()):
        if date_index[i] in field_sampling:    
            df = sfr_df[sfr_df['time'] == i]
            Ini_cond = (df.iloc[0]['Qin'], 0., 300.)
            df.loc[:, 'Flow'], df.loc[:, 'Rn'], df.loc[:, 'EC'] = modflow_model.Calculate_Rn_from_SFR_with_simple_model(df, Ini_cond)
            df['Cumulative Length km'] = df['Cumulative Length'] / 1000.
            df_list += [df[['Cumulative Length', 'Flow', 'Rn', 'EC']]]
            df.plot(x='Cumulative Length km', y='Rn', style='-', ax=ax, label=date_index[i].date())
    ax.set_ylabel("Radon (mBq/L)")
    ax.set_xlabel("River chainage (m)")

    fig2 = plt.figure(figsize=(12, 12))
    plt.subplots_adjust(hspace = 0.1)
    n_subplots = len(field_sampling)
    counter = 0
    for i in range(0,33): #sfr_df['time'].max()):
        if date_index[i] in field_sampling:    
            counter += 1
            df = sfr_df[sfr_df['time'] == i]
            df['Qaquifer_adj'] = -df['Qaquifer'] / (df['dx'])# * df['width'])
            df['Cumulative Length km'] = df['Cumulative Length'] / 1000.
            df_list += [df[['Cumulative Length', 'Flow', 'Rn', 'EC']]]
            df.plot(x='Cumulative Length', y='Rn', style='-', ax=ax, label=date_index[i].date())
            ax2 = fig2.add_subplot(n_subplots, 1, counter)
            #df.plot(x='Cumulative Length', style='o',  y=['Flow', 'Qin'], ax=ax2, label=date_index[i].date())
            df[df['Qaquifer_adj'] != 0.].rolling(3).mean().dropna().plot(x='Cumulative Length km', y='Qaquifer_adj', style='-', ax=ax2, color=palette[counter - 1])#, label=date_index[i].date())  
            ax2.set_ylabel("SW-GW fluxes (m$^3$/m/day)")
            ax2.set_xticks(np.arange(0,160,20))
            ax2.set_xticklabels([])
            ax2.set_ylim(-3., 3.)
            ax2.set_xlabel('')
            ax2.axhline(0., color='black', linestyle='--')
            ax2.axvline(73.4, color='black', linestyle='-', lw=2.)
    
    ax2.set_xticklabels(np.arange(0,160,20))
    ax2.set_xlabel('Distance from lake Eppalock (km)')
    fig2.tight_layout(rect=[0.02,0.05,1,0.99])
    
    import matplotlib.patches as mpatches

    for Camp in labels:
        handles.append(mpatches.Patch(color=palette[labels.index(Camp)], alpha = 0.5))
    fig2.legend(handles = handles, labels = labels,  ncol = 4, 
                         prop = {'size': 12}, frameon = True, numpoints=1,
                         loc = 'center', bbox_to_anchor=(0.52, 0.02))

    
def custom_boxplot(bp, boxcol='blue', boxfillcol='white', lw=1, 
                   whiskcol='black', capcol='black', mediancol='red'):
    ## change outline color, fill color and linewidth of the boxes
    for box in bp['boxes']:
        # change outline color
        box.set( color=boxcol, linewidth=lw)
        # change fill color
        box.set( facecolor = boxfillcol )
    
    ## change color and linewidth of the whiskers
    for whisker in bp['whiskers']:
        whisker.set(color=whiskcol, linewidth=lw, linestyle='-')
    
    ## change color and linewidth of the caps
    for cap in bp['caps']:
        cap.set(color=capcol, linewidth=lw)
    
    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set(color=mediancol, linewidth=lw)
    
    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='+', color='black', alpha=0.5)     