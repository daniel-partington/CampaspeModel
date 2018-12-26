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
        
def load_pest_file_and_jacobian_into_Schur(model_folder, forecasts=None):
    jco = os.path.join(model_folder, "pest_emu.jco")
    pst_file = jco.replace(".jco", ".pst")
    unc_file = jco.replace(".jco", ".unc")
    unix2dos(pst_file)
    # SETUP PST OBJECT
    pst = pyemu.Pst(pst_file)             

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
           lambda x:((x - df.loc["base", :]) / df.loc["base", :]) * 100., axis=1)             
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

def bar_chart(df_ob, ticklabels, title=True, ylabel=True, ax=None, 
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
                                      adjust=1.0, unc_type='', vlim=(0, 100),
                                      flow_low_quantile=0.25, flow_high_quantile=0.75,
                                      flow_type_colours=None,
                                      title_prefix="Percent reduction in uncertainty in SW-GW exchange: ",
                                      cbar_label='[%]'):
    
    fig = plt.figure(figsize=(2*9/2., 3*5/2.))
     
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
        if index < 4:
            ax.set_xticklabels("")
        if index in [1, 3, 5]:
            ax.set_yticklabels("")    
        else:            
            ax.set_yticklabels(heat.index, rotation=0)    
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

    plot_stream_reaches(m, ax1, x_offset=8000)

    df = df_active.loc[novel_data_sets_alias[0]]
    heat = create_spatiotemporal_unc_array(df, months, reaches) 
    heat = heat.astype(np.float64)
    heat = heat * adjust
    heat = heat.transpose()
    heat = heat.reindex(heat.index.tolist()[::-1])
    sns.heatmap(heat, vmin=vlim[0], vmax=vlim[1], cmap='viridis', ax=ax2, cbar=False, yticklabels=1)
    ax2.set_yticklabels(heat.index, rotation=0)    
    plt.title("{}{}".format(title_prefix, novel_data_sets_alias[0]))
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

def get_cost(row):
    ind = row.name
    ind = ind.split('_')
    cost = 0.
    for letter in ind:
        cost += obs_costs[letter]
    return cost

def get_cost_unc(row):
    ind = row.name
    ind = ind.split('_')
    cost_unc = 0.
    for letter in ind:
        cost_unc += obs_costs_unc[letter]
    return cost_unc

def all_combos(la, obs_groups_types, obs_groups_types_map, obs_groups_types_abbrev_map):
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
                                           ob_interest, sfr_df, v_bounds=(0, 100)):
    months = ['June 2016', 'July 2016', 'August 2016',
                 'September 2016', 'October 2016', 'November 2016', 'December 2016', 
                 'January 2017', 'February 2017', 'March 2017', 'April 2017', 'May 2017']

    fig = plt.figure(figsize=(7, len(select_months) * 3.5))

    for index, select_month in enumerate(select_months):
        month_num = 32 - 11 + select_month
        exchange = sfr_df[sfr_df['time'] == month_num].reset_index(
                    range(sfr_df.shape[0]))['Qaquifer'] / (sfr_orig['rchlen'] * sfr_seg['width2']).tolist()
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
        ax.text(-26, 15, "{}{}".format(months[select_month], extra_text[index]), 
                rotation=90, fontdict=font, horizontalalignment='center',
                verticalalignment='center')
        #fig.colorbar(ax2)
        ax.set_ylim((-2, 60))
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
                    range(sfr_df.shape[0]))['Qaquifer'] / (sfr_orig['rchlen'] * sfr_seg['width2']).tolist()
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
        sfr_df['Qaquifer_adj'] = sfr_df['Qaquifer'] / (sfr_df['rchlen'] * sfr_df['width'])
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
        ax.set_ylabel('Exchange (m/d)')
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
    cbar_ax2 = fig.add_axes([0.85, 0.663, 0.03, 0.286])
    fig.colorbar(ax3, cax=cbar_ax2)    
    cbar_ax2.get_yaxis().labelpad = 16
    cbar_ax2.set_ylabel('Reach Exchange (m${^3}$/d)', rotation=270)
    cbar_ax2.text(0.1, 0.85, 'Losing', color='white', #horizontalalignment='center',
                  transform=cbar_ax2.transAxes, rotation=270)
    cbar_ax2.text(0.1, 0.3, 'Gaining', color='white', #horizontalalignment='center',
                  transform=cbar_ax2.transAxes, rotation=270)
    
    plt.subplots_adjust(wspace=0, hspace=0.1, bottom=0.06, top=0.95, left=0.14, right=0.82)
    plt.savefig(p_j(save_folder, 'Instream_pot_obs_nrf_a_stacked.png'), dpi=300)        

def plot_SWGW_exchange(sfr_df, ax=None, show_gauges=False, fontsize=8,
                       inflow_data=None, colour_dict=None, date_index=None,
                       plot_only=None, adjust_left=0.14,  adjust_right=0.97, linewidth=0.8, bbox_to_anchor=(1.05, 1),
                       alpha=0.5):
    
    sfr_df.loc[:, 'Cumulative Length (km)'] = sfr_df['Cumulative Length'] / 1000.
    sfr_df['Qaquifer_adj'] = sfr_df['Qaquifer'] / (sfr_df['rchlen'] * sfr_df['width'])
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
    ax.set_ylabel('SW-GW exchange [m/d]', fontsize=fontsize)
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

def plot_ground_pot_obs_predunc_stream(obtype, dfs, hgus, m, axes, vlim=None):
    plot_titles = ["{} {}".format(x, obtype) for x in ['Shallow', 'Deep']]
    for index in range(2):
        hgu = hgus[index]
        sim_locs_bool = zone2D_info[0][hgu[0]] | zone2D_info[0][hgu[1]]
        ax = axes[index]
        plot_stream_reaches_basic(m, ax)
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
    
def plot_stream_reaches(m, ax, new_fig=False, x_offset=5000., c2=['red', 'blue']):
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

def plot_stream_reaches_basic(m, ax, new_fig=False, bounds=None, zones=None):
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
    
if __name__ == '__main__':

    # Setup all of the requisite data for Data Worth Analysis and pretty plotting
    p_j = os.path.join

    model_folder = r"C:\Workspace\part0075\MDB modelling\testbox\PEST1000\HPC_Feb2018"
    save_folder = p_j(model_folder, r"original")
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
    
    la = load_pest_file_and_jacobian_into_Schur(model_folder)
    forecasts_df = forecast_names(la.pst.observation_data)
    forecasts_df.loc[:, 'var'] = [la.posterior_forecast[row] for row in forecasts_df.index]
    forecasts_df.loc[:, 'unc'] = forecasts_df['var'].apply(np.sqrt)

    # Adjustment of weights based on PEST run:
    
#    observation_data = la.pst.observation_data        
#    obscov_adjust = la.obscov.to_dataframe()
#    
#    for key in new_weights_dict:
#        print key
#        observation_data.loc[observation_data['obgnme'] == key, 'weight'] = new_weights_dict[key]
#
#    for key in new_weights_perc_dict:
#        print key
#        observation_data.loc[observation_data['obgnme'] == key, 'weight'] = \
#            new_weights_perc_dict[key] * observation_data.loc[observation_data['obgnme'] == key, 'obsval']
#
#    for ob in obscov_adjust.columns:
#        obscov_adjust.loc[ob, ob] = observation_data.loc[ob, 'weight']
#
#    la.obscov = obscov_adjust
#    
#    la.obscov.to_dataframe()
    
    #la.pst.adjust_weights_resfile(resfile='pest.rei.11')
    #la.adjust_obscov_resfile(resfile='pest_emu.res')

    # Some usefule definitions of what is old, new and potential data    
    old_data = ['head1', 'head3', 'head5', 'head6', 'stage', 'gflow', 'gstrec']
    new_field_data = ['fstrec', 'radon', 'c14'] #ignore ['fflow, 'fdepth'] 
    potential_data = ['st_sim', 'fl_sim', 'rn_sim', 'ec_sim', 'shshal', 'shdeep', 'c14shal', 'c14deep']
    potential_data_stream = ['st_sim', 'fl_sim', 'rn_sim', 'ec_sim']
    potential_data_ground = ['shshal', 'shdeep', 'c14shal', 'c14deep']

    df_worth, df_perc, df_worth_added, df_perc_add = compare_adding_and_subtraction_of_obs_groups(la)
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
                       "all" 
                       ]
    
    obs_groups_al = [
                  ['head1', 'head3', 'head5', 'head6'], 
                  ['stage', 'fdepth'],
                  ['gflow', 'fflow'], 
                  ['c14'], 
                  ['radon'], 
                  ['gstrec', 'fstrec'], 
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'fdepth', 'gflow', 'fflow', 'c14', 'radon', 'gstrec', 'fstrec']]

    combo_obs_dict_al = create_combos(la, obs_groups_name_al, obs_groups_al)

    df_worth_add_combo_al, df_worth_add_combo_unc_al, df_unc_perc_add_combo_al = process_combos(la, combo_obs_dict_al, obs_groups_name_al)
    ob_interest_al = 'nrf_a'
    #box_plot_for_the_punters(630., df_worth_add_combo_unc_al[df_worth_add_combo_unc_al.index != 'base'][ob_interest_al] * 365. / 1000000.)

    # Just the bar charts
    #fig = plt.figure(figsize=(10,4))
    #ax = fig.add_subplot(1, 1, 1)
    #ticklabels_al = ("Head", "Stage", "Flow", "$^{14}$C", "$^{222}$Rn", 
    #                 "EC", "All data")
    #bar_chart(df_unc_perc_add_combo_al['nrf_a'], ticklabels_al, ax=ax)
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
    #box_plot_for_the_punters(630., df_worth_add_combo_unc[df_worth_add_combo_unc.index != 'base'][ob_interest] * 365. / 1000000.)


    #fig = plt.figure(figsize=(4.92,3))
    #ticklabels_c = ("Head", "+ Stage", "+ Flow", "Hydraulic \n+ $^{14}$C",
    #                "Hydraulic \n+ $^{222}$Rn", "Hydraulic \n+ EC", "All data")
    #bar_chart(df_unc_perc_add_combo[ob_interest], ticklabels_c,
    #          title=False, ylabel=False, save=False)
    #plt.savefig(p_j(save_folder, 'Predictive uncert cumulative hydraulic vs chem.png'), dpi=300)

    # Combined plot of individual and cumulative 
    fig = plt.figure(figsize=(4.92,5))
    ax = fig.add_subplot(2, 1, 1)
    ticklabels_al = ("Head", "Stage", "Flow", "$^{14}$C", "$^{222}$Rn", 
                     "EC", "All data")
    arxe = bar_chart(df_unc_perc_add_combo_al[ob_interest], ticklabels_al, ax=ax, 
                     save=False)
    arxe.yaxis.set_label_coords(-0.13,-0.1)
    xticklabels1 = arxe.get_xticklabels()
    arxe.set_xticklabels(xticklabels1, rotation=45)    
    #arxe.text(0.05, 0.9, 'a.', horizontalalignment='center',
    #             verticalalignment='center', transform=ax.transAxes, fontsize=12)    

    ax = fig.add_subplot(2, 1, 2)
    ticklabels_c = ("Head", "+ Stage", "+ Flow", "Hydraulic \n+ $^{14}$C",
                    "Hydraulic \n+ $^{222}$Rn", "Hydraulic \n+ EC", "All data")
    arxe2 = bar_chart(df_unc_perc_add_combo[ob_interest], ticklabels_c, ax=ax, 
                      title=False, ylabel=False, save=False)
    xticklabels = arxe2.get_xticklabels()
    arxe2.set_xticklabels(xticklabels, rotation=45)    
    #arxe2.text(0.05, 0.9, 'b.', horizontalalignment='center',
    #             verticalalignment='center', transform=ax.transAxes, fontsize=12)    
    fig.subplots_adjust(bottom=0.2, left=0.18, right=0.95, top=0.96, hspace=0.45)
    plt.savefig(p_j(save_folder, 'Predictive uncert individual vs cumulative w hydraulic vs chem.png'), dpi=300)


    # Looking at all predictions together:
    #df_unc_perc_add_combo.apply(pd.Series.describe, axis=1)
    #df_unc_perc_add_combo_al.apply(pd.Series.describe, axis=1)
    #df_unc_perc_add_combo_al['nrf_a']
    #df_unc_perc_add_combo['nrf_a']
    #df_unc_perc_add_combo_al.transpose().plot(kind='box')

    # Combined plot of individual and cumualtive 
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
        
    df_noseason_al.loc[losing.index].plot(kind='hist')
    gaining.plot(kind='hist')
    neutral.plot(kind='hist')
                            
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

#    df_rrf_al = df_unc_perc_add_combo_al[[col for col in df_unc_perc_add_combo_al.columns if 'rrf_m' in col]]
#    df_rrf = df_unc_perc_add_combo[[col for col in df_unc_perc_add_combo.columns if 'rrf_m' in col]]
#    combined_alone_and_combo_boxplots_for_select_preds(df_rrf_al, df_rrf, 'monthly_reach')
#    df_rrfa_al = df_unc_perc_add_combo_al[[col for col in df_unc_perc_add_combo_al.columns if 'rrf_a' in col]]
#    df_rrfa = df_unc_perc_add_combo[[col for col in df_unc_perc_add_combo.columns if 'rrf_a' in col]]
#    combined_alone_and_combo_boxplots_for_select_preds(df_rrfa_al, df_rrfa, 'annual_reach')
#    df_nrfm_al = df_unc_perc_add_combo_al[[col for col in df_unc_perc_add_combo_al.columns if 'nrf_m' in col]]
#    df_nrfm = df_unc_perc_add_combo[[col for col in df_unc_perc_add_combo.columns if 'nrf_m' in col]]
#    combined_alone_and_combo_boxplots_for_select_preds(df_nrfm_al, df_nrfm, 'monthly_net')

    #_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
    # Uncertatinty for all the forecast observations 
    df = la.obscov.to_dataframe()
    obscov = pd.Series(np.diag(df), index=[df.index, df.columns])

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

    #_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
    # ANALYSE SAME COMBOS ACROSS VARIOUS SPATIOTEMPORAL RANGES
    
    #df_worth_add_combo.columns = [forecasts_df['unique_name'] for x in df_worth_add_combo.columns]
    #df_worth_add_combo_unc.columns = forecasts_df['unique_name'].tolist()
    #df_unc_perc_add_combo.columns = forecasts_df['unique_name'].tolist()
    
    
    columns_name = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 'nrf']

    rows_name_orig = ['Annual', 'June', 'July', 'August',
                 'September', 'October', 'November', 'December', 'January', 'February', 'March', 'April', 'May', 
                     'Quarter1', 'Quarter2', 'Quarter3', 'Quarter4', ]

    rows_name = ['Annual', 'Quarter1', 'Quarter2', 'Quarter3', 'Quarter4', 'June', 'July', 'August',
                 'September', 'October', 'November', 'December', 'January', 'February', 'March', 'April', 'May']

    rows_name2 = ['Annual', 'Jun, 16', 'Jul, 16', 'Aug, 16',
                 'Sep, 16', 'Oct, 16', 'Nov, 16', 'Dec, 16', 'Jan, 17', 'Feb, 17', 'Mar, 17', 'Apr, 17', 'May, 17']

    columns = {}

    novel_data_sets = ['h', 'h_s'    , 'h_s_f'      , 'h_s_f_c14', 'h_s_f_radon', 'h_s_f_stream_ec', 'all']
    novel_data_sets_alias = ['Head', '+ Stage', '+ Flow', 'Hydraulic + $^{14}$C' , 'Hydraulic + $^{222}$Rn', 'Hydraulic + EC',  'All data']

    novel_data_set_to_alias = {key:val for key, val in zip(novel_data_sets, novel_data_sets_alias)}

    month_and_annual = [x for x in df_worth_add_combo_unc.columns if 's' not in x]
    df_worth_add_combo_unc_filter = df_worth_add_combo_unc[month_and_annual] 
    df_worth_add_combo_unc_filter = df_worth_add_combo_unc_filter[df_worth_add_combo_unc_filter.index != 'base']
    df_worth_add_combo_unc_filter.index = [novel_data_set_to_alias[x] for x in df_worth_add_combo_unc_filter.index]
    df_perc_add_combo_unc_filter = df_unc_perc_add_combo[month_and_annual] 
    df_perc_add_combo_unc_filter.index = [novel_data_set_to_alias[x] for x in df_perc_add_combo_unc_filter.index]
    months = rows_name2[1:]
    reaches = 11

    ###########################################################################
    #
    # 
    #
    ###########################################################################
    red_fac = 0.7
    fig = plt.figure(figsize=(9., 8 * red_fac))  
    gs = gridspec.GridSpec(2, 3,
                       width_ratios=[0.85, 3, 0.1],
                       height_ratios=[4, 1.5]
                       )
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    plot_stream_reaches(m, ax1, x_offset=8000, c2=['teal', 'black'])

    heat = pd.DataFrame(columns=["r{}".format(x) for x in range(1, reaches + 1)] + ['nrf'], index=['Annual'] + months)
    heat.loc['Annual', 'nrf'] = sfr_df_reaches_relevant.groupby(by='time').mean()['Qaquifer_lin2'].mean()
    for i in range(reaches):
        heat.loc['Annual', "r{}".format(i + 1)] = \
             sfr_df_reaches_relevant[sfr_df_reaches_relevant['reach_dw'] == i]['Qaquifer_lin2'].mean()
    # end for 
    for i, m2 in enumerate(months):
        heat.loc[m2, 'nrf'] = \
             sfr_df_reaches_relevant[sfr_df_reaches_relevant['time'] == 21 + i]['Qaquifer_lin2'].mean()
    # end for
    for i in range(reaches):
        for j, m2 in enumerate(months):
            sfr_df_reaches_relevant_reach = \
                 sfr_df_reaches_relevant[sfr_df_reaches_relevant['reach_dw'] == i]
            heat.loc[m2, "r{}".format(i + 1)] = \
                 sfr_df_reaches_relevant_reach[sfr_df_reaches_relevant_reach['time'] == 21 + j]['Qaquifer_lin2'].mean()
        # end for
    # end for
    heat = heat.astype(np.float64)
    heat = heat.transpose()
    heat = heat.reindex(heat.index.tolist()[::-1])

    biggest = max(abs(heat.values.min()), abs(heat.values.max()))
    vlim = (-biggest, biggest)

    sns.heatmap(heat, vmin=vlim[0], vmax=vlim[1], annot=True, annot_kws={"size": 8}, cmap='seismic', fmt='.2f', ax=ax2, cbar=False, yticklabels=1)
    ax2.set_yticklabels(heat.index, rotation=0)    
    ax2.tick_params(direction='out')
    plt.title("SW-GW exchange fluxes")
    ax3 = plt.subplot(gs[2])
    ax3.text(-0.7, 0.2, 'Gaining', color='blue', horizontalalignment='center',
             verticalalignment='center', transform=ax3.transAxes, fontsize=12, rotation=90)
    ax3.text(-0.7, 0.8, 'Losing', color='red', horizontalalignment='center',
             verticalalignment='center', transform=ax3.transAxes, fontsize=12, rotation=90)
    norm = mpl.colors.Normalize(vmin=vlim[0], vmax=vlim[1])
    cb1 = mpl.colorbar.ColorbarBase(ax3, cmap='seismic',
                                    norm=norm,
                                    orientation='vertical')
    cb1.set_label('[m/d]')   
    cb1.outline.set_visible(False)
    ax4 = plt.subplot(gs[4])
    
    plot_inflow_w_average(ax4, Eppalock_inflow, date_range, '30-06-2016', 
                          rows_name2, lower=flow_low_quantile, upper=flow_high_quantile,
                          colour_dict=flow_type_colours)
    fig.subplots_adjust(bottom=0.02, left=0.08, right=0.91, top=0.94, hspace=0.40, wspace=0.21)    
    plt.savefig(p_j(save_folder, 'Spatiotemporal_exchange_heatmap.png'), dpi=300)    
    
    
    
    #plot_spatiotemporal_individual(novel_data_sets_alias, df_worth_add_combo_unc_filter, 
    #                              adjust=365. * 1000. / 1000000000. , cbar_text='[Gl/yr]', unc_type='',
    #                              vlim=(0., 10.))
        
    #plot_spatiotemporal_individual(novel_data_sets_alias, df_perc_add_combo_unc_filter, 
    #                              adjust=1.0, cbar_text='[%]', unc_type='perc_red',
    #                              vlim=(0., 100.))

    plot_spatiotemporal_unc_top_plus6(novel_data_sets_alias, 
                                      df_perc_add_combo_unc_filter, 
                                      months, reaches, unc_type='',
                                      flow_low_quantile=flow_low_quantile,
                                      flow_high_quantile=flow_high_quantile,
                                      flow_type_colours=flow_type_colours
                                      )
    
    plot_spatiotemporal_unc_top_plus6(novel_data_sets_alias, 
                                      df_worth_add_combo_unc_filter,
                                      months, reaches,
                                      adjust=365. * 1000. / 1000000000., unc_type='unc',
                                      vlim=(0., 10.))
    
    #_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
    # Getting the average of SW-GW exchanges:

    plot_SWGW_exchange(sfr_df, show_gauges=True, 
                       inflow_data=Eppalock_inflow_monthly,
                       colour_dict=flow_type_colours,#{'low':'red', 'high':'blue', 'regular':'orange'},
                       date_index=date_index,
                       #plot_only='low',
                       adjust_left=0.1,
                       adjust_right=0.83,
                       linewidth=2,
                       bbox_to_anchor=(1.0, 1),
                       alpha=0.5)

    
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
    
    obs_costs = {'base':0,
                 'h':100,
                 'sh':100,
                 'dh':100,
                 'f':100,
                 's':40,
                 'e':10,
                 'r':170,
                 'c':577
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

#    #_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
#    # ANALYSE ALL COMBOS IN THE CONTEXT OF COST
#    
#    # Let's look at ALL combos for DATA additions
#
#    obs_groups_types_p = ['head', 'stage', 
#                        'flow', 'c14', 'radon', 
#                        'strec'] 
#    
#    obs_groups_types_map_p = {'head': ['head1', 'head3', 'head5', 'head6', 'shshal', 'shdeep'], 
#                            'stage': ['stage', 'st_sim'], 
#                            'flow': ['gflow', 'fl_sim'], 
#                            'c14': ['c14','c14shal', 'c14deep'], 
#                            'radon': ['radon', 'rn_sim'], 
#                            'strec': ['gstrec', 'ec_sim']}
#    
#    obs_groups_types_abbrev_map_p = {'head': 'h', 
#                                    'stage': 's', 
#                                    'flow': 'f', 
#                                    'c14': 'c', 
#                                    'radon': 'r', 
#                                    'strec': 'e'}
#    
#    obs_combos_df_p = all_combos(la, obs_groups_types_p, obs_groups_types_map_p, obs_groups_types_abbrev_map_p) 
#        
#    obs_map_p = {'h':'Head',
#                 'f':'Flow',
#                 's':'Stage',
#                 'e':'EC',
#                 'r':'$^{222}$Rn',
#                 'c':'$^{14}$C'
#                 }
#                 
#    combo_df_p = obs_combos_df_p
#    combo_df_p = combo_df_p[combo_df_p.index != 'base']
#    combo_df_p = combo_df_p.apply(np.sqrt)
#    combo_df_p.loc[:, 'cost'] = combo_df_p.apply(get_cost, axis=1)
#    combo_df_p.loc[:, 'cost_unc'] = combo_df_p.apply(get_cost_unc, axis=1)
#    
#    combo_df_lite_p = combo_df_p
#    combo_cost_unc_plot(combo_df_lite_p, 'nrf_a',
#                        fname='Bar_unc_vs_cost_{}_potential',
#                        ax_text='h=head \ns=stage \nf=flow \ne=ec \nr=radon \nc=$^{14}$C')
#    
#    # Let's look at clustering the combos by uncertainty and cost using
#    # unsupervised machine learning via K-means clustering
#    clusters = 10
#    combo_df_pred_p, combo_df_grouped_p = cluster_combos_kmeans(combo_df_p, obs_map_p, clusters=clusters, normalise=False)
#    scatter_unc_cost(combo_df_pred_p, combo_df_grouped_p, clusters=clusters, method='sklearn',
#                     title="Data worth (potential data) for SW-GW exchange in the Campaspe River",
#                     append_text="_potential", xlim=(2500, 24000))
#
#    #combo_df_pred2, combo_df_grouped2 = cluster_combos_dbscan(combo_df, eps=5000, min_samples=2)
#    #scatter_unc_cost(combo_df_pred2, combo_df_grouped2, method='DBSCAN')    
#
#    #_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
#    # ANALYSE ALL COMBOS IN THE CONTEXT OF COST
#    
#    # Let's look at ALL combos for DATA additions
#
#    obs_groups_types_h = ['shallow head', 'deep head', 'stage', 
#                        'flow', 'c14', 'radon', 
#                        'strec'] 
#    
#    obs_groups_types_map_h = {'shallow head': ['head1', 'head3'], 
#                            'deep head':  ['head5', 'head6'], 
#                            'stage': ['stage'], 
#                            'flow': ['gflow'], 
#                            'c14': ['c14'], 
#                            'radon': ['radon'], 
#                            'strec': ['gstrec']}
#    
#    obs_groups_types_abbrev_map_h = {'shallow head': 'sh', 
#                                    'deep head': 'dh',
#                                    'stage': 's', 
#                                    'flow': 'f', 
#                                    'c14': 'c', 
#                                    'radon': 'r', 
#                                    'strec': 'e'}
#    
#    obs_combos_df_h = all_combos(la, obs_groups_types_h, obs_groups_types_map_h, obs_groups_types_abbrev_map_h) 
#    
#    obs_map_h = {'sh':'Shallow Head',
#                 'dh':'Deep Head',
#                 'f':'Flow',
#                 's':'Stage',
#                 'e':'EC',
#                 'r':'$^{222}$Rn',
#                 'c':'$^{14}$C'
#                 }
#                 
#    combo_df_h = obs_combos_df_h
#    combo_df_h = combo_df_h[combo_df_h.index != 'base']
#    combo_df_h = combo_df_h.apply(np.sqrt)
#    combo_df_h.loc[:, 'cost'] = combo_df_h.apply(get_cost, axis=1)
#    combo_df_h.loc[:, 'cost_unc'] = combo_df_h.apply(get_cost_unc, axis=1)
#    
#    combo_df_lite_h = combo_df_h
#    combo_cost_unc_plot(combo_df_lite_h, 'nrf_a',
#                        fname='Bar_unc_vs_cost_{}_head_split',
#                        ax_text='sh=shallow head \ndh=deep head \ns=stage \nf=flow \ne=ec \nr=radon \nc=$^{14}$C')
#    
#    # Let's look at clustering the combos by uncertainty and cost using
#    # unsupervised machine learning via K-means clustering
#    clusters = 14
#    combo_df_pred_h, combo_df_grouped_h = cluster_combos_kmeans(combo_df_h, obs_map_h, clusters=clusters, normalise=True)
#    scatter_unc_cost(combo_df_pred_h, combo_df_grouped_h, clusters=clusters, method='sklearn',
#                     title="Data worth (existing data) for SW-GW exchange in the Campaspe River",
#                     append_text="_head_split", xlim=(5000, 25000))

    #_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
    # Potential Observations

    obs_groups_name2 = []
    for ob_gp in obs_groups_name:
        obs_groups_name2 += [ob_gp, ob_gp + '_new']

    obs_groups2 = [
                  ['head1', 'head3', 'head5', 'head6'], 
                  ['head1', 'head3', 'head5', 'head6', 'shshal', 'shdeep'], 
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'fdepth'],
                  ['head1', 'head3', 'head5', 'head6', 'shshal', 'shdeep', 'stage', 'fdepth', 'st_sim'],
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'fdepth', 'gflow', 'fflow'], 
                  ['head1', 'head3', 'head5', 'head6', 'shshal', 'shdeep', 'stage', 'fdepth', 'st_sim', 'gflow', 'fflow', 'fl_sim'], 
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'fdepth', 'gflow', 'fflow', 'c14'], 
                  ['head1', 'head3', 'head5', 'head6', 'shshal', 'shdeep', 'stage', 'fdepth', 'st_sim', 'gflow', 'fflow', 'fl_sim', 'c14', 'c14shal', 'c14deep'], 
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'fdepth', 'gflow', 'fflow', 'radon'], 
                  ['head1', 'head3', 'head5', 'head6', 'shshal', 'shdeep', 'stage', 'fdepth', 'st_sim', 'gflow', 'fflow', 'fl_sim', 'radon', 'rn_sim'], 
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'fdepth', 'gflow', 'fflow', 'gstrec', 'fstrec'], 
                  ['head1', 'head3', 'head5', 'head6', 'shshal', 'shdeep', 'stage', 'fdepth', 'st_sim', 'gflow', 'fflow', 'fl_sim', 'gstrec', 'fstrec', 'ec_sim'], 
                  ['head1', 'head3', 'head5', 'head6', 'stage', 'fdepth', 'gflow', 'fflow', 'c14', 'radon', 'gstrec', 'fstrec'],
                  ['head1', 'head3', 'head5', 'head6', 'shshal', 'shdeep', 'stage', 'fdepth', 'st_sim', 'gflow', 'fflow', 'fl_sim', 'c14', 'c14shal', 'c14deep', 'radon', 'rn_sim', 'gstrec', 'fstrec', 'ec_sim']]

    obs_groups3 = [
                  ['head1', 'head3', 'shshal'],
                  ['head5', 'head6'],
                  ['head5', 'head6', 'shdeep'],
                  ['c14'],
                  ['c14', 'c14shal', 'c14deep'],
                  ['stage'],
                  ['stage', 'st_sim'],
                  ['gflow'],
                  ['gflow', 'fl_sim'],
                  ['radon'],
                  ['radon', 'rn_sim'],
                  ['gstrec', 'fstrec'],
                  ['gstrec', 'fstrec', 'ec_sim'],
                  old_data + new_field_data,
                  old_data + new_field_data + potential_data]

    combo_obs_dict2 = create_combos(la, obs_groups_name2, obs_groups2)
    df_worth_add_combo2, df_worth_add_combo_unc2, df_unc_perc_add_combo2 = process_combos(la, combo_obs_dict2, obs_groups_name2)
    
    pst_df = la.pst.observation_data
 
    fig = plt.figure(figsize=(4.92,3))
   
    df_ob = df_unc_perc_add_combo2['nrf_a']        
    existing = df_ob[~df_ob.index.str.contains('new')].tolist() 
    new = df_ob[df_ob.index.str.contains('new')].tolist()
    
    ind = np.arange(len(existing))  # the x locations for the groups
    width = 0.35  # the width of the bars

    ax = fig.add_subplot(111)
    ax.text(0.2, 1.05, 'Hydraulic', horizontalalignment='center',
         verticalalignment='center', transform=ax.transAxes, fontsize=10)
    ax.text(0.64, 1.05, 'Chemical', horizontalalignment='center',
         verticalalignment='center', transform=ax.transAxes, fontsize=10)
    #ax.text(0.93, 1.05, 'Combined', horizontalalignment='center',
    #     verticalalignment='center', transform=ax.transAxes, fontsize=10)
    ax.fill_between([2.68, 5.65], [0, 0], [100, 100], color='lightgrey', alpha=0.6) #, interpolate=True)
    ax.plot([5.65, 5.65], [0, 100], color='grey', linestyle='--')
    rects1 = ax.bar(ind - width/2, existing, width, 
                    color='IndianRed', label='Existing', edgecolor='none')
    rects2 = ax.bar(ind + width/2, new, width, 
                    color='SkyBlue', label='Potential', edgecolor='none')
    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('% reduction in predictive\n uncertainty')
    ax.set_ylim(0, 100)        
    ax.set_xlim(0 - width, len(existing) - width)        
    #ax.set_title('Observation types reduction in uncertainty')
    ax.set_xticks([x + width/2 for x in ind])
    comp_ticklabels = ("Head", "+ Stage", "+ Flow", "Hydraulic \n+ $^{14}$C",
                        "Hydraulic \n+ $^{222}$Rn", "Hydraulic \n+ EC", "All data")
    ax.set_xticklabels(comp_ticklabels, rotation=45)

    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    #ax.legend(loc='top center', fontsize=12)  
    ax.text(0.03, 0.37, 'Existing', transform=ax.transAxes, fontsize=10, rotation=90)#, color='white')
    ax.text(0.084, 0.40, 'Potential', transform=ax.transAxes, fontsize=10, rotation=90)

    fig.subplots_adjust(bottom=0.28, left=0.18, right=0.95, top=0.92, hspace=0.45)
    plt.savefig(p_j(save_folder, 'Comp_existing_potential_combos_nrf_a.png'), dpi=300)

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
    comp_ticklabels = ("Head", "+ Stage", "+ Flow", "Hydraulic \n+ $^{14}$C",
                        "Hydraulic \n+ $^{222}$Rn", "Hydraulic \n+ EC", "All data")
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

    novel_data_set_to_alias2 = {'all_new': 'All data',
                                 'h_new': 'Head',
                                 'h_s_new': '+ Stage',
                                 'h_s_f_new': '+ Flow',
                                 'h_s_f_c14_new': 'Hydraulic + $^{14}$C',
                                 'h_s_f_radon_new': 'Hydraulic + $^{222}$Rn',
                                 'h_s_f_stream_ec_new': 'Hydraulic + EC'}

    novel_data_set_to_alias3 = {'all_diff': 'All data',
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
    df_perc_add_combo_unc_filter2 = df_unc_perc_add_combo2[month_and_annual] 
    df_perc_add_combo_unc_filter2 = df_perc_add_combo_unc_filter2.loc[[col for col in df_perc_add_combo_unc_filter2.index if 'new' in col], :]    
    df_perc_add_combo_unc_filter2.index = [novel_data_set_to_alias2[x] for x in df_perc_add_combo_unc_filter2.index]

    plot_spatiotemporal_unc_top_plus6(novel_data_sets_alias, 
                                      df_perc_add_combo_unc_filter2, 
                                      months, reaches, unc_type='pot',
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
                                      months, reaches, unc_type='pot_diff',
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

    if os.path.exists('All_potential_data_stream.csv') and not force_recalc:
        stream_pot_obs = pd.read_csv('All_potential_data_stream.csv', index_col=0)
    else:
        stream_pot_obs = la.get_added_obs_importance(
            {ob:[ob] for ob in pst_df[pst_df['obgnme'].isin(potential_data_stream)].index.tolist()})
        stream_pot_obs.to_csv('All_potential_data_stream.csv')
    
    stream_pot_obs_unc = stream_pot_obs[~stream_pot_obs.index.isin(['base'])].apply(
        lambda x:(1 - x / stream_pot_obs.loc["base", :]) * 100., axis=1)

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
    ax.set_ylabel('Exisiting data')
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
    mega_plot_river_predunc_red_month_axes(select_months, stream_pot_obs_unc, pst_df, 
                                           potential_data_stream, 'ob3454', sfr_df)    
   
    if os.path.exists('All_potential_data_ground.csv') and not force_recalc:
        ground_pot_obs = pd.read_csv('All_potential_data_ground.csv', index_col=0)
    else:
        ground_pot_obs = la.get_added_obs_importance(
            {ob:[ob] for ob in pst_df[pst_df['obgnme'].isin(potential_data_ground)].index.tolist()})
        ground_pot_obs.to_csv('All_potential_data_ground.csv')
    # end if
    
    ground_pot_obs_unc = ground_pot_obs[~ground_pot_obs.index.isin(['base'])].apply(
        lambda x:(1 - x / ground_pot_obs.loc["base", :]) * 100., axis=1)

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
            plot_stream_reaches_basic(m, ax, zones=hgu)            
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
        plot_stream_reaches_basic(m, ax, bounds=m.model_boundary)
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
        plot_stream_reaches_basic(m, ax, bounds=m.model_boundary)
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

    #existing_obs_groups = ['head1', 'head3', 'head5', 'head6', 'stage', 'fdepth', 'gflow', 'fflow', 'c14', 'radon', 'gstrec', 'fstrec']
    #next_most = la.next_most_important_added_obs(
    #    forecast=pst_df[pst_df['obgnme'] == 'nrf_a']['obsnme'].tolist()[0], 
    #    base_obslist=pst_df[pst_df['obgnme'].isin(existing_obs_groups)]['obsnme'].tolist())
