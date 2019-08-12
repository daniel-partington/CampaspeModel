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
    os.chdir(model_folder)    

    if not os.path.exists(save_folder):
        os.mkdir(save_folder)
    # end if

    # Some useful bits from the simulation outputs

    field_sampling = [datetime.datetime(2016,0o3,31),
                  datetime.datetime(2016,12,31),
                  datetime.datetime(2017,0o4,30)]

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

    months = rows_name2[1:]
    reaches = 11

    ###########################################################################
    #
    # 
    #
    ###########################################################################
    red_fac = 0.7
    col_to_use = 'Qaquifer_lin'
    fig = plt.figure(figsize=(9., 8 * red_fac))  
    gs = gridspec.GridSpec(2, 3,
                       width_ratios=[0.85, 3, 0.1],
                       height_ratios=[4, 1.5]
                       )
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    plot_stream_reaches(m, ax1, zone2D_info, x_offset=8000, c2=['teal', 'black'])

    heat = pd.DataFrame(columns=["r{}".format(x) for x in range(1, reaches + 1)] + ['nrf'], index=['Annual'] + months)
    heat.loc['Annual', 'nrf'] = sfr_df_reaches_relevant.groupby(by='time').mean()[col_to_use].mean()
    for i in range(reaches):
        heat.loc['Annual', "r{}".format(i + 1)] = \
             sfr_df_reaches_relevant[sfr_df_reaches_relevant['reach_dw'] == i][col_to_use].mean()
    # end for 
    for i, m2 in enumerate(months):
        heat.loc[m2, 'nrf'] = \
             sfr_df_reaches_relevant[sfr_df_reaches_relevant['time'] == 21 + i][col_to_use].mean()
    # end for
    for i in range(reaches):
        for j, m2 in enumerate(months):
            sfr_df_reaches_relevant_reach = \
                 sfr_df_reaches_relevant[sfr_df_reaches_relevant['reach_dw'] == i]
            heat.loc[m2, "r{}".format(i + 1)] = \
                 sfr_df_reaches_relevant_reach[sfr_df_reaches_relevant_reach['time'] == 21 + j][col_to_use].mean()
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
    if col_to_use == 'Qaquifer_lin2':
        cb1.set_label('[m/d]')   
    if col_to_use == 'Qaquifer_lin':
        cb1.set_label('[m$^2$/d]')   

    cb1.outline.set_visible(False)
    ax4 = plt.subplot(gs[4])
    
    plot_inflow_w_average(ax4, Eppalock_inflow, date_range, '30-06-2016', 
                          rows_name2, lower=flow_low_quantile, upper=flow_high_quantile,
                          colour_dict=flow_type_colours)
    fig.subplots_adjust(bottom=0.02, left=0.08, right=0.91, top=0.94, hspace=0.40, wspace=0.21)    
    plt.savefig(p_j(save_folder, 'Spatiotemporal_exchange_heatmap.png'), dpi=300)    
    

    #_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
    # Getting the average of SW-GW exchanges:

    plot_SWGW_exchange(sfr_df, save_folder, show_gauges=True, 
                       inflow_data=Eppalock_inflow_monthly,
                       colour_dict=flow_type_colours,
                       date_index=date_index,
                       adjust_left=0.1,
                       adjust_right=0.83,
                       linewidth=2,
                       bbox_to_anchor=(1.0, 1),
                       alpha=0.5)

