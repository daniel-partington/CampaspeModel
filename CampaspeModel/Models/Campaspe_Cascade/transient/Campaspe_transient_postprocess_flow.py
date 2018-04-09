import sys
import os
sys.path.append('C:\Workspace\part0075\GIT_REPOS')

import numpy as np
import flopy.utils.binaryfile as bf

from HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface

# Configuration Loader
from HydroModelBuilder.Utilities.Config.ConfigLoader import ConfigLoader


def run(model_folder, data_folder, mf_exe, param_file="", verbose=True):
    
    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder, "02_transient_flow_packaged.pkl"))
    name = MM.GW_build.keys()[0]
    m = MM.GW_build[name]
    
    # Load in the new parameters based on parameters.txt or dictionary of new parameters
 

    if verbose:
        print "************************************************************************"
        print " Build and run MODFLOW model "
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ## Currently using flopyInterface directly rather than running from the ModelManager ...

    modflow_model = flopyInterface.ModflowModel(m, data_folder=data_folder)

    modflow_model.executable = mf_exe_folder

    modflow_model.buildMODFLOW(transport=True, write=False)

    #modflow_model.compareAllObs()
    
    #modflow_model.viewHeadsByZone()
   
    #plots = {}
    
    #plots['wb'] = modflow_model.waterBalanceTS(plot=True)
       
    #modflow_model.viewHeads()
    #modflow_model.viewHeads2()
       
    #return
    #modflow_model.sfr
    
    
    from flopy.utils import sfroutputfile
    sfr = sfroutputfile.SfrFile(os.path.join(modflow_model.data_folder, modflow_model.sfr.file_name[0] + '.out'))
    
    date_index = m.model_time.t['dateindex']
    sfr_df = sfr.get_dataframe()
    fn = lambda x: x.iloc[0]
    sfr_df_group = sfr_df.groupby(by='time').agg({'Qaquifer':np.sum,'Qin':fn, 'Qovr':np.sum, 'Qprecip':np.sum, 'Qet':np.sum})
    #sfr_df_group.plot(x=date_index[1:], y=['Qaquifer', 'Qin', 'Qovr', 'Qprecip', 'Qet'])
    sfr_info = m.river_mapping['Campaspe']
    cum_len = sfr_info['Cumulative Length'].tolist()
    sfr_df.loc[:, 'Cumulative Length'] = cum_len * (sfr_df['time'].max() + 1)

#    modflow_model.    
    
#    for i in range(sfr_df['time'].max()):
#        ax = sfr_df[sfr_df['time'] == i].plot(x='Cumulative Length', y=[x for x in sfr_df.columns if 'Q' in x])    
#        ax.set_ylim(-80000, 400000)

    num_reaches = m.pilot_points['Campaspe'].num_points #4
    known_points = m.pilot_points['Campaspe'].points
    
    # Hyporheic zone depth         
    hz_depth_vals = [m.parameters.param['hz_dpth{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    R_depth_HZ = np.interp(sfr_info['Cumulative Length'].tolist(), known_points, hz_depth_vals)

    df_size = sfr_info.shape[0]
    t_steps = sfr_df['time'].max() + 1
    # Hyporheic zone porosity
    sfr_df.loc[:, 'HZ_poro'] = [m.parameters.param['hz_poro']['PARVAL1']] * df_size * t_steps
    # Hyporheic zone production of radon
    sfr_df.loc[:, 'HZ_Prod_Rate'] = [m.parameters.param['hz_prod']['PARVAL1']] * df_size * t_steps
    # Hyporheic zone residence time
    sfr_df.loc[:, 'HZ_RTime'] = [m.parameters.param['hz_rt']['PARVAL1']] * df_size * t_steps
    # Depth of the hyporheic zone
    sfr_df.loc[:, 'R_depth_HZ'] = R_depth_HZ.tolist() * t_steps              
    # Gas transfer velocity
    sfr_df.loc[:, 'GTV'] = [m.parameters.param['gtv']['PARVAL1']] * df_size * t_steps
    # Groundwater radon concentration
    sfr_df.loc[:, 'GW_Rn_conc'] = [m.parameters.param['gw_conc']['PARVAL1']] * df_size * t_steps
    # Groundwater EC
    sfr_df.loc[:, 'GW_EC'] = [10.] * df_size * t_steps # ARBITRARY!!!!
    # EC of the inflowing tributary water if present
    sfr_df.loc[:, 'Trib_EC'] = [10.] * df_size * t_steps # ARBITRARY!!!! 
    # Radon concentration of inflowing tributary water if present
    sfr_df.loc[:, 'Trib_Rn'] = [10.] * df_size * t_steps # ARBITRARY!!!!
    # Reach lengths
    sfr_df.loc[:, 'dx'] = sfr_info['rchlen'].tolist() * t_steps 

    df_list = []
    # March 2016, Dec 2016, May 2017
    import datetime
    import matplotlib.pyplot as plt
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    field_sampling = [datetime.datetime(2016,03,31),
                      datetime.datetime(2016,12,31),
                      datetime.datetime(2017,04,30)]
    for i in range(0,33): #sfr_df['time'].max()):
        if date_index[i] in field_sampling:    
            df = sfr_df[sfr_df['time'] == i]
            #df.to_csv(r"C:\Workspace\part0075\badRn.csv")
            Ini_cond = (df.iloc[0]['Qin'], 0., 300.)
            df.loc[:, 'Flow'], df.loc[:, 'Rn'], df.loc[:, 'EC'] = modflow_model.Calculate_Rn_from_SFR_with_simple_model(df, Ini_cond)
            df['Qaquifer_adj'] = df['Qaquifer'] / (df['dx'] * df['width'])
            df_list += [df[['Cumulative Length', 'Flow', 'Rn', 'EC']]]
            df.plot(x='Cumulative Length', y='Rn', style='-', ax=ax, label=date_index[i].date())  
            df.plot(x='Cumulative Length', style='o',  y=['Flow', 'Qin'], ax=ax2, label=date_index[i].date())  
            df[df['Qaquifer_adj'] != 0.].plot(x='Cumulative Length', style='o',  y=['Qaquifer_adj'], ax=ax2, secondary_y=True, label=date_index[i].date())  
            
            #df.plot(x='Cumulative Length', secondary_y=True, style='o',  y=['Flow', 'Qaquifer'], ax=ax)  
            #ax.set_title("Radon (mBq/L)")
            #ax.set_ylim(0, 0.2)
    ax.set_ylabel("Radon (mBq/L)")
    ax.set_xlabel("River chainage (m)")


#    import matplotlib.animation
#    
#    x = df_list
#    fig, ax = plt.subplots()
#    line = ax.plot(x=df_list[0]['Cumulative Length'].tolist(), y=df_list[0]['Rn'].tolist())
#    #df_list[0].plot(x='Cumulative Length', y='Flow', secondary_y=True, ax=ax)
#
#    def update(num, x, y, line):
#        #ax.cla()
#        line.set_data(x=df_list[n]['Cumulative Length'].tolist()[:num], y=df_list[n]['Rn'].tolist()[:num])
#        #df_list[n].plot(x='Cumulative Length', y='Rn', ax=ax)
#        #df_list[n].plot(x='Cumulative Length', y='Flow', secondary_y=True, ax=ax)
#        return line
#    
#    ani = matplotlib.animation.FuncAnimation(fig, update, range(len(df_list)), blit=True)
#    
#    ani.save(r"C:\Workspace\part0075\test.gif")
#    
#    plt.show()        

#
#    import matplotlib.pyplot as plt
#    import matplotlib.animation as animation
#    
#    #npdata = numpy.random.randint(100, size=(5,6,10))
#    plotlays, plotcols = [x for x in sfr_df.columns if 'Q' in x], "bgroyi"
#    
#    fig = plt.figure()
#    ax = plt.axes(xlim=(cum_len[0], cum_len[-1]), ylim=(-40000, 140000))
#    timetext = ax.text(0.5,50,'')
#    
#    lines = []
#    for index,lay in enumerate(plotlays):
#        lobj = ax.plot([],[],lw=2,color=plotcols[index])[0]
#        lines.append(lobj)
#    
#    def init():
#        for line in lines:
#            line.set_data([],[])
#        return lines
#    
#    def animate(i):
#        timetext.set_text(i)
#        x = cum_len
#        for lnum,line in enumerate(lines):
#            line.set_data(x, sfr_df[:,plotlays[lnum]-1,i])
#        return tuple(lines) + (timetext,)
#    
#    anim = animation.FuncAnimation(fig, animate, init_func=init,
#                                   frames=range(sfr_df['time'].max() + 1), interval=100, blit=True)
#    
#    plt.show()

    
    return sfr_df
    return modflow_model, plots

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
