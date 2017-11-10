"""Update parameters and run transient flow model for Campaspe."""

import os
import sys

import numpy as np
import flopy
import flopy.utils.binaryfile as bf

#sys.path.append('C:\Workspace\part0075\GIT_REPOS')
from HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface

## Configuration Loader
#from HydroModelBuilder.Utilities.Config.ConfigLoader import CONFIG

# Configuration Loader
from HydroModelBuilder.Utilities.Config.ConfigLoader import ConfigLoader

def run(model_folder, data_folder, mt_exe_folder, param_file=None, verbose=True):
    """Model Runner."""

    # MM is short for model manager
    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder, "02_transient_flow_packaged.pkl"))

    name = MM.GW_build.keys()[0]
    m = MM.GW_build[name]

    # Load in the new parameters based on parameters.txt or dictionary of new parameters

    if param_file:
        MM.GW_build[name].updateModelParameters(os.path.join(data_folder, 'parameters.txt'), verbose=verbose)

    if verbose:
        print "************************************************************************"
        print " Instantiate MODFLOW model "

    ###########################################################################
    ###########################################################################
    ###########################################################################
    # Currently using flopyInterface directly rather than running from the ModelManager ...
    modflow_model = flopyInterface.ModflowModel(MM.GW_build[name], 
                                                data_folder=os.path.join(data_folder, "model_" + m.name))
    modflow_model.buildMODFLOW(transport=True, write=False, verbose=False)

    if verbose:
        print "************************************************************************"
        print " Set initial conc from ss solution "

    species = ['EC', 'C14']

    path=os.path.join(data_folder, "model_01_steady_state")
    concobj = bf.UcnFile(os.path.join(path, 'MT3D001.UCN'))
    times = concobj.get_times()        
    C14_init = concobj.get_data(totim=times[-1])
    EC_init = np.ones_like(C14_init)
    EC_init_simple = 2000.
    EC_init = EC_init * EC_init_simple


    ibound = modflow_model.model_data.model_mesh3D[1]        
    ibound[ibound == -1] = 0

    # This needs to be automatically generated with the map_raster2mesh routine ...
    zone_map = {1: 'qa', 2: 'utb', 3: 'utqa', 4: 'utam', 
                5: 'utaf', 6: 'lta', 7: 'bse'}

    prsity = np.zeros_like(ibound, dtype=float) 
    for zone in zone_map:          
        prsity[ibound == zone] = m.parameters.param['por_{}'\
                                  .format(zone_map[zone])]['PARVAL1']
    # Dicts to take care of type of transport sim ... multi not working at the moment
    sconc_dict         = {'C14':C14_init, 'EC':500.}
    rech_conc_dict     = {'C14':100.,     'EC':EC_init_simple}
    ghb_conc_dict      = {'C14':0.,       'EC':EC_init_simple}
    riv_conc_dict      = {'C14':100.,     'EC':EC_init_simple}
    precip_conc_dict   = {'C14':100.,     'EC':0.0}
    str_conc_init_dict = {'C14':100.,     'EC':0.0}

    for specimen in species:
        if verbose:
            print "************************************************************************"
            print " Instantiate MT3D model "
            
        mt = flopy.mt3d.Mt3dms(modelname=modflow_model.name + \
                                        '_transport_{}'.format(specimen),
                               ftlfilename='mt3d_link.ftl',
                               ftlfree = True,
                               modflowmodel=modflow_model.mf, 
                               model_ws=path, #modflow_model.data_folder, 
                               exe_name=mt_exe_folder) #'MT3D-USGS_64.exe')
        
        if verbose:
            print "************************************************************************"
            print " Setup MT3D packages "
    
        if verbose:
            print("Add the BTN package to the model")

        flopy.mt3d.Mt3dBtn(mt, DRYCell=True, OmitDryBud=False, icbund=ibound, 
                           prsity=prsity, ncomp=1, mcomp=1, 
                           cinact=-9.9E1, thkmin=-1.0E-6, 
                           ifmtcn=5, ifmtnp=0, ifmtrf=0, ifmtdp=0, nprs=0, 
                           timprs=None, savucn=1, nprobs=0,  laycon=1,
                           chkmas=1, nprmas=1, dt0=10000.0, ttsmax=100000.0,
                           sconc=sconc_dict[specimen], 
                           species_names=[specimen])
    
        if verbose:
            print("Add the ADV package to the model")
    
        flopy.mt3d.Mt3dAdv(mt, mixelm=0, percel=1, 
                           mxpart=250000, nadvfd=1, itrack=3,
                           wd=0.5, dceps=1.0E-4, nplane=0, npl=5,
                           nph=8, npmin=1, npmax=16, nlsink=0, 
                           npsink=8, dchmoc=0.01)
    
        if verbose:
            print("Add the DSP package to the model")
    
        al = m.parameters.param['disp']['PARVAL1']
        flopy.mt3d.Mt3dDsp(mt, multiDiff=True, al=al, trpt=0.1, trpv=0.1, 
                           dmcoef=0.0)
    
        if specimen == 'C14':
            if verbose:
                print("Add the RCT package to the model")
            flopy.mt3d.Mt3dRct(mt, isothm=1, ireact=1, igetsc=0, 
                               rc1=np.log(2)/(5730.0*365.0))
    
        if verbose:
            print("Add the GCG package to the model")
        
        flopy.mt3d.Mt3dGcg(mt, mxiter=1000, iter1=100, isolve=1, 
                           ncrs=0, accl=1, cclose=1.0E-4, iprgcg=0)
    
        ssm_data = {}
        itype = flopy.mt3d.Mt3dSsm.itype_dict()
        for per in range(modflow_model.nper):
            ssm_data[per] = []
        ibound = modflow_model.model_data.model_mesh3D[1]        
        ibound[ibound == -1] = 0
        
        crch = {}
        crch2 = {}
        for per in range(modflow_model.nper):
            crch[per] = []
            crch2[per] = []
        # End for
    
        bc = modflow_model.model_data.boundaries.bc
        for boundary in bc:
            bc_boundary = bc[boundary]
            bc_type = bc_boundary['bc_type']
            bc_array = bc_boundary['bc_array']
    
            if bc_type == 'recharge':
                for key in bc_array.keys():
                    crch[key] = np.ones_like(bc_array[key])               
                    crch[key] = crch[key] * rech_conc_dict[specimen]
                    crch[key][ibound[0]==0] = 0.0
                
            if bc_type == 'wells':
                for key in bc_array.keys():
                    for well in bc_array[key]:
                        ssm_data[key].append((well[0], well[1], well[2], 
                                              0.0, itype['WEL']))
                max_sp = max(bc_array.keys())
                if max_sp < modflow_model.nper - 1:
                    for per in range(max_sp + 1, modflow_model.nper):
                        for well in bc_array[max_sp]:
                            ssm_data[per].append((well[0], well[1], well[2], 
                                                  0.0, itype['WEL']))
    
            if bc_type == 'drain':
                for key in bc_array.keys():
                    for drain in bc_array[key]:
                        ssm_data[key].append((drain[0], drain[1], drain[2], 
                                              0.0, itype['DRN']))
                max_sp = max(bc_array.keys())
                if max_sp < modflow_model.nper - 1:
                    for per in range(max_sp + 1, modflow_model.nper):
                        for drain in bc_array[max_sp]:
                            ssm_data[per].append((drain[0], drain[1], drain[2], 
                                                  0.0, itype['DRN']))
    
            if bc_type == 'general head':
                for key in bc_array.keys():
                    for ghb in bc_array[key]:
                        ssm_data[key].append((ghb[0], ghb[1], ghb[2], 
                                              ghb_conc_dict[specimen], itype['GHB']))

                max_sp = max(bc_array.keys())
                if max_sp < modflow_model.nper - 1:
                    for per in range(max_sp + 1, modflow_model.nper):
                        for ghb in bc_array[max_sp]:
                            ssm_data[per].append((ghb[0], ghb[1], ghb[2], 
                                                  ghb_conc_dict[specimen], itype['GHB']))
    
        river_exists = False
        river_flow_exists = False
    
        for boundary in bc:
            bc_boundary = bc[boundary]
            bc_type = bc_boundary['bc_type']
            bc_array = bc_boundary['bc_array']
            if bc_type == 'river':
                river_exists = True
            if bc_type == 'river_flow':
                river_flow_exists = True
                stream_flow_bc = bc_array
        
        if river_exists:
            river = {}
            river[0] = []
            for boundary in bc:
                bc_boundary = bc[boundary]
                bc_type = bc_boundary['bc_type']
                bc_array = bc_boundary['bc_array']
                if bc_type == 'river':
                    time_key = bc_array.keys()[0]
                    river[0] += bc_array[time_key]
                if bc_type == 'channel':
                    time_key = bc_array.keys()[0]
                    river[0] += bc_array[time_key]
            
            for key in river.keys():
                for riv in river[key]:
                    ssm_data[key].append((riv[0], riv[1], riv[2], 
                                          riv_conc_dict[specimen], itype['RIV'])) 
                max_sp = max(river.keys())
                if max_sp < modflow_model.nper - 1:
                    for per in range(max_sp + 1, modflow_model.nper):
                        for riv in river[max_sp]:
                            ssm_data[per].append((riv[0], riv[1], riv[2], 
                                                  riv_conc_dict[specimen], itype['RIV']))
    
        if river_flow_exists:
            if verbose:
                print("Add the SFT package to the model")
    
            nsfinit = stream_flow_bc[0].shape[0]
            seg_len = np.unique(stream_flow_bc[0]['iseg'], return_counts=True)
            obs_sf = np.cumsum(seg_len[1])
            obs_sf = obs_sf.tolist()
            Lake_Eppalock_EC = m.boundaries.bc['Eppalock_EC']['bc_array']
            sf_stress_period_data = {}
            for per in range(modflow_model.nper):
                sf_stress_period_data[per] = []
    
            for per in range(m.model_time.t['steps']):
                if specimen == 'EC':
                    sf_stress_period_data[per].append((0, 0, Lake_Eppalock_EC.iloc[per]['Mean']))
                elif specimen == 'C14':
                    sf_stress_period_data[per].append((0, 0, 100.))
                # end if
                for seg in range(nsfinit):
                    # For precipitation
                    sf_stress_period_data[per].append((seg, 1, precip_conc_dict[specimen]))
                    # For evaporation
                    # ... forged aboud it ... can't even handle! But it would be:
                    # sf_stress_period_data[per].append((seg, 5, 0.0))
                
            gage_output = None
            dispsf = m.parameters.param['sfdisp']['PARVAL1']        
            
            flopy.mt3d.Mt3dSft(mt, 
                               nsfinit=nsfinit, # Number of simulated stream reaches in SFR 
                               mxsfbc=nsfinit*3, # Maximum number of stream boundary conditions
                               icbcsf=81, # Integer directing to write reach-by-reach conc info to unit specified by integer 
                               ioutobs=82, # Unit number for output for concs at specified gage locations
                               ietsfr=1, #Specifies that mass is not removed with ET
                               wimp=0.5, # Stream solver time weighting factor (between 0.0 and 1.0)
                               wups=1.0, # Space weighting factor employed in the stream network solver (between 0 and 1)
                               isfsolv=1, # This is the default and only supported value at this stage
                               cclosesf=1.0E-6, # This is the closure criterion for the SFT solver 
                               mxitersf=10, # Maximum number of iterations for the SFT solver
                               crntsf=1.0,  # Courant constraint specific to SFT
                               iprtxmd=0, # flag to print SFT solution info to standard output file (0 = print nothing)
                               coldsf=str_conc_init_dict[specimen], # Initial concentration in the stream network (can also be specified as an array for each reach)
                               #coldsf2=np.array([EC_init_simple] * nsfinit), # Initial concentration in the stream network (can also be specified as an array for each reach)
                               #coldsf= 0.0, #np.array([EC_init_simple] * nsfinit), # Initial concentration in the stream network (can also be specified as an array for each reach)
                               dispsf=dispsf, # Dispersivity in each stream reach
                               #dispsf2=dispsf, # Dispersivity in each stream reach
                               nobssf=len(obs_sf), #, 
                               obs_sf=obs_sf, 
                               sf_stress_period_data=sf_stress_period_data, # Dict for each stress period with each dict containing ISFNBC, ISFBCTYP, (CBCSF(n), n=1, NCOMP)
                               filenames=gage_output)
    
                   
        ssm_data_temp = {}
        for key in ssm_data.keys():
            if len(ssm_data[key]) > 0:
                ssm_data_temp[key] = ssm_data[key]
            #end if
        #end for
        ssm_data = ssm_data_temp
    
        if verbose:
            print("Add the SSM package to the model")
    
        flopy.mt3d.Mt3dSsm(mt, 
                           stress_period_data=ssm_data, 
                           crch=crch, 
                           #crch2=crch2
                           )
        
        mt.write_input()
        
        success, buff = mt.run_model(silent=True)
        
        post = flopyInterface.MT3DPostProcess(modflow_model, 
                                              mt_name=mt.name)
        
        post.writeObservations(specimen)

    # end for

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
    hz_depth_vals = [m.parameters.param['hz_dpth{}'.format(x)]['PARVAL1'] for \
                                        x in range(num_reaches)] 
    R_depth_HZ = np.interp(sfr_info['Cumulative Length'].tolist(), 
                           known_points, 
                           hz_depth_vals)

    df_size = sfr_info.shape[0]
    t_steps = sfr_df['time'].max() + 1
    # Hyporheic zone porosity
    sfr_df.loc[:, 'HZ_poro'] = [m.parameters.param['hz_poro']['PARVAL1']] \
                               * df_size * t_steps
    # Hyporheic zone production of radon
    sfr_df.loc[:, 'HZ_Prod_Rate'] = [m.parameters.param['hz_prod']['PARVAL1']] \
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
    sfr_df.loc[:, 'GW_Rn_conc'] = [m.parameters.param['gw_conc']['PARVAL1']] \
                                  * df_size * t_steps
    # Groundwater EC
    sfr_df.loc[:, 'GW_EC'] = [10.] * df_size * t_steps # ARBITRARY!!!!
    # EC of the inflowing tributary water if present
    sfr_df.loc[:, 'Trib_EC'] = [10.] * df_size * t_steps # ARBITRARY!!!! 
    # Radon concentration of inflowing tributary water if present
    sfr_df.loc[:, 'Trib_Rn'] = [10.] * df_size * t_steps # ARBITRARY!!!!
    # Reach lengths
    sfr_df.loc[:, 'dx'] = sfr_info['rchlen'].tolist() * t_steps 

    Radon_obs = m.observations.obs_group['radon']
    Radon_obs_ts = Radon_obs['time_series']
    intervals_of_interest = Radon_obs_ts['interval'].unique()
    radon_df_dict = {}
    for i in intervals_of_interest:
        df = sfr_df[sfr_df['time'] == i]
        Ini_cond = (df.iloc[0]['Qin'], 0., 300.)
        df.loc[:, 'Flow'], df.loc[:, 'Rn'], df.loc[:, 'EC'] = \
              modflow_model.Calculate_Rn_from_SFR_with_simple_model(df, 
                                                                    Ini_cond)
        radon_df_dict[i] = df

    obs_set = 'radon'
    obs_df = m.observations.obs_group[obs_set]['time_series']
    obs_df = obs_df[obs_df['active'] == True]
    sfr_location = m.observations.obs_group[obs_set]['locations']['seg_loc']
    with open(os.path.join(modflow_model.data_folder, 'observations_radon.txt')
              , 'w') as f:
        for observation in obs_df.index:
            interval = int(obs_df['interval'].loc[observation])
            name = obs_df['name'].loc[observation]
            seg = sfr_location.loc[name]
            col_of_interest = 'Rn'
            df_radon = radon_df_dict[interval] 
            sim_obs = df_radon[df_radon['segment'] == seg] \
                              [col_of_interest]
            f.write('%f\n' % sim_obs)                

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
              
    #post.compareAllObs()
    #post.viewConcsByZone()
    
#    #
#    concobj = bf.UcnFile(modflow_model.data_folder + 'MT3D001.UCN')
#    arry = concobj.get_alldata()[32]
#    import matplotlib.pyplot as plt
#    #vmin, vmax = 0.0, 100.0
#    for i in range(7):
#        plt.figure()
#        #plt.imshow(arry[i], vmin=vmin, vmax=vmax, interpolation='none')
#        plt.imshow(arry[i], interpolation='none')
#        plt.colorbar()


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
        CONFIG = ConfigLoader('../../config/model_config.json')\
                        .set_environment("02_transient_flow")
        model_config = CONFIG.model_config
        model_folder = model_config['model_folder'] + model_config['grid_resolution'] + os.path.sep
        data_folder = model_config['data_folder']
        mt_exe_folder = model_config['mt_exe_folder']
        param_file = model_config['param_file']

    if param_file:
        run = run(model_folder, data_folder, mt_exe_folder, param_file=param_file, verbose=verbose)
    else:
        run = run(model_folder, data_folder, mt_exe_folder, verbose=verbose)
