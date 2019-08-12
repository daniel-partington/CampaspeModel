"""Update parameters and run Steady State transport model for Campaspe."""

import os
import sys

import numpy as np
import flopy

#sys.path.append('C:\Workspace\part0075\GIT_REPOS')
from HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface

# Configuration Loader
from HydroModelBuilder.Utilities.Config.ConfigLoader import ConfigLoader



def run(model_folder, data_folder, mt_exe_folder, param_file=None, verbose=True):
    """Model Runner."""
    # MM is short for model manager
    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder, r"01_steady_state_packaged.pkl"))

    name = list(MM.GW_build.keys())[0]
    m = MM.GW_build[name]

    # Load in the new parameters based on parameters.txt or dictionary of new parameters

    if param_file:
        MM.GW_build[name].updateModelParameters(os.path.join(data_folder, 'parameters.txt'), verbose=False)

    if verbose:
        print("************************************************************************")
        print(" Instantiate MODFLOW model ")

    ###########################################################################
    ###########################################################################
    ###########################################################################
    # Currently using flopyInterface directly rather than running from the ModelManager ...
    modflow_model = flopyInterface.ModflowModel(MM.GW_build[name], data_folder=os.path.join(data_folder, "model_" + m.name))

    modflow_model.nper = 1  # This is the number of stress periods which is set to 1 here
    modflow_model.perlen = 40000 * 365  # 
    modflow_model.nstp = 1  # This is the number of sub-steps to do in each stress period
    modflow_model.steady = True # This is to tell FloPy that is a steady model ... or not
    
    modflow_model.buildMODFLOW(transport=True, verbose=False, write=False) #, write=True)

    mt = flopy.mt3d.Mt3dms(modelname=modflow_model.name + '_transport', 
                           ftlfilename='mt3d_link.ftl', 
                           ftlfree=True,
                           modflowmodel=modflow_model.mf, 
                           model_ws=modflow_model.data_folder, 
                           exe_name=mt_exe_folder) #'MT3D-USGS_64.exe')

    #Add the BTN package to the model
    ibound = np.copy(modflow_model.model_data.model_mesh3D[1])        
    ibound[ibound == -1] = 0
    
    # This needs to be automatically generated with the map_raster2mesh routine ...
    zone_map = {1: 'qa', 2: 'utb', 3: 'utqa', 4: 'utam', 5: 'utaf', 6: 'lta', 7: 'bse'}

    prsity = np.zeros_like(ibound, dtype=float) 
    for zone in zone_map:          
        prsity[ibound == zone] = m.parameters.param['por_{}'.format(zone_map[zone])]['PARVAL1']
    
    flopy.mt3d.Mt3dBtn(mt, DRYCell=True, icbund=ibound, prsity=prsity,
                       ncomp=1, mcomp=1, cinact=-9.9E1, thkmin=-1.0E-6, ifmtcn=5, 
                       ifmtnp=0, ifmtrf=0, ifmtdp=0, nprs=0, 
                       timprs=None, savucn=1, nprobs=0,  laycon=1,
                       chkmas=1, nprmas=1, dt0=100000.0, ttsmax=100000.0)

    #Add the ADV package to the model
    flopy.mt3d.Mt3dAdv(mt, mixelm=0, percel=1, 
                       mxpart=250000, nadvfd=1, itrack=3,
                       wd=0.5, dceps=1.0E-4, nplane=0, npl=5,
                       nph=8, npmin=1, npmax=16, nlsink=0, 
                       npsink=8, dchmoc=0.01)

    #Add the DSP package to the model
    al = m.parameters.param['disp']['PARVAL1']
    flopy.mt3d.Mt3dDsp(mt, multiDiff=True, al=al, trpt=0.1, trpv=0.1, 
                       dmcoef=0.0)

    #Add the RCT package to the model
    flopy.mt3d.Mt3dRct(mt, isothm=1, ireact=1, igetsc=0, 
                       rc1=np.log(2)/(5730.0*365.0))

    #Add the GCG package to the model
    flopy.mt3d.Mt3dGcg(mt, mxiter=1000, iter1=100, isolve=1, 
                       ncrs=0, accl=1, cclose=1.0E-4, iprgcg=0)

    ssm_data = {}
    itype = flopy.mt3d.Mt3dSsm.itype_dict()
    for per in range(modflow_model.nper):
        ssm_data[per] = []
    ibound = modflow_model.model_data.model_mesh3D[1]        
    ibound[ibound == -1] = 0
    
    crch = {}
    for per in range(modflow_model.nper):
        crch[per] = []

    bcs = modflow_model.model_data.boundaries.bc

    for boundary in bcs:
        if bcs[boundary]['bc_type'] == 'recharge':
            for key in list(bcs[boundary]['bc_array'].keys()):
                crch[key] = np.ones_like(bcs[boundary]['bc_array'][key])               
                crch[key] = crch[key] * 100.0
                crch[key][ibound[0]==0] = 0.0
            
        if bcs[boundary]['bc_type'] == 'wells':
            for key in list(bcs[boundary]['bc_array'].keys()):
                for well in bcs[boundary]['bc_array'][key]:
                    ssm_data[key].append((well[0], well[1], well[2], 100.0, itype['WEL']))
                    
        if bcs[boundary]['bc_type'] == 'general head':
            for key in list(bcs[boundary]['bc_array'].keys()):
                for ghb in bcs[boundary]['bc_array'][key]:
                    ssm_data[key].append((ghb[0], ghb[1], ghb[2], 0.0, itype['GHB']))

    river_exists = False
    river_flow_exists = False
    for boundary in bcs:
        if bcs[boundary]['bc_type'] == 'river':
            river_exists = True
        if bcs[boundary]['bc_type'] == 'river_flow':
            river_flow_exists = True
            stream_flow_bc = bcs[boundary]['bc_array']
    
    if river_exists:
        river = {}
        river[0] = []
        for boundary in bcs:
            if bcs[boundary]['bc_type'] == 'river':
                time_key = list(bcs[boundary]['bc_array'].keys())[0]
                river[0] += bcs[boundary]['bc_array'][time_key]
            if bcs[boundary]['bc_type'] == 'channel':
                time_key = list(bcs[boundary]['bc_array'].keys())[0]
                river[0] += bcs[boundary]['bc_array'][time_key]
        
        for key in list(river.keys()):
            for riv in river[key]:
                ssm_data[key].append((riv[0], riv[1], riv[2], 100.0, itype['RIV'])) 

    if river_flow_exists:
        nsfinit = stream_flow_bc[0].shape[0]
        seg_len = np.unique(stream_flow_bc[0]['iseg'], return_counts=True)
        obs_sf = np.cumsum(seg_len[1])
        obs_sf = obs_sf.tolist()
        sf_stress_period_data = {0: [0, 0, 100.]}
        gage_output = None
        
        flopy.mt3d.Mt3dSft(mt, 
                           nsfinit=nsfinit, # Number of simulated stream reaches in SFR 
                           mxsfbc=nsfinit, # Maximum number of stream boundary conditions
                           icbcsf=81, # Integer directing to write reach-by-reach conc info to unit specified by integer 
                           ioutobs=82, # Unit number for output for concs at specified gage locations
                           ietsfr=0, #Specifies that mass is not removed with ET
                           wimp=1.0, # Stream solver time weighting factor (between 0.0 and 1.0)
                           wups=1.0, # Space weighting factor employed in the stream network solver (between 0 and 1)
                           isfsolv=1, # This is the default and only supported value at this stage
                           cclosesf=1.0E-6, # This is the closure criterion for the SFT solver 
                           mxitersf=10, # Maximum number of iterations for the SFT solver
                           crntsf=1.0,  # Courant constraint specific to SFT
                           iprtxmd=0, # flag to print SFT solution info to standard output file (0 = print nothing)
                           coldsf=100.0, # Initial concentration in the stream network (can also be specified as an array for each reach)
                           dispsf=100.0, # Dispersivity in each stream reach
                           nobssf=len(obs_sf), #, 
                           obs_sf=obs_sf, 
                           sf_stress_period_data=sf_stress_period_data, # Dict for each stress period with each dict containing ISFNBC, ISFBCTYP, (CBCSF(n), n=1, NCOMP)
                           filenames=gage_output)

               
    flopy.mt3d.Mt3dSsm(mt, stress_period_data=ssm_data, crch=crch)
    
    mt.write_input()
    
    success, buff = mt.run_model(silent=True)

#    import pandas as pd
#    sfr_transport = pd.read_csv(os.path.join(data_folder, "model_" + m.name,"01_steady_state_transport"), delim_whitespace=True, skiprows=1)
#    for sfrnode in [1, 20, 120]:
#        ax = sfr_transport[sfr_transport['SFR-NODE'] == sfrnode].plot(x='TIME', y=['SFR-CONCENTRATION'])
#        ax.set_title("C14")

    
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
                        .set_environment("01_steady_state")
        model_config = CONFIG.model_config
        model_folder = model_config['model_folder'] + model_config['grid_resolution'] + os.path.sep
        data_folder = model_config['data_folder']
        mt_exe_folder = model_config['mt_exe_folder']
        param_file = model_config['param_file']

    if param_file:
        run(model_folder, data_folder, mt_exe_folder, param_file=param_file, verbose=verbose)
    else:
        run(model_folder, data_folder, mt_exe_folder, verbose=verbose)
