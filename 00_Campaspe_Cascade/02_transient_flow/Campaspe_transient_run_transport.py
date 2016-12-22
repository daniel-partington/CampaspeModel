"""Update parameters and run Steady State model for Campaspe."""

import os
import sys

import numpy as np
import flopy
import flopy.utils.binaryfile as bf

sys.path.append('C:\Workspace\part0075\GIT_REPOS')
from HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface

# Configuration Loader
from HydroModelBuilder.Utilities.Config.ConfigLoader import CONFIG
# import flopy.utils.binaryfile as bf


# MM is short for model manager


def run(model_folder, data_folder, mf_exe_folder, param_file=None):
    """Model Runner."""

    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder, "02_transient_flow_packaged.pkl"))

    name = MM.GW_build.keys()[0]

    # Load in the new parameters based on parameters.txt or dictionary of new parameters

    #if param_file:
    #    MM.GW_build[name].updateModelParameters(os.path.join(data_folder, 'parameters.txt'))

    print "************************************************************************"
    print " Instantiate MODFLOW model "

    ###########################################################################
    ###########################################################################
    ###########################################################################
    # Currently using flopyInterface directly rather than running from the ModelManager ...
    modflow_model = flopyInterface.ModflowModel(MM.GW_build[name], data_folder=data_folder)

    modflow_model.executable = mf_exe_folder

    modflow_model.buildMODFLOW()

    #print modflow_model.model_data.observations.obs_group.keys()
    #sys.exit()
    print "************************************************************************"
    print " Instantiate MT3D model "
    
    mt = flopy.mt3d.Mt3dms(modelname=modflow_model.name+'_transport', ftlfilename='mt3d_link.ftl', 
                           modflowmodel=modflow_model.mf, model_ws=modflow_model.data_folder, 
                           exe_name='MT3D-USGS_64.exe')

    
    print "************************************************************************"
    print " Set initial conc from ss solution "

    path=os.path.join(data_folder,"model_01_steady_state\\")
    concobj = bf.UcnFile(path + 'MT3D001.UCN')
    times = concobj.get_times()        
    conc_init = concobj.get_data(totim=times[-1])
    
    print("Add the BTN package to the model")
    ibound = modflow_model.model_data.model_mesh3D[1]        
    ibound[ibound == -1] = 0
    prsity = MM.GW_build[name].parameters.param['porosity']['PARVAL1']
    flopy.mt3d.Mt3dBtn(mt, optional_args='DRYCELL', icbund=ibound, prsity=prsity
                       sconc=conc_init, ncomp=1, mcomp=1, 
                       cinact=-9.9E1, thkmin=-1.0E-6, ifmtcn=5, 
                       ifmtnp=0, ifmtrf=0, ifmtdp=0, nprs=0, 
                       timprs=None, savucn=1, nprobs=0, 
                       chkmas=1, nprmas=1, dt0=10000.0, ttsmax=100000.0)

    print("Add the ADV package to the model")
    flopy.mt3d.Mt3dAdv(mt, mixelm=0, percel=1, 
                       mxpart=250000, nadvfd=1, itrack=3,
                       wd=0.5, dceps=1.0E-4, nplane=0, npl=5,
                       nph=8, npmin=1, npmax=16, nlsink=0, 
                       npsink=8, dchmoc=0.01)

    print("Add the DSP package to the model")
    al = MM.GW_build[name].parameters.param['disp']['PARVAL1']
    flopy.mt3d.Mt3dDsp(mt, multiDiff=True, al=al, trpt=0.1, trpv=0.1, dmcoef=0.0)

    print("Add the RCT package to the model")
    flopy.mt3d.Mt3dRct(mt, isothm=1, ireact=1, igetsc=0, rc1=np.log(2)/(5730.0*365.0))

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
    for per in range(modflow_model.nper):
        crch[per] = []
    # End for

    river = {}
    river[0] = []

    bc = modflow_model.model_data.boundaries.bc
    for boundary in bc:
        bc_boundary = bc[boundary]
        bc_type = bc_boundary['bc_type']
        bc_array = bc_boundary['bc_array']

        if bc_type == 'recharge':
            for key in bc_array.keys():
                crch[key] = np.ones_like(bc_array[key])               
                crch[key] = crch[key] * 100.0
                crch[key][ibound[0]==0] = 0.0
        
        if (bc_type == 'river') or (bc_type == 'channel'):
            time_key = bc_array.keys()[0]
            river[0] += bc_array[time_key]
        # End if
            
        if bc_type == 'wells':
            for key in bc_array.keys():
                for well in bc_array[key]:
                    ssm_data[key].append((well[0], well[1], well[2], 100.0, itype['WEL']))
                    
        if bc_type == 'drain':
            for key in bc_array.keys():
                for drain in bc_array[key]:
                    ssm_data[key].append((drain[0], drain[1], drain[2], 100.0, itype['DRN']))

        if bc_type == 'general head':
            for key in bc_array.keys():
                for ghb in bc_array[key]:
                    ssm_data[key].append((ghb[0], ghb[1], ghb[2], 0.0, itype['GHB']))

    if len(river) > 0:
        for key in river.keys():
            for riv in river[key]:
                ssm_data[key].append((riv[0], riv[1], riv[2], 100.0, itype['RIV']))
            # End for
        # End for
    # End if

    ssm_data_temp = {}
    for key in ssm_data.keys():
        if len(ssm_data[key]) > 0:
            ssm_data_temp[key] = ssm_data[key]
        #end if
    #end for
    ssm_data = ssm_data_temp
    

    print("Add the SSM package to the model")
    flopy.mt3d.Mt3dSsm(mt, stress_period_data=ssm_data, crch=crch)
    
    print(" Add LMT package to the MODFLOW model to allow linking with MT3DMS")
    flopy.modflow.ModflowLmt(modflow_model.mf)

    mt.write_input()
    
    success, buff = mt.run_model()
    
    post = flopyInterface.MT3DPostProcess(modflow_model)
    
    post.writeObservations()
    #post.compareAllObs()
    #post.viewConcsByZone()
    
#    #
#    concobj = bf.UcnFile(modflow_model.data_folder + 'MT3D001.UCN')
#    arry = concobj.get_alldata()[680]
#    import matplotlib.pyplot as plt
#    vmin, vmax = 0.0, 100.0
#    for i in range(7):
#        plt.figure()
#        plt.imshow(arry[i], vmin=vmin, vmax=vmax, interpolation='none')


if __name__ == "__main__":
    args = sys.argv
    if len(args) > 1:
        model_folder = sys.argv[1]
        data_folder = sys.argv[2]
        mf_exe_folder = sys.argv[3]
        if len(args) > 4:
            param_file = sys.argv[4]
    else:
        model_config = CONFIG.model_config
        model_folder = model_config['model_folder'] + model_config['grid_resolution'] + os.path.sep
        data_folder = model_config['data_folder']
        mf_exe_folder = model_config['mf_exe_folder']
        param_file = model_config['param_file']

    if param_file:
        run = run(model_folder, data_folder, mf_exe_folder, param_file=param_file)
    else:
        run = run(model_folder, data_folder, mf_exe_folder)
