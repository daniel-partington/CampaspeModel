import sys
sys.path.append('C:\Workspace\part0075\GIT_REPOS')
import os
import numpy as np
from osgeo import osr

from HydroModelBuilder.HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface
from HydroModelBuilder.HydroModelBuilder.GISInterface.GDALInterface.GDALInterface import GDALInterface

import flopy.utils.binaryfile as bf
# MM is short for model manager

def run(model_folder, data_folder, mf_exe_folder, param_file=None):

    # Define basic model parameters:
    Proj_CS = osr.SpatialReference()
    Proj_CS.ImportFromEPSG(28355) # 28355 is the code for gda94 mga zone 55; http://spatialreference.org/ref/epsg/gda94-mga-zone-55/
    
    Interface = GDALInterface()
    Interface.projected_coordinate_system = Proj_CS 
    Interface.pcs_EPSG = "EPSG:28355"
    
    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder, r"01_steady_state_packaged.pkl"))
    
    name = MM.GW_build.keys()[0]
    
    # Load in the new parameters based on parameters.txt or dictionary of new parameters

    MM.GW_build[name].GISInterface = Interface
    MM.GW_build[name].updateGISinterface()
    
    if param_file:
        MM.GW_build[name].updateModelParameters(os.path.join(data_folder, 'parameters.txt'))
    
    print "************************************************************************"
    print " Updating HGU parameters "
    
    # This needs to be automatically generated from with the map_raster2mesh routine ...
    zone_map = {1:'qa', 2:'utb', 3:'utqa', 4:'utam', 5:'utaf', 6:'lta', 7:'bse'}
    
    Zone = MM.GW_build[name].model_mesh3D[1].astype(float)
    Kh = MM.GW_build[name].model_mesh3D[1].astype(float)
    Kv = MM.GW_build[name].model_mesh3D[1].astype(float)
    Sy = MM.GW_build[name].model_mesh3D[1].astype(float)
    SS = MM.GW_build[name].model_mesh3D[1].astype(float)
    for key in zone_map.keys():
        Kh[Zone == key] = MM.GW_build[name].parameters.param['Kh_' + zone_map[key]]['PARVAL1']
        Kv[Zone == key] = MM.GW_build[name].parameters.param['Kv_' + zone_map[key]]['PARVAL1']
        Sy[Zone == key] = MM.GW_build[name].parameters.param['Sy_' + zone_map[key]]['PARVAL1']
        SS[Zone == key] = MM.GW_build[name].parameters.param['SS_' + zone_map[key]]['PARVAL1']
    
    MM.GW_build[name].properties.assign_model_properties('Kh', Kh)
    MM.GW_build[name].properties.assign_model_properties('Kv', Kv)
    MM.GW_build[name].properties.assign_model_properties('Sy', Sy)
    MM.GW_build[name].properties.assign_model_properties('SS', SS)
    
    print "************************************************************************"
    print " Updating river parameters "
    
    mapped_river = MM.GW_build[name].load_obj(os.path.join(model_folder, r"Campaspe_Riv_model.shp_mapped.pkl"))
    
    simple_river = []
    riv_width_avg = 10.0 #m
    riv_bed_thickness = 0.10 #m
    for riv_cell in mapped_river: #MM.GW_build[name].polyline_mapped['Campaspe_Riv_model.shp']:
        row = riv_cell[0][0]
        col = riv_cell[0][1]
        if MM.GW_build[name].model_mesh3D[1][0][row][col] == -1:
            continue
        #print test_model.model_mesh3D
        stage = MM.GW_build[name].model_mesh3D[0][0][row][col] - 0.01
        bed = MM.GW_build[name].model_mesh3D[0][0][row][col] - 0.1 - MM.GW_build[name].parameters.param['bed_depress']['PARVAL1']
        cond = riv_cell[1] * riv_width_avg * MM.GW_build[name].parameters.param['Kv_riv']['PARVAL1'] / riv_bed_thickness
        simple_river += [[0, row, col, stage, cond, bed]]
    
    riv = {}
    riv[0] = simple_river
    #MM.GW_build[name].boundaries.create_model_boundary_condition('Campaspe River', 'river', bc_static=True)
    MM.GW_build[name].boundaries.assign_boundary_array('Campaspe River', riv)
    
    mapped_river = MM.GW_build[name].load_obj(os.path.join(model_folder, r"River_Murray_model.shp_mapped.pkl"))
    
    simple_river = []
    riv_width_avg = 10.0 #m
    riv_bed_thickness = 0.10 #m
    for riv_cell in mapped_river: #MM.GW_build[name].polyline_mapped['Campaspe_Riv_model.shp']:
        row = riv_cell[0][0]
        col = riv_cell[0][1]
        if MM.GW_build[name].model_mesh3D[1][0][row][col] == -1:
            continue
        #print test_model.model_mesh3D
        stage = MM.GW_build[name].model_mesh3D[0][0][row][col] - 0.01
        bed = MM.GW_build[name].model_mesh3D[0][0][row][col] -0.1 - MM.GW_build[name].parameters.param['RMstage']['PARVAL1']
        cond = riv_cell[1] * riv_width_avg * MM.GW_build[name].parameters.param['Kv_RM']['PARVAL1'] / riv_bed_thickness
        simple_river += [[0, row, col, stage, cond, bed]]
    
    riv = {}
    riv[0] = simple_river
    #MM.GW_build[name].boundaries.create_model_boundary_condition('Campaspe River', 'river', bc_static=True)
    MM.GW_build[name].boundaries.assign_boundary_array('Murray River', riv)
    
    
    print "************************************************************************"
    print " Updating recharge boundary "

    # Adjust rainfall to recharge using 10% magic number
    interp_rain = np.copy(MM.GW_build[name].boundaries.bc['Rainfall']['bc_array'])

    for i in [1,2,3,7]:
        interp_rain[MM.GW_build[name].model_mesh3D[1][0]==i] = interp_rain[MM.GW_build[name].model_mesh3D[1][0]==i] * MM.GW_build[name].parameters.param['rch_red_'+zone_map[i]]['PARVAL1']
    for i in [4,5,6]:
        interp_rain[MM.GW_build[name].model_mesh3D[1][0]==i] = interp_rain[MM.GW_build[name].model_mesh3D[1][0]==i] * 0.
    
    rch = {}
    rch[0] = interp_rain
    
    MM.GW_build[name].boundaries.assign_boundary_array('Rain_reduced', rch)
    
    print "************************************************************************"
    print " Updating Murray River GHB boundary"
    
    mapped_river = MM.GW_build[name].load_obj(os.path.join(model_folder, r"River_Murray_model.shp_mapped.pkl"))
    
    MurrayGHB = []
    for MurrayGHB_cell in mapped_river:
        row = MurrayGHB_cell[0][0]
        col = MurrayGHB_cell[0][1]
        #print MM.GW_build[name].model_mesh3D
        for lay in range(MM.GW_build[name].model_mesh3D[1].shape[0]):    
            if MM.GW_build[name].model_mesh3D[1][0][row][col] == -1:
                continue
            MurrayGHBstage = MM.GW_build[name].model_mesh3D[0][0][row][col] + MM.GW_build[name].parameters.param['MGHB_stage']['PARVAL1']
            dx = MM.GW_build[name].gridHeight
            dz = MM.GW_build[name].model_mesh3D[0][lay][row][col] - MM.GW_build[name].model_mesh3D[0][lay+1][row][col]
            MGHBconductance = dx * dz * MM.GW_build[name].parameters.param['MGHBcond']['PARVAL1']
            MurrayGHB += [[lay, row, col, MurrayGHBstage, MGHBconductance]]
    
    ghb = {}
    ghb[0] = MurrayGHB
    
    print "************************************************************************"
    print " Updating Western GHB boundary"
    
    mapped_west = MM.GW_build[name].load_obj(os.path.join(model_folder, r"western_head_model.shp_mapped.pkl"))
    
    WestGHB = []
    for WestGHB_cell in mapped_west:
        row = WestGHB_cell[0][0]
        col = WestGHB_cell[0][1]
        for lay in range(MM.GW_build[name].model_mesh3D[1].shape[0]):    
            if MM.GW_build[name].model_mesh3D[1][0][row][col] == -1:
                continue
            WestGHBstage = MM.GW_build[name].model_mesh3D[0][0][row][col] + MM.GW_build[name].parameters.param['WGHB_stage']['PARVAL1']
            dx = MM.GW_build[name].gridHeight
            dz = MM.GW_build[name].model_mesh3D[0][lay][row][col] - MM.GW_build[name].model_mesh3D[0][lay+1][row][col]
            WGHBconductance = dx * dz * MM.GW_build[name].parameters.param['WGHBcond']['PARVAL1']
            WestGHB += [[lay, row, col, WestGHBstage, WGHBconductance]]
    
    ghb[0] += WestGHB
    
    print "************************************************************************"
    print " Updating Eastern GHB boundary"
    
    mapped_east = MM.GW_build[name].load_obj(os.path.join(model_folder, r"eastern_head_model.shp_mapped.pkl"))
    
    EastGHB = []
    for EastGHB_cell in mapped_east:
        row = EastGHB_cell[0][0]
        col = EastGHB_cell[0][1]
        for lay in range(MM.GW_build[name].model_mesh3D[1].shape[0]):    
            if MM.GW_build[name].model_mesh3D[1][0][row][col] == -1:
                continue
            EastGHBstage = MM.GW_build[name].model_mesh3D[0][0][row][col] + MM.GW_build[name].parameters.param['EGHB_stage']['PARVAL1']
            dx = MM.GW_build[name].gridHeight
            dz = MM.GW_build[name].model_mesh3D[0][lay][row][col] - MM.GW_build[name].model_mesh3D[0][lay+1][row][col]
            EGHBconductance = dx * dz * MM.GW_build[name].parameters.param['EGHBcond']['PARVAL1']
            EastGHB += [[lay, row, col, EastGHBstage, EGHBconductance]]
    
    ghb[0] += EastGHB
    
    print "************************************************************************"
    print " Updating GHB boundary"
    
    MM.GW_build[name].boundaries.assign_boundary_array('GHB', ghb)

#    print "************************************************************************"
#    print " Set initial head "
#
#    path=os.path.join(data_folder)
#    fname="01_steady_stateMainProcess"
#    headobj = bf.HeadFile(path + fname +'.hds')
#    times = headobj.get_times()        
#    head = headobj.get_data(totim=times[-1])
#    
#    MM.GW_build[name].initial_conditions.set_as_initial_condition("Head", head)

    import pandas as pd
    
    C14_wells_info_file = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Chemistry\C14_bore_depth.csv"
    df_C14_info = pd.read_csv(C14_wells_info_file)    
    df_C14_info = df_C14_info.dropna()
    df_C14_info = df_C14_info.set_index('Bore_id')    
    
    print "************************************************************************"
    print " Build and run MODFLOW model "
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ## Currently using flopyInterface directly rather than running from the ModelManager ...
    modflow_model = flopyInterface.ModflowModel(MM.GW_build[name], data_folder=data_folder)
    
    modflow_model.executable = mf_exe_folder
    
    modflow_model.runMODFLOW()
    
    #print " Return the stream-aquifer exchange for reaches as list "
    #
    #SWGWexchange = [1]
    #
    #print " Return the average depth to the GW table in areas as list "
    
    #AvgDepthToGWTable = 1   
    #DepthToGWTable = [1]
    
    #modflow_model.writeObservations()
    
    #modflow_model.viewHeadsByZone()
    
    modflow_model.checkCovergence()
    
    #modflow_model.viewHeads()
    #modflow_model.viewHeads2()    
    
    import flopy    
    
    mf = modflow_model.mf    
    model_path = os.path.join(data_folder, 'model_01_steady_state')
    
    mpp = flopy.modpath.Modpath('wells', exe_name='mp6x64', modflowmodel=mf, model_ws=model_path)
    mppbas = flopy.modpath.ModpathBas(mpp, hnoflo=mf.bas6.hnoflo, hdry=mf.upw.hdry, 
                                     ibound=mf.bas6.ibound.array, prsity=0.2)
    #sim = mpp.create_mpsim(trackdir='backward', simtype='endpoint', packages='WEL')
    #mpp.write_input()
    success, buff = mpp.run_model()
    
#    endpoints_file = r"C:\Workspace\part0075\MDB modelling\testbox\PEST2\model_01_steady_state\wells.mpend"
    # load the endpoint data
    endfile = os.path.join(model_path, 'wells.mpend')
    endfile = os.path.join(model_path, mpp.sim.endpoint_file)
    endobj = flopy.utils.EndpointFile(endfile)
    ept = endobj.get_alldata()

    # # of wells
#    final_time = {}
#    for i in range(len(well_name.keys())):
#        final_time[well_name[i]] = np.mean([p[5] for p in ept if p[1] == i])
#        print final_time[well_name[i]]
        
    # load the pathline data
    #pthfile = os.path.join(model_path, 'wells.mppth')
    #pthfile = os.path.join(model_path, pathline_file)
    #pthobj = flopy.utils.PathlineFile(pthfile)
    #plines = pthobj.get_alldata()
   
    
    # Do some post-processing
    
    import matplotlib.pyplot as plt
    import flopy.utils.binaryfile as bf
    plt.subplot(1,1,1,aspect='equal')
#    hds = bf.HeadFile(modelname + '.hds')
#    times = hds.get_times()
#    head = hds.get_data(totim=times[0])
#    levels = np.arange(1,10,1)
#    extent = (delr/2., Lx - delr/2., Ly - delc/2., delc/2.)
    modelmap = flopy.plot.ModelMap(model=mf)
    #quiver = modelmap.plot_discharge(frf, fff, head=head)
    
    
    modelmap.plot_ibound()
    modelmap.plot_bc('WEL')
    linecollection = modelmap.plot_grid(color='lightgray')
    #CS = plt.contour(head[0, :, :], levels=levels, extent=extent)
    #plt.clabel(CS)
    
    for d in mf.wel.stress_period_data[0]:
        modelmap.plot_endpoint(ept, direction='starting', selection_direction='ending') #, selection=(d[0], d[1], d[2]), zorder=100)
    
    # construct maximum travel time to plot (200 years - MODFLOW time unit is seconds)
    #travel_time_max = 200. * 365.25 * 24. * 60. * 60. 
    #ctt = '<={}'.format(travel_time_max)
    
    # plot the pathlines
    #modelmap.plot_pathline(plines, layer='all', colors='red', travel_time=ctt)
#    modelmap.plot_pathline(plines, layer=4, colors='red', alpha=0.5, travel_time=ctt)
    
    #plt.colorbar()
    plt.show()        

if __name__ == "__main__":
    args = sys.argv
    if len(args) > 1:
        model_folder = sys.argv[1]
        data_folder = sys.argv[2]
        mf_exe_folder = sys.argv[3]
        if len(args) > 4:
            param_file = sys.argv[4]
    else:
        grid_resolution = '5000'
        model_folder = r"C:\Workspace\part0075\MDB modelling\testbox\01_steady_state\structured_model_grid_" + grid_resolution + r"m\\"
        data_folder = r"C:\Workspace\part0075\MDB modelling\testbox\PEST5000\master\\"    
        mf_exe_folder = r"C:\Workspace\part0075\GIT_REPOS\CampaspeModel\MODFLOW-NWT_64.exe"
        param_file = r"C:\Workspace\part0075\MDB modelling\testbox\PEST5000\master\parameters.txt"
    
    if param_file:
        run(model_folder, data_folder, mf_exe_folder, param_file=param_file)
    else:
        run(model_folder, data_folder, mf_exe_folder)