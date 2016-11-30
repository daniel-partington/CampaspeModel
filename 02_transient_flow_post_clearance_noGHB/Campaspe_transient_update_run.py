import sys
import os
sys.path.append('C:\Workspace\part0075\GIT_REPOS')

import flopy.utils.binaryfile as bf

from HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface
# MM is short for model manager

def run(model_folder, data_folder, mf_exe, param_file=None):
    
    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder, "02_transient_flow_packaged.pkl"))
    
    name = MM.GW_build.keys()[0]
    # modify data folder
    #MM.GW_build['Campaspe'].data_folder
    #modify output folder
    #MM.GW_build['Campaspe'].out_data_folder
    
    # Load in the new parameters based on parameters.txt or dictionary of new parameters
    
    #MM.updateParameters('Campaspe', 'parameters.txt')
    
    if param_file:
        with open(param_file, 'r') as f:
            text = f.readlines()
            # Remove header    
            text = text[1:]
            # Read in parameters and replace values in parameters class for param
            for line in text:
                param_name, value = line.strip('\n').split('\t')
                value = value.lstrip()
                MM.GW_build[name].parameters.param[param_name]['PARVAL1'] = float(value)

   
    
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
    
    mapped_river = MM.GW_build[name].load_obj(os.path.join(model_folder,"Campaspe_Riv_model.shp_mapped.pkl"))
    print mapped_river
    
    simple_river = []
    riv_width_avg = 10.0 #m
    riv_bed_thickness = 0.10 #m
    for riv_cell in mapped_river: #MM.GW_build[name].polyline_mapped['Campaspe_Riv_model.shp']:
        row = riv_cell[0][0]
        col = riv_cell[0][1]
        if MM.GW_build[name].model_mesh3D[1][0][row][col] == -1:
            continue
        #print test_model.model_mesh3D
        stage = MM.GW_build[name].model_mesh3D[0][0][row][col]
        bed = MM.GW_build[name].model_mesh3D[0][0][row][col] - MM.GW_build[name].parameters.param['bed_depress']['PARVAL1']
        cond = riv_cell[1] * riv_width_avg * MM.GW_build[name].parameters.param['Kv_riv']['PARVAL1'] / riv_bed_thickness
        simple_river += [[0, row, col, stage, cond, bed]]
    
    riv = {}
    riv[0] = simple_river
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
        stage = MM.GW_build[name].model_mesh3D[0][0][row][col]
        bed = MM.GW_build[name].model_mesh3D[0][0][row][col] - MM.GW_build[name].parameters.param['RMstage']['PARVAL1']
        cond = riv_cell[1] * riv_width_avg * MM.GW_build[name].parameters.param['Kv_RM']['PARVAL1'] / riv_bed_thickness
        simple_river += [[0, row, col, stage, cond, bed]]
    
    riv = {}
    riv[0] = simple_river
    #MM.GW_build[name].boundaries.create_model_boundary_condition('Campaspe River', 'river', bc_static=True)
    MM.GW_build[name].boundaries.assign_boundary_array('Murray River', riv)
    
    
    print "************************************************************************"
    print " Updating recharge boundary "

    interp_rain = MM.GW_build[name].boundaries.bc['Rainfall']['bc_array']
    # Adjust rainfall to recharge using rainfall reduction
    for i in [1,2,3,7]:
        for key in interp_rain.keys():
            interp_rain[key][MM.GW_build[name].model_mesh3D[1][0]==i] = interp_rain[key][MM.GW_build[name].model_mesh3D[1][0]==i] * MM.GW_build[name].parameters.param['rch_red_'+zone_map[i]]['PARVAL1']

    for i in [4,5,6,]:
        for key in interp_rain.keys():
            interp_rain[key][MM.GW_build[name].model_mesh3D[1][0]==i] = interp_rain[key][MM.GW_build[name].model_mesh3D[1][0]==i] * 0.
    
    rch = {}
    for key in MM.GW_build[name].boundaries.bc['Rainfall']['bc_array'].keys():
        rch[key] = interp_rain[key] 
    
    MM.GW_build[name].boundaries.assign_boundary_array('Rain_reduced', rch)
    
    #print " Include irrigation in the recharge array"
    
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
    
    print "************************************************************************"
    print " Updating Drains boundary"
    
    mapped_drains = MM.GW_build[name].load_obj(os.path.join(model_folder,"Drain_Clip_model.shp_mapped.pkl"))
    
    simple_drain = []
    drain_width_avg = 3.0 #m
    drain_bed_thickness = 0.10 #m
    for drain_cell in mapped_drains:
        row = drain_cell[0][0]
        col = drain_cell[0][1]
        if MM.GW_build[name].model_mesh3D[1][0][row][col] == -1:
            continue
        #print MM.GW_build[name].model_mesh3D
        drain_bed = MM.GW_build[name].model_mesh3D[0][0][row][col] - MM.GW_build[name].parameters.param['drain_drop']['PARVAL1']
        drain_cond = drain_cell[1] * drain_width_avg * MM.GW_build[name].parameters.param['Kv_drain']['PARVAL1'] / drain_bed_thickness
        simple_drain += [[0, row, col, drain_bed, drain_cond]]
    
    drain = {}
    drain[0] = simple_drain
    
    MM.GW_build[name].boundaries.assign_boundary_array('Drain', drain)
    
    print "************************************************************************"
    print " Updating Channels boundary"
    
    simple_channel = []
    channel_width_avg = 10.0 #m
    channel_bed_thickness = 0.10 #m
    for channel_cell in MM.GW_build[name].polyline_mapped['Channel_Clip_model.shp']:
        row = channel_cell[0][0]
        col = channel_cell[0][1]
        if MM.GW_build[name].model_mesh3D[1][0][row][col] == -1:
            continue
        #print MM.GW_build[name].model_mesh3D
        channel_stage = MM.GW_build[name].model_mesh3D[0][0][row][col]
        channel_bed = MM.GW_build[name].model_mesh3D[0][0][row][col] - MM.GW_build[name].parameters.param['chan_drop']['PARVAL1']
        channel_cond = channel_cell[1] * channel_width_avg * MM.GW_build[name].parameters.param['Kv_chan']['PARVAL1'] / channel_bed_thickness
        simple_channel += [[0, row, col, channel_stage, channel_cond, channel_bed]]
    
    channel = {}
    channel[0] = simple_channel
    
    MM.GW_build[name].boundaries.assign_boundary_array('Channel', channel)

    print "************************************************************************"
    print " Updating pumping boundary"

    
    print "************************************************************************"
    print " Set initial head "

    path=os.path.join(data_folder,"model_01_steady_state\\")
    fname="01_steady_stateMainProcess"
    headobj = bf.HeadFile(path + fname +'.hds')
    times = headobj.get_times()        
    head = headobj.get_data(totim=times[-1])
    
    MM.GW_build[name].initial_conditions.set_as_initial_condition("Head", head)
    
    print "************************************************************************"
    print " Build and run MODFLOW model "
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ## Currently using flopyInterface directly rather than running from the ModelManager ...

    modflow_model = flopyInterface.ModflowModel(MM.GW_build[name], data_folder=data_folder)
    
    modflow_model.executable = mf_exe
    
    modflow_model.runMODFLOW()
    
    #print " Return the stream-aquifer exchange for reaches as list "
    #
    #SWGWexchange = [1]
    #
    #print " Return the average depth to the GW table in areas as list "
    
    #AvgDepthToGWTable = 1   
    #DepthToGWTable = [1]

    ss_converge = modflow_model.checkCovergence(path=os.path.join(data_folder,"model_01_steady_state") , name="01_steady_stateMainProcess")
    tr_converge = modflow_model.checkCovergence()
    
    if ss_converge:
        if tr_converge:
            modflow_model.writeObservations()

    
    #modflow_model.compareAllObs()
    
    #modflow_model.viewHeadsByZone2(nper='all')
 
    #ts = MM.GW_build[name].observations.obs_group['head']['time_series']
    #wells_of_interest = ['79234', '62589']    
    #wells = {}
    #for well in wells_of_interest:    
    #    wells[well] = ts[ts['name'] == well]
    #    wells[well] = ts[ts['name'] == well]

    #print modflow_model.getObservation(wells[wells_of_interest[1]]['obs_map'].tolist()[0], 0, 'head')    
    #obs = []
    #sim = []
    #time = []
    #import matplotlib.pyplot as plt
    #for well in wells_of_interest:
    #    for interval in wells[well]['interval']:
    #        obs += [wells[well][wells[well]['interval'] == interval]['value'].tolist()[0]]
    #        time += [wells[well][wells[well]['interval'] == interval]['datetime'].tolist()[0]]
    #        sim += [modflow_model.getObservation(wells_of_interest[1], interval, 'head')[0]]   
    
    #    plt.figure()
    #    plt.plot(time, obs)
    #    plt.plot(time, sim)
   
    #modflow_model.viewHeads()
    #modflow_model.viewHeads2()    
    #modflow_model.getRiverFlux('Campaspe River')
    #return modflow_model

if __name__ == "__main__":
    args = sys.argv
    if len(args) > 1:
        model_folder = sys.argv[1]
        data_folder = sys.argv[2]
        mf_exe = sys.argv[3]
        if len(args) > 4:
            param_file = sys.argv[4]
    else:
        grid_resolution = '5000'
        model_folder = r"C:\Workspace\part0075\MDB modelling\testbox\02_transient_flow\structured_model_grid_" + grid_resolution + r"m\\"
        data_folder = r"C:\Workspace\part0075\MDB modelling\testbox\PEST5000\master\\"    
        mf_exe = r"C:\Workspace\part0075\GIT_REPOS\CampaspeModel\MODFLOW-NWT_64.exe"
        param_file = r"C:\Workspace\part0075\MDB modelling\testbox\PEST5000\master\parameters.txt"
    
    if param_file:
        modflow_model = run(model_folder, data_folder, mf_exe, param_file=param_file)
    else:
        modflow_model = run(model_folder, data_folder, mf_exe)