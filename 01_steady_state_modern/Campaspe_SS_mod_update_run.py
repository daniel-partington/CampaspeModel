import sys
sys.path.append('C:\Workspace\part0075\GIT_REPOS')
import os
import numpy as np

from HydroModelBuilder.HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface
import flopy.utils.binaryfile as bf
# MM is short for model manager

def run(model_folder, data_folder, mf_exe_folder, param_file=None):
    
    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder, r"01_steady_state_packaged.pkl"))
    
    name = MM.GW_build.keys()[0]
  
    # Load in the new parameters based on parameters.txt or dictionary of new parameters
    
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

    camp_river = []
    riv_width_avg = 10.0 #m
    riv_bed_thickness = 0.10 #m
    for riv_cell in mapped_river:
        row = riv_cell[0][0]
        col = riv_cell[0][1]
        if MM.GW_build[name].model_mesh3D[1][0][row][col] == -1:
            continue
        #print SS_model.model_mesh3D
        stage = MM.GW_build[name].model_mesh3D[0][0][row][col]
        bed = MM.GW_build[name].model_mesh3D[0][0][row][col] - MM.GW_build[name].parameters.param['bed_depress']['PARVAL1']
        cond = riv_cell[1] * riv_width_avg * MM.GW_build[name].parameters.param['Kv_riv']['PARVAL1'] / riv_bed_thickness
        camp_river += [[0, row, col, stage, cond, bed]]
    
    riv = {}
    riv[0] = camp_river
    
    # Create list of campaspe river cells and make sure other riv boundaries are not in cells e.g. with the channels
    camp_river_locs = [(x[1],x[2]) for x in camp_river]

    MM.GW_build[name].boundaries.assign_boundary_array('Campaspe River', riv)

    # Now update the Murray River
    mapped_river = MM.GW_build[name].load_obj(os.path.join(model_folder, r"River_Murray_model.shp_mapped.pkl"))
    
    simple_river = []
    riv_width_avg = 10.0 #m
    riv_bed_thickness = 0.10 #m
    for riv_cell in mapped_river:
        row = riv_cell[0][0]
        col = riv_cell[0][1]
        if MM.GW_build[name].model_mesh3D[1][0][row][col] == -1:
            continue
    
        # Don't assign boundary if in same location as Campaspe river
        if (row, col) in camp_river_locs:
            continue
        #print SS_model.model_mesh3D
        stage = MM.GW_build[name].model_mesh3D[0][0][row][col]
        bed = MM.GW_build[name].model_mesh3D[0][0][row][col] - MM.GW_build[name].parameters.param['RMstage']['PARVAL1']
        cond = riv_cell[1] * riv_width_avg * MM.GW_build[name].parameters.param['Kv_RM']['PARVAL1'] / riv_bed_thickness
        simple_river += [[0, row, col, stage, cond, bed]]
    
    riv = {}
    riv[0] = simple_river

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
    
    #print " Include irrigation in the recharge array"
    
    #print MM.GW_build[name].observations.obs_group['head'].keys()
    #print MM.GW_build[name].observations.obs_group['head']['mapped_observations']
    
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
    print " Updating Drains BC"
    
    
    simple_drain = []
    drain_width_avg = 3.0 #m
    drain_bed_thickness = 0.10 #m
    for drain_cell in MM.GW_build[name].polyline_mapped['Drain_Clip_model.shp']:
        row = drain_cell[0][0]
        col = drain_cell[0][1]
        if MM.GW_build[name].model_mesh3D[1][0][row][col] == -1:
            continue
        #print SS_model.model_mesh3D
        drain_bed = MM.GW_build[name].model_mesh3D[0][0][row][col] - MM.GW_build[name].parameters.param['drain_drop']['PARVAL1']
        drain_cond = drain_cell[1] * drain_width_avg * MM.GW_build[name].parameters.param['Kv_drain']['PARVAL1'] / drain_bed_thickness
        simple_drain += [[0, row, col, drain_bed, drain_cond]]
    
    drain = {}
    drain[0] = simple_drain
   
    MM.GW_build[name].boundaries.assign_boundary_array('Drain', drain)
    
    
    print "************************************************************************"
    print " Updating Channels BC"
    
    simple_channel = []
    channel_width_avg = 10.0 #m
    channel_bed_thickness = 0.10 #m
    for channel_cell in MM.GW_build[name].polyline_mapped['Channel_Clip_model.shp']:
        row = channel_cell[0][0]
        col = channel_cell[0][1]
        if MM.GW_build[name].model_mesh3D[1][0][row][col] == -1:
            continue
        # Don't assign boundary if in same location as Campaspe river
        if (row, col) in camp_river_locs:
            continue
        #print SS_model.model_mesh3D
        channel_stage = MM.GW_build[name].model_mesh3D[0][0][row][col]
        channel_bed = MM.GW_build[name].model_mesh3D[0][0][row][col] - MM.GW_build[name].parameters.param['chan_drop']['PARVAL1']
        channel_cond = channel_cell[1] * channel_width_avg * MM.GW_build[name].parameters.param['Kv_chan']['PARVAL1'] / channel_bed_thickness
        simple_channel += [[0, row, col, channel_stage, channel_cond, channel_bed]]
    
    channel = {}
    channel[0] = simple_channel
    
    MM.GW_build[name].boundaries.assign_boundary_array('Channel', channel)


    print "************************************************************************"
    print " Set initial head "

    path=os.path.join(data_folder)
    fname="Initial"
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
    
    riv_exch = modflow_model.getRiverFlux('Campaspe River')
    for key in riv_exch.keys():
        print 'Campaspe River net flux: ' + str(round(sum([x[0] for x in riv_exch[key]]))) + ' m3/d'    

    riv_exch = modflow_model.getRiverFlux('Murray River')
    for key in riv_exch.keys():
        print 'Murray River net flux: ' + str(round(sum([x[0] for x in riv_exch[key]]))) + ' m3/d'    

    #modflow_model.viewHeads()
    modflow_model.viewHeads2()    

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
        model_folder = r"C:\Workspace\part0075\MDB modelling\testbox\01_steady_state_modern\structured_model_grid_" + grid_resolution + r"m\\"
        data_folder = r"C:\Workspace\part0075\MDB modelling\testbox\PEST_SS_MOD5000\\"    
        mf_exe_folder = r"C:\Workspace\part0075\GIT_REPOS\CampaspeModel\MODFLOW-NWT_64.exe"
        param_file = r"C:\Workspace\part0075\MDB modelling\testbox\PEST_SS_MOD5000\parameters.txt"
    
    if param_file:
        run(model_folder, data_folder, mf_exe_folder, param_file=param_file)
    else:
        run(model_folder, data_folder, mf_exe_folder)