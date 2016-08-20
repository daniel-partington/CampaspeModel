import sys
sys.path.append('C:\Workspace\part0075\GIT_REPOS')
import os

from HydroModelBuilder.HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface
# MM is short for model manager

def run(model_folder, data_folder, mf_exe_folder, param_file=None):
    
    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder, r"01_steady_state_packaged.pkl"))
    
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
                if param_name in MM.GW_build[name].parameters.param.keys():        
                    MM.GW_build[name].parameters.param[param_name]['PARVAL1'] = float(value)
           
    
    print "************************************************************************"
    print " Updating HGU parameters "
    
    # This needs to be automatically generated from with the map_raster2mesh routine ...
    zone_map = {1:'qa', 2:'utb', 3:'utqa', 4:'utam', 5:'utaf', 6:'lta', 7:'bse'}
    
    Kh = MM.GW_build[name].model_mesh3D[1].astype(float)
    Kv = MM.GW_build[name].model_mesh3D[1].astype(float)
    Sy = MM.GW_build[name].model_mesh3D[1].astype(float)
    SS = MM.GW_build[name].model_mesh3D[1].astype(float)
    for key in zone_map.keys():
        Kh[Kh == key] = MM.GW_build[name].parameters.param['Kh_' + zone_map[key]]['PARVAL1']
        Kv[Kv == key] = MM.GW_build[name].parameters.param['Kv_' + zone_map[key]]['PARVAL1']
        Sy[Sy == key] = MM.GW_build[name].parameters.param['Sy_' + zone_map[key]]['PARVAL1']
        SS[SS == key] = MM.GW_build[name].parameters.param['SS_' + zone_map[key]]['PARVAL1']
    
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
        stage = MM.GW_build[name].model_mesh3D[0][0][row][col]
        bed = MM.GW_build[name].model_mesh3D[0][0][row][col] - MM.GW_build[name].parameters.param['bed_depress']['PARVAL1']
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
    
    interp_rain = MM.GW_build[name].boundaries.bc['Rainfall']['bc_array'] * MM.GW_build[name].parameters.param['magic_rain']['PARVAL1']
    
    rch = {}
    rch[0] = interp_rain
    
    #MM.GW_build[name].boundaries.create_model_boundary_condition('Rain_reduced', 'recharge', bc_static=True)
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
            MurrayGHBstage = (MM.GW_build[name].model_mesh3D[0][lay+1][row][col] + MM.GW_build[name].model_mesh3D[0][lay][row][col])/2. + MM.GW_build[name].parameters.param['MGHB_stage']['PARVAL1']
            dx = MM.GW_build[name].model_mesh3D[0][0][0][1] - MM.GW_build[name].model_mesh3D[0][0][0][0]
            MGHBconductance = dx * MM.GW_build[name].parameters.param['MGHBcond']['PARVAL1']
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
            WestGHBstage = (MM.GW_build[name].model_mesh3D[0][lay+1][row][col] + MM.GW_build[name].model_mesh3D[0][lay][row][col])/2. + MM.GW_build[name].parameters.param['WGHB_stage']['PARVAL1']
            dx = MM.GW_build[name].model_mesh3D[0][0][0][1] - MM.GW_build[name].model_mesh3D[0][0][0][0]
            WGHBconductance = dx * MM.GW_build[name].parameters.param['WGHBcond']['PARVAL1']
            WestGHB += [[lay, row, col, WestGHBstage, WGHBconductance]]
    
    ghb[0] += WestGHB
    
    print "************************************************************************"
    print " Updating Western GHB boundary"
    
    mapped_east = MM.GW_build[name].load_obj(os.path.join(model_folder, r"eastern_head_model.shp_mapped.pkl"))
    
    EastGHB = []
    for EastGHB_cell in mapped_east:
        row = EastGHB_cell[0][0]
        col = EastGHB_cell[0][1]
        for lay in range(MM.GW_build[name].model_mesh3D[1].shape[0]):    
            if MM.GW_build[name].model_mesh3D[1][0][row][col] == -1:
                continue
            EastGHBstage = (MM.GW_build[name].model_mesh3D[0][lay+1][row][col] + MM.GW_build[name].model_mesh3D[0][lay][row][col])/2. + MM.GW_build[name].parameters.param['EGHB_stage']['PARVAL1']
            dx = MM.GW_build[name].model_mesh3D[0][0][0][1] - MM.GW_build[name].model_mesh3D[0][0][0][0]
            EGHBconductance = dx * MM.GW_build[name].parameters.param['EGHBcond']['PARVAL1']
            EastGHB += [[lay, row, col, EastGHBstage, EGHBconductance]]
    
    ghb[0] += EastGHB
    
    print "************************************************************************"
    print " Updating GHB boundary"
    
    MM.GW_build[name].boundaries.assign_boundary_array('GHB', ghb)
    
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
    
    
    
    #modflow_model.viewHeads()
    #modflow_model.viewHeads2()    

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
        data_folder = r"C:\Workspace\part0075\MDB modelling\testbox\PEST\\"    
        mf_exe_folder = r"C:\Workspace\part0075\GIT_REPOS\CampaspeModel\MODFLOW-NWT_64.exe"
        param_file = r"C:\Workspace\part0075\MDB modelling\testbox\PEST\parameters.txt"
    
    if param_file:
        run(model_folder, data_folder, mf_exe_folder, param_file=param_file)
    else:
        run(model_folder, data_folder, mf_exe_folder)