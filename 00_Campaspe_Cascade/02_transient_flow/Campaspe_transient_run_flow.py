import sys
import os
sys.path.append('C:\Workspace\part0075\GIT_REPOS')

import numpy as np
import flopy.utils.binaryfile as bf

from HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface

# Configuration Loader
from HydroModelBuilder.Utilities.Config.ConfigLoader import ConfigLoader


def run(model_folder, data_folder, mf_exe, param_file=None, verbose=True):
    
    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder, "02_transient_flow_packaged.pkl"))
    
    name = MM.GW_build.keys()[0]
    m = MM.GW_build[name]
    
    # Load in the new parameters based on parameters.txt or dictionary of new parameters
 
    if param_file:
        m.updateModelParameters(os.path.join(data_folder, 'parameters.txt'), verbose=verbose)
    
    if verbose:
        print "************************************************************************"
        print " Updating HGU parameters "
    
    # This needs to be automatically generated from with the map_raster2mesh routine ...
    zone_map = {1:'qa', 2:'utb', 3:'utqa', 4:'utam', 5:'utaf', 6:'lta', 7:'bse'}
    
    Zone = m.model_mesh3D[1].astype(float)
    #Kh = m.model_mesh3D[1].astype(float)
    Kv = m.model_mesh3D[1].astype(float)
    Sy = m.model_mesh3D[1].astype(float)
    SS = m.model_mesh3D[1].astype(float)

    #points_values_dict = {}
    for index, key in enumerate(zone_map.keys()):
        # if key ==7:
        #    continue
        #Kh[Zone == key] = m.parameters.param['Kh_' + zone_map[key]]['PARVAL1']
        #for index2, param in enumerate(m.parameters.param_set['Kh_' + zone_map[key]]):
        #    if index2 == 0:
        #        points_values_dict[index] = [m.parameters.param[param]['PARVAL1']]
        #    else: 
        #        points_values_dict[index] += [m.parameters.param[param]['PARVAL1']]
            
        Kv[Zone == key] = m.parameters.param['Kv_' + zone_map[key]]['PARVAL1']
        Sy[Zone == key] = m.parameters.param['Sy_' + zone_map[key]]['PARVAL1']
        SS[Zone == key] = m.parameters.param['SS_' + zone_map[key]]['PARVAL1']

    Kh = m.load_array(os.path.join(data_folder, 'hk_val_array.npy'))    

    m.properties.assign_model_properties('Kh', Kh)
    m.properties.assign_model_properties('Kv', Kv)
    m.properties.assign_model_properties('Sy', Sy)
    m.properties.assign_model_properties('SS', SS)

    if verbose:
        print "************************************************************************"
        print " Updating river parameters "
    
    riv_width_avg = 10.0 #m
    riv_bed_thickness = 0.10 #m
    
    cond = []
    for index, riv_cell in enumerate(m.polyline_mapped['Campaspe_Riv_model.shp']):
        row = riv_cell[0][0]
        col = riv_cell[0][1]
        if m.model_mesh3D[1][0][row][col] == -1:
            continue
        cond += [riv_cell[1] * riv_width_avg * m.parameters.param['Kv_riv']['PARVAL1'] / riv_bed_thickness]

    riv = m.boundaries.bc['Campaspe River']['bc_array'].copy()
    for key in riv.keys():
        riv[key] = [[x[0], x[1], x[2], x[3], cond[ind], x[5]] for ind, x in enumerate(riv[key])]
            
    m.boundaries.assign_boundary_array('Campaspe River', riv)
    
    if verbose:
        print "************************************************************************"
        print " Updating Murray River boundary"
    
    simple_river = []
    riv_width_avg = 10.0 #m
    riv_bed_thickness = 0.10 #m
#    for riv_cell in mapped_river: #m.polyline_mapped['Campaspe_Riv_model.shp']:
#        row = riv_cell[0][0]
#        col = riv_cell[0][1]
#        if m.model_mesh3D[1][0][row][col] == -1:
#            continue
#        stage = m.model_mesh3D[0][0][row][col]
#        bed = m.model_mesh3D[0][0][row][col] - m.parameters.param['RMstage']['PARVAL1']
#        cond = riv_cell[1] * riv_width_avg * m.parameters.param['Kv_RM']['PARVAL1'] / riv_bed_thickness
#        simple_river += [[0, row, col, stage, cond, bed]]

    riv_cells = [[x[1], x[2], x[3]] for x in m.boundaries.bc['Murray River']['bc_array'][0]]
    for riv_cell in riv_cells:
        row, col, stage = riv_cell
        bed = m.model_mesh3D[0][0][row][col] - m.parameters.param['RMstage']['PARVAL1']
        if bed < m.model_mesh3D[0][1][row][col]:
            bed = m.model_mesh3D[0][0][row][col] + 0.01
        #end if
        cond = riv_cell[1] * riv_width_avg * m.parameters.param['Kv_RM']['PARVAL1'] / riv_bed_thickness
        simple_river += [[0, row, col, stage, cond, bed]]

    
    riv = {}
    riv[0] = simple_river
    #m.boundaries.create_model_boundary_condition('Campaspe River', 'river', bc_static=True)
    m.boundaries.assign_boundary_array('Murray River', riv)
    
    
    if verbose:
        print "************************************************************************"
        print " Updating recharge boundary "

    interp_rain = m.boundaries.bc['Rainfall']['bc_array']
    # Adjust rainfall to recharge using rainfall reduction
    for i in [1,2,3,7]:
        for key in interp_rain.keys():
            interp_rain[key][m.model_mesh3D[1][0]==i] = interp_rain[key][m.model_mesh3D[1][0]==i] * m.parameters.param['rch_red_'+zone_map[i]]['PARVAL1']

    for i in [4,5,6,]:
        for key in interp_rain.keys():
            interp_rain[key][m.model_mesh3D[1][0]==i] = interp_rain[key][m.model_mesh3D[1][0]==i] * 0.
    
    rch = {}
    for key in m.boundaries.bc['Rainfall']['bc_array'].keys():
        rch[key] = interp_rain[key] 
    
    m.boundaries.assign_boundary_array('Rain_reduced', rch)
    
    #print " Include irrigation in the recharge array"
    
    if verbose:
        print "************************************************************************"
        print " Updating Murray River GHB boundary"
    
    MurrayGHB = []
#    for MurrayGHB_cell in mapped_river:
#        row = MurrayGHB_cell[0][0]
#        col = MurrayGHB_cell[0][1]
#        #print m.model_mesh3D
#        for lay in range(m.model_mesh3D[1].shape[0]):    
#            if m.model_mesh3D[1][0][row][col] == -1:
#                continue
#            MurrayGHBstage = m.model_mesh3D[0][0][row][col] + m.parameters.param['MGHB_stage']['PARVAL1']
#            dx = m.gridHeight
#            dz = m.model_mesh3D[0][lay][row][col] - m.model_mesh3D[0][lay+1][row][col]
#            MGHBconductance = dx * dz * m.parameters.param['MGHBcond']['PARVAL1']
#            MurrayGHB += [[lay, row, col, MurrayGHBstage, MGHBconductance]]

    MurrayGHB_cells = [[x[0], x[1], x[2]] for x in m.boundaries.bc['GHB']['bc_array'][0]]
    for MurrayGHB_cell in MurrayGHB_cells:
        lay, row, col = MurrayGHB_cell
        MurrayGHBstage = m.model_mesh3D[0][0][row][
            col] + m.parameters.param['MGHB_stage']['PARVAL1']
        if MurrayGHBstage < m.model_mesh3D[0][lay + 1][row][col]:
            continue
        dx = m.gridHeight
        dz = m.model_mesh3D[0][lay][row][col] - \
            m.model_mesh3D[0][lay + 1][row][col]
        MGHBconductance = dx * dz * m.parameters.param['MGHBcond']['PARVAL1']
        MurrayGHB += [[lay, row, col, MurrayGHBstage, MGHBconductance]]
    
    ghb = {}
    ghb[0] = MurrayGHB

    
#    print "************************************************************************"
#    print " Updating Western GHB boundary"
#    
#    mapped_west = m.load_obj(os.path.join(model_folder, r"western_head_model.shp_mapped.pkl"))
#    
#    WestGHB = []
#    for WestGHB_cell in mapped_west:
#        row = WestGHB_cell[0][0]
#        col = WestGHB_cell[0][1]
#        for lay in range(m.model_mesh3D[1].shape[0]):    
#            if m.model_mesh3D[1][0][row][col] == -1:
#                continue
#            WestGHBstage = m.model_mesh3D[0][0][row][col] + m.parameters.param['WGHB_stage']['PARVAL1']
#            dx = m.gridHeight
#            dz = m.model_mesh3D[0][lay][row][col] - m.model_mesh3D[0][lay+1][row][col]
#            WGHBconductance = dx * dz * m.parameters.param['WGHBcond']['PARVAL1']
#            WestGHB += [[lay, row, col, WestGHBstage, WGHBconductance]]
#    
#    ghb[0] += WestGHB
#    
#    print "************************************************************************"
#    print " Updating Eastern GHB boundary"
#    
#    mapped_east = m.load_obj(os.path.join(model_folder, r"eastern_head_model.shp_mapped.pkl"))
#    
#    EastGHB = []
#    for EastGHB_cell in mapped_east:
#        row = EastGHB_cell[0][0]
#        col = EastGHB_cell[0][1]
#        for lay in range(m.model_mesh3D[1].shape[0]):    
#            if m.model_mesh3D[1][0][row][col] == -1:
#                continue
#            EastGHBstage = m.model_mesh3D[0][0][row][col] + m.parameters.param['EGHB_stage']['PARVAL1']
#            dx = m.gridHeight
#            dz = m.model_mesh3D[0][lay][row][col] - m.model_mesh3D[0][lay+1][row][col]
#            EGHBconductance = dx * dz * m.parameters.param['EGHBcond']['PARVAL1']
#            EastGHB += [[lay, row, col, EastGHBstage, EGHBconductance]]
#    
#    ghb[0] += EastGHB
    
    if verbose:
        print "************************************************************************"
        print " Updating GHB boundary"
    
    m.boundaries.assign_boundary_array('GHB', ghb)
    
    if verbose:
        print "************************************************************************"
        print " Updating Drains boundary"
    
    mapped_drains = m.polyline_mapped['Drain_Clip_model.shp']
    
    simple_drain = []
    drain_width_avg = 3.0 #m
    drain_bed_thickness = 0.10 #m
    for drain_cell in mapped_drains:
        row = drain_cell[0][0]
        col = drain_cell[0][1]
        if m.model_mesh3D[1][0][row][col] == -1:
            continue
        #print m.model_mesh3D
        drain_bed = m.model_mesh3D[0][0][row][col] - m.parameters.param['drain_drop']['PARVAL1']
        drain_cond = drain_cell[1] * drain_width_avg * m.parameters.param['Kv_drain']['PARVAL1'] / drain_bed_thickness
        simple_drain += [[0, row, col, drain_bed, drain_cond]]
    
    drain = {}
    drain[0] = simple_drain
    
    m.boundaries.assign_boundary_array('Drain', drain)
    
    if verbose:
        print "************************************************************************"
        print " Updating Channels boundary"
    
    simple_channel = []
    channel_width_avg = 10.0 #m
    channel_bed_thickness = 0.10 #m
    for channel_cell in m.polyline_mapped['Channel_Clip_model.shp']:
        row = channel_cell[0][0]
        col = channel_cell[0][1]
        if m.model_mesh3D[1][0][row][col] == -1:
            continue
        #print m.model_mesh3D
        channel_stage = m.model_mesh3D[0][0][row][col]
        channel_bed = m.model_mesh3D[0][0][row][col] - m.parameters.param['chan_drop']['PARVAL1']
        channel_cond = channel_cell[1] * channel_width_avg * m.parameters.param['Kv_chan']['PARVAL1'] / channel_bed_thickness
        simple_channel += [[0, row, col, channel_stage, channel_cond, channel_bed]]
    
    channel = {}
    channel[0] = simple_channel
    
    m.boundaries.assign_boundary_array('Channel', channel)

    if verbose:
        print "************************************************************************"
        print " Updating pumping boundary"

    #pumpy = m.boundaries.bc['licenced_wells']['bc_array']
    #wel = {key: [[b[0], b[1], b[2], b[3] * pumping] for b in a] for key, a in pumpy.iteritems()}

    #m.boundaries.assign_boundary_array('licenced_wells', wel)
    
    if verbose:
        print "************************************************************************"
        print " Set initial head "

    path=os.path.join(data_folder,"model_01_steady_state")
    fname="01_steady_state"
    headobj = bf.HeadFile(os.path.join(path, fname + '.hds'))
    times = headobj.get_times()        
    head = headobj.get_data(totim=times[-1])
    
    m.initial_conditions.set_as_initial_condition("Head", head)
    
    if verbose:
        print "************************************************************************"
        print " Build and run MODFLOW model "
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ## Currently using flopyInterface directly rather than running from the ModelManager ...

    modflow_model = flopyInterface.ModflowModel(m, data_folder=data_folder)

    modflow_model.executable = mf_exe_folder

    modflow_model.buildMODFLOW(transport=True, write=True)

    #modflow_model.checkMODFLOW()  # Note this is slow af, so use only in setting up model

    modflow_model.runMODFLOW(silent=True)

    converge = modflow_model.checkCovergence()

    if converge:
        if verbose:
            print('model converged')
            #break
    
        modflow_model.writeObservations()

    Campaspe_riv_flux = modflow_model.getRiverFlux('Campaspe River')
    net_riv_flux =  np.sum(np.array([x[0] for x in Campaspe_riv_flux[max(Campaspe_riv_flux.keys())]]))
    
    with open(modflow_model.data_folder + os.path.sep + 'observations_net_riv_flux.txt', 'w') as f:
        f.write('%f\n' %net_riv_flux)  
        
    #modflow_model.compareAllObs()
        
    #if ss_converge:
    #    if tr_converge:
    #        modflow_model.writeObservations()


    return modflow_model

if __name__ == "__main__":

    # Get general model config information
    CONFIG = ConfigLoader('../../config/model_config.json')\
                    .set_environment("02_transient_flow")

    verbose = False
                    
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
        run = run(model_folder, data_folder, mf_exe_folder, param_file=param_file, verbose=verbose)
    else:
        run = run(model_folder, data_folder, mf_exe_folder, verbose=verbose)
