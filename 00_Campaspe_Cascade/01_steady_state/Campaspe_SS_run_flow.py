"""Update parameters and run Steady State flow model for Campaspe."""

import os
import sys
import shutil

import numpy as np

from HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface

# Configuration Loader
from HydroModelBuilder.Utilities.Config.ConfigLoader import ConfigLoader

def run(model_folder, data_folder, mf_exe_folder, param_file=None, verbose=True):
    """Model Runner."""

    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder, r"01_steady_state_packaged.pkl"))

    name = MM.GW_build.keys()[0]
    m = MM.GW_build[name]

    # Load in the new parameters based on parameters.txt or dictionary of new parameters

    if param_file is not None:
        m.updateModelParameters(os.path.join(data_folder, 'parameters.txt'), verbose=verbose)
    
    if verbose:
        print "************************************************************************"
        print " Updating HGU parameters "

    # This needs to be automatically generated from with the map_raster2mesh routine ...
    zone_map = {1: 'qa', 2: 'utb', 3: 'utqa', 4: 'utam', 5: 'utaf', 6: 'lta', 7: 'bse'}

    Zone = m.model_mesh3D[1].astype(float)
    Kh = m.model_mesh3D[1].astype(float)
    Kv = m.model_mesh3D[1].astype(float)
    Sy = m.model_mesh3D[1].astype(float)
    SS = m.model_mesh3D[1].astype(float)
    
    points_values_dict = {}
    for index, key in enumerate(zone_map.keys()):
        # if key ==7:
        #    continue
#        Kh[Zone == key] = m.parameters.param['Kh_' + zone_map[key]]['PARVAL1']
        for index2, param in enumerate(m.parameters.param_set['Kh_' + zone_map[key]]):
            if index2 == 0:
                points_values_dict[index] = [m.parameters.param[param]['PARVAL1']]
            else: 
                points_values_dict[index] += [m.parameters.param[param]['PARVAL1']]
            
        Kv[Zone == key] = m.parameters.param['Kv_' + zone_map[key]]['PARVAL1']
        Sy[Zone == key] = m.parameters.param['Sy_' + zone_map[key]]['PARVAL1']
        SS[Zone == key] = m.parameters.param['SS_' + zone_map[key]]['PARVAL1']

    hk = m.pilot_points['hk']
    zones = len(zone_map.keys())
    hk.output_directory = os.path.join(data_folder, 'pilot_points')
    hk.update_pilot_points_files_by_zones(zones, points_values_dict)
    import time
    time.sleep(3)
    hk.run_pyfac2real_by_zones(zones) 
    hk.save_mesh3D_array(filename=os.path.join(data_folder, 'hk_val_array'))
    Kh = hk.val_array
    m.save_array(os.path.join(data_folder, 'Kh'), Kh)
    
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
        print " Updating River Murray"
    
    mapped_river = m.polyline_mapped['River_Murray_model.shp']

    simple_river = []
    riv_width_avg = 10.0  # m
    riv_bed_thickness = 0.10  # m
    for riv_cell in mapped_river:  # m.polyline_mapped['Campaspe_Riv_model.shp']:
        row = riv_cell[0][0]
        col = riv_cell[0][1]
        if m.model_mesh3D[1][0][row][col] == -1:
            continue
        stage = m.model_mesh3D[0][0][row][col] - 0.01
        bed = m.model_mesh3D[0][0][row][col] - 0.1 - \
            m.parameters.param['RMstage']['PARVAL1']
        cond = riv_cell[1] * riv_width_avg * \
            m.parameters.param['Kv_RM']['PARVAL1'] / riv_bed_thickness
        simple_river += [[0, row, col, stage, cond, bed]]

    riv = {}
    riv[0] = simple_river
    m.boundaries.assign_boundary_array('Murray River', riv)
    
    if verbose:
        print "************************************************************************"
        print " Updating recharge boundary "

    # Adjust rainfall to recharge using 10% magic number
    interp_rain = np.copy(m.boundaries.bc['Rainfall']['bc_array'])

    for i in [1, 2, 3, 7]:
        interp_rain[m.model_mesh3D[1][0] == i] = interp_rain[m.model_mesh3D[
            1][0] == i] * m.parameters.param['SSrch_' + zone_map[i]]['PARVAL1']

    for i in [4, 5, 6, ]:
        interp_rain[m.model_mesh3D[1][0] == i] = interp_rain[
            m.model_mesh3D[1][0] == i] * 0.

    rch = {}
    rch[0] = interp_rain

    # m.boundaries.create_model_boundary_condition('Rain_reduced',
    # 'recharge', bc_static=True)
    m.boundaries.assign_boundary_array('Rain_reduced', rch)

    # print " Include irrigation in the recharge array"

    # print m.observations.obs_group['head'].keys()
    # print m.observations.obs_group['head']['mapped_observations']

    if verbose:
        print "************************************************************************"
        print " Updating River Murray GHB boundary"

    MurrayGHB = []
#    for MurrayGHB_cell in mapped_river:
#        row = MurrayGHB_cell[0][0]
#        col = MurrayGHB_cell[0][1]
#        # print m.model_mesh3D
#        for lay in range(m.model_mesh3D[1].shape[0]):
#            
#            if m.model_mesh3D[1][0][row][col] == -1:
#                continue
#            if lay == 0:
#                continue
#            
#            MurrayGHBstage = m.model_mesh3D[0][0][row][
#                col] + m.parameters.param['MGHB_stage']['PARVAL1']
#            
#            if MurrayGHBstage < m.model_mesh3D[0][lay + 1][row][col]:
#                continue
#
#            dx = m.gridHeight
#            dz = m.model_mesh3D[0][lay][row][col] - \
#                m.model_mesh3D[0][lay + 1][row][col]
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

    if verbose:
        print "************************************************************************"
        print " Updating GHB boundary"

    m.boundaries.assign_boundary_array('GHB', ghb)


    """
    Model is very sensitive to recharge conditions and in particular low recharge
    rates. In order to make sure the model converges, we can test if it failed
    after execution by calling checkConvergence fomr the flopyInterface object:
        modflow_model
    If this is the case we can start with  a higher recharge value and then 
    run, and if we hit convergence we can save the final heads and use as an 
    initial condition for a lower recharge rate which will hopefully allow the 
    model to run.
    """
    retries = 1
    
    for i in range(retries):
    
    
        if verbose:
            print "************************************************************************"
            print " Build and run MODFLOW model "
    
        #init_head = [400., 100., 100., 50.]
        #head = np.full(m.model_mesh3D[1].shape, init_head[i])

        #m.initial_conditions.set_as_initial_condition("Head", head)

        # Override temporal aspects of model build:
        modflow_model = flopyInterface.ModflowModel(m, data_folder=data_folder)

        if verbose:
            print "************************************************************************"
            print " Set initial head "
    
        if os.path.exists(os.path.join(data_folder, '01_steady_state.hds')):
            if verbose:
                print(" Trying last head solution as initial condition")

            filename = os.path.join(data_folder, '01_steady_state.hds')
            head = modflow_model.getFinalHeads(str(filename))
            modflow_model.strt = head
        else:
            if verbose:
                print(" No initial head file found")

                
        modflow_model.nper = 1  # This is the number of stress periods which is set to 1 here
        modflow_model.perlen = 40000 * 365  # This is the period of time which is set to 40000 yrs
        modflow_model.nstp = 1  # This is the number of sub-steps to do in each stress period
        modflow_model.steady = True #True #True # This is to tell FloPy that is a transient model
                
        modflow_model.executable = mf_exe_folder
    
        modflow_model.buildMODFLOW(transport=True)
    
        #modflow_model.checkMODFLOW()
    
        modflow_model.runMODFLOW(silent=True)
    
        converge = modflow_model.checkCovergence()

        if converge:
            break
        
        if verbose:
            print(" Trying new initial conditions")

    if not converge:
        retries = 5 #5
        if verbose:
            print(" Modifying recharge")

        for i in range(retries):
            
            rech_steps = [0.1, 0.08, 0.05, 0.02, 0.01]
            if verbose:
                print(" Trying recharge reduction of: {}".format(rech_steps[i]))
            interp_rain = np.copy(m.boundaries.bc['Rainfall']['bc_array'])
        
            for k in [1, 2, 3, 7]:
                interp_rain[m.model_mesh3D[1][0] == k] = interp_rain[m.model_mesh3D[
                    1][0] == k] * rech_steps[i]
        
            for k in [4, 5, 6, ]:
                interp_rain[m.model_mesh3D[1][0] == k] = interp_rain[
                    m.model_mesh3D[1][0] == k] * 0.
        
            rch = {}
            rch[0] = interp_rain
        
            m.boundaries.assign_boundary_array('Rain_reduced', rch)
    
            if i > 0:
                filename = os.path.join(data_folder, 'init{}'.format(rech_steps[i-1])) + '.hds'
                head = modflow_model.getFinalHeads(str(filename))
                modflow_model.strt = head

            #modflow_model = flopyInterface.ModflowModel(m, data_folder=data_folder)
            modflow_model.nper = 1  # This is the number of stress periods which is set to 1 here
            modflow_model.perlen = 40000 * 365  # This is the period of time which is set to 40000 yrs
            modflow_model.nstp = 100 #1 # This is the number of sub-steps to do in each stress period
            modflow_model.steady = False #True # This is to tell FloPy that is a transient model
                
            modflow_model.buildMODFLOW(transport=True)
        
            #modflow_model.checkMODFLOW()
        
            modflow_model.runMODFLOW(silent=True)
        
            converge = modflow_model.checkCovergence()
    
            if not converge:
                if verbose:
                    print("Absolute failure to find solution")
                #end if
                break
            else:
                if verbose:
                    print("Testing steady state solution from transient solution")
                #end if
                shutil.copy(os.path.join(modflow_model.data_folder, modflow_model.name + '.hds'), 
                            os.path.join(data_folder, 'init{}.hds'.format(rech_steps[i])))
                
#        if converge:
                interp_rain = np.copy(m.boundaries.bc['Rainfall']['bc_array'])
            
                for k in [1, 2, 3, 7]:
                    interp_rain[m.model_mesh3D[1][0] == k] = interp_rain[m.model_mesh3D[
                        1][0] == k] * m.parameters.param['SSrch_' + zone_map[k]]['PARVAL1']
            
                for k in [4, 5, 6, ]:
                    interp_rain[m.model_mesh3D[1][0] == k] = interp_rain[
                        m.model_mesh3D[1][0] == k] * 0.
            
                rch = {}
                rch[0] = interp_rain
            
                m.boundaries.assign_boundary_array('Rain_reduced', rch)
        
                modflow_model.nper = 1  # This is the number of stress periods which is set to 1 here
                modflow_model.perlen = 40000 * 365  # This is the period of time which is set to 40000 yrs
                modflow_model.nstp = 1  # This is the number of sub-steps to do in each stress period
                modflow_model.steady = True #True # This is to tell FloPy that is a transient model
        
                filename = os.path.join(data_folder, 'init{}.hds'.format(rech_steps[i]))
                head = modflow_model.getFinalHeads(str(filename))
                modflow_model.strt = head
                
                modflow_model.buildMODFLOW(transport=True)
                modflow_model.runMODFLOW(silent=True)
                converge = modflow_model.checkCovergence()            
                if converge:
                    shutil.copy(os.path.join(modflow_model.data_folder, modflow_model.name + '.hds'), 
                                data_folder)
                    break
                
    #modflow_model.waterBalance()
    #modflow_model.viewGHB()
    #modflow_model.viewHeads()
    #modflow_model.viewHeadLayer(figsize=(20,10))

    #riv_exch = modflow_model.getRiverFlux('Campaspe River')
    #for key in riv_exch.keys():
    #    print 'Campaspe River net flux: ' + str(round(sum([x[0] for x in riv_exch[key]]))) + ' m3/d'

    #riv_exch = modflow_model.getRiverFlux('Murray River')
    #for key in riv_exch.keys():
    #    print 'Murray River net flux: ' + str(round(sum([x[0] for x in riv_exch[key]]))) + ' m3/d'
  

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
        # Get general model config information
        CONFIG = ConfigLoader('../../config/model_config.json')\
                        .set_environment("01_steady_state")
        model_config = CONFIG.model_config
        model_folder = model_config['model_folder'] + model_config['grid_resolution'] + os.path.sep
        data_folder = model_config['data_folder']
        mf_exe_folder = model_config['mf_exe_folder']
        param_file = model_config['param_file']

    if param_file:
        run(model_folder, data_folder, mf_exe_folder, param_file=param_file, verbose=verbose)
    else:
        run(model_folder, data_folder, mf_exe_folder, verbose=verbose)
