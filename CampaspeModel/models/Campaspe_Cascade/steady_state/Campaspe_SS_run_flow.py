"""Update parameters and run Steady State flow model for Campaspe."""

import os
import sys
import shutil
import time

import numpy as np

from HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface

# Configuration Loader
from HydroModelBuilder.Utilities.Config.ConfigLoader import ConfigLoader

def run(model_folder, data_folder, mf_exe_folder, param_file="", verbose=False,
        plots=False):
    """Model Runner."""

    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder, r"01_steady_state_packaged.pkl"))

    name = MM.GW_build.keys()[0]
    m = MM.GW_build[name]

    # Load in the new parameters based on parameters.txt or dictionary of new parameters

    if param_file != "":
        m.updateModelParameters(os.path.join(data_folder, 'parameters.txt'), verbose=verbose)
    
    if verbose:
        print "************************************************************************"
        print " Updating HGU parameters "

    # This needs to be automatically generated with the map_raster2mesh routine ...
    zone_map = {1: 'qa', 2: 'utb', 3: 'utqa', 4: 'utam', 5: 'utaf', 6: 'lta', 7: 'bse'}

    default_array = m.model_mesh3D[1].astype(float)
    Zone = np.copy(default_array)
    Kh = np.copy(default_array)
    Kv = np.copy(default_array)
    Sy = np.copy(default_array)
    SS = np.copy(default_array)
    
    def create_pp_points_dict(zone_map, Zone, prop_array, prop_name, m):
        points_values_dict = {}
        for index, key in enumerate(zone_map.keys()):
            for index2, param in enumerate(m.parameters.param_set[prop_name + zone_map[key]]):
                if index2 == 0:
                    points_values_dict[index] = [m.parameters.param[param]['PARVAL1']]
                else: 
                    points_values_dict[index] += [m.parameters.param[param]['PARVAL1']]
        return points_values_dict    

    def update_pilot_points(zone_map, Zone, prop_array, par_name, prop_name, prop_folder, m, prop_array_fname):
        points_values_dict = create_pp_points_dict(zone_map, Zone, prop_array, prop_name, m)
        p = m.pilot_points[par_name]
        zones = len(zone_map.keys())
        p.output_directory = os.path.join(data_folder, prop_folder)
        p.update_pilot_points_files_by_zones(zones, points_values_dict)
        time.sleep(3)
        p.run_pyfac2real_by_zones(zones) 
        #p.save_mesh3D_array(filename=os.path.join(data_folder, prop_array_fname))
        return p.val_array

    Kh = update_pilot_points(zone_map, Zone, Kh, 'hk', 'kh', 'hk_pilot_points',
                             m, 'hk_val_array')  
    m.save_array(os.path.join(data_folder, 'Kh'), Kh)

    print("Erroneous K pilot cells: {}".format(len(Kh[Kh > 10000.])))
    Kh[Kh > 10000.] = 25.
    Kv = Kh * 0.1
    m.save_array(os.path.join(data_folder, 'Kv'), Kv)

    Sy = update_pilot_points(zone_map, Zone, Sy, 'sy', 'sy', 'sy_pilot_points',
                             m, 'sy_val_array')
    print("Erroneous Sy pilot cells: {}".format(len(Sy[Sy > 0.5])))
    Sy[Sy > 0.5] = 0.25
    m.save_array(os.path.join(data_folder, 'Sy'), Sy)
    
    SS = update_pilot_points(zone_map, Zone, SS, 'ss', 'ss', 'ss_pilot_points',
                             m, 'ss_val_array')
    print("Erroneous Ss pilot cells: {}".format(len(SS[SS > 0.01])))
    SS[SS > 0.01] = 1E-5
    m.save_array(os.path.join(data_folder, 'SS'), SS)


    m.properties.update_model_properties('Kh', Kh)
    m.properties.update_model_properties('Kv', Kv)
    m.properties.update_model_properties('Sy', Sy)
    m.properties.update_model_properties('SS', SS)



    if verbose:
        print "************************************************************************"
        print " Updating river parameters "

    reach_df = m.mf_sfr_df['Campaspe']['reach_df'] 
    segment_data = m.mf_sfr_df['Campaspe']['seg_df']

    num_reaches = m.pilot_points['Campaspe'].num_points #4
    known_points = m.pilot_points['Campaspe'].points
    # Create reach data
    river_seg = m.river_mapping['Campaspe']
    
    strcond_val = [m.parameters.param['kv_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #@ TODO: Check this next part to make sure that zero conductance reaches 
    #@ that are a result of params won't be ignored!!
    river_seg.loc[river_seg['strhc1'] != 0.0, 'strhc1'] = np.interp(
            river_seg[river_seg['strhc1'] != 0.0]['Cumulative Length'].tolist(), 
            known_points, strcond_val)
    
    strthick_val = [m.parameters.param['bedthck{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    river_seg['strthick'] = np.interp(river_seg['Cumulative Length'].tolist(), known_points, strthick_val)
    
    
    reach_df = river_seg[['k','i','j','iseg','ireach','rchlen','strtop','slope','strthick','strhc1']]
    reach_data = reach_df.to_records(index=False)
    
    # Set the roughness for the channel
    roughch_val = [m.parameters.param['mn_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    roughch = np.interp(river_seg['Cumulative Length'].tolist(), 
                                    known_points, roughch_val)
    # Set the roughness for the banks
    roughbk_val = [m.parameters.param['mn_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    roughbk = np.interp(river_seg['Cumulative Length'].tolist(), 
                                    known_points, roughbk_val)
    river_seg['roughch'] = roughch
    river_seg['roughbk'] = roughbk
    
    
    width1_val = [m.parameters.param['rivwdth{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    width1 = np.interp(river_seg['Cumulative Length'].tolist(), 
                                    known_points, width1_val)
    segment_data['width2'] = segment_data['width1'] = width1
    
    segment_data1 = segment_data.to_records(index=False)
    seg_dict = {0: segment_data1}
    
    if verbose:
        print "************************************************************************"
        print " Updating Campaspe river boundary"
    
    m.boundaries.update_boundary_array('Campaspe River', [reach_data, seg_dict])

    if verbose:
        print "************************************************************************"
        print " Updating River Murray"

    mriver_seg = m.river_mapping['Murray']
    mriver_seg['strhc1'] = m.parameters.param['kv_rm']['PARVAL1']                      
    
    simple_river = []
    for row in mriver_seg.iterrows():
        row = row[1]
        simple_river += [[row['k'], row['i'], row['j'], row['stage'], \
                          row['strhc1'] * row['rchlen'] * row['width1'], row['strtop']]]
    
    riv = {}
    riv[0] = simple_river
    m.boundaries.update_boundary_array('Murray River', riv)
    
    if verbose:
        print "************************************************************************"
        print " Updating recharge boundary "

    # Adjust rainfall to recharge using zoned rainfall reduction parameters
    interp_rain = np.copy(m.boundaries.bc['Rainfall']['bc_array'])
    recharge_zone_array = m.boundaries.bc['Rain_reduced']['zonal_array']
    rch_zone_dict = m.boundaries.bc['Rain_reduced']['zonal_dict']

    rch_zones = len(rch_zone_dict.keys())

    par_rech_vals = [m.parameters.param['ssrch{}'.format(i)]['PARVAL1'] \
                     for i in range(rch_zones - 1)]

    def update_recharge(vals):
        for i in range(rch_zones - 1):
            interp_rain[recharge_zone_array == rch_zone_dict[i+1]] = \
                interp_rain[recharge_zone_array == rch_zone_dict[i+1]] * \
                vals[i]
        
        interp_rain[recharge_zone_array==rch_zone_dict[0]] = interp_rain[recharge_zone_array == rch_zone_dict[0]] * 0.0
        interp_rain[m.model_mesh3D[1][0] == -1] = 0.
        return interp_rain

    rch = {}
    rch[0] = update_recharge(par_rech_vals)

    m.boundaries.update_boundary_array('Rain_reduced', rch)


    if verbose:
        print "************************************************************************"
        print " Updating River Murray GHB boundary"

    MurrayGHB = []

    MurrayGHB_cells = [[x[0], x[1], x[2], x[3]] for x in m.boundaries.bc['GHB']['bc_array'][0]]
    for MurrayGHB_cell in MurrayGHB_cells:
        lay, row, col = MurrayGHB_cell[:3]
        MurrayGHBstage = MurrayGHB_cell[3] #m.parameters.param['mghb_stage']['PARVAL1']
        dx = m.gridHeight
        dz = m.model_mesh3D[0][lay][row][col] - \
            m.model_mesh3D[0][lay + 1][row][col]
        MGHBconductance = dx * dz * m.parameters.param['mghbk']['PARVAL1'] #/ 10000.
        MurrayGHB += [[lay, row, col, MurrayGHBstage, MGHBconductance]]

    ghb = {}
    ghb[0] = MurrayGHB

    if verbose:
        print "************************************************************************"
        print " Updating GHB boundary"

    m.boundaries.update_boundary_array('GHB', ghb)


    if verbose:
        print "************************************************************************"
        print " Check for boundary condition updating"
        m.generate_update_report()


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
    

        # Override temporal aspects of model build:
        modflow_model = flopyInterface.ModflowModel(m, data_folder=os.path.join(data_folder, "model_" + m.name))

        if verbose:
            print "************************************************************************"
            print " Set initial head "
    
        if os.path.exists(os.path.join(data_folder, '01_steady_state.hds')):
            if verbose:
                print(" Trying last head solution as initial condition")

            filename = os.path.join(data_folder, '01_steady_state.hds')
            head = modflow_model.get_final_heads(str(filename))
            modflow_model.strt = head
        else:
            if verbose:
                print(" No initial head file found")

                
        modflow_model.nper = 1  # This is the number of stress periods which is set to 1 here
        modflow_model.perlen = 40000 * 365  # This is the period of time which is set to 40000 yrs
        modflow_model.nstp = 1  # This is the number of sub-steps to do in each stress period
        modflow_model.steady = True #True #True # This is to tell FloPy that is a transient model
                
        modflow_model.executable = mf_exe_folder
    
        modflow_model.buildMODFLOW(transport=True, verbose=verbose, check=True)

        converge = modflow_model.runMODFLOW(silent=~verbose)

       # return
    
        if converge:
            break
        
        if verbose:
            print(" Trying new initial conditions")

    if not converge:
        retries = 1 #5
        if verbose:
            print(" Modifying recharge")

        for i in range(retries):
            
            rech_steps = [0.01] # [0.1, 0.08, 0.05, 0.02, 0.01]
            if verbose:
                print(" Trying recharge reduction of: {}".format(rech_steps[i]))
            interp_rain = np.copy(m.boundaries.bc['Rainfall']['bc_array'])
        
            rch = {}
            rch[0] = update_recharge([rech_steps[i]] * (rch_zones - 1))
        
            m.boundaries.update_boundary_array('Rain_reduced', rch)
    
            if i > 0:
                filename = os.path.join(data_folder, 'init{}'.format(rech_steps[i-1])) + '.hds'
                head = modflow_model.get_final_heads(str(filename))
                modflow_model.strt = head

            modflow_model.nper = 1  # This is the number of stress periods which is set to 1 here
            modflow_model.perlen = 40000 * 365  # This is the period of time which is set to 40000 yrs
            modflow_model.nstp = 100 # This is the number of sub-steps to do in each stress period
            modflow_model.steady = False # This is to tell FloPy that is a transient model
            modflow_model.headtol = 1 #1E-3 # Initially relax head convergence criteria to get convergence
            modflow_model.fluxtol = 1E6 # 
                
            modflow_model.buildMODFLOW(transport=True)
        
            converge = modflow_model.runMODFLOW(silent=False)
        
            print("Convergence?: {}".format(converge))
        
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
                
                interp_rain = np.copy(m.boundaries.bc['Rainfall']['bc_array'])
            
                rch = {}
                rch[0] = update_recharge(par_rech_vals)
            
                m.boundaries.update_boundary_array('Rain_reduced', rch)
        
                modflow_model.nper = 1  # This is the number of stress periods which is set to 1 here
                modflow_model.perlen = 40000 * 365  # This is the period of time which is set to 40000 yrs
                modflow_model.nstp = 1  # This is the number of sub-steps to do in each stress period
                modflow_model.steady = True # This is to tell FloPy this is a transient model

                filename = os.path.join(data_folder, 'init{}.hds'.format(rech_steps[i]))
                head = modflow_model.get_final_heads(str(filename))
                modflow_model.strt = head

                for converge_criterion in [1E-4]:
                    print("Testing steady state solution with convergence criteria set to {}".format(converge_criterion))
                    modflow_model.headtol = converge_criterion # Initially relax head convergence criteria to get convergence
                    modflow_model.buildMODFLOW(transport=True)
                    converge = modflow_model.runMODFLOW(silent=True)
                    #break
                    if not converge:
                        if verbose:
                            print("Absolute failure to find solution")
                        #end if
                        break
                    filename = os.path.join(modflow_model.data_folder, modflow_model.name + '.hds')
                    head = modflow_model.get_final_heads(str(filename))
                    modflow_model.strt = head

                if converge:
                    shutil.copy(os.path.join(modflow_model.data_folder, modflow_model.name + '.hds'), 
                                data_folder)
                    break
                
    if converge:
        if plots:
            modflow_model.waterBalance(1)
            modflow_model.viewGHB()
            modflow_model.viewHeads()
            modflow_model.viewHeads2()
            #modflow_model.viewHeadsByZone()
            #modflow_model.viewHeadLayer(figsize=(20,10))
    
    return modflow_model

if __name__ == "__main__":

    verbose = False
    plots = False
                    
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
                        .set_environment("01_steady_state")
        model_config = CONFIG.model_config
        model_folder = model_config['model_folder'] + model_config['grid_resolution']
        data_folder = model_config['data_folder']
        mf_exe_folder = model_config['mf_exe_folder']
        param_file = model_config['param_file']

    if param_file:
        run = run(model_folder, data_folder, mf_exe_folder, param_file=param_file, verbose=verbose, plots=verbose)
    else:
        run = run(model_folder, data_folder, mf_exe_folder, verbose=verbose, plots=verbose)
