import sys
import os
import time

import numpy as np
import flopy.utils.binaryfile as bf

from HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface

# Configuration Loader
from HydroModelBuilder.Utilities.Config.ConfigLoader import ConfigLoader


def run(model_folder, data_folder, mf_exe, param_file="", verbose=True):

    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder, "03_recharge_flow_packaged.pkl"))
    name = list(MM.GW_build.keys())[0]
    m = MM.GW_build[name]
    
    # Load in the new parameters based on parameters.txt or dictionary of new parameters
 
    #if param_file != "":
    #    m.updateModelParameters(os.path.join(data_folder, 'parameters.txt'), verbose=verbose)
    
    if verbose:
        print("************************************************************************")
        print(" Updating HGU parameters ")
    
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
        zones = len(list(zone_map.keys()))
        p.output_directory = os.path.join(data_folder, prop_folder)
        p.update_pilot_points_files_by_zones(zones, points_values_dict)
        p.run_pyfac2real_by_zones(zones) 
        #p.save_mesh3D_array(filename=os.path.join(data_folder, prop_array_fname))
        return p.val_array

    Kh = update_pilot_points(zone_map, Zone, Kh, 'hk', 'kh_', 'hk_pilot_points',
                             m, 'hk_val_array')  
    m.save_array(os.path.join(data_folder, 'Kh'), Kh)

    print(("Erroneous pilot cells: {}".format(len(Kh[Kh > 10000.]))))
    Kh[Kh > 10000.] = 25.
    Kv = Kh * 0.1
    m.save_array(os.path.join(data_folder, 'Kv'), Kv)

    Sy = update_pilot_points(zone_map, Zone, Sy, 'sy', 'sy_', 'sy_pilot_points',
                             m, 'sy_val_array')
    Sy[Sy > 1] = 0.25
    m.save_array(os.path.join(data_folder, 'Sy'), Sy)
    
    SS = update_pilot_points(zone_map, Zone, SS, 'ss', 'ss_', 'ss_pilot_points',
                             m, 'ss_val_array')
    SS[SS > 1] = 1E-5
    m.save_array(os.path.join(data_folder, 'SS'), SS)

    m.properties.update_model_properties('Kh', Kh)
    m.properties.update_model_properties('Kv', Kv)
    m.properties.update_model_properties('Sy', Sy)
    m.properties.update_model_properties('SS', SS)
    
    if verbose:
        print("************************************************************************")
        print(" Updating river parameters ")
    
#    reach_df = m.mf_sfr_df['Campaspe'].reach_df 
#    segment_data = m.mf_sfr_df['Campaspe'].seg_df
#
#    num_reaches = m.pilot_points['Campaspe'].num_points #4
#    known_points = m.pilot_points['Campaspe'].points
#    # Create reach data
#    river_seg = m.river_mapping['Campaspe']
#    
#    strcond_val = [m.parameters.param['kv_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
#    river_seg.loc[river_seg['strhc1'] != 0.0, 'strhc1'] = np.interp(
#            river_seg[river_seg['strhc1'] != 0.0]['Cumulative Length'].tolist(), 
#            known_points, strcond_val)
#    
#    strthick_val = [m.parameters.param['bedthck{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
#    river_seg['strthick'] = np.interp(river_seg['Cumulative Length'].tolist(), known_points, strthick_val)
#    
#    reach_df = river_seg[['k','i','j','iseg','ireach','rchlen','strtop','slope','strthick','strhc1']]
#    reach_data = reach_df.to_records(index=False)
#    
#    # Set the roughness for the channel
#    roughch_val = [m.parameters.param['mn_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
#    roughch = np.interp(river_seg['Cumulative Length'].tolist(), 
#                                    known_points, roughch_val)
#    # Set the roughness for the banks
#    roughbk_val = [m.parameters.param['mn_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
#    roughbk = np.interp(river_seg['Cumulative Length'].tolist(), 
#                                    known_points, roughbk_val)
#    river_seg['roughch'] = roughch
#    river_seg['roughbk'] = roughbk
#    
#    width1_val = [m.parameters.param['rivwdth{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
#    width1 = np.interp(river_seg['Cumulative Length'].tolist(), 
#                                    known_points, width1_val)
#
#    segment_data1 = {}
#    for key in segment_data.keys():
#        segment_data[key]['width2'] = segment_data[key]['width1'] = width1
#        segment_data1[key] = segment_data[key].to_records(index=False)
#    
#    seg_dict = segment_data1    

    if verbose:
        print("************************************************************************")
        print(" Updating Campaspe river boundary")
    
#    m.boundaries.update_boundary_array('Campaspe River', [reach_data, seg_dict])
    
    if verbose:
        print("************************************************************************")
        print(" Updating Murray River boundary")
    
#    mriver_seg = m.river_mapping['Murray']
#    mriver_seg['strhc1'] = m.parameters.param['kv_rm']['PARVAL1']                      
#    
#    simple_river = []
#    for row in mriver_seg.iterrows():
#        row = row[1]
#        simple_river += [[row['k'], row['i'], row['j'], row['stage'], \
#                          row['strhc1'] * row['rchlen'] * row['width1'], row['strtop']]]
#    
#    riv = {}
#    riv[0] = simple_river
#    m.boundaries.update_boundary_array('Murray River', riv)
    
    if verbose:
        print("************************************************************************")
        print(" Updating recharge boundary ")

    # Adjust rainfall to recharge using zoned rainfall reduction parameters
    # Need to make copies of all rainfall arrays
    interp_rain = m.boundaries.bc['Rainfall']['bc_array']
    for key in list(interp_rain.keys()):
        interp_rain[key] = np.copy(np.zeros_like(interp_rain[key]))

    recharge_zone_array = m.boundaries.bc['Rain_reduced']['zonal_array']
    rch_zone_dict = m.boundaries.bc['Rain_reduced']['zonal_dict']

    par_rech_file = os.path.join(data_folder, 'parameters_rech.txt')
    
    rech_dict = {}
    with open(par_rech_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            key, val = line.split()
            rech_dict[key] = float(val) 

    import matplotlib.pyplot as plt
    rech_zone_array_cp = np.copy(recharge_zone_array)
    rech_zone_array_cp[m.model_mesh3D[1][0] == -1] = 0
    rech_zone_array_cp[rech_zone_array_cp < 0] = 0
    rech_zone_array_cp_mask = np.ma.masked_where(rech_zone_array_cp == 0, rech_zone_array_cp)
    unique_vals = np.unique(rech_zone_array_cp)[1:]
    import matplotlib as mpl
    colors = mpl.cm.jet(unique_vals)
    cmap = mpl.colors.ListedColormap(colors)
    ax = plt.imshow(rech_zone_array_cp_mask, cmap=cmap)
    cax = plt.colorbar(ax, ticks=unique_vals)
    cax.set_ticklabels(unique_vals)
    
    zone_coverage = {}
    def update_recharge(rech_dict):
        for key in list(interp_rain.keys()):
            for key2 in rech_dict:
                zone_coverage[key2] = len(recharge_zone_array[(recharge_zone_array == float(key2)) & (m.model_mesh3D[1][0] != -1)])
                interp_rain[key][recharge_zone_array == float(key2)] = \
                    rech_dict[key2] / (1000. * 365. * 86400.)
            interp_rain[key][recharge_zone_array == rch_zone_dict[0]] = \
                interp_rain[key][recharge_zone_array == rch_zone_dict[0]] * 0.0
            interp_rain[key][m.model_mesh3D[1][0] == -1] = 0.

        return interp_rain

    interp_rain = update_recharge(rech_dict)
    rch = interp_rain
    m.boundaries.update_boundary_array('Rain_reduced', rch)

    if verbose:
        print("************************************************************************")
        print(" Updating Murray River GHB boundary")
    
#    MurrayGHB = []
#
#    MurrayGHB_cells = [[x[0], x[1], x[2], x[3]] for x in m.boundaries.bc['GHB']['bc_array'][0]]
#    for MurrayGHB_cell in MurrayGHB_cells:
#        lay, row, col = MurrayGHB_cell[:3]
#        MurrayGHBstage = MurrayGHB_cell[3] #m.parameters.param['mghb_stage']['PARVAL1']
#        dx = m.gridHeight
#        dz = m.model_mesh3D[0][lay][row][col] - \
#            m.model_mesh3D[0][lay + 1][row][col]
#        MGHBconductance = dx * dz * m.parameters.param['mghbk']['PARVAL1']
#        MurrayGHB += [[lay, row, col, MurrayGHBstage, MGHBconductance]]
#        
#
#    ghb = {}
#    ghb[0] = MurrayGHB

    if verbose:
        print("************************************************************************")
        print(" Updating GHB boundary")

#    m.boundaries.update_boundary_array('GHB', ghb)

    if verbose:
        print("************************************************************************")
        print(" Updating Drains boundary")
    
#    mapped_drains = m.polyline_mapped['Drain_Clip_model.shp']
#    
#    simple_drain = []
#    drain_width_avg = 3.0 #m
#    drain_bed_thickness = 0.10 #m
#    for drain_cell in mapped_drains:
#        row = drain_cell[0][0]
#        col = drain_cell[0][1]
#        if m.model_mesh3D[1][0][row][col] == -1:
#            continue
#        #print m.model_mesh3D
#        drain_bed = m.model_mesh3D[0][0][row][col] - m.parameters.param['drain_drop']['PARVAL1']
#        drain_cond = drain_cell[1] * drain_width_avg * m.parameters.param['kv_drain']['PARVAL1'] / drain_bed_thickness
#        simple_drain += [[0, row, col, drain_bed, drain_cond]]
#    
#    drain = {}
#    drain[0] = simple_drain
#    
#    m.boundaries.assign_boundary_array('Drain', drain)
    
    if verbose:
        print("************************************************************************")
        print(" Updating Channels boundary")

    if verbose:
        print("************************************************************************")
        print(" Updating pumping boundary")

    if verbose:
        print("************************************************************************")
        print(" Check for boundary condition updating")
        m.generate_update_report()
    
    if verbose:
        print("************************************************************************")
        print(" Set initial head ")

    path=os.path.join(data_folder)
    fname="01_steady_state"
    head = flopyInterface.get_previous_conditions(os.path.join(path, fname+'.hds'))
    
    m.initial_conditions.set_as_initial_condition("Head", head)
    
    if verbose:
        print("************************************************************************")
        print(" Build and run MODFLOW model ")
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ## Currently using flopyInterface directly rather than running from the ModelManager ...

    modflow_model = flopyInterface.ModflowModel(m, data_folder=data_folder)

    modflow_model.executable = mf_exe_folder

    modflow_model.headtol = 1E-5
    modflow_model.fluxtol = 1E-1
    modflow_model.steady = False

    modflow_model.buildMODFLOW(transport=False, write=True, verbose=True, check=True)

    converge = modflow_model.runMODFLOW(silent=True)

    if converge:
        if verbose:
            print('model converged')
            #break
    
        modflow_model.write_observations()

    wat_bal_df = modflow_model.waterBalance(1, plot=False)
    rech_all = wat_bal_df.loc['RECHARGE_pos']['Flux m^3/d']
    rech_cells = len(rch[0][(rch[0] > 0) & (m.model_mesh3D[1][0] != -1)]) 
    rech_area = rech_cells * 1000 * 1000
    for key in list(zone_coverage.keys()):
        print((key, zone_coverage[key] / float(rech_cells) * 100))
        
    rech = rech_all / rech_area * 1000 * 365.
    #print("Recharge = {} mm/yr".format(rech))
    with open(os.path.join(data_folder, "model_03_recharge_flow", "recharge.txt"), 'w') as f:
        f.write("{}".format(rech))
    
    modflow_model.compareAllObs_metrics(to_file=True)

    modflow_model.compareAllObs()

    return modflow_model

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
            param_file = ""
    else:
        # Get general model config information
        CONFIG = ConfigLoader('../../config/model_config.json')\
                        .set_environment("03_recharge_flow")
        model_config = CONFIG.model_config
        model_folder = model_config['model_folder'] + model_config['grid_resolution'] + os.path.sep
        data_folder = model_config['data_folder']
        mf_exe_folder = model_config['mf_exe_folder']
        param_file = model_config['param_file']
        

    if param_file:
        run = run(model_folder, data_folder, mf_exe_folder, param_file=param_file, verbose=verbose)
    else:
        run = run(model_folder, data_folder, mf_exe_folder, verbose=verbose)
