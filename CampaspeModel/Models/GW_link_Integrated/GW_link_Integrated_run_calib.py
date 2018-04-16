"""
Run Campaspe GW model using inputs from Farm model and SW model and
then return SW/GW exchanges, avg depth to GW, depth to GW at ecology sites and
head at trigger bores.
"""

import cPickle as pickle
import os
import sys
import warnings
import time

import flopy.utils.binaryfile as bf
import numpy as np

from HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface
# Configuration Loader
from HydroModelBuilder.Utilities.Config.ConfigLoader import ConfigLoader


# TODO: Set the stream gauges, ecology bores, policy bores at the start in some
# other class or in here but so they are available in the run function.


def process_line(line):
    return [x.strip() for x in line.split(':')[1].strip().split(',')]
# end process_line


def run(model_folder, data_folder, mf_exe_folder, farm_zones=None, param_file=None, riv_stages=None,
        rainfall_irrigation=None, pumping=None, verbose=True, MM=None, recharge_amt=0.03, is_steady=False, cache={}):
    """
    GW Model Runner.

    :param model_folder: str, path to model folder
    :param data_folder: str, path to data
    :param mf_exe_folder: str, path to MODFLOW executable
    :param farm_zones: list[str], farm zone IDs to map ground water level values to (in order)
    :param param_file: str, path to parameter file
    :param riv_stages: np.recarray, gauge numbers and stage
    :param rainfall_irrigation: np.ndarray, array representing rainfall and irrigation input.
                                Must match the model extent.
    :param pumping: float, daily pumping amount in m^3/day (ML/day to m^3/day => ML * 1000)

    :returns: tuple[np.recarray], four elements
              xchange: exchange for each gauge by gauge ID
              avg_gw_depth: average depth for each zone
              ecol_depth_to_gw: average gw depth at each cell that contains an ecological bore of interest.
              trigger_head: Policy trigger well heads
    """
    p_j = os.path.join

    # DEBUG: Setting constant high pumping rate to see if gw levels go down
    # pumping = 1.0

    if MM is None:
        MM = GWModelManager()
        MM.load_GW_model(p_j(model_folder, "GW_link_Integrated_packaged.pkl"))
    # End if

    # Complimentary models requirements, i.e. bore and gauge data that should be
    # referenceable to this model for parsing specific outputs and receiving inputs:
    if not hasattr(MM, 'ext_linkage_bores') or MM.ext_linkage_bores is None:

        # Read in bores that relate to external models
        model_linking = p_j(data_folder, "model_linking.csv")

        # Need to check given data folder and its parent directory
        # Dirty hack, I know :(
        if not os.path.exists(model_linking):
            model_linking = p_j(data_folder, "..", "model_linking.csv")
            if not os.path.exists(model_linking):
                raise IOError("Could not find bore linkages information (`model_linking.csv`)")
            # End if
        # End if

        with open(model_linking, 'r') as f:
            lines = f.readlines()

            for line in lines:
                if line.split(':')[0] == 'Ecology':
                    Ecology_bores = process_line(line)
                elif line.split(':')[0] == 'Policy':
                    Policy_bores = process_line(line)
                elif line.split(':')[0] == 'SW_stream_gauges':
                    Stream_gauges = process_line(line)
            # End for
        # End with

        MM.ext_linkage_bores = {
            "Ecology_bores": Ecology_bores,
            "Policy_bores": Policy_bores,
            "Stream_gauges": Stream_gauges
        }

        cache['Ecology_bores'] = Ecology_bores
        cache['Policy_bores'] = Policy_bores
        cache['Stream_gauges'] = Stream_gauges
    else:
        Ecology_bores = cache['Ecology_bores']
        Policy_bores = cache['Policy_bores']
        Stream_gauges = cache['Stream_gauges']
    # End if

    name = MM.GW_build.keys()[0]
    this_model = MM.GW_build[name]
    mesh = this_model.model_mesh3D
    mesh_0, mesh_1 = mesh[0], mesh[1]

    # Load in the new parameters based on parameters.txt or dictionary of new parameters
    if param_file:
        this_model.updateModelParameters(p_j(data_folder, 'parameters.txt'), verbose=verbose)
    # End if

    model_params = this_model.parameters.param
    model_boundaries = this_model.boundaries
    model_boundaries_bc = model_boundaries.bc

    if verbose:
        print "************************************************************************"
        print " Updating HGU parameters "
    
    # This needs to be automatically generated with the map_raster2mesh routine ...
    zone_map = {1: 'qa', 2: 'utb', 3: 'utqa', 4: 'utam', 5: 'utaf', 6: 'lta', 7: 'bse'}

    default_array = this_model.model_mesh3D[1].astype(float)
    zone = np.copy(default_array)
    kh = np.copy(default_array)
    kv = np.copy(default_array)
    sy = np.copy(default_array)
    ss = np.copy(default_array)
    
    alt_k_vals = {'khutqa': 10.,
                  'khutb' : 10., #1.,
                  'khqa'  : 44.355,
                  'khutam': 10., #0.1,
                  'khutaf': 60.,
                  'khlta' : 60.,
                  'khbse' : 10.}
    
    k_factor = 1.05
    ghb_k_factor = 3.
    for key in alt_k_vals:
        alt_k_vals[key] = alt_k_vals[key] * k_factor
                  
                  
    alt_ss_vals = {'ssutqa': 1E-5 ,
                  'ssutb' : 1E-5,
                  'ssqa'  : 1E-5,
                  'ssutam': 1E-5,
                  'ssutaf': 1E-5,
                  'sslta' : 1E-5,
                  'ssbse' : 1E-5}

    alt_sy_vals = {'syutqa': 0.25,
                  'syutb' : 0.25,
                  'syqa'  : 0.25,
                  'syutam': 0.25,
                  'syutaf': 0.25,
                  'sylta' : 0.25,
                  'sybse' : 0.25}
                  
                  
    def create_pp_points_dict(zone_map, zone, prop_array, prop_name, m):
        points_values_dict = {}
        for index, key in enumerate(zone_map.keys()):
            for index2, param in enumerate(m.parameters.param_set[prop_name + zone_map[key]]):
                if index2 == 0:
                    if prop_name == 'ss':
                        points_values_dict[index] = [alt_ss_vals[prop_name + zone_map[key]]]
                    elif prop_name == 'sy':
                        points_values_dict[index] = [alt_sy_vals[prop_name + zone_map[key]]]
                    elif prop_name == 'kh':
                        points_values_dict[index] = [alt_k_vals[prop_name + zone_map[key]]]
                    else:
                        points_values_dict[index] = [m.parameters.param[param]['PARVAL1']]
                else: 
                    if prop_name == 'ss':
                        points_values_dict[index] += [alt_ss_vals[prop_name + zone_map[key]]]
                    elif prop_name == 'sy':
                        points_values_dict[index] += [alt_sy_vals[prop_name + zone_map[key]]]
                    elif prop_name == 'kh':
                        points_values_dict[index] += [alt_k_vals[prop_name + zone_map[key]]]
                    else:
                        points_values_dict[index] += [m.parameters.param[param]['PARVAL1']]
        return points_values_dict    
        
    def update_pilot_points(zone_map, zone, prop_array, par_name, prop_name, prop_folder, m, prop_array_fname):
        points_values_dict = create_pp_points_dict(zone_map, zone, prop_array, prop_name, m)
        p = m.pilot_points[par_name]
        zones = len(zone_map.keys())
        p.output_directory = os.path.join(model_folder, prop_folder)
        p.update_pilot_points_files_by_zones(zones, points_values_dict)
        time.sleep(3)
        p.run_pyfac2real_by_zones(zones) 
        #p.save_mesh3D_array(filename=os.path.join(data_folder, prop_array_fname))
        return p.val_array

    kh = update_pilot_points(zone_map, zone, kh, 'hk', 'kh', 'hk_pilot_points',
                             this_model, 'hk_val_array')  
    #this_model.save_array(os.path.join(data_folder, 'Kh'), kh)
    if verbose:
        print("Erroneous K pilot cells: {}".format(len(kh[kh > 10000.])))
    kh[kh > 10000.] = 25.
    kv = kh * 0.1
    #this_model.save_array(os.path.join(data_folder, 'Kv'), kv)

    sy = update_pilot_points(zone_map, zone, sy, 'sy', 'sy', 'sy_pilot_points',
                             this_model, 'sy_val_array')
    if verbose:
        print("Erroneous Sy pilot cells: {}".format(len(sy[sy > 0.5])))
    sy[sy > 0.5] = 0.5
    #this_model.save_array(os.path.join(data_folder, 'Sy'), sy)
    
    ss = update_pilot_points(zone_map, zone, ss, 'ss', 'ss', 'ss_pilot_points',
                             this_model, 'ss_val_array')
    if verbose:
        print("Erroneous Ss pilot cells: {}".format(len(ss[ss > 0.01])))
    ss[ss > 0.01] = 1E-5
    #this_model.save_array(os.path.join(data_folder, 'SS'), ss)
    
    this_model.properties.update_model_properties('Kh', kh)
    this_model.properties.update_model_properties('Kv', kv)
    this_model.properties.update_model_properties('Sy', sy)
    this_model.properties.update_model_properties('SS', ss)    

    if verbose:
        print "************************************************************************"
        print " Updating Campaspe River boundary"

    river_seg = this_model.river_mapping['Campaspe']
    num_reaches = this_model.pilot_points['Campaspe'].num_points
    known_points = this_model.pilot_points['Campaspe'].points

    strcond_val = [model_params['kv_riv{}'.format(x)]['PARVAL1'] for x in xrange(num_reaches)]
    river_seg['strhc1'] = np.interp(river_seg['Cumulative Length'].values.tolist(), known_points, strcond_val) * 0.05
    river_seg['multi'] = river_seg['strhc1'] * river_seg['rchlen'] * river_seg['width1']
    river_bc = this_model.boundaries.bc['Campaspe River']['bc_array']

    if not is_steady:
        for key in river_bc:
            riv_stress_period = river_bc[key]
            for idx, rc in enumerate(riv_stress_period):
                river_bc[key][idx] = [rc[0], rc[1], rc[2], rc[3], river_seg['multi'].tolist()[idx], rc[5]] 
    else:
        mean_stage = {}
        for key in river_bc:
            riv_stress_period = river_bc[key]
            for idx, rc in enumerate(riv_stress_period):
                if key == 0:
                    mean_stage[idx] = rc[3] 
                else:
                    mean_stage[idx] = mean_stage[idx] * (1. - 1. / (key + 1.)) + rc[3] / (key + 1.)
        
        for key in river_bc:
            riv_stress_period = river_bc[key]
            for idx, rc in enumerate(riv_stress_period):
                river_bc[key][idx] = [rc[0], rc[1], rc[2], mean_stage[idx], river_seg['multi'].tolist()[idx], rc[5]] 

     #simple_river = river_seg[['k', 'i', 'j', 'stage', 'multi', 'strtop']].values.tolist()

    model_boundaries.update_boundary_array('Campaspe River', river_bc)

    if verbose:
        print "************************************************************************"
        print " Updating Murray River boundary"

    # Updating Murray River
    mriver_seg = this_model.river_mapping['Murray']
    mriver_seg['strhc1'] = model_params['kv_rm']['PARVAL1']
    mriver_seg['multi'] = mriver_seg['strhc1'] * mriver_seg['rchlen'] * mriver_seg['width1']

    mriver_bc = this_model.boundaries.bc['Murray River']['bc_array']
    for key in mriver_bc:
        mriv_stress_period = mriver_bc[key]
        for idx, rc in enumerate(mriv_stress_period):
            mriver_bc[key][idx] = [rc[0], rc[1], rc[2], rc[3], mriver_seg['multi'].tolist()[idx], rc[5]] 
    #msimple_river = mriver_seg[['k', 'i', 'j', 'stage', 'multi', 'strtop']].values.tolist()

    model_boundaries.update_boundary_array('Murray River', mriver_bc)

    if verbose:
        print "************************************************************************"
        print " Updating recharge boundary "

    # Adjust rainfall to recharge using zoned rainfall reduction parameters
    # Need to make copies of all rainfall arrays
    interp_rain = model_boundaries_bc['Rainfall']['bc_array']
    for key in interp_rain.keys():
        interp_rain[key] = np.copy(interp_rain[key])
    
    if is_steady:
        interp_rain[0] = np.mean([interp_rain[x] for x in interp_rain], axis=0)
        
    recharge_zone_array = model_boundaries_bc['Rain_reduced']['zonal_array']
    rch_zone_dict = model_boundaries_bc['Rain_reduced']['zonal_dict']

    rch_zones = len(rch_zone_dict.keys())

    par_rech_vals = [model_params['rchred{}'.format(i)]['PARVAL1'] \
                     for i in xrange(rch_zones - 1)]

    def update_recharge(vals):
        for key in interp_rain.keys():
            for i in xrange(rch_zones - 1):
#                interp_rain[key][recharge_zone_array == rch_zone_dict[i + 1]] = \
#                    interp_rain[key][recharge_zone_array == rch_zone_dict[i + 1]] * \
#                    vals[i]
                if i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]:
                    interp_rain[key][recharge_zone_array == rch_zone_dict[i + 1]] = \
                        interp_rain[key][recharge_zone_array == rch_zone_dict[i + 1]] * \
                        vals[i] * 0.
                else:
                    interp_rain[key][recharge_zone_array == rch_zone_dict[i + 1]] = \
                        interp_rain[key][recharge_zone_array == rch_zone_dict[i + 1]] * \
                        vals[i] #* 0.005

            interp_rain[key][recharge_zone_array == rch_zone_dict[0]] = \
                interp_rain[key][recharge_zone_array == rch_zone_dict[0]] * 0.0
            # Killing recharge on inactive cells
            interp_rain[key][this_model.model_mesh3D[1][0] == -1] = 0.
            # Killing recharge across the outcropping bedrock              
            interp_rain[key][this_model.model_mesh3D[1][0] == 7] = 0.
        return interp_rain

    interp_rain = update_recharge(par_rech_vals)
    rch = interp_rain

    
    model_boundaries.update_boundary_array('Rain_reduced', rch)
    #model_boundaries.assign_boundary_array('Rain_reduced', {0: interp_rain})

    pumpy = model_boundaries_bc['licenced_wells']['bc_array']
    wel = {key: [[b[0], b[1], b[2], b[3] * pumping] for b in a] for key, a in pumpy.iteritems()}
    if is_steady:
        wel[0] = wel[wel.keys()[-1]]
           
    model_boundaries.assign_boundary_array('licenced_wells', wel)

    MurrayGHB = []
    MurrayGHB_cells = [[x[0], x[1], x[2], x[3]] for x in model_boundaries_bc['GHB']['bc_array'][0]]
    for MurrayGHB_cell in MurrayGHB_cells:
        lay, row, col = MurrayGHB_cell[:3]
        MurrayGHBstage = MurrayGHB_cell[3]
        dx = this_model.gridHeight
        dz = mesh_0[lay][row][col] - mesh_0[lay + 1][row][col]
        MGHBconductance = dx * dz * model_params['mghbk']['PARVAL1'] * ghb_k_factor # / 10000.
        MurrayGHB += [[lay, row, col, MurrayGHBstage, MGHBconductance]]
    # End for

    ghb = {}
    ghb[0] = MurrayGHB    
    
    model_boundaries.assign_boundary_array('GHB', {0: MurrayGHB})
    
    if not is_steady:
        fname = 'model_{}'.format(name)
        try:
            headobj = bf.HeadFile(p_j(data_folder, fname, name + '.hds'))
            times = headobj.get_times()
            head = headobj.get_data(totim=times[-1])
            this_model.initial_conditions.set_as_initial_condition("Head", head)
        except IndexError:
            raise IndexError(
                "Corrupted MODFLOW hds file - check, replace, or clear {}".format(
                    p_j(data_folder, fname, name + '.hds')))
    #    except IOError:
    #        warnings.warn("MODFLOW hds file not found. Recreating head state from existing values.")
    #        head = np.stack([this_model.model_mesh3D[0][0] for i in range(7)], axis=0)
    #        this_model.initial_conditions.set_as_initial_condition("Head", head)
        # End try

    # Currently using flopyInterface directly rather than running from the ModelManager ...
    modflow_model = flopyInterface.ModflowModel(this_model, data_folder=p_j(data_folder, "model_{}".format(name)))

    # Override temporal aspects of model build:
    modflow_model.steady = is_steady  # This is to tell FloPy that is a transient model
    modflow_model.executable = mf_exe_folder

    modflow_model.buildMODFLOW()
    modflow_model.runMODFLOW()

    if not is_steady:
        #modflow_model.viewHeadsByZone2(1, head_name='policy_bores')
        modflow_model.compareAllObs('policy_bores') 
    
#    try:
#        swgw_exchanges = cache['swgw_exchanges']
#        avg_depth_to_gw = cache['avg_depth_to_gw']
#        ecol_depth_to_gw = cache['ecol_depth_to_gw']
#        trigger_heads = cache['trigger_heads']
#
#        river_reach_cells = cache['river_reach_cells']
#        river_reach_ecol = cache['river_reach_ecol']
#        farm_map = cache['farm_map']
#        farm_map_dict = cache['farm_map_dict']
#        geo_mask = cache['geo_mask']
#    except KeyError:
#        # SW-GW exchanges:
#        swgw_exchanges = np.recarray((1,), dtype=[(gauge, np.float) for gauge in Stream_gauges])
#        avg_depth_to_gw = np.recarray((1,), dtype=[(farm_zone, np.float) for farm_zone in farm_zones])
#        ecol_depth_to_gw = np.recarray((1,), dtype=[(bore, np.float) for bore in Ecology_bores])
#        trigger_heads = np.recarray((1, ), dtype=[(bore, np.float) for bore in Policy_bores])
#
#        cache['swgw_exchanges'] = swgw_exchanges
#        cache['avg_depth_to_gw'] = avg_depth_to_gw
#        cache['ecol_depth_to_gw'] = ecol_depth_to_gw
#        cache['trigger_heads'] = trigger_heads
#
#        river_reach_cells = river_seg[['gauge_id', 'k', 'j', 'i', 'amalg_riv_points']]
#        river_reach_cells.set_value(0, 'gauge_id', 'none')
#
#        river_reach_cells.loc[river_reach_cells['gauge_id'] == 'none', 'gauge_id'] = np.nan
#        river_reach_cells = river_reach_cells.bfill()
#        river_reach_cells['cell'] = river_reach_cells.loc[river_reach_cells['gauge_id'].isin(Stream_gauges),
#                                                          ['k', 'i', 'j']].values.tolist()
#
#        reach_cells = {}
#        for gauge in Stream_gauges:
#            reach_cells[gauge] = river_reach_cells.loc[river_reach_cells['gauge_id'] == gauge, 'cell'].values
#        # End for
#
#        river_reach_cells = reach_cells
#        cache['river_reach_cells'] = reach_cells
#
#        river_reach_ecol = river_seg.loc[:, ['gauge_id', 'k', 'j', 'i']]
#        river_reach_ecol = river_reach_ecol[river_reach_ecol['gauge_id'] != 'none']
#        ecol_gauge_id = river_reach_ecol['gauge_id']
#        river_reach_ecol['cell'] = river_reach_ecol.loc[:, ['k', 'i', 'j']].values.tolist()
#
#        eco_cells = {}
#        for ind, ecol_bore in enumerate(Ecology_bores):
#            eco_cells[Stream_gauges[ind]] = river_reach_ecol.loc[ecol_gauge_id == Stream_gauges[ind], 'cell'].values[0]
#        # End for
#
#        river_reach_ecol = eco_cells
#        cache['river_reach_ecol'] = eco_cells
#
#        # Mask all cells that are either Coonambidgal or Shepparton formation
#        # A mask could also be constructed for farm areas by using the mapped farms
#        # from the groundwater model builder object
#        tgt_mesh = modflow_model.model_data.model_mesh3D[1][0]
#        geo_mask = (tgt_mesh == 3) | (tgt_mesh == 1)
#        farm_map, farm_map_dict = this_model.polygons_mapped['farm_v1_prj_model.shp']
#        cache['farm_map'] = farm_map
#        cache['farm_map_dict'] = farm_map_dict
#        cache['geo_mask'] = geo_mask
#    # End try
#
#    for gauge in Stream_gauges:
#        swgw_exchanges[gauge] = modflow_model.getRivFluxNodes(river_reach_cells[gauge])
#    # End for
#
#    # Average depth to GW table:
#    # Mask all cells that are either Coonambidgal or Shepparton formation
#    # A mask could also be constructed for farm areas by using the mapped farms
#    # from the groundwater model builder object
#    for farm_zone in farm_zones:
#        mask = geo_mask & (farm_map == farm_map_dict[int(farm_zone)])
#        avg_depth_to_gw[farm_zone] = modflow_model.getAverageDepthToGW(mask=mask)
#    # End for
#
#    """
#    Ecology heads of importance
#    River gauges of importance for ecology: 406201, 406202, 406207, 406218, 406265
#    Corresponding GW bores nearest to:      83003,  89586,  82999,  5662,   44828
#
#    The bores are being ignored for nw as they are too far away from the stream
#    gauges. Instead the simulated head in the cell with the stream gauge is
#    used which represents the average head over the 5 km x 5 km cell.
#    """
#
#    heads = modflow_model.getHeads()
#    for ind, ecol_bore in enumerate(Ecology_bores):
#        _i, _h, _j = river_reach_ecol[Stream_gauges[ind]]
#        ecol_depth_to_gw[ecol_bore] = mesh_0[_i, _h, _j] - heads[_i, _h, _j]
#    # End for
#
#    # TODO: Check that all of the wells listed were mapped to the model mesh and
#    # are available for inspection
#
#    """
#    Trigger heads of importance for policy
#
#    1. The reference trigger bore for the Elmore-Rochester/Echuca/Bamawn zone
#       was selected to represent the interactions between the river and
#       groundwater extractions. So we'd need to make sure we have baseflow
#       represented in the river between Lake Eppalock  (which I think is an
#       input to the ecology model so already covered).
#
#    2. The reference trigger bore for the Barnadown zone was selected to
#       represent the gradient of groundwater flow between the Campaspe and the
#       Murray. So can we have the gradient of flow as an indicator in the
#       integrated model as well (if it isn't already)?
#    """
#    for trigger_bore in Policy_bores:
#        # NOTE: This returns the head in mAHD
#        trigger_heads[trigger_bore] = modflow_model.getObservation(trigger_bore, 0, 'head')[0]
#    # End for

    # TODO: Check that all of the wells listed were mapped to the model mesh and
    # are available for inspection
    #return swgw_exchanges, avg_depth_to_gw, ecol_depth_to_gw, trigger_heads
    return modflow_model
# End run()


def main():
    print("Running from: " + os.getcwd())
    CONFIG = ConfigLoader(os.path.join(os.path.dirname(os.path.dirname(__file__)).replace('models',''),
                                       "config", "model_config.json"))\
        .set_environment("GW_link_Integrated")

    def load_obj(filename):
        if filename[-4:] == '.pkl':
            with open(filename, 'r') as f:
                return pickle.load(f)
        else:
            print('File type not recognised as "pkl"')
        # End if

    # Example river level data (to be inputted from SW Model)
    fname = "initial_river_levels.pkl"
    riv_stages = load_obj(os.path.join(CONFIG.settings['data_folder'], fname))
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
        data_folder = os.path.join(model_config['data_folder'], 'hindcast')
        mf_exe_folder = model_config['mf_exe_folder']
        param_file = model_config['param_file']
    # End if

    model_folder = model_folder.replace("structured_model_grid_5000m", "hindcast/structured_model_grid_5000m")
    
    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder, r"GW_link_Integrated_packaged.pkl"))

    # SS first
    run_params = {
        "model_folder": model_folder,
        "data_folder": data_folder,
        "mf_exe_folder": mf_exe_folder,
        "farm_zones": ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'],
        "param_file": param_file if param_file else None,
        "riv_stages": riv_stages,
        "rainfall_irrigation": None,
        "pumping": 0.0,  # m^3/day
        "MM": MM,
        "verbose": False,
        "is_steady": True
    }

    model_results_ss = run(**run_params)

    ## TR second
    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder, r"GW_link_Integrated_packaged.pkl"))
    run_params = {
        "model_folder": model_folder,
        "data_folder": data_folder,
        "mf_exe_folder": mf_exe_folder,
        "farm_zones": ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'],
        "param_file": param_file if param_file else None,
        "riv_stages": riv_stages,
        "rainfall_irrigation": None,
        "pumping": 1.0,  # m^3/day
        "MM": MM,
        "verbose": False,
        "is_steady": False
    }

    model_results = run(**run_params)
    
#    swgw_exchanges, avg_depth_to_gw, ecol_depth_to_gw, trigger_heads = run(**run_params)
#    print("swgw_exchanges", swgw_exchanges)
#    print("avg_depth_to_gw", avg_depth_to_gw)
#    print("ecol_depth_to_gw", ecol_depth_to_gw)
#    print("trigger_heads", trigger_heads)
    return model_results_ss, model_results

if __name__ == "__main__":
    model_results_ss, model_results = main()
