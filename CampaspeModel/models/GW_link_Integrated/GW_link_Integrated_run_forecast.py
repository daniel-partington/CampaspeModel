"""
Run Campaspe GW model using inputs from Farm model and SW model and
then return SW/GW exchanges, avg depth to GW, depth to GW at ecology sites and
head at trigger bores.
"""

import copy
try:
    import cPickle as pickle
except ImportError:
    import pickle


import os
from os.path import join as p_j
import sys
import warnings

import flopy.utils.binaryfile as bf
import numpy as np
import pandas as pd

from .common_run_funcs import *
from HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface
# Configuration Loader
from HydroModelBuilder.Utilities.Config.ConfigLoader import ConfigLoader

if sys.version_info[0] < 3:
    range = xrange


# TODO: Set the stream gauges, ecology bores, policy bores at the start in some
# other class or in here but so they are available in the run function.
def run(model_folder, data_folder, mf_exe_folder, farm_zones=None, param_file=None, riv_stages=None,
        rainfall_irrigation=None, pumping=None, verbose=True, MM=None, is_steady=False, cache={}):
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

    :returns: tuple[object], five elements
              xchange: np.recarray, exchange for each gauge by gauge ID
              avg_depth_to_gw: np.recarray, average depth for each zone (by Zone ID)
              ecol_depth_to_gw: np.recarray, average gw depth at each cell that contains an ecological bore of interest.
              trigger_heads: np.recarray, Policy trigger well heads
              modflow_model: gw_model object
    """
    if MM is None:
        MM = GWModelManager()
        MM.load_GW_model(p_j(model_folder, "GW_link_Integrated_packaged.pkl"))
    # End if

    # Complimentary models requirements, i.e. bore and gauge data that should be
    # referenceable to this model for parsing specific outputs and receiving inputs:
    if not hasattr(MM, 'ext_linkage_bores') or MM.ext_linkage_bores is None:

        # Read in bores that relate to external models
        this_file_loc = os.path.dirname(os.path.realpath(__file__))
        model_linking = "model_linking.csv"
        model_linking = p_j(this_file_loc, model_linking)

        with open(model_linking, 'r') as f:
            lines = f.readlines()

            for line in lines:
                if line.split(':')[0] == 'Ecology':
                    ecology_bores = process_line(line)
                elif line.split(':')[0] == 'Policy':
                    policy_bores = process_line(line)
                elif line.split(':')[0] == 'SW_stream_gauges':
                    stream_gauges = process_line(line)
            # End for
        # End with

        MM.ext_linkage_bores = {
            "Ecology_bores": ecology_bores,
            "Policy_bores": policy_bores,
            "Stream_gauges": stream_gauges
        }

        cache['Ecology_bores'] = ecology_bores
        cache['Policy_bores'] = policy_bores
        cache['Stream_gauges'] = stream_gauges
    else:
        ecology_bores = cache['Ecology_bores']
        policy_bores = cache['Policy_bores']
        stream_gauges = cache['Stream_gauges']
    # End if

    name = MM.name
    this_model = MM.GW_build[name]

    mesh = this_model.model_mesh3D
    mesh_0, mesh_1 = mesh[0], mesh[1]

    # Load in the new parameters based on parameters.txt or dictionary of new parameters
    if param_file:
        if 'parameters.txt' not in param_file:
            param_file = p_j(param_file, 'parameters.txt')
        # End if
        this_model.updateModelParameters(param_file, verbose=verbose)
    # End if

    model_params = this_model.parameters.param
    model_boundaries = this_model.boundaries
    model_boundaries_bc = model_boundaries.bc

    campaspe_rivermap = this_model.river_mapping['Campaspe'].copy()
    campaspe_rivermap['stage_from_gauge'] = campaspe_rivermap['stage']
    campaspe_rivermap['multi'] = np.nan
    river_seg = campaspe_rivermap.to_records()

    campaspe_pp = this_model.pilot_points['Campaspe']
    num_reaches = campaspe_pp.num_points
    known_points = campaspe_pp.points

    strcond_val = [model_params['kv_riv{}'.format(
        x)]['PARVAL1'] for x in range(num_reaches)]
    river_seg['strhc1'] = np.interp(
        river_seg['Cumulative Length'].tolist(), known_points, strcond_val)
    campaspe_stage = river_seg[river_seg['gauge_id'] != 'none']

    # common_gauges = campaspe_stage[np.isin(campaspe_stage['gauge_id'], common_gauges)].tolist().intersection(riv_stages.dtype.names)
    gauges = set(np.unique(campaspe_stage['gauge_id']))
    common_gauges = gauges.intersection(riv_stages.dtype.names)
    common_gauges_idx = campaspe_stage[np.isin(
        campaspe_stage['gauge_id'], list(common_gauges))]['index']

    iseg_idx = river_seg[np.isin(
        river_seg['iseg'], campaspe_stage['iseg'])]['index']
    ######

    tmp = riv_stages[list(common_gauges)][0].tolist()
    for idx, row in enumerate(common_gauges_idx):
        campaspe_stage[campaspe_stage['index'] ==
                       row]['stage_from_gauge'] = tmp[idx]
    # End for

    # Need data in reverse order
    tmp = np.sort(campaspe_stage['stage_from_gauge'])[::-1]
    for idx, row in enumerate(iseg_idx):
        river_seg['stage_from_gauge'][row] = tmp[idx]
    # End for

    # Convert stages to depths
    campaspe_stage['depth'] = campaspe_stage['stage'] - \
        campaspe_stage['gauge_zero']

    river_seg['depth'] = np.nan

    matching_idx = np.isin(river_seg['index'], campaspe_stage['index'])
    river_seg['depth'][matching_idx] = campaspe_stage['depth']

    # interpolate depths

    # Numpy approach - much faster but different results
    # nan_idx = np.isnan(river_seg['depth'])
    # not_nans = np.logical_not(nan_idx)
    # good_data = river_seg['depth'][not_nans]

    # nan_idx = np.where(nan_idx)[0]
    # interpolated = np.interp(
    #     nan_idx, not_nans.nonzero()[0], good_data)
    # river_seg['depth'][nan_idx] = interpolated

    # Using Pandas to interpolate
    campaspe_rivermap.loc[:, 'depth'] = campaspe_rivermap.set_index(
                                            campaspe_rivermap['Cumulative Length']
                                        )['depth'].interpolate(method='values', 
                                                               limit_direction='both').values
    river_seg['depth'] = campaspe_rivermap['depth'].bfill()

    # Original interpolation method
    # river_seg.loc[:, 'depth'] = river_seg.set_index(river_seg['Cumulative Length'])[
    #     'depth'].interpolate(method='values', limit_direction='both').tolist()
    # river_seg.loc[:, 'depth'] = river_seg['depth'].bfill()

    # Calculate stages using bed elevation and depth
    river_seg['stage'] = river_seg['strtop'] + river_seg['depth']
    river_seg['multi'] = river_seg['strhc1'] * river_seg['rchlen'] * river_seg['width1']

    simple_river = river_seg[['k', 'i', 'j',
                              'stage', 'multi', 'strtop']].tolist()
    model_boundaries.update_boundary_array('Campaspe River', {0: simple_river})

    if verbose:
        print("************************************************************************")
        print(" Updating HGU parameters ")

    # This needs to be automatically generated with the map_raster2mesh routine ...
    # zone_map = {1: 'qa', 2: 'utb', 3: 'utqa', 4: 'utam', 5: 'utaf', 6: 'lta', 7: 'bse'}

    if verbose:
        print("************************************************************************")
        print(" Updating Murray River boundary")

    # Updating Murray River
    mriver_seg = this_model.river_mapping['Murray']
    mriver_seg['strhc1'] = model_params['kv_rm']['PARVAL1']
    mriver_seg['multi'] = mriver_seg['strhc1'] * \
        mriver_seg['rchlen'] * mriver_seg['width1']
    msimple_river = mriver_seg[['k', 'i', 'j',
                                'stage', 'multi', 'strtop']].values.tolist()

    model_boundaries.update_boundary_array('Murray River', {0: msimple_river})

    if verbose:
        print("************************************************************************")
        print(" Updating recharge boundary ")

    # Adjust rainfall to recharge using zoned rainfall reduction parameters
    # Need to make copies of all rainfall arrays
    interp_rain = model_boundaries_bc['Rainfall']['bc_array']
    if rainfall_irrigation is not None:
        for key in list(interp_rain.keys()):
            interp_rain[key] = np.copy(rainfall_irrigation)
            # if key > 0:
            #     interp_rain[key] = np.copy(interp_rain[key])
            # else:
            #     # first time period (integrated model is run on a single time period)
            #     interp_rain[key] = np.copy(rainfall_irrigation)
            # # End if
        # End for
    else:
        warnings.warn("Rainfall+Irrigation input ignored by GW model")
        for key in list(interp_rain.keys()):
            interp_rain[key] = np.copy(interp_rain[key])
    # End if

    recharge_zone_array = model_boundaries_bc['Rain_reduced']['zonal_array']
    rch_zone_dict = model_boundaries_bc['Rain_reduced']['zonal_dict']
    rch_zones = len(list(rch_zone_dict.keys()))

    par_rech_vals = [model_params['rchred{}'.format(i)]['PARLBND']  # ['PARVAL1']
                     for i in range(rch_zones - 1)]

    rch = update_recharge(par_rech_vals, interp_rain,
                          rch_zones, recharge_zone_array, rch_zone_dict, mesh_1)

    model_boundaries.update_boundary_array('Rain_reduced', rch)

    if verbose:
        print("************************************************************************")
        print(" Updating pumping boundary ")

    try:
        # summed_vals = cache['summed_pump']
        normalized_well = cache['normalized_pumpy']
    except KeyError:
        pumpy = model_boundaries_bc['licenced_wells']['bc_array']

        # Need a copy of the list of lists
        well_boundary = copy.deepcopy(pumpy)

        PUMPING_VAL_POS = 3  # zero-based index
        # sum the 4th item in the list of lists
        summed_vals = {key: sum(zip(*vals)[PUMPING_VAL_POS]) 
                       for key, vals in well_boundary.items()}
        summed_vals = float(sum(summed_vals.values()))

        normalized_well = {key: [[b[0], b[1], b[2], b[3] / summed_vals] for b in a]
                           for key, a in well_boundary.items()}

        # cache['summed_pump'] = summed_vals
        cache['normalized_pumpy'] = normalized_well
    # End try

    ################
    # pumpy = model_boundaries_bc['licenced_wells']['bc_array']

    # PUMPING_VAL_POS = 3  # zero-based index
    # # sum the 4th item in the list of lists
    # summed_vals = {key: sum(zip(*vals)[PUMPING_VAL_POS]) 
    #                for key, vals in pumpy.items()}
    # summed_vals = float(sum(summed_vals.values()))

    # normalized_well = {key: [[b[0], b[1], b[2], b[3] / summed_vals] for b in a] 
    #                    for key, a in pumpy.items()}
    # cache['normalized_pumpy'] = normalized_well

    # well dict element structure: [ active_layer, row, col, -time[1][pump] ]
    # DEBUG: ORIGINAL
    # wel = {key: [[b[0], b[1], b[2], b[3] * pumping] for b in a]
    #        for key, a in pumpy.iteritems()}
    # #############

    """
    What should be done is to get the bc array for pumping first,
    cycle through and get the sum of pumping and then normalise it all so it sums to 1,
    then you can go through applying your pumping value such that the model pumping is as you want it to be specified.
    """
    wel = {key: [[b[0], b[1], b[2], b[3] * pumping] for b in a]
           for key, a in normalized_well.items()}

    model_boundaries.assign_boundary_array('licenced_wells', wel)

    if verbose:
        print("************************************************************************")
        print(" Updating GHB boundary ")

    kh = this_model.properties.properties['Kh']

    MurrayGHB = []
    MurrayGHB_cells = [[x[0], x[1], x[2], x[3]]
                       for x in model_boundaries_bc['GHB']['bc_array'][0]]
    dx = this_model.gridHeight
    for MurrayGHB_cell in MurrayGHB_cells:
        lay, row, col = MurrayGHB_cell[:3]
        MurrayGHBstage = MurrayGHB_cell[3]
        dz = mesh_0[lay][row][col] - mesh_0[lay + 1][row][col]
        # model_params['mghbk']['PARVAL1']
        MGHBconductance = dx * dz * kh[lay][row][col]
        MurrayGHB += [[lay, row, col, MurrayGHBstage, MGHBconductance]]
    # End for

    model_boundaries.assign_boundary_array('GHB', {0: MurrayGHB})

    if verbose:
        print("************************************************************************")
        print(" Updating heads boundary ")

    fname = 'model_{}'.format(name)
    try:
        head = flopyInterface.get_previous_conditions(
            p_j(data_folder, fname, name + '.hds')
        )

        this_model.initial_conditions.set_as_initial_condition("Head", head)
    except IndexError:
        raise IndexError(
            "Corrupted MODFLOW hds file - check, replace, or clear {}".format(
                p_j(data_folder, fname, name + '.hds')))
    except IOError:
        warnings.warn(
            "MODFLOW hds file not found. Recreating head state from existing values.")
        head = np.stack([this_model.model_mesh3D[0][0]
                         for i in range(7)], axis=0)
        this_model.initial_conditions.set_as_initial_condition("Head", head)
    # End try

    # Currently using flopyInterface directly rather than running from the ModelManager ...
    modflow_model = flopyInterface.ModflowModel(this_model,
                                                data_folder=p_j(
                                                    data_folder, "model_{}".format(name)),
                                                executable=mf_exe_folder,
                                                nstp=1,
                                                nper=1,
                                                perlen=1,
                                                steady=is_steady)

    # Override temporal aspects of model build:
    # modflow_model.steady = is_steady  # This is to tell FloPy that is a transient model
    # modflow_model.executable = mf_exe_folder

    # modflow_model.nstp = 1
    # modflow_model.nper = 1
    # modflow_model.perlen = 1

    ######### REMOVE THIS WHEN READY ###############
    ### Collecting recharge volume data

    # np.save('recharge_zone_array.npy', recharge_zone_array)

    # with open('rch_zone_dict.pkl', 'wb') as fp:
    #     pickle.dump(rch_zone_dict, fp)

    # with open('par_rech_vals.pkl', 'wb') as fp:
    #     pickle.dump(par_rech_vals, fp)

    modflow_model.buildMODFLOW()
    modflow_model.runMODFLOW()

    rain_red = np.ones_like(recharge_zone_array, dtype=np.float64)
    for key in rch_zone_dict:
        if key == 0:
            rain_red[recharge_zone_array == rch_zone_dict[key]] = 0.0
        else:
            rain_red[recharge_zone_array ==
                        rch_zone_dict[key]] = par_rech_vals[key - 1]
        # End if
    # End for

    bas = modflow_model.bas  # get_package('BAS6')
    rain_red[bas.ibound.array[0] == 0] = 0.
    rain_red_w_area = rain_red * (5000**2)

    # recharge vol in m^3
    recharge_volume = np.sum(rain_red_w_area * rainfall_irrigation)

    ##################

    # SW-GW exchanges:
    swgw_exchanges = np.recarray(
        (1,), dtype=[(gauge, np.float) for gauge in stream_gauges])
    avg_depth_to_gw = np.recarray(
        (1,), dtype=[(farm_zone, np.float) for farm_zone in farm_zones])
    ecol_depth_to_gw = np.recarray(
        (1,), dtype=[(bore, np.float) for bore in ecology_bores])
    trigger_heads = np.recarray(
        (1, ), dtype=[(bore, np.float) for bore in policy_bores])
    rech_vol = np.recarray(
        (1, ), dtype=[(bore, np.float) for bore in ['river']])
    
    rech_vol['river'] = recharge_volume

    river_reach_cells = river_seg[['gauge_id', 'k', 'j', 'i', 
                                   'amalg_riv_points']]
    # Clear out any gauges not required
    # river_reach_cells.set_value(0, 'gauge_id', 'none')
    river_reach_cells['gauge_id'][0] = 'none'
    river_reach_gauges = river_reach_cells['gauge_id']

    river_reach_cells['gauge_id'][river_reach_gauges == 'none'] = np.nan
    river_reach_cells['gauge_id'][~np.isin(river_reach_gauges, stream_gauges)] = np.nan

    # backward fill nans
    river_reach_cells['gauge_id'] = pd.Series(river_reach_cells['gauge_id']).bfill().values


    # mask = np.isnan(river_reach_cells['gauge_id'])
    # idx = np.where(~mask, np.arange(mask.shape[1]), mask.shape[1] - 1)
    # idx = np.minimum.accumulate(idx[:, ::-1], axis=1)[:, ::-1]
    # river_reach_cells = river_reach_cells[np.arange(idx.shape[0])[:, None], idx]
    # # river_reach_cells = river_reach_cells.bfill()
    river_reach_gauges = river_reach_cells['gauge_id']

    river_cells = river_reach_cells[[
        'k', 'i', 'j']][np.isin(river_reach_gauges, stream_gauges)]
    # river_reach_cells['cell'] = river_reach_cells.loc[river_reach_gauges.isin(stream_gauges),
    #                                                   ['k', 'i', 'j']].values.tolist()

    reach_cells = {}
    for gauge in stream_gauges:
        reach_cells[gauge] = river_cells[np.where(river_reach_gauges == gauge)[0]]
    # End for

    river_reach_cells = reach_cells

    river_reach_ecol = river_seg[['gauge_id', 'k', 'j', 'i']]
    river_reach_ecol = river_reach_ecol[river_reach_ecol['gauge_id'] != 'none']
    ecol_gauge_id = river_reach_ecol['gauge_id']
    river_reach_ecol_cells = river_reach_ecol[['k', 'i', 'j']]

    eco_cells = {}
    for ind, ecol_bore in enumerate(ecology_bores):
        eco_cells[stream_gauges[ind]] = river_reach_ecol_cells[np.where(ecol_gauge_id ==
                                                                        stream_gauges[ind])[0]][0]
    # End for

    river_reach_ecol = eco_cells

    # Mask all cells that are either Coonambidgal or Shepparton formation
    # A mask could also be constructed for farm areas by using the mapped farms
    # from the groundwater model builder object
    tgt_mesh = modflow_model.model_data.model_mesh3D[1][0]
    geo_mask = (tgt_mesh == 3) | (tgt_mesh == 1)
    farm_map, farm_map_dict = this_model.polygons_mapped['farm_v1_prj_model.shp']

    ###############

    for gauge in stream_gauges:
        swgw_exchanges[gauge] = modflow_model.getRivFluxNodes(
            river_reach_cells[gauge])
    # End for

    # Average depth to GW table:
    # Mask all cells that are either Coonambidgal or Shepparton formation
    # A mask could also be constructed for farm areas by using the mapped farms
    # from the groundwater model builder object
    for farm_zone in farm_zones:
        mask = geo_mask & (farm_map == farm_map_dict[int(farm_zone)])
        avg_depth_to_gw[farm_zone] = modflow_model.get_average_depth_to_GW(
            mask=mask)
    # End for

    """
    Ecology heads of importance
    River gauges of importance for ecology: 406201, 406202, 406207, 406218, 406265
    Corresponding GW bores nearest to:      83003,  89586,  82999,  5662,   44828

    The bores are being ignored for now as they are too far away from the stream
    gauges. Instead the simulated head in the cell with the stream gauge is
    used which represents the average head over the 5 km x 5 km cell.
    """
    heads = modflow_model.get_heads()
    for ind, ecol_bore in enumerate(ecology_bores):
        _i, _h, _j = river_reach_ecol[stream_gauges[ind]]
        ecol_depth_to_gw[ecol_bore] = mesh_0[_i, _h, _j] - heads[_i, _h, _j]
    # End for

    # TODO: Check that all of the wells listed were mapped to the model mesh and
    # are available for inspection

    """
    Trigger heads of importance for policy

    1. The reference trigger bore for the Elmore-Rochester/Echuca/Bamawn zone
       was selected to represent the interactions between the river and
       groundwater extractions. So we'd need to make sure we have baseflow
       represented in the river between Lake Eppalock  (which I think is an
       input to the ecology model so already covered).

    2. The reference trigger bore for the Barnadown zone was selected to
       represent the gradient of groundwater flow between the Campaspe and the
       Murray. So can we have the gradient of flow as an indicator in the
       integrated model as well (if it isn't already)?
    """
    for trigger_bore in policy_bores:
        # NOTE: This returns the head in mAHD
        trigger_heads[trigger_bore] = modflow_model.get_observation(
            trigger_bore, 0, 'policy_bores')[0]
    # End for

    modflow_model.cleanup()

    # TODO: Check that all of the wells listed were mapped to the model mesh and
    # are available for inspection
    return swgw_exchanges, avg_depth_to_gw, ecol_depth_to_gw, trigger_heads, rech_vol, modflow_model

# End run()


def main():
    print("Running from: " + os.getcwd())

    CONFIG = load_model_config()

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
        model_folder = model_config['model_folder'] + \
            model_config['grid_resolution'] + os.path.sep
        data_folder = os.path.join(model_config['data_folder'], 'forecast')
        mf_exe_folder = model_config['mf_exe_folder']
        param_file = model_config['param_file']
    # End if

    model_folder = model_folder.replace(
        "structured_model_grid_5000m", "forecast/structured_model_grid_5000m")
    MM = GWModelManager()
    MM.load_GW_model(os.path.join(
        model_folder, r"GW_link_Integrated_packaged.pkl"))

    this_model = MM.GW_build[MM.name]
    update_campaspe_pilot_points(this_model, model_folder, use_alt_vals=True)

    model_name = this_model.name
    model_dir = 'model_{}'.format(model_name)
    heads_file_loc = os.path.join(data_folder.replace(
        'forecast', ''), '{}.hds'.format('forecast_initial'))
    dst_loc = os.path.join(data_folder, model_dir, '{}.hds'.format(model_name))

    from shutil import copyfile
    print('Copying {} to {}'.format(heads_file_loc, dst_loc))
    copyfile(heads_file_loc, dst_loc)

    for i in range(300):
        run_params = {
            "model_folder": model_folder,
            "data_folder": data_folder,
            "mf_exe_folder": mf_exe_folder,
            "farm_zones": ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'],
            "param_file": param_file if param_file else None,
            "riv_stages": riv_stages,
            "rainfall_irrigation": None,
            "pumping": 1.01,  # m^3/day
            "MM": MM,
            "verbose": False,
            "is_steady": False
        }

        swgw_exchanges, avg_depth_to_gw, ecol_depth_to_gw, trigger_heads, modflow_model = run(
            **run_params)
        # print("swgw_exchanges", swgw_exchanges)
        # print("avg_depth_to_gw", avg_depth_to_gw)
        # print("ecol_depth_to_gw", ecol_depth_to_gw)
        # print("trigger_heads", trigger_heads)
    # end for
    return modflow_model


if __name__ == "__main__":
    modflow_model = main()
