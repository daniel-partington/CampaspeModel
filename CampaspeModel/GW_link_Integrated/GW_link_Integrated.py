"""
Run Campaspe GW model using inputs from Farm model and SW model and
then return SW/GW exchanges, avg depth to GW, depth to GW at ecology sites and
head at trigger bores.
"""

import cPickle as pickle
import os
import sys
import warnings

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
    :param farm_zones: list, string IDs of farm zones to map ground water level values to (in order)
    :param param_file: str, path to parameter file
    :param riv_stages: np recarray, gauge numbers and stage
    :param rainfall_irrigation: np array, array representing rainfall and irrigation input.
                                Must match the model extent.
    :param pumping: float, daily pumping amount in m/day

    :returns: tuple, four elements
                xchange: numpy recarray, exchange for each gauge by gauge ID
                avg_gw_depth: numpy recarray, Average depth for each zone
                ecol_depth_to_gw: numpy recarray, TODO
                trigger_head: numpy recarray, Trigger well heads
    """
    p_j = os.path.join

    if MM is None:
        MM = GWModelManager()
        MM.load_GW_model(p_j(model_folder, "GW_link_Integrated_packaged.pkl"))
    # End if

    # Complimentary models requirements, i.e. bore and gauge data that should be
    # referenceable to this model for parsing specific outputs and receiving inputs:
    if hasattr(MM, 'ext_linkage_bores'):
        Ecology_bores = cache['Ecology_bores']
        Policy_bores = cache['Policy_bores']
        Stream_gauges = cache['Stream_gauges']
    elif not hasattr(MM, 'ext_linkage_bores') or MM.ext_linkage_bores is None:

        # Read in bores that relate to external models
        model_linking = p_j(data_folder, 'model_linking.csv')
        if not os.path.isfile(model_linking):
            # assume model run with generated subdirectory
            model_linking = p_j(data_folder, '..', 'model_linking.csv')
            if not os.path.isfile(model_linking):
                raise IOError("Model linkage file could not be found")
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

    river_seg = this_model.river_mapping['Campaspe']
    num_reaches = this_model.pilot_points['Campaspe'].num_points
    known_points = this_model.pilot_points['Campaspe'].points

    strcond_val = [model_params['kv_riv{}'.format(x)]['PARVAL1'] for x in xrange(num_reaches)]
    river_seg['strhc1'] = np.interp(river_seg['Cumulative Length'].values.tolist(), known_points, strcond_val)
    river_seg['stage_from_gauge'] = river_seg['stage']
    Campaspe_stage = river_seg.loc[river_seg['gauge_id'] != 'none', :]

    try:
        common_gauges = cache['common_gauges']
        common_gauges_loc = cache['common_gauges_loc']
        iseg_loc = cache['iseg_loc']
    except KeyError:
        # Get common gauges, and fill specific values
        common_gauges = list(set(Campaspe_stage['gauge_id'].values).intersection(riv_stages.dtype.names))
        common_gauges_loc = Campaspe_stage.loc[Campaspe_stage.gauge_id.isin(common_gauges), 'gauge_id']
        iseg_loc = river_seg.loc[river_seg.iseg.isin(Campaspe_stage.iseg), 'stage_from_gauge']

        cache['common_gauges'] = common_gauges
        cache['common_gauges_loc'] = common_gauges_loc
        cache['iseg_loc'] = iseg_loc
    # End try

    tmp = riv_stages[common_gauges][0].tolist()
    for idx, row in enumerate(common_gauges_loc.index):
        Campaspe_stage.set_value(row, 'stage_from_gauge', tmp[idx])
    # End for

    tmp = sorted(Campaspe_stage['stage_from_gauge'].values.tolist(), reverse=True)
    for idx, row in enumerate(iseg_loc.index):
        river_seg.set_value(row, 'stage_from_gauge', tmp[idx])
    # End for

    river_seg['stage'] = \
        river_seg.set_index(river_seg['Cumulative Length'])['stage_from_gauge']. \
        interpolate(method='values', limit_direction='both').values.tolist()

    river_seg['multi'] = river_seg['strhc1'] * river_seg['rchlen'] * river_seg['width1']
    simple_river = river_seg[['k', 'i', 'j', 'stage', 'multi', 'strtop']].values.tolist()

    model_boundaries.update_boundary_array('Campaspe River', {0: simple_river})

    # Updating Murray River
    mriver_seg = this_model.river_mapping['Murray']
    mriver_seg['strhc1'] = model_params['kv_rm']['PARVAL1']
    mriver_seg['multi'] = mriver_seg['strhc1'] * mriver_seg['rchlen'] * mriver_seg['width1']
    msimple_river = mriver_seg[['k', 'i', 'j', 'stage', 'multi', 'strtop']].values.tolist()

    model_boundaries.update_boundary_array('Murray River', {0: msimple_river})

    if rainfall_irrigation is not None:
        interp_rain = np.copy(rainfall_irrigation)
    else:
        warnings.warn("Rainfall+Irrigation input ignored by GW model")
        interp_rain = np.copy(model_boundaries_bc['Rainfall']['bc_array'])
    # End if

    # integers (1 - 7) reflect number of layers
    temp_lst = [4, 5, 6]
    for i in xrange(1, 8):
        match = mesh_1[0] == i

        if i in temp_lst:  # (4 - 6)
            interp_rain[match] = 0
            continue
        # End if

        # Adjust rainfall to recharge using a magic number (defaults to 0.03 -> 3%)
        interp_rain[match] = interp_rain[match] * recharge_amt
    # End for

    model_boundaries.assign_boundary_array('Rain_reduced', {0: interp_rain})

    pumpy = model_boundaries_bc['licenced_wells']['bc_array']
    wel = {key: [[b[0], b[1], b[2], b[3] * pumping] for b in a] for key, a in pumpy.iteritems()}

    model_boundaries.assign_boundary_array('licenced_wells', wel)

    MurrayGHB = []
    MurrayGHB_cells = [[x[0], x[1], x[2], x[3]] for x in model_boundaries_bc['GHB']['bc_array'][0]]
    for MurrayGHB_cell in MurrayGHB_cells:
        lay, row, col = MurrayGHB_cell[:3]
        MurrayGHBstage = MurrayGHB_cell[3]
        dx = this_model.gridHeight
        dz = mesh_0[lay][row][col] - mesh_0[lay + 1][row][col]
        MGHBconductance = dx * dz * model_params['mghbk']['PARVAL1']  # / 10000.
        MurrayGHB += [[lay, row, col, MurrayGHBstage, MGHBconductance]]
    # End for

    model_boundaries.assign_boundary_array('GHB', {0: MurrayGHB})
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
    except IOError:
        if verbose:
            print "Using initial head conditions"
        # End if
    # End try

    # Currently using flopyInterface directly rather than running from the ModelManager ...
    modflow_model = flopyInterface.ModflowModel(this_model, data_folder=data_folder)

    # Override temporal aspects of model build:
    modflow_model.nper = 1  # This is the number of stress periods which is set to 1 here
    modflow_model.perlen = 1  # This is the period of time which is set to 1 day here
    modflow_model.nstp = 1  # This is the number of sub-steps to do in each stress period
    modflow_model.steady = is_steady  # This is to tell FloPy that is a transient model
    modflow_model.executable = mf_exe_folder

    modflow_model.buildMODFLOW()
    modflow_model.runMODFLOW()

    try:
        swgw_exchanges = cache['swgw_exchanges']
        avg_depth_to_gw = cache['avg_depth_to_gw']
        ecol_depth_to_gw = cache['ecol_depth_to_gw']
        trigger_heads = cache['trigger_heads']

        river_reach_cells = cache['river_reach_cells']
        river_reach_ecol = cache['river_reach_ecol']
        farm_map = cache['farm_map']
        farm_map_dict = cache['farm_map_dict']
        geo_mask = cache['geo_mask']
    except KeyError:
        # SW-GW exchanges:
        swgw_exchanges = np.recarray((1,), dtype=[(gauge, np.float) for gauge in Stream_gauges])
        avg_depth_to_gw = np.recarray((1,), dtype=[(farm_zone, np.float) for farm_zone in farm_zones])
        ecol_depth_to_gw = np.recarray((1,), dtype=[(bore, np.float) for bore in Ecology_bores])
        trigger_heads = np.recarray((1, ), dtype=[(bore, np.float) for bore in Policy_bores])

        cache['swgw_exchanges'] = swgw_exchanges
        cache['avg_depth_to_gw'] = avg_depth_to_gw
        cache['ecol_depth_to_gw'] = ecol_depth_to_gw
        cache['trigger_heads'] = trigger_heads

        river_reach_cells = river_seg[['gauge_id', 'k', 'j', 'i', 'amalg_riv_points']]
        river_reach_cells.set_value(0, 'gauge_id', 'none')
        river_reach_cells.loc[river_reach_cells['gauge_id'] == 'none', 'gauge_id'] = np.nan
        river_reach_cells = river_reach_cells.bfill()
        river_reach_cells['cell'] = river_reach_cells.loc[river_reach_cells['gauge_id'].isin(Stream_gauges), [
            'k', 'i', 'j']].values.tolist()

        reach_cells = {}
        for gauge in Stream_gauges:
            reach_cells[gauge] = river_reach_cells.loc[river_reach_cells['gauge_id'] == gauge, 'cell'].values
        # End for

        river_reach_cells = reach_cells
        cache['river_reach_cells'] = reach_cells

        river_reach_ecol = river_seg.loc[:, ['gauge_id', 'k', 'j', 'i']]
        river_reach_ecol = river_reach_ecol[river_reach_ecol['gauge_id'] != 'none']
        river_reach_ecol['cell'] = river_reach_ecol.loc[:, ['k', 'i', 'j']].values.tolist()

        eco_cells = {}
        for ind, ecol_bore in enumerate(Ecology_bores):
            eco_cells[Stream_gauges[ind]] = river_reach_ecol.loc[river_reach_ecol['gauge_id']
                                                                 == Stream_gauges[ind], 'cell'].values[0]
        # End for

        river_reach_ecol = eco_cells
        cache['river_reach_ecol'] = eco_cells

        # Mask all cells that are either Coonambidgal or Shepparton formation
        # A mask could also be constructed for farm areas by using the mapped farms
        # from the groundwater model builder object
        tgt_mesh = modflow_model.model_data.model_mesh3D[1][0]
        geo_mask = (tgt_mesh == 3) | (tgt_mesh == 1)
        farm_map, farm_map_dict = this_model.polygons_mapped['farm_v1_prj_model.shp']
        cache['farm_map'] = farm_map
        cache['farm_map_dict'] = farm_map_dict
        cache['geo_mask'] = geo_mask
    # End try

    for gauge in Stream_gauges:
        swgw_exchanges[gauge] = modflow_model.getRivFluxNodes(river_reach_cells[gauge])
    # End for

    # Average depth to GW table:
    # Mask all cells that are either Coonambidgal or Shepparton formation
    # A mask could also be constructed for farm areas by using the mapped farms
    # from the groundwater model builder object
    for farm_zone in farm_zones:
        mask = geo_mask & (farm_map == farm_map_dict[int(farm_zone)])
        avg_depth_to_gw[farm_zone] = modflow_model.getAverageDepthToGW(mask=mask)
    # End for

    """
    Ecology heads of importance
    River gauges of importance for ecology: 406201, 406202, 406207, 406218, 406265
    Corresponding GW bores nearest to:      83003,  89586,  82999,  5662,   44828

    The bores are being ignored for nw as they are too far away from the stream
    gauges. Instead the simulated head in the cell with the stream gauge is
    used which represents the average head over the 5 km x 5 km cell.
    """

    heads = modflow_model.getHeads()
    for ind, ecol_bore in enumerate(Ecology_bores):
        _i, _h, _j = river_reach_ecol[Stream_gauges[ind]]
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

    for trigger_bore in Policy_bores:
        # NOTE: This returns the head in mAHD
        trigger_heads[trigger_bore] = modflow_model.getObservation(trigger_bore, 0, 'head')[0]
    # End for

    # TODO: Check that all of the wells listed were mapped to the model mesh and
    # are available for inspection
    return swgw_exchanges, avg_depth_to_gw, ecol_depth_to_gw, trigger_heads
# End run()


def main():
    print("Running from: " + os.getcwd())
    CONFIG = ConfigLoader(os.path.join(os.path.dirname(os.path.dirname(__file__)),
                                       "config", "model_config.json"))\
        .set_environment("GW_link_Integrated")

    def load_obj(filename):
        if filename[-4:] == '.pkl':
            try:
                with open(filename, 'r') as f:
                    f = pickle.load(f)
            except ImportError:
                with open(filename, 'rb') as f:
                    f = pickle.load(f)
            # End if
        elif filename.endswith(".npy"):
            f = np.load(filename)
        else:
            print 'File type not recognised as "pkl"'
        # End if

        return f

    # Example river level data (to be inputted from SW Model)
    try:
        fname = "initial_river_levels.pkl"
        riv_stages = load_obj(os.path.join(CONFIG.settings['data_folder'], fname))
    except (IOError):
        fname = "initial_river_levels.npy"
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
        data_folder = model_config['data_folder']
        mf_exe_folder = model_config['mf_exe_folder']
        param_file = model_config['param_file']
    # End if

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
        "pumping": 0.1,
        "MM": MM,
        "verbose": False,
        "is_steady": True
    }

    run(**run_params)
    run_params.update({"is_steady": False})
    for i in xrange(100):
        run(**run_params)
    # End for
    swgw_exchanges, avg_depth_to_gw, ecol_depth_to_gw, trigger_heads = run(**run_params)
    print "swgw_exchanges", swgw_exchanges
    print "avg_depth_to_gw", avg_depth_to_gw
    print "ecol_depth_to_gw", ecol_depth_to_gw
    print "trigger_heads", trigger_heads


if __name__ == "__main__":
    main()
