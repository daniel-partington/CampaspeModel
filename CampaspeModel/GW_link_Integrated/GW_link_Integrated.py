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


def loadObj(data_folder, model_name, filename):
    """
    Interface to Model Manager object loader.

    Attempts to load model object from alternate source when file cannot be found.

    :param model_folder: str, Folder where the model is
    :param model_name: str, Name of the model to load object from
    :param filename: str, Filename of picked object.
                     Attempts to load from shapefile with same name on exception.
    """
    filename_no_ext = os.path.splitext(filename)[0].split('.')[0]

    try:
        model_obj = MM.GW_build[model_name].load_obj(os.path.join(
            data_folder, filename))
    except IOError:
        model_obj = MM.GW_build[model_name].polyline_mapped[filename_no_ext + ".shp"]
    # End try

    return model_obj
# End loadObj()


def closest_node(node, nodes):
    nodes = np.asarray(nodes)
    dist_2 = np.sum((nodes - node)**2, axis=1)
    return np.argmin(dist_2)
# End closest_node()


def process_line(line):
    processed = [x.strip() for x in line.split(':')[1].strip().split(',')]
    return processed
# end process_line


def run(model_folder, data_folder, mf_exe_folder, farm_zones=None, param_file=None, riv_stages=None,
        rainfall_irrigation=None, pumping=None, verbose=True, MM=None, recharge_amt=0.03, is_steady=False):
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

    warnings.warn("This function uses hardcoded values for Farm Zones and SW Gauges")

    p_j = os.path.join

    if MM is None:
        MM = GWModelManager()
        MM.load_GW_model(p_j(model_folder, "GW_link_Integrated_packaged.pkl"))
    # End if

    # Complimentary models requirements, i.e. bore and gauge data that should be
    # referenceable to this model for parsing specific outputs and receiving inputs:
    if not hasattr(MM, 'ext_linkage_bores') or MM.ext_linkage_bores is None:

        # Read in bores that relate to external models
        model_linking = p_j(data_folder, "model_linking.csv")
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
    else:
        Ecology_bores = MM.ext_linkage_bores["Ecology_bores"]
        Policy_bores = MM.ext_linkage_bores["Policy_bores"]
        Stream_gauges = MM.ext_linkage_bores["Stream_gauges"]
    # End if

    name = MM.GW_build.keys()[0]
    this_model = MM.GW_build[name]
    mesh = this_model.model_mesh3D
    mesh_0 = mesh[0]
    mesh_1 = mesh[1]

    # Load in the new parameters based on parameters.txt or dictionary of new parameters
    if param_file:
        this_model.updateModelParameters(p_j(data_folder, 'parameters.txt'), verbose=verbose)
    # End if

    if verbose:
        print "************************************************************************"
        print " Updating river parameters "

    river_seg = this_model.river_mapping['Campaspe']

    num_reaches = this_model.pilot_points['Campaspe'].num_points  # 4
    known_points = this_model.pilot_points['Campaspe'].points

    strcond_val = [this_model.parameters.param['kv_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)]
    river_seg['strhc1'] = np.interp(river_seg['Cumulative Length'].tolist(),
                                    known_points, strcond_val)

    river_seg.loc[:, 'stage_from_gauge'] = np.nan

    gauge_val_dict = {}
    for gauge_id in riv_stages.dtype.names:
        if gauge_id in Stream_gauges:
            val = riv_stages[gauge_id]
            gauge_val_dict[gauge_id] = val

    Campaspe_stage = river_seg[river_seg['gauge_id'] != 'none']
    new_riv_stages = []
    for row in Campaspe_stage.iterrows():
        ind = row[1]['gauge_id']
        try:
            print gauge_val_dict[str(ind)]
            new_riv_stages += [gauge_val_dict[str(ind)]]
        except:
            print("No value for: {}".format(ind))
            print river_seg[river_seg['gauge_id'] == ind]['stage']
            new_riv_stages += river_seg[river_seg['gauge_id'] == ind]['stage'].tolist()

    Campaspe_stage.loc[:, 'stage_from_gauge'] = new_riv_stages

    river_seg.loc[river_seg['iseg'].isin(Campaspe_stage['iseg'].tolist()),
                  'stage_from_gauge'] = \
        sorted(Campaspe_stage['stage_from_gauge'].tolist(),
               reverse=True)
#
    river_seg['stage_from_gauge'] = \
        river_seg.set_index(river_seg['Cumulative Length'])['stage_from_gauge']. \
        interpolate(method='values', limit_direction='both').tolist()

    #river_seg['stage_from_gauge'] = river_seg['stage_from_gauge'].bfill()

    river_seg['stage'] = river_seg['stage_from_gauge']

    simple_river = []

    for row in river_seg.iterrows():
        row = row[1]
        simple_river += [[row['k'], row['i'], row['j'], row['stage'],
                          row['strhc1'] * row['rchlen'] * row['width1'], row['strtop']]]
    # end for

    this_model.boundaries.update_boundary_array('Campaspe River', {0: simple_river})

    # Updating Murray River

    mriver_seg = this_model.river_mapping['Murray']
    mriver_seg['strhc1'] = this_model.parameters.param['kv_rm']['PARVAL1']

    msimple_river = []

    for row in mriver_seg.iterrows():
        row = row[1]
        msimple_river += [[row['k'], row['i'], row['j'], row['stage'],
                           row['strhc1'] * row['rchlen'] * row['width1'], row['strtop']]]
    # end for

    this_model.boundaries.update_boundary_array('Murray River', {0: msimple_river})

    if verbose:
        print "************************************************************************"
        print " Updating recharge boundary "

    if rainfall_irrigation is not None:
        interp_rain = np.copy(rainfall_irrigation)
    else:
        warnings.warn("Rainfall+Irrigation input ignored by GW model")
        interp_rain = np.copy(this_model.boundaries.bc['Rainfall']['bc_array'])
    # End if

    # integers reflect number of layers
    for i in [1, 2, 3, 4, 5, 6, 7]:
        match = mesh_1[0] == i

        if i in [4, 5, 6]:
            interp_rain[match] = 0
            continue
        # End if

        # Adjust rainfall to recharge using a magic number (defaults to 0.03 -> 3%)
        interp_rain[match] = interp_rain[match] * recharge_amt
    # End for

    this_model.boundaries.assign_boundary_array('Rain_reduced', {0: interp_rain})

    if verbose:
        print "************************************************************************"
        print " Updating pumping boundary"

    pumpy = this_model.boundaries.bc['licenced_wells']['bc_array']
    wel = {key: [[b[0], b[1], b[2], b[3] * pumping] for b in a] for key, a in pumpy.iteritems()}

    this_model.boundaries.assign_boundary_array('licenced_wells', wel)

    if verbose:
        print "************************************************************************"
        print " Updating Murray River GHB boundary"

    MurrayGHB = []

    MurrayGHB_cells = [[x[0], x[1], x[2], x[3]] for x in this_model.boundaries.bc['GHB']['bc_array'][0]]
    for MurrayGHB_cell in MurrayGHB_cells:
        lay, row, col = MurrayGHB_cell[:3]
        MurrayGHBstage = MurrayGHB_cell[3]  # m.parameters.param['mghb_stage']['PARVAL1']
        dx = this_model.gridHeight
        dz = this_model.model_mesh3D[0][lay][row][col] - \
            this_model.model_mesh3D[0][lay + 1][row][col]
        MGHBconductance = dx * dz * this_model.parameters.param['mghbk']['PARVAL1']  # / 10000.
        MurrayGHB += [[lay, row, col, MurrayGHBstage, MGHBconductance]]
    # End for
    ghb = {}
    ghb[0] = MurrayGHB

    if verbose:
        print "************************************************************************"
        print " Updating GHB boundary"

    this_model.boundaries.assign_boundary_array('GHB', {0: MurrayGHB})

    if verbose:
        print "************************************************************************"
        print " Setting initial head "

    fname = 'model_{}'.format(name)
    if os.path.exists(p_j(data_folder, fname, name + '.hds')):
        try:
            headobj = bf.HeadFile(p_j(data_folder, fname, name + '.hds'))
            times = headobj.get_times()
            head = headobj.get_data(totim=times[-1])
            this_model.initial_conditions.set_as_initial_condition("Head", head)
        except IndexError as e:
            raise IndexError(
                "Corrupted MODFLOW hds file - check, replace, or remove {}".format(
                    p_j(data_folder, fname, name + '.hds')))
    else:
        if verbose:
            print "Using initial head conditions"
        # End if
    # End if

    if verbose:
        print "************************************************************************"
        print " Build and run MODFLOW model "

    ###########################################################################
    ###########################################################################
    ###########################################################################
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

    if is_steady is False:
        modflow_model.checkConvergence()

    # SW-GW exchanges:
    swgw_exchanges = np.recarray((1,), dtype=[(str(gauge), np.float) for gauge
                                              in Stream_gauges])

#    for gauge in Stream_gauges:
#        swgw_exchanges[gauge] = modflow_model.getRivFluxNodes(riv_reach_nodes[gauge])
#
#    if verbose:
#        print("Upstream of weir", swgw_exchanges[
#              sw_stream_gauges[0]] + swgw_exchanges[sw_stream_gauges[1]])
#        print("Downstream of weir", swgw_exchanges[
#              sw_stream_gauges[2]] + swgw_exchanges[sw_stream_gauges[3]])

    # Average depth to GW table:
    avg_depth_to_gw = np.recarray((1,), dtype=[(farm_zone, np.float) for farm_zone in farm_zones])
    tgt_mesh = modflow_model.model_data.model_mesh3D[1][0]

    # Mask all cells that are either Coonambidgal or Shepparton formation
    # A mask could also be constructed for farm areas by using the mapped farms
    # from the groundwater model builder object
    geo_mask = (tgt_mesh == 3) | (tgt_mesh == 1)
    farm_map, farm_map_dict = this_model.polygons_mapped['farm_v1_prj_model.shp']
    for farm_zone in farm_zones:
        mask = geo_mask & (farm_map == farm_map_dict[int(farm_zone)])
        avg_depth_to_gw[farm_zone] = modflow_model.getAverageDepthToGW(mask=mask)
    # End for

    """
    Ecology heads of importance
    River gauges of importance for ecology: 406201, 406202, 406207, 406218, 406265
    Corresponding GW bores nearest to:      83003,  89586,  82999,  5662,   44828
    """
    ecol_depth_to_gw_bores = Ecology_bores
    ecol_depth_to_gw = np.recarray((1,), dtype=[(bore, np.float)
                                                for bore in ecol_depth_to_gw_bores])
    # to set
    for ecol_bore in ecol_depth_to_gw_bores:
        # NOTE: This returns the depth to groundwater below the surface in metres
        ecol_depth_to_gw[ecol_bore] = modflow_model.getObservation(ecol_bore, 0, 'head')[1]
    # end for

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

    trigger_head_bores = Policy_bores
    trigger_heads = np.recarray((1, ), dtype=[(bore, np.float) for bore in trigger_head_bores])
    for trigger_bore in trigger_head_bores:
        # NOTE: This returns the head in mAHD
        trigger_heads[trigger_bore] = modflow_model.getObservation(trigger_bore, 0, 'head')[0]
    # end for

    # TODO: Check that all of the wells listed were mapped to the model mesh and
    # are available for inspection
    return swgw_exchanges, avg_depth_to_gw, ecol_depth_to_gw, trigger_heads


if __name__ == "__main__":

    print("Running from: " + os.getcwd())
    CONFIG = ConfigLoader(os.path.join(os.path.dirname(os.path.dirname(__file__)),
                                       "config", "model_config.json"))\
        .set_environment("GW_link_Integrated")

    def load_obj(filename):
        if filename[-4:] == '.pkl':
            with open(filename, 'r') as f:
                return pickle.load(f)
        else:
            print 'File type not recognised as "pkl"'
        # End if

    # Example river level data (to be inputted from SW Model)
    # folder = r"C:\Workspace\part0075\GIT_REPOS\CampaspeModel\testbox\integrated\data"
    fname = r"dev_river_levels.pkl"
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
        "pumping": 10.0,  # {'5': 10},
        "MM": MM,
        "verbose": False,
        "is_steady": False
    }

    swgw_exchanges, avg_depth_to_gw, ecol_depth_to_gw, trigger_heads = run(**run_params)
    print "swgw_exchanges", swgw_exchanges
    print "avg_depth_to_gw", avg_depth_to_gw
    print "ecol_depth_to_gw", ecol_depth_to_gw
    print "trigger_heads", trigger_heads
