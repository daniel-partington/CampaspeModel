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


def run(model_folder, data_folder, mf_exe_folder, farm_zones=None, param_file=None, riv_stages=None,
        rainfall_irrigation=None, pumping=None, verbose=True, MM=None, is_steady=False):
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

    # ******************************************************************************
    # ******************************************************************************
    #
    # Complimentary models requirements, i.e. bore and gauge data that should be
    # referenceable to this model for parsing specific outputs and receiving inputs:
    if not hasattr(MM, 'ext_linkage_bores') or MM.ext_linkage_bores is None:
        # Read in bores that relate to external models
        model_linking = os.path.join(data_folder, "model_linking.csv")
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

    # ******************************************************************************
    # ******************************************************************************

    if MM is None:
        MM = GWModelManager()
        MM.load_GW_model(os.path.join(model_folder, "GW_link_Integrated_packaged.pkl"))
    # End if

    name = MM.GW_build.keys()[0]

    this_model = MM.GW_build[name]
    mesh = this_model.model_mesh3D
    mesh_0 = mesh[0]
    mesh_1 = mesh[1]

    # Load in the new parameters based on parameters.txt or dictionary of new parameters
    if param_file:
        this_model.updateModelParameters(os.path.join(data_folder, 'parameters.txt'),
                                         verbose=verbose)

    # This needs to be automatically generated from with the map_raster2mesh routine ...
    # zone_map = {1: 'qa', 2: 'utb', 3: 'utqa', 4: 'utam', 5: 'utaf', 6: 'lta', 7: 'bse'}

    if verbose:
        print "************************************************************************"
        print " Updating river parameters "

    # loadObj(data_folder, name, r"Campaspe_Riv_model.shp_mapped.pkl")
    Campaspe_river = this_model.polyline_mapped['Campaspe_Riv_model.shp']

    simple_river = []
    riv_width_avg = 10.0  # m
    riv_bed_thickness = 0.10  # m

    # Need to account for ~8m drop after Campaspe weir.

    Campaspe_river_gauges = this_model.points_mapped['processed_river_sites_stage_clipped.shp']

    filter_gauges = [riv_gauge for riv_gauge in Campaspe_river_gauges
                     if str(riv_gauge[1][0]) in Stream_gauges]

    simple_river = []
    riv_width_avg = 10.0  # m
    riv_bed_thickness = 0.10  # m

    # Map river from high to low
    new_riv = Campaspe_river
    for index, riv_cell in enumerate(Campaspe_river):
        row = riv_cell[0][0]
        col = riv_cell[0][1]
        new_riv[index] += [mesh_0[0][row][col]]

    new_riv = sorted(new_riv, key=lambda x: (x[0][1]), reverse=False)
    new_riv = sorted(new_riv, key=lambda x: (x[0][0]), reverse=True)

    stages = np.full((len(new_riv)), np.nan, dtype=np.float64)
    beds = np.full((len(new_riv)), np.nan, dtype=np.float64)

    # Identify cells that correspond to river gauges
    riv_gauge_logical = np.full((len(new_riv)), False, dtype=np.bool)

    # To account for fact that river shapefile and gauges shapefile are not perfect
    # we get the closest river cell to the gauge cell

    # Define river gauges at start of river cell
    new_riv_cells = [x[0] for x in new_riv]
    filter_gauge_loc = [new_riv_cells[x]
                        for x in [closest_node(x[0], new_riv_cells) for x in filter_gauges]]

    for index, riv in enumerate(new_riv):
        # Create logical array to identify those which are gauges and those which are not
        if riv[0] in filter_gauge_loc:
            riv_gauge_logical[index] = True
            gauge_ind = [i for i, x in enumerate(filter_gauge_loc) if x == riv[0]]

            stages[index] = filter_gauges[gauge_ind[0]][1][0]
            beds[index] = stages[index] - 1.0

        # Add chainage to new_riv array:
        if index == 0:
            new_riv[index] += [0.0]
        else:
            new_riv[index] += [new_riv[index - 1][3] + new_riv[index - 1][1]]

    # River x in terms of chainage:
    river_x = np.array([x[3] for x in new_riv])
    river_x_unknown = river_x[~riv_gauge_logical]
    river_x_known = river_x[riv_gauge_logical]

    # Now interpolate know values of stage and bed to unknown river locations:
    stages[~riv_gauge_logical] = np.interp(
        river_x_unknown, river_x_known, stages[riv_gauge_logical])
    beds[~riv_gauge_logical] = np.interp(river_x_unknown, river_x_known, beds[riv_gauge_logical])

    # Need to create list of relevant riv nodes for each reach by
    # splitting the riv nodes list up

    riv_reach_nodes = {}
    for index, gauge in enumerate(Stream_gauges):
        if index == 0:
            split_loc = new_riv_cells.index(filter_gauge_loc[index])
            riv_reach_nodes[gauge] = new_riv_cells[0:split_loc + 1]
        else:
            split_loc_upstream = new_riv_cells.index(filter_gauge_loc[index - 1])
            split_loc_downstream = new_riv_cells.index(filter_gauge_loc[index])
            riv_reach_nodes[gauge] = new_riv_cells[split_loc_upstream:split_loc_downstream + 1]
        # End if
    # End for

    for index, riv_cell in enumerate(Campaspe_river):
        row = riv_cell[0][0]
        col = riv_cell[0][1]

        if mesh_1[0][row][col] == -1:
            continue
        stage_temp = stages[index]
        if stage_temp < mesh_0[1][row][col]:
            stage = mesh_0[1][row][col] + 0.01
        else:
            stage = stage_temp
        # End if

        bed_temp = beds[index]
        if bed_temp < mesh_0[1][row][col]:
            bed = mesh_0[1][row][col]
        else:
            bed = bed_temp
        # End if

        cond = riv_cell[1] * riv_width_avg * \
            this_model.parameters.param['Kv_riv']['PARVAL1'] / riv_bed_thickness
        simple_river += [[0, row, col, stage, cond, bed]]

    riv = {0: simple_river}

    this_model.boundaries.assign_boundary_array('Campaspe River', riv)

    # Updating Murray River

    # loadObj(data_folder, name, r"River_Murray_model.shp_mapped.pkl")
    mapped_river = this_model.polyline_mapped['River_Murray_model.shp']

    simple_river = []
    riv_width_avg = 10.0  # m
    riv_bed_thickness = 0.10  # m
    for riv_cell in mapped_river:
        cell = riv_cell[0]
        row = cell[0]
        col = cell[1]

        if mesh_1[0][row][col] == -1:
            continue
        # End if

        stage = mesh_0[0][row][col] - 0.01
        bed = mesh_0[0][row][col] - 0.1 - \
            this_model.parameters.param['RMstage']['PARVAL1']
        cond = riv_cell[1] * riv_width_avg * \
            this_model.parameters.param['Kv_RM']['PARVAL1'] / riv_bed_thickness
        simple_river += [[0, row, col, stage, cond, bed]]

    riv = {0: simple_river}

    this_model.boundaries.assign_boundary_array('Murray River', riv)

    if verbose:
        print "************************************************************************"
        print " Updating recharge boundary "

    # Adjust rainfall to recharge using 10% magic number
    if rainfall_irrigation is not None:
        interp_rain = rainfall_irrigation
    else:
        warnings.warn("Rainfall+Irrigation input currently ignored by GW model")
        interp_rain = np.copy(this_model.boundaries.bc['Rainfall']['bc_array'])

        for i in [1, 2, 3, 7]:
            match = interp_rain[mesh_1[0] == i]
            # match = interp_rain[mesh[
            #     1][0] == i] * 0.1
            interp_rain[mesh_1[0] == i] = match * 0.05

        for i in [4, 5, 6]:
            interp_rain[mesh_1[0] == i] = 0
    # End if

    rch = {0: interp_rain}
    this_model.boundaries.assign_boundary_array('Rain_reduced', rch)

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
    for MurrayGHB_cell in mapped_river:
        row = MurrayGHB_cell[0][0]
        col = MurrayGHB_cell[0][1]

        MMparams = this_model.parameters.param
        for lay in xrange(mesh_1.shape[0]):
            if mesh_1[0][row][col] == -1:
                continue

            MurrayGHBstage = mesh_0[0][row][col] + MMparams['MGHB_stage']['PARVAL1']
            dx = this_model.gridHeight
            dz = mesh_0[lay][row][col] - mesh_0[lay + 1][row][col]
            MGHBconductance = dx * dz * MMparams['MGHBcond']['PARVAL1']
            MurrayGHB += [[lay, row, col, MurrayGHBstage, MGHBconductance]]

    ghb = {0: MurrayGHB}

    if verbose:
        print "************************************************************************"
        print " Updating GHB boundary"

    this_model.boundaries.assign_boundary_array('GHB', ghb)

    if verbose:
        print "************************************************************************"
        print " Setting initial head "

    fname = 'model_{}'.format(name)
    if os.path.exists(os.path.join(data_folder, fname, name) + '.hds'):
        headobj = bf.HeadFile(os.path.join(data_folder, fname, name) + '.hds')
        times = headobj.get_times()
        head = headobj.get_data(totim=times[-1])
        this_model.initial_conditions.set_as_initial_condition("Head", head)
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
        modflow_model.checkCovergence()

    # modflow_model.waterBalance()

#    print("Campaspe River flux NET: ", np.array([x[0] for x in modflow_model.getRiverFlux('Campaspe River')[0]]).sum())
#    print("Campaspe River flux +'ve: ", np.array([x[0] for x in modflow_model.getRiverFlux('Campaspe River')[0] if x[0] > 0.0]).sum())
#    print("Campaspe River flux -'ve: ", np.array([x[0] for x in modflow_model.getRiverFlux('Campaspe River')[0] if x[0] < 0.0]).sum())
#
#    camp_riv_flux = modflow_model.getRiverFlux('Campaspe River')[0]
#    import matplotlib.pyplot as plt
#
#    plt.figure()
#
#    points = [x[1] for x in camp_riv_flux]
#
#    cmap = plt.get_cmap('spectral')
#
#    x = [ x[2] for x in points ]
#    y = [ y[1] for y in points ]
#    vals = [v[0] for v in camp_riv_flux]
#    ax = plt.scatter(x, y, c=vals, cmap=cmap, vmin=-1000, vmax= 1000)
#
#    plt.colorbar(ax)

    # modflow_model.viewHeads2()

    # modflow_model.waterBalance(plot=False)
    # modflow_model.viewHeads2()

    """
    SW-GW exchanges:
    """
    swgw_exchanges = np.recarray((1,), dtype=[(str(gauge), np.float) for gauge
                                              in Stream_gauges])

    for gauge in Stream_gauges:
        swgw_exchanges[gauge] = modflow_model.getRiverFluxNodes(
            riv_reach_nodes[gauge])

    if verbose:
        print("Upstream of weir", swgw_exchanges[
              sw_stream_gauges[0]] + swgw_exchanges[sw_stream_gauges[1]])
        print("Downstream of weir", swgw_exchanges[
              sw_stream_gauges[2]] + swgw_exchanges[sw_stream_gauges[3]])

    # Average depth to GW table:
    avg_depth_to_gw = np.recarray(
        (1,), dtype=[(farm_zone, np.float) for farm_zone in farm_zones])

    tgt_mesh = modflow_model.model_data.model_mesh3D[1][0]

    # Mask all cells that are either Coonambidgal or Shepparton formation
    geo_mask = (tgt_mesh == 3) | (tgt_mesh == 1)
    # A mask could also be constructed for farm areas by using the mapped farms
    # from the groundwater model builder object
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
        # ecol_depth_to_gw[ecol_bore] = modflow_model.getObservation(ecol_bore, 0, 'head')[1]
        ecol_depth_to_gw[ecol_bore] = np.random.rand()
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
    # to set
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
    # print "swgw_exchanges", swgw_exchanges
    # print "avg_depth_to_gw", avg_depth_to_gw
    # print "ecol_depth_to_gw", ecol_depth_to_gw
    # print "trigger_heads", trigger_heads
