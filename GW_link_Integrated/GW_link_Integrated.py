"""
Run Campaspe GW model using inputs from Farm model and SW model and
then return SW/GW exchanges, avg depth to GW, depth to GW at ecology sites and
head at trigger bores.
"""

import os
import sys
import warnings

import numpy as np

import flopy.utils.binaryfile as bf

from HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface

# TODO: Set the stream gauges, ecology bores, policy bores at the start in some
# other class or in here but so they are available in the run function.

# Configuration Loader
from HydroModelBuilder.Utilities.Config.ConfigLoader import ConfigLoader


def run(model_folder, data_folder, mf_exe_folder, param_file=None, riv_stages=None,
        rainfall_irrigation=None, pumping=None):
    """
    GW Model Runner.

    :param model_folder: str, path to model folder
    :param data_folder: str, path to data
    :param mf_exe_folder: str, path to MODFLOW executable
    :param param_file: str, path to parameter file
    :param riv_stages: np recarray, gauge numbers and stage
    :param rainfall_irrigation: np array, array representing rainfall and irrigation input.
                                Must match the model extent.
    :param pumping: dict, {Farm Zone ID: daily pumping amount (m/day)}

    """

    warnings.warn("This function uses hardcoded values for Farm Zones and SW Gauges")

    def loadObj(data_folder, model_name, filename):
        """
        Interface to Model Manager object loader.

        Attempts to load model object from alternate source when file cannot be found.

        :param model_folder: Folder where the model is
        :param model_name: Name of the model to load object from
        :param filename: Filename of picked object.
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

    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder, r"GW_link_Integrated_packaged.pkl"))

    name = MM.GW_build.keys()[0]

    # Load in the new parameters based on parameters.txt or dictionary of new parameters
    if param_file:
        MM.GW_build[name].updateModelParameters(os.path.join(data_folder, 'parameters.txt'))

    # This needs to be automatically generated from with the map_raster2mesh routine ...
#    zone_map = {1: 'qa', 2: 'utb', 3: 'utqa', 4: 'utam', 5: 'utaf', 6: 'lta', 7: 'bse'}

    print "************************************************************************"
    print " Updating river parameters "

    # loadObj(data_folder, name, r"Campaspe_Riv_model.shp_mapped.pkl")
    Campaspe_river = MM.GW_build[name].polyline_mapped['Campaspe_Riv_model.shp']

    simple_river = []
    riv_width_avg = 10.0  # m
    riv_bed_thickness = 0.10  # m

    #sw_stream_gauges = [406214, 406219, 406201, 406224, 406218, 406202, 406265]

    sw_stream_gauges = ['406201', '406218', '406202', '406265']

    # Need to account for ~8m drop after Campaspe weir.

    Campaspe_river_gauges = MM.GW_build[name].points_mapped[
        'processed_river_sites_stage_clipped.shp']

    filter_gauges = []
    for riv_gauge in Campaspe_river_gauges:
        # if riv_gauge[1][0] in use_gauges:
        if str(riv_gauge[1][0]) in sw_stream_gauges:
            filter_gauges += [riv_gauge]

    simple_river = []
    riv_width_avg = 10.0  # m
    riv_bed_thickness = 0.10  # m

    # Map river from high to low
    new_riv = Campaspe_river
    for index, riv_cell in enumerate(Campaspe_river):
        row = riv_cell[0][0]
        col = riv_cell[0][1]
        new_riv[index] += [MM.GW_build[name].model_mesh3D[0][0][row][col]]

    new_riv = sorted(new_riv, key=lambda x: (x[0][1]), reverse=False)
    new_riv = sorted(new_riv, key=lambda x: (x[0][0]), reverse=True)

    stages = np.full((len(new_riv)), np.nan, dtype=np.float64)
    beds = np.full((len(new_riv)), np.nan, dtype=np.float64)

    # Identify cells that correspond to river gauges
    riv_gauge_logical = np.full((len(new_riv)), False, dtype=np.bool)

    # To account for fact that river shapefile and gauges shapefile are not perfect
    # we get the closest river cell to the gauge cell

    def closest_node(node, nodes):
        nodes = np.asarray(nodes)
        dist_2 = np.sum((nodes - node)**2, axis=1)
        return np.argmin(dist_2)

    # Define river gauges at start of river cell
    new_riv_cells = [x[0] for x in new_riv]
    filter_gauge_loc = [new_riv_cells[x]
                        for x in [closest_node(x[0], new_riv_cells) for x in filter_gauges]]

    #river_stage_file = os.path.join(data_folder, r"dev_river_levels_recarray.pkl")
    #river_stage_data = MM.GW_build[name].load_obj(river_stage_file)

    for index, riv in enumerate(new_riv):
        # Create logical array to identify those which are gauges and those which are not
        if riv[0] in filter_gauge_loc:
            riv_gauge_logical[index] = True
            gauge_ind = [i for i, x in enumerate(filter_gauge_loc) if x == riv[0]]
#            location =  str(filter_gauges[gauge_ind[0]][1][0])
            stages[index] = riv_stages[riv_stages["g_id"] ==
                                       str(filter_gauges[gauge_ind[0]][1][0])][0][1]
#            stages[index] = river_stage_data["Mean stage (m)"].loc[river_stage_data["Site ID"] == filter_gauges[gauge_ind[0]][1][0]]
            # river_stage_data["Mean stage (m)"].loc[river_stage_data["Site ID"]== ??]
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

    # 1st order the gauges from upstream to downstream:
    ordered_gauges = ['406201', '406218', '406202', '406265']

    riv_reach_nodes = {}
    for index, gauge in enumerate(ordered_gauges):
        if index == 0:
            split_loc = new_riv_cells.index(filter_gauge_loc[index])
            riv_reach_nodes[gauge] = new_riv_cells[0:split_loc + 1]
        else:
            split_loc_upstream = new_riv_cells.index(filter_gauge_loc[index - 1])
            split_loc_downstream = new_riv_cells.index(filter_gauge_loc[index])
            riv_reach_nodes[gauge] = new_riv_cells[split_loc_upstream:split_loc_downstream + 1]

    #adjusted = 0
    #kept = 0
    #riv_vis = np.zeros_like(MM.GW_build[name].model_mesh3D[0][0])
    #riv_vis_on = np.zeros_like(MM.GW_build[name].model_mesh3D[0][0])

    for index, riv_cell in enumerate(Campaspe_river):
        row = riv_cell[0][0]
        col = riv_cell[0][1]
        if MM.GW_build[name].model_mesh3D[1][0][row][col] == -1:
            continue
        stage_temp = stages[index]
        if stage_temp < MM.GW_build[name].model_mesh3D[0][1][row][col]:
            stage = MM.GW_build[name].model_mesh3D[0][1][row][col] + 0.01
            #print("Adjusting stage by: ", stage-stage_temp)
            #adjusted += 1
            #riv_vis_on[row][col] = -1
            #riv_vis_on[row][col] = stage
        else:
            stage = stage_temp
            #print("Keeping stage at: ", row, col)
            #kept += 1
            #riv_vis_on[row][col] = 1
            #riv_vis[row][col] = stage
        bed_temp = beds[index]
        if bed_temp < MM.GW_build[name].model_mesh3D[0][1][row][col]:
            bed = MM.GW_build[name].model_mesh3D[0][1][row][col]
        else:
            bed = bed_temp

        cond = riv_cell[1] * riv_width_avg * \
            MM.GW_build[name].parameters.param['Kv_riv']['PARVAL1'] / riv_bed_thickness
        simple_river += [[0, row, col, stage, cond, bed]]

    #print("For river cells, kept: ", kept)
    #print("For river cells, adjusted: ", adjusted)

    #import matplotlib.pyplot as plt
    # plt.imshow(riv_vis_on)

    riv = {}
    riv[0] = simple_river
    # MM.GW_build[name].boundaries.create_model_boundary_condition('Campaspe
    # River', 'river', bc_static=True)

    MM.GW_build[name].boundaries.assign_boundary_array('Campaspe River', riv)
    print("Model Boundary")
    print(MM.GW_build[name].model_boundary)

    # Updating Murray River

    # loadObj(data_folder, name, r"River_Murray_model.shp_mapped.pkl")
    mapped_river = MM.GW_build[name].polyline_mapped['River_Murray_model.shp']

    simple_river = []
    riv_width_avg = 10.0  # m
    riv_bed_thickness = 0.10  # m
    for riv_cell in mapped_river:  # MM.GW_build[name].polyline_mapped['Campaspe_Riv_model.shp']:
        row = riv_cell[0][0]
        col = riv_cell[0][1]
        if MM.GW_build[name].model_mesh3D[1][0][row][col] == -1:
            continue
        # print test_model.model_mesh3D
        stage = MM.GW_build[name].model_mesh3D[0][0][row][col] - 0.01
        bed = MM.GW_build[name].model_mesh3D[0][0][row][col] - 0.1 - \
            MM.GW_build[name].parameters.param['RMstage']['PARVAL1']
        cond = riv_cell[1] * riv_width_avg * \
            MM.GW_build[name].parameters.param['Kv_RM']['PARVAL1'] / riv_bed_thickness
        simple_river += [[0, row, col, stage, cond, bed]]

    riv = {}
    riv[0] = simple_river
    # MM.GW_build[name].boundaries.create_model_boundary_condition('Campaspe
    # River', 'river', bc_static=True)
    MM.GW_build[name].boundaries.assign_boundary_array('Murray River', riv)

    print "************************************************************************"
    print " Updating recharge boundary "

    # Adjust rainfall to recharge using 10% magic number
    interp_rain = np.copy(MM.GW_build[name].boundaries.bc['Rainfall']['bc_array'])

    # TODO: Replace interp_rain with input from farm model, i.e. rainfall_irrigation

    for i in [1, 2, 3, 7]:
        interp_rain[MM.GW_build[name].model_mesh3D[1][0] == i] = interp_rain[MM.GW_build[name].model_mesh3D[
            1][0] == i] * 0.1  # MM.GW_build[name].parameters.param['rch_red_' + zone_map[i]]['PARVAL1']

    for i in [4, 5, 6, ]:
        interp_rain[MM.GW_build[name].model_mesh3D[1][0] == i] = interp_rain[
            MM.GW_build[name].model_mesh3D[1][0] == i] * 0.

    rch = {}
    rch[0] = interp_rain

    MM.GW_build[name].boundaries.assign_boundary_array('Rain_reduced', rch)

    print "************************************************************************"
    print " Updating pumping boundary"

    pumpy = MM.GW_build[name].boundaries.bc['licenced_wells']['bc_array']
    wel = {key: [[b[0], b[1], b[2], b[3] * pumping] for b in a] for key, a in pumpy.iteritems()}

    MM.GW_build[name].boundaries.assign_boundary_array('licenced_wells', wel)

    print "************************************************************************"
    print " Updating Murray River GHB boundary"

    MurrayGHB = []
    for MurrayGHB_cell in mapped_river:
        row = MurrayGHB_cell[0][0]
        col = MurrayGHB_cell[0][1]
        # print MM.GW_build[name].model_mesh3D
        for lay in range(MM.GW_build[name].model_mesh3D[1].shape[0]):
            if MM.GW_build[name].model_mesh3D[1][0][row][col] == -1:
                continue
            MurrayGHBstage = MM.GW_build[name].model_mesh3D[0][0][row][
                col] + MM.GW_build[name].parameters.param['MGHB_stage']['PARVAL1']
            dx = MM.GW_build[name].gridHeight
            dz = MM.GW_build[name].model_mesh3D[0][lay][row][col] - \
                MM.GW_build[name].model_mesh3D[0][lay + 1][row][col]
            MGHBconductance = dx * dz * MM.GW_build[name].parameters.param['MGHBcond']['PARVAL1']
            MurrayGHB += [[lay, row, col, MurrayGHBstage, MGHBconductance]]

    ghb = {}
    ghb[0] = MurrayGHB

    print "************************************************************************"
    print " Updating GHB boundary"

    MM.GW_build[name].boundaries.assign_boundary_array('GHB', ghb)

    print "************************************************************************"
    print " Set initial head "

    # TODO: Update head based on last iteration rather than initial head
    fname = "initial"
    headobj = bf.HeadFile(os.path.join(data_folder, fname) + '.hds')

    times = headobj.get_times()
    head = headobj.get_data(totim=times[-1])

    MM.GW_build[name].initial_conditions.set_as_initial_condition("Head", head)

    print "************************************************************************"
    print " Build and run MODFLOW model "

    ###########################################################################
    ###########################################################################
    ###########################################################################
    # Currently using flopyInterface directly rather than running from the ModelManager ...
    modflow_model = flopyInterface.ModflowModel(MM.GW_build[name], data_folder=data_folder)

    # Override temporal aspects of model build:
    modflow_model.nper = 1  # This is the number of stress periods which is set to 1 here
    modflow_model.perlen = 1  # This is the period of time which is set to 1 day here
    modflow_model.nstp = 1  # This is the number of sub-steps to do in each stress period
    modflow_model.steady = False  # This is to tell FloPy that is a transient model

    modflow_model.executable = mf_exe_folder

    modflow_model.buildMODFLOW()

    # modflow_model.runMODFLOW()

    modflow_model.checkCovergence()

    modflow_model.waterBalance()
    # modflow_model.viewHeads2()

    """
    SW-GW exchanges:
    """

    swgw_exchanges = np.recarray((1,), dtype=[(str(gauge), np.float) for gauge
                                              in sw_stream_gauges])

    for gauge in sw_stream_gauges:
        swgw_exchanges[gauge] = modflow_model.getRiverFluxNodes(
            riv_reach_nodes[gauge])

    """
    Average depth to GW table:

    """

    farm_zones = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']
    avg_depth_to_gw = np.recarray(
        (1,), dtype=[(str(farm_zone), np.float) for farm_zone in farm_zones])

    for farm_zone in farm_zones:
        # The mask below just chooses all cells that are either Coonambidgal or
        # Shepparton formation. A mask could also be constructed for farm areas
        # by using the
        mask = (modflow_model.model_data.model_mesh3D[1][0] == 3) | (
            modflow_model.model_data.model_mesh3D[1][0] == 1)
        avg_depth_to_gw[str(farm_zone)] = modflow_model.getAverageDepthToGW(mask=mask)

    """
    Ecology heads of importance

    River gauges of importance for ecology: 406201, 406202, 406207, 406218, 406265
    Corresponding GW bores nearest to:      83003,  89586,  82999,  5662,   44828

    """

    ecol_depth_to_gw_bores = ['83003', '89586', '82999', '5662', '44828']
    ecol_depth_to_gw = np.recarray((1,), dtype=[(bore, np.float)
                                                for bore in ecol_depth_to_gw_bores])
    # to set
    for ecol_bore in ecol_depth_to_gw_bores:
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

    #trigger_head_bores = ['79234', '62589']
    # Removed bore 79234 as it sits above Lake Eppalock which is not included
    # in the groudnwater model.

    trigger_head_bores = ['62589']
    trigger_heads = np.recarray((1, ), dtype=[(bore, np.float) for bore in trigger_head_bores])
    # to set
    for trigger_bore in trigger_head_bores:
        trigger_heads[trigger_bore] = np.random.rand()
    # end for

    # TODO: Check that all of the wells listed were mapped to the model mesh and
    # are available for inspection

    return swgw_exchanges, avg_depth_to_gw, ecol_depth_to_gw, trigger_heads

if __name__ == "__main__":

    CONFIG = ConfigLoader(os.path.join('../config', "model_config.json"))\
        .set_environment("GW_link_Integrated")

    print CONFIG.model_config

    import pickle

    def load_obj(filename):
        if filename[-4:] == '.pkl':
            with open(filename, 'rb') as f:
                return pickle.load(f)
        else:
            print 'File type not recognised as "pkl"'
        # end if

    #folder = r"C:\Workspace\part0075\GIT_REPOS\CampaspeModel\testbox\integrated\data"
    folder = r"C:/UserData/takuyai/ownCloud/CampaspeModel/testbox/integrated/data"
    fname = r"dev_river_levels_recarray.pkl"

    riv_stages = load_obj(os.path.join(folder, fname))

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

    run_params = {
        "model_folder": model_folder,
        "data_folder": data_folder,
        "mf_exe_folder": mf_exe_folder,
        "param_file": param_file if param_file else None,
        "riv_stages": riv_stages,
        "rainfall_irrigation": None,
        "pumping": None,
    }

    result = run(**run_params)
