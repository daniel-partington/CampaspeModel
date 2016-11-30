"""
Run Campaspe GW model using inputs from Farm model and SW model and 
then return SW/GW exchanges, avg depth to GW, depth to GW at ecology sites and
head at trigger bores.
"""

import os
import sys
import numpy as np

import flopy.utils.binaryfile as bf

from HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface

# Configuration Loader
from HydroModelBuilder.Utilities.Config.ConfigLoader import CONFIG

# TODO: Set the stream gauges, ecology bores, policy bores at the start in some
# other class or in here but so they are available in the run function.


def run(model_folder, data_folder, mf_exe_folder, param_file=None, riv_stages=None,
        rainfall_irrigation=None, pumping=None):
    """
    GW Model Runner

    :param riv_stages: np rec array fo gauge number and stage


    """

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
    MM.load_GW_model(os.path.join(model_folder, r"01_steady_state_packaged.pkl"))

    name = MM.GW_build.keys()[0]

    # Load in the new parameters based on parameters.txt or dictionary of new parameters
    if param_file:
        MM.GW_build[name].updateModelParameters(os.path.join(data_folder, 'parameters.txt'))

    # This needs to be automatically generated from with the map_raster2mesh routine ...
#    zone_map = {1: 'qa', 2: 'utb', 3: 'utqa', 4: 'utam', 5: 'utaf', 6: 'lta', 7: 'bse'}

    print "************************************************************************"
    print " Updating river parameters "

    Campaspe_river = loadObj(data_folder, name, r"Campaspe_Riv_model.shp_mapped.pkl")

    simple_river = []
    riv_width_avg = 10.0  # m
    riv_bed_thickness = 0.10  # m
    for riv_cell in Campaspe_river:  # MM.GW_build[name].polyline_mapped['Campaspe_Riv_model.shp']:
        row = riv_cell[0][0]
        col = riv_cell[0][1]
        if MM.GW_build[name].model_mesh3D[1][0][row][col] == -1:
            continue

    # TODO: Update stage data with interpolated values based on riv_stages passed in to function

        stage = MM.GW_build[name].model_mesh3D[0][0][row][col] - 0.01
        bed = MM.GW_build[name].model_mesh3D[0][0][row][col] - 0.1 - \
            MM.GW_build[name].parameters.param['bed_depress']['PARVAL1']
        cond = riv_cell[1] * riv_width_avg * \
            MM.GW_build[name].parameters.param['Kv_riv']['PARVAL1'] / riv_bed_thickness
        simple_river += [[0, row, col, stage, cond, bed]]

    riv = {}
    riv[0] = simple_river
    # MM.GW_build[name].boundaries.create_model_boundary_condition('Campaspe
    # River', 'river', bc_static=True)

    MM.GW_build[name].boundaries.assign_boundary_array('Campaspe River', riv)

    mapped_river = loadObj(data_folder, name, r"River_Murray_model.shp_mapped.pkl")

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

    modflow_model.runMODFLOW()

    modflow_model.checkCovergence()

    """
    SW-GW exchanges:
    """

    sw_stream_gauges = [406214, 406219, 406201, 406224, 406218, 406202, 406265]
    swgw_exchanges = np.recarray((1,), dtype=[(str(gauge), np.float) for gauge in sw_stream_gauges])

    # TODO: Generate segment subsets of river nodes list for each gauging station
    # and use the following function already written to grab the flux for the
    # chosen nodes
    #riv_exch = modflow_model.getRiverFlux('Campaspe River')

    """
    Average depth to GW table:

    """

    farm_zones = [1]
    avg_depth_to_gw = np.recarray(
        (1,), dtype=[(str(farm_zone), np.float) for farm_zone in farm_zones])

    # TODO: Get average depth to GW using:
    # 1. Active cells heads in the first layer of the model
    # 2. Surface elevation in at the top of those active cells:
    # 3. Average the difference between the two arrays

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
       groundwater extractions. So we’d need to make sure we have baseflow 
       represented in the river between Lake Eppalock  (which I think is an 
       input to the ecology model so already covered).
   
    2. The reference trigger bore for the Barnadown zone was selected to 
       represent the gradient of groundwater flow between the Campaspe and the 
       Murray. So can we have the gradient of flow as an indicator in the 
       integrated model as well (if it isn’t already)?
    """

    trigger_head_bores = ['79234', '62589']
    trigger_heads = np.recarray((1, ), dtype=[(bore, np.float) for bore in trigger_head_bores])
    # to set
    for trigger_bore in trigger_head_bores:
        trigger_heads[trigger_bore] = np.random.rand()
    # end for

    # TODO: Check that all of the wells listed were mapped to the model mesh and
    # are available for inspection

    return swgw_exchanges, avg_depth_to_gw, ecol_depth_to_gw, trigger_heads

if __name__ == "__main__":
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
        run(model_folder, data_folder, mf_exe_folder, param_file=param_file)
    else:
        run(model_folder, data_folder, mf_exe_folder)
