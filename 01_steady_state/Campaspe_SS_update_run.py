"""Update parameters and run Steady State model for Campaspe."""

import os
import sys
import numpy as np

from HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface

# Configuration Loader
from HydroModelBuilder.Utilities.Config.ConfigLoader import CONFIG
# import flopy.utils.binaryfile as bf

sys.path.append('C:\Workspace\part0075\GIT_REPOS')

# MM is short for model manager


def run(model_folder, data_folder, mf_exe_folder, param_file=None):
    """Model Runner."""
    def loadObj(model_folder, model_name, filename):
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
                model_folder, filename))
        except IOError:
            model_obj = MM.GW_build[model_name].polyline_mapped[filename_no_ext + ".shp"]
        # End try

        return model_obj
    # End loadObj()

    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder, r"01_steady_state_packaged.pkl"))

    name = MM.GW_build.keys()[0]

    # Load in the new parameters based on parameters.txt or dictionary of new parameters

    #if param_file:
    #    MM.GW_build[name].updateModelParameters(os.path.join(data_folder, 'parameters.txt'))
    
    print "************************************************************************"
    print " Updating HGU parameters "

    # This needs to be automatically generated from with the map_raster2mesh routine ...
    zone_map = {1: 'qa', 2: 'utb', 3: 'utqa', 4: 'utam', 5: 'utaf', 6: 'lta', 7: 'bse'}

    Zone = MM.GW_build[name].model_mesh3D[1].astype(float)
    Kh = MM.GW_build[name].model_mesh3D[1].astype(float)
    Kv = MM.GW_build[name].model_mesh3D[1].astype(float)
    Sy = MM.GW_build[name].model_mesh3D[1].astype(float)
    SS = MM.GW_build[name].model_mesh3D[1].astype(float)
    for key in zone_map.keys():
        # if key ==7:
        #    continue
        Kh[Zone == key] = MM.GW_build[name].parameters.param['Kh_' + zone_map[key]]['PARVAL1']
        Kv[Zone == key] = MM.GW_build[name].parameters.param['Kv_' + zone_map[key]]['PARVAL1']
        Sy[Zone == key] = MM.GW_build[name].parameters.param['Sy_' + zone_map[key]]['PARVAL1']
        SS[Zone == key] = MM.GW_build[name].parameters.param['SS_' + zone_map[key]]['PARVAL1']

    MM.GW_build[name].properties.assign_model_properties('Kh', Kh)
    MM.GW_build[name].properties.assign_model_properties('Kv', Kv)
    MM.GW_build[name].properties.assign_model_properties('Sy', Sy)
    MM.GW_build[name].properties.assign_model_properties('SS', SS)

    print "************************************************************************"
    print " Updating river parameters "

    mapped_river = loadObj(model_folder, name, r"Campaspe_Riv_model.shp_mapped.pkl")

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
            MM.GW_build[name].parameters.param['bed_depress']['PARVAL1']
        cond = riv_cell[1] * riv_width_avg * \
            MM.GW_build[name].parameters.param['Kv_riv']['PARVAL1'] / riv_bed_thickness
        simple_river += [[0, row, col, stage, cond, bed]]

    riv = {}
    riv[0] = simple_river
    # MM.GW_build[name].boundaries.create_model_boundary_condition('Campaspe
    # River', 'river', bc_static=True)
    MM.GW_build[name].boundaries.assign_boundary_array('Campaspe River', riv)

    mapped_river = loadObj(model_folder, name, r"River_Murray_model.shp_mapped.pkl")

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

    for i in [1, 2, 3, 7]:
        interp_rain[MM.GW_build[name].model_mesh3D[1][0] == i] = interp_rain[MM.GW_build[name].model_mesh3D[
            1][0] == i] * MM.GW_build[name].parameters.param['rch_red_' + zone_map[i]]['PARVAL1']

    for i in [4, 5, 6, ]:
        interp_rain[MM.GW_build[name].model_mesh3D[1][0] == i] = interp_rain[
            MM.GW_build[name].model_mesh3D[1][0] == i] * 0.

    rch = {}
    rch[0] = interp_rain

    # MM.GW_build[name].boundaries.create_model_boundary_condition('Rain_reduced',
    # 'recharge', bc_static=True)
    MM.GW_build[name].boundaries.assign_boundary_array('Rain_reduced', rch)

    # print " Include irrigation in the recharge array"

    # print MM.GW_build[name].observations.obs_group['head'].keys()
    # print MM.GW_build[name].observations.obs_group['head']['mapped_observations']

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
    print " Updating Western GHB boundary"

    mapped_west = loadObj(model_folder, name, r"western_head_model.shp_mapped.pkl")

    WestGHB = []
    for WestGHB_cell in mapped_west:
        row = WestGHB_cell[0][0]
        col = WestGHB_cell[0][1]
        for lay in range(MM.GW_build[name].model_mesh3D[1].shape[0]):
            if MM.GW_build[name].model_mesh3D[1][0][row][col] == -1:
                continue
            WestGHBstage = MM.GW_build[name].model_mesh3D[0][0][row][
                col] + MM.GW_build[name].parameters.param['WGHB_stage']['PARVAL1']
            dx = MM.GW_build[name].gridHeight
            dz = MM.GW_build[name].model_mesh3D[0][lay][row][col] - \
                MM.GW_build[name].model_mesh3D[0][lay + 1][row][col]
            WGHBconductance = dx * dz * MM.GW_build[name].parameters.param['WGHBcond']['PARVAL1']
            WestGHB += [[lay, row, col, WestGHBstage, WGHBconductance]]

    ghb[0] += WestGHB

    print "************************************************************************"
    print " Updating Eastern GHB boundary"

    mapped_east = loadObj(model_folder, name, r"eastern_head_model.shp_mapped.pkl")

    EastGHB = []
    for EastGHB_cell in mapped_east:
        row = EastGHB_cell[0][0]
        col = EastGHB_cell[0][1]
        for lay in range(MM.GW_build[name].model_mesh3D[1].shape[0]):
            if MM.GW_build[name].model_mesh3D[1][0][row][col] == -1:
                continue
            EastGHBstage = MM.GW_build[name].model_mesh3D[0][0][row][
                col] + MM.GW_build[name].parameters.param['EGHB_stage']['PARVAL1']
            dx = MM.GW_build[name].gridHeight
            dz = MM.GW_build[name].model_mesh3D[0][lay][row][col] - \
                MM.GW_build[name].model_mesh3D[0][lay + 1][row][col]
            EGHBconductance = dx * dz * MM.GW_build[name].parameters.param['EGHBcond']['PARVAL1']
            EastGHB += [[lay, row, col, EastGHBstage, EGHBconductance]]

    ghb[0] += EastGHB

    print "************************************************************************"
    print " Updating GHB boundary"

    MM.GW_build[name].boundaries.assign_boundary_array('GHB', ghb)

#    print "************************************************************************"
#    print " Set initial head "
#
#    path=os.path.join(data_folder)
#    fname="01_steady_stateMainProcess"
#    headobj = bf.HeadFile(path + fname +'.hds')
#    times = headobj.get_times()
#    head = headobj.get_data(totim=times[-1])
#
#    MM.GW_build[name].initial_conditions.set_as_initial_condition("Head", head)

    print "************************************************************************"
    print " Build and run MODFLOW model "

    ###########################################################################
    ###########################################################################
    ###########################################################################
    # Currently using flopyInterface directly rather than running from the ModelManager ...
    modflow_model = flopyInterface.ModflowModel(MM.GW_build[name], data_folder=data_folder)

    modflow_model.executable = mf_exe_folder

    modflow_model.buildMODFLOW()

    modflow_model.checkMODFLOW()

    modflow_model.runMODFLOW()

    # print " Return the stream-aquifer exchange for reaches as list "
    #
    # SWGWexchange = [1]
    #
    # print " Return the average depth to the GW table in areas as list "

    # AvgDepthToGWTable = 1
    # DepthToGWTable = [1]

    modflow_model.checkCovergence()

   
    modflow_model.writeObservations()

    # modflow_model.viewHeadsByZone()

    riv_exch = modflow_model.getRiverFlux('Campaspe River')
    for key in riv_exch.keys():
        print 'Campaspe River net flux: ' + str(round(sum([x[0] for x in riv_exch[key]]))) + ' m3/d'

    riv_exch = modflow_model.getRiverFlux('Murray River')
    for key in riv_exch.keys():
        print 'Murray River net flux: ' + str(round(sum([x[0] for x in riv_exch[key]]))) + ' m3/d'

    # ts = MM.GW_build[name].observations.obs_group['head']['time_series']
    # wells_of_interest = ['79234', '62589']
    # wells = {}
    # for well in wells_of_interest:
    #    wells[well] = ts[ts['name'] == well]
    #    wells[well] = ts[ts['name'] == well]

    # print modflow_model.getObservation(wells[wells_of_interest[1]]['obs_map'].tolist()[0], 0, 'head')
    # print modflow_model.getObservation(wells_of_interest[1], 0, 'head')

    #modflow_model.viewHeads()
    
    modflow_model.viewHeadsByZone()

    #modflow_model.viewHeads2()

    modflow_model.waterBalance()
    
    #modflow_model.buildMT3D()
    
    #modflow_model.runMT3D()
    
    #modflow_model.viewConcsByZone()
   

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
