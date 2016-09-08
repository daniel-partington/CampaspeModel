from HydroModelBuilder.HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface
# MM is short for model manager


def run_GW_model(testMM, river_stages, rainfall, data_folder):

    # Load in the new parameters based on parameters.txt or dictionary of new parameters
    #testMM.updateParameters('Campaspe', 'paramaters.txt')

    print "************************************************************************"
    print " Updating HGU parameters "

    print "************************************************************************"
    print " Updating recharge boundary "

    old_val = testMM.GW_build['Campaspe'].parameters.param['magic_rain']['PARVAL1']
    new_val = 0.05

    rch = {}
    rch[0] = testMM.GW_build['Campaspe'].boundaries.bc[
        'Rain_reduced']['bc_array'][0] / old_val * new_val

    testMM.GW_build['Campaspe'].boundaries.create_model_boundary_condition(
        'Rain_reduced', 'recharge', bc_static=True)
    testMM.GW_build['Campaspe'].boundaries.assign_boundary_array('Rain_reduced', rch)

    # print " Include irrigation in the recharge array"

    print " Updating river conditions "

    print "************************************************************************"
    print " Build and run MODFLOW model "

    ###########################################################################
    ###########################################################################
    ###########################################################################
    # Currently using flopyInterface directly rather than running from the ModelManager ...
    modflow_model = flopyInterface.ModflowModel(
        testMM.GW_build['Campaspe'], data_folder=data_folder)

    modflow_model.runMODFLOW()

    print " Return the stream-aquifer exchange for reaches as list "

    SWGWexchange = [1]

    print " Return the average depth to the GW table in areas as list "

    AvgDepthToGWTable = 1
    DepthToGWTable = [1]

    modflow_model.viewHeads()
    modflow_model.viewHeads2()

    return SWGWexchange, AvgDepthToGWTable, DepthToGWTable


if __name__ == "__main__":
    # To do: add irrigation inputs to modify recharge boundary condition
    river_stages = 1
    rainfall = 1

    testMM = GWModelManager()
    # testMM.load_GW_model(r"C:\Workspace\part0075\GIT_REPOS\CampaspeModel\Campaspe_packaged.pkl")
    testMM.load_GW_model(
        r"C:/UserData/takuyai/ownCloud/CampaspeModel/CampaspeIntegrated/Campaspe_packaged.pkl")

    data_folder = r"C:/UserData/takuyai/ownCloud/CampaspeModel"

    # modify data folder
    testMM.GW_build['Campaspe'].data_folder = data_folder
    # modify output folder
    # testMM.GW_build['Campaspe'].out_data_folder

    run_GW_model(testMM, river_stages, rainfall, data_folder)
