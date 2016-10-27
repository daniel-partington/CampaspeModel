"""Update parameters and run Steady State model for Campaspe."""

import os
import sys

from HydroModelBuilder.HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface

# Configuration Loader
from HydroModelBuilder.HydroModelBuilder.Utilities.Config.ConfigLoader import CONFIG
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
    print " Instantiate MODFLOW model "

    ###########################################################################
    ###########################################################################
    ###########################################################################
    # Currently using flopyInterface directly rather than running from the ModelManager ...
    modflow_model = flopyInterface.ModflowModel(MM.GW_build[name], data_folder=data_folder)

    modflow_model.executable = mf_exe_folder

    modflow_model.buildMODFLOW()

    modflow_model.buildMT3D()
    
    modflow_model.runMT3D()
    
    modflow_model.viewConcsByZone()
    
    return modflow_model

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
        play = run(model_folder, data_folder, mf_exe_folder, param_file=param_file)
    else:
        run(model_folder, data_folder, mf_exe_folder)
