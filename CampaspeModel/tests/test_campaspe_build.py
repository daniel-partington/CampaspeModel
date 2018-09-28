import subprocess
import os
import sys

import argparse
import datetime

import cPickle as pickle

from osgeo import osr

from CampaspeModel.build_common import campaspe_data, campaspe_mesh
from HydroModelBuilder.GISInterface.GDALInterface.GDALInterface import \
    GDALInterface
from HydroModelBuilder.GWModelBuilder import GWModelBuilder
from HydroModelBuilder.Utilities.Config.ConfigLoader import ConfigLoader


def build_setup(forecast_run):

    try:
        tr_model = getattr(build_setup, 'tr_model')
    except AttributeError:

        p_j = os.path.join
        CONFIG = ConfigLoader(p_j('..', 'config', 'model_config.json')).set_environment("GW_link_Integrated")

        model_mode = 'forecast' if forecast_run else 'hindcast'

        # Define basic model parameters:
        proj_cs = osr.SpatialReference()

        # 28355 is the code for gda94 mga zone 55; http://spatialreference.org/ref/epsg/gda94-mga-zone-55/
        EPSG_CODE = 28355
        proj_cs.ImportFromEPSG(EPSG_CODE)

        interface = GDALInterface()
        interface.projected_coordinate_system = proj_cs
        interface.pcs_EPSG = "EPSG:{}".format(EPSG_CODE)

        get_conf_set = CONFIG.get_setting
        model_config = CONFIG.model_config

        temp_data_path = get_conf_set(['model_build', 'temp_data'])
        input_data_path = get_conf_set(['model_build', 'input_data'])
        campaspe_data_folder = p_j(temp_data_path, "Campaspe_data")
        build_setup.river_path = p_j(input_data_path, "Waterways")
        build_setup.sw_data_path = p_j(temp_data_path, "Campaspe_data/SW/All_streamflow_Campaspe_catchment/Updated")
        build_setup.hu_raster_path = p_j(campaspe_data_folder, "ESRI_GRID_raw", "ESRI_GRID")
        build_setup.campaspe_data_folder = campaspe_data_folder

        build_setup.bore_levels_file = "bore_levels"
        build_setup.bore_info_file = "bore_info"
        model_build_input_path = get_conf_set(['model_build', 'input_data'])

        build_setup.data_folder = model_config['data_folder']
        build_setup.mf_exe_folder = model_config['mf_exe_folder']
        build_setup.param_file = model_config['param_file']
        build_setup.climate_path = p_j(get_conf_set(['model_build', 'campaspe_data']), 'Climate')

        model_params = {
            "name": "GW_link_Integrated",
            "data_folder": model_build_input_path,
            "campaspe_data": get_conf_set(['model_build', 'campaspe_data']),
            "model_data_folder": model_config['data_folder'] + model_mode,
            "out_data_folder": get_conf_set(['model_build', 'data_build']),
            "GISInterface": interface,
            "model_type": "Modflow",
            "mesh_type": "structured"
        }
        tr_model = GWModelBuilder(**model_params)
        build_setup.tr_model = tr_model
    # End try

    return tr_model
# End build_setup()

def build_props(tr_model, campaspe_data_folder):
    custom_data = \
        campaspe_data.process_custom_scripts_and_spatial_data(tr_model,
                                                              campaspe_data_folder,
                                                              verbose=True,
                                                              GW_link_Integrated=True)

    return custom_data

    # Returned elements and corresponding variables used in build script
    # hgu_props = custom_data['HGU_props']
    # rain_gauges = custom_data['rain_gauges']
    # long_term_historic_weather = custom_data['long_term_historic_weather']
    # recharge_zones = custom_data['recharge_zones']
    # recharge_info = custom_data['recharge_zone_info_detailed']
    # surface_raster_high_res = custom_data['surface_raster_high_res']
    # surface_raster_high_res_GSA = custom_data['surface_raster_high_res_GSA']
    # river_gauges = custom_data['river_gauges']
    # campaspe_river_poly_file = custom_data['Campaspe_river_poly_file']
    # murray_river_poly_file = custom_data['Murray_river_poly_file']
    # river_stage_data = custom_data['river_stage_data']
    # river_flow_data = custom_data['river_flow_data']  # Never used!
    # campaspe = custom_data['Campaspe']
    # campaspe_field_elevations = custom_data['Campaspe_field_elevations']
    # campaspe_relevant = custom_data['Campaspe_relevant']
    # bores_shpfile = custom_data['bores_shpfile']
    # final_bores = custom_data['final_bores']
    # pumping_data = custom_data['pumping_data']
    # pumps_points = custom_data['pumps_points']
    # bore_data_info = custom_data['bore_data_info']
    # bore_data_levels = custom_data['bore_data_levels']
    # farms_poly = custom_data['farms_poly']
# End build_props()

def test_build_mesh():
    tr_model = build_setup(forecast_run=True)
    custom_data = build_props(tr_model, build_setup.campaspe_data_folder)
    res = 5000

    hgu, hu_raster_files_reproj = \
        campaspe_mesh.build_mesh_and_set_properties(tr_model,
                                                    build_setup.hu_raster_path,
                                                    custom_data['HGU_props'],
                                                    resolution=int(res)
                                                    )

    with open('data/build_test/hgu_dump.pkl', 'rb') as f_hgu:
        hgu_data = pickle.load(f_hgu)

    with open('data/build_test/hgu_raster_files_reproj_dump.pkl', 'rb') as f_rpj:
        reproj = pickle.load(f_rpj)

    assert hgu == hgu_data, "Mesh build mismatch! Something has changed!"
    assert hu_raster_files_reproj == reproj, "HU Raster file collection mismatch! Something has changed!"
# End test_build_mesh()


if __name__ == '__main__':
    tr_model = build_setup(forecast_run=True)

    test_build_mesh()
