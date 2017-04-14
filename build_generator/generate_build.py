import os
from script_pieces import *


def generate_build_script(script_name, *args):

    script = ""
    for i in args:
        script += i
    # End for

    if not script_name.endswith(".py"):
        script_name = "{}.py".format(script_name)
    # End if

    if not os.path.isdir("gen_builds"):
        os.mkdir("gen_builds")
    # End if

    with open("gen_builds/{}".format(script_name), 'w') as output:
        output.write(script)
    # End with

if __name__ == '__main__':
    verbose = False

    project = "GW_link_Integrated"

    # TODO: MAKE CLIMATE PATH GENERIC
    generate_build_script("test_build.py",
                          s_imports(),
                          s_config(project),
                          s_basic_model_params(28355),
                          s_gen_modelbuilder(prj_name=project,
                                             model_type="Modflow",
                                             mesh_type="structured",
                                             climate_path="C:/Workspace/part0075/MDB modelling/Campaspe_data/Climate/"
                                             ),
                          s_set_boundaries("GW_model_area.shp", 20000),
                          s_weather_stations(['Kyneton', 'Elmore', 'Rochester', 'Echuca'],
                                             "rain_processed",
                                             "Rain_gauges.shp"),
                          s_bore_data(),
                          s_hydrogeo_properties(),
                          s_include_c14(),
                          s_process_river_stations(),
                          s_load_gw_boundary_shp(),
                          s_generate_mesh(),
                          s_interp_rainfall_to_grid(),
                          s_create_recharge_boundary(),
                          s_map_bores_to_grid(),
                          s_create_c14_obs_well(),
                          s_map_pumping_wells_to_grid(),
                          s_create_pumping_boundary(),
                          s_map_river_to_grid(),
                          s_create_river_boundary(),
                          s_map_river_to_grid2(),
                          s_create_river_boundary2(),
                          s_create_GHB_boundary(),
                          s_collate_obs(),
                          s_package_model()
                          )
