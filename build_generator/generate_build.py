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
                          s_load_bores(),
                          s_load_pumping_wells(),
                          s_hydrogeo_properties(),
                          s_load_carbon14_data(),
                          s_process_river_stations(),
                          s_load_river_shp(),
                          s_load_gw_boundary_shp(),
                          s_build_mesh(),
                          s_assign_mesh_properties(),
                          s_interp_rainfall_to_grid(),
                          s_create_recharge_boundary(),
                          s_map_bores_to_grid(),
                          s_set_initial_heads(),
                          s_the_rest()
                          )
