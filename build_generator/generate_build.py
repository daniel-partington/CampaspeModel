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

    # Map river variable name to its shapefile
    riv_to_shp_map = {"Campaspe": "Campaspe_Riv.shp",
                      "Murray": "River_Murray.shp"}

    gw_boundary_shps = {"WGWbound": "western_head.shp",
                        "EGWbound": "eastern_head.shp"}

    # TODO: MAKE CLIMATE PATH GENERIC
    generate_build_script("test_build.py",
                          s_imports(),
                          s_config(project),
                          s_basic_model_params(28355),
                          s_gen_modelbuilder(prj_name=project,
                                             model_type="Modflow",
                                             mesh_type="structured",
                                             study_area="campaspe"
                                             ),
                          s_model_linkage(),
                          s_set_boundaries("GW_model_area.shp", 20000),
                          s_weather_stations(['Kyneton', 'Elmore', 'Rochester', 'Echuca'],
                                             "rain_processed",
                                             "Rain_gauges.shp"),
                          s_bore_data("Groundwater licence information for Dan Partington bc301115.xlsx",
                                      ),
                          s_hydrogeo_properties("Hydrogeologic_variables.xlsx"),
                          # s_include_c14("C14_locs.xlsx"),
                          s_process_river_stations(riv_to_shp_map),
                          s_load_gw_boundary_shp(gw_boundary_shps),
                          s_generate_mesh(),
                          s_interp_rainfall_to_grid(),
                          s_create_recharge_boundary(),
                          s_map_bores_to_grid(),
                          # s_create_c14_obs_well(),
                          s_map_pumping_wells_to_grid(),
                          s_create_pumping_boundary(),
                          # ['406201', '406203', '406218', '406202', '406265']
                          s_map_river_to_grid(),
                          s_create_river_boundary(),
                          s_map_river_to_grid2(),
                          s_create_river_boundary2(),
                          s_create_GHB_boundary(),
                          s_collate_obs(),
                          s_package_model()
                          )
