import argparse
import datetime
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from osgeo import osr

from CampaspeModel.build_common import campaspe_data, campaspe_mesh, rivers
from CampaspeModel.build_common.drains import prepare_drain_data_for_model
from CampaspeModel.build_common.groundwater_boundary import \
    prepare_ghb_boundary_from_murray_data
from CampaspeModel.build_common.initial_head_from_bores import \
    generate_initial_head_from_bores
from CampaspeModel.build_common.pumping import prepare_pumping_data_for_model
from CampaspeModel.build_common.rainfall_recharge import \
    prepare_transient_rainfall_data_for_model
from CampaspeModel.build_utils.multifrequency_resampling import (resample_obs_time_series_to_model_data_index,
                                                                 resample_to_model_data_index)
from HydroModelBuilder.GISInterface.GDALInterface.GDALInterface import \
    GDALInterface
from HydroModelBuilder.GWModelBuilder import GWModelBuilder
from HydroModelBuilder.Utilities.Config.ConfigLoader import ConfigLoader


def parse_string_to_bool(opt):
    if opt.lower().startswith('t'):
        val = True
    elif opt.lower().startswith('f'):
        val = False
    else:
        raise ValueError('Invalid option, must be true or false')
    # End if

    return val
# End parse_string_to_bool()


def process_line(line):
    processed = [x.strip() for x in line.split(':')[1].strip().split(',')]
    return processed


def main():

    parser = argparse.ArgumentParser("Campaspe Model Build Processor")
    parser.add_argument("--forecast",
                        type=parse_string_to_bool,
                        default=True,
                        help="Run the build process to generate a forecasting model")
    parser.add_argument("--verbose", type=parse_string_to_bool, default=False, help="Output all status messages")
    args = parser.parse_args()

    VERBOSE = bool(args.verbose)
    forecast_run = args.forecast

    p_j = os.path.join
    CONFIG = ConfigLoader(p_j('..', '..', 'config', 'model_config.json')).set_environment("GW_link_Integrated")

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
    data_folder = model_config['data_folder']
    mf_exe_folder = model_config['mf_exe_folder']
    param_file = model_config['param_file']
    climate_path = p_j(get_conf_set(['model_build', 'campaspe_data']), 'Climate')

    temp_data_path = get_conf_set(['model_build', 'temp_data'])
    input_data_path = get_conf_set(['model_build', 'input_data'])
    river_path = p_j(input_data_path, "Waterways")
    sw_data_path = p_j(temp_data_path, "Campaspe_data/SW/All_streamflow_Campaspe_catchment/Updated")
    campaspe_data_folder = p_j(temp_data_path, "Campaspe_data")
    hu_raster_path = p_j(campaspe_data_folder, "ESRI_GRID_raw", "ESRI_GRID")

    bore_levels_file = "bore_levels"
    bore_info_file = "bore_info"
    model_build_input_path = get_conf_set(['model_build', 'input_data'])

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

    custom_data = \
        campaspe_data.process_custom_scripts_and_spatial_data(tr_model,
                                                              campaspe_data_folder,
                                                              verbose=True,
                                                              GW_link_Integrated=True)

    hgu_props = custom_data['HGU_props']
    rain_gauges = custom_data['rain_gauges']
    long_term_historic_weather = custom_data['long_term_historic_weather']
    recharge_zones = custom_data['recharge_zones']
    recharge_info = custom_data['recharge_zone_info_detailed']
    surface_raster_high_res = custom_data['surface_raster_high_res']
    surface_raster_high_res_GSA = custom_data['surface_raster_high_res_GSA']
    river_gauges = custom_data['river_gauges']
    campaspe_river_poly_file = custom_data['Campaspe_river_poly_file']
    murray_river_poly_file = custom_data['Murray_river_poly_file']
    river_stage_data = custom_data['river_stage_data']
    river_flow_data = custom_data['river_flow_data']  # Never used!
    campaspe = custom_data['Campaspe']
    campaspe_field_elevations = custom_data['Campaspe_field_elevations']
    campaspe_relevant = custom_data['Campaspe_relevant']
    bores_shpfile = custom_data['bores_shpfile']
    final_bores = custom_data['final_bores']
    pumping_data = custom_data['pumping_data']
    pumps_points = custom_data['pumps_points']
    bore_data_info = custom_data['bore_data_info']
    bore_data_levels = custom_data['bore_data_levels']
    farms_poly = custom_data['farms_poly']

    # ******************************************************************************
    # Complementary models requirements, i.e. bore and gauge data that should be
    # referenceable to this model for parsing specific outputs and receiving inputs:

    if VERBOSE:
        print("************************************************************************")
        print("Attempting to map these bores...")
    # End if

    model_linking = os.path.join(data_folder, "model_linking.csv")
    with open(model_linking, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.split(':')[0] == 'Ecology':
                ecology_bores = process_line(line)
                print(('Ecology: {}'.format(ecology_bores)))
            elif line.split(':')[0] == 'Policy':
                policy_bores = process_line(line)
                print(('Policy: {}'.format(policy_bores)))
            elif line.split(':')[0] == 'SW_stream_gauges':
                stream_gauges = process_line(line)
                print(('SW: {}'.format(stream_gauges)))
            # End if
        # End for
    # End with

    print(("""
    Attempting to match the following bores/gauges
    Ecology: {}
    Policy: {}
    Stream: {}
    If the above indicated bore/gauges are incorrect, modify the model_linking.csv file in:
    {}

    If the location/data for the indicated bores is found not to be appropriate, the closest 
    'best' bore will be used instead.
    """.format(ecology_bores, policy_bores, stream_gauges, data_folder)))
    print(("=" * 30))

    if VERBOSE:
        print('########################################################################')
        print('########################################################################')
        print('## Model specific building ')
        print('########################################################################')
        print('########################################################################')

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #@@@ CONSTRUCTION OF TIME PERIODS FOR MODEL @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    start_time_interest = datetime.date(2000, 1, 1)

    if forecast_run:
        start, end = datetime.date(2014, 1, 1), datetime.date(2015, 1, 1)
        frequencies = ['A']
        date_index = pd.date_range(start=start, end=end, freq=frequencies[0])
        date_group = [start, end]
        tr_model.model_time.set_temporal_components(steady_state=False,
                                                    start_time=start,
                                                    end_time=end,
                                                    date_index=date_index)
        # OVERRIDE model temporal components:
        tr_model.model_time.t['steps'] = 1
        tr_model.model_time.t['intervals'] = [end - start]

    else:
        # 1980 - 2000 Annual
        # 2000 - 2015 Monthly
        frequencies = ['A', 'M']

        start, end = datetime.date(1980, 1, 1), datetime.date(2015, 1, 1)
        date_index_start = pd.date_range(start=start, end=start_time_interest, freq=frequencies[0])
        date_index = date_index_start[:-1].append(pd.date_range(start=start_time_interest, end=end, freq=frequencies[1]))
        date_group = [start, start_time_interest, end]

        tr_model.model_time.set_temporal_components(steady_state=False, start_time=start,
                                                    end_time=end, date_index=date_index)

    res = model_config['grid_resolution']
    res = res[:-1] if res.endswith('m') else res

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    # Define the grid width and grid height for the model mesh which is stored as a multipolygon
    # shapefile GDAL object
    if VERBOSE:
        print("************************************************************************")
        print((" Defining structured mesh at {res}x{res}".format(res=res)))

    hgu, hu_raster_files_reproj = \
        campaspe_mesh.build_mesh_and_set_properties(tr_model,
                                                    hu_raster_path,
                                                    hgu_props,
                                                    resolution=int(res)
                                                    )

    if VERBOSE:
        print("************************************************************************")
        print(" Interpolating rainfall data to grid and time steps")

    interp_rain, interp_et, recharge_zone_array, rch_zone_dict = \
        prepare_transient_rainfall_data_for_model(tr_model,
                                                  recharge_zones,
                                                  recharge_info,
                                                  long_term_historic_weather,
                                                  date_index,
                                                  frequencies,
                                                  date_group,
                                                  start,
                                                  end,
                                                  rain_gauges)

    if VERBOSE:
        print("************************************************************************")
        print(" Creating recharge boundary ")

    tr_model.boundaries.create_model_boundary_condition('Rain_reduced', 'recharge', bc_static=False)
    tr_model.boundaries.associate_zonal_array_and_dict('Rain_reduced', recharge_zone_array, rch_zone_dict)
    tr_model.boundaries.assign_boundary_array('Rain_reduced', interp_rain)

    if VERBOSE:
        print("************************************************************************")
        print(" Mapping farm areas to grid")

    tr_model.map_polygon_to_grid(farms_poly, feature_name="ZoneID")

    if VERBOSE:
        print("************************************************************************")
        print(" Mapping bores to grid ")

    tr_model.map_points_to_grid(bores_shpfile, feature_id='HydroCode')

    # clip bores by those within the farm area:
    # import subprocess
    # from osgeo import ogr
    # def clip_points_data(MBO, filename, clipping_poly, path=None):
    #    """Read in point data from shapefile, e.g. observation bore locations, rainfall gauges
    #
    #    :param filename: filename for the point shapefile that is to be read in.
    #    :param path: Path of the files. (Default value = None)
    #
    #    :returns: GDAL file pointer
    #    """
    #    base = os.path.splitext(os.path.basename(filename))[0]
    #    new_f = base + "_farm.shp"
    #    new_file = os.path.join(MBO.out_data_folder, new_f)
    #    print(new_file)
    #    driver = ogr.GetDriverByName("ESRI Shapefile")
    #
    #    if os.path.isfile(new_file):
    #        print 'Using previously generated file: ' + new_file
    #    else:
    #        # Clip first using boundary polygon
    #        # target_srs = MBO.GISInterface.pcs_EPSG
    #
    ##        fn = os.path.join(MBO.out_data_folder, base + '_reproj.shp')
    ##        command = 'ogr2ogr -t_srs "' + target_srs + '" "' + fn + '" "' + filename + '"'
    # print command
    # try:
    ##            print(subprocess.check_output(command, shell=True))
    # except subprocess.CalledProcessError as e:
    ##            print("stdout output on error:\n" + e.output)
    #
    #        command = 'ogr2ogr -clipsrc "' + clipping_poly + '" "' + new_file + '" "' + filename + '" -f "ESRI Shapefile"'
    #        print command
    #        try:
    #            print(subprocess.check_output(command, shell=True))
    #        except subprocess.CalledProcessError as e:
    #            print("stdout output on error:\n" + e.output)
    #
    #    # End if
    #
    #    ds = driver.Open(new_file, 0)
    #    MBO.GISInterface._test_osgeo_load(ds, new_file)
    #    return ds
    # End clip_points_data()
    #
    #farm_bores_shpfile = clip_points_data(tr_model, bores_shpfile.name, farms_poly.name)

    # ****** ^^^^ The above is not working at the moment due to a bug that occurs when
    # polygon shapefiles generated through py GDAL are returned from a separate module in
    # which case the .shx and .prj are missing and the files are corrrupted. They need
    # to be closed and then loaded again. An issue already exists for this and in the
    # future I will remove the passing of GDAL objects around the classes to avoid this
    # error and instead just make a loading and closing class for polygonal shape files

    farm_bores_shpfile = tr_model.read_points_data(os.path.join(temp_data_path, 'Campaspe_data/SW/Farm/Farm_bores.shp'))
    tr_model.map_points_to_grid(farm_bores_shpfile, feature_id='HydroCode')

    farm_bores = []
    for bores in tr_model.points_mapped[farm_bores_shpfile.name.split('\\')[-1]]:
        farm_bores += bores[1]

    bores_more_filter = []
    bores_more_filter_policy = []
    bores_above_surface = []
    bores_below_top_of_bedrock = []
    bores_in_top_layer = []
    mesh3D_0 = tr_model.model_mesh3D[0]

    for bores in tr_model.points_mapped["Farm_bores_clipped.shp"]:
        row, col = bores[0][0], bores[0][1]
        for bore in bores[1]:
            try:
                bore_depth = bore_data_info.loc[bore, 'depth']
            except Exception as e:
                if bore in policy_bores:
                    print(e)
                    print(('Policy bore not in info: ', bore))
                continue
            # End try

            if bore_depth > mesh3D_0[0][row][col]:
                bores_above_surface += [bore]
                if bore in policy_bores:
                    print(('Policy bore above surf: ', bore))
                continue
            if bore_depth > mesh3D_0[1][row][col]:
                bores_in_top_layer += [bore]
            if bore_depth <= mesh3D_0[-2][row][col]:
                bores_below_top_of_bedrock += [bore]
                if bore in policy_bores:
                    print(('Policy bore in bedrock: {}'.format(bore)))
                    print(('Bore depth is at: {}'.format(bore_depth)))
                    print(('Bedrock top is at: {}'.format(mesh3D_0[-2][row][col])))
                    print('Using above cell in Deep Lead for head by redefining bore_depth')
                    final_bores.set_value(final_bores[final_bores['HydroCode'] == bore].index,
                                          'depth', mesh3D_0[-2][row][col] + 1.0)
                    bores_more_filter_policy += [bore]
                    bores_more_filter += [bore]
                continue
                # sys.exit('Halting model build due to bore not being mapped')
            if bore in policy_bores:
                bores_more_filter_policy += [bore]
            bores_more_filter += [bore]

    print(('Bores above the surface: {}'.format(len(bores_above_surface))))
    print(('Bores below top of bedrock: {}'.format(len(bores_below_top_of_bedrock))))
    print(('Final bores within aquifers: {}'.format(len(bores_more_filter))))

    final_bores = final_bores[final_bores["HydroCode"].isin(bores_more_filter)]
    bore_points = [[final_bores.loc[x, "Easting"], final_bores.loc[x, "Northing"]] for x in final_bores.index]

    bore_points_3d = final_bores[["HydroCode", "Easting", "Northing", "depth"]]
    bore_points_3d = bore_points_3d.set_index("HydroCode")


    def generate_bore_observations_for_model(bores_obs_time_series, name):
        # Modify into standard format for the GWModelBuilder class
        bores_obs_time_series = bores_obs_time_series.rename(columns={'HydroCode': 'name',
                                                                      'bore_date': 'datetime',
                                                                      'result': 'value'})

        bores_obs_time_series['datetime'] = pd.to_datetime(bores_obs_time_series['datetime'])

        # Kill all values where head observation is 0.0 indicating bad data that was missed
        # by quality control codes in pre-processing function getBoreData
        bores_obs_time_series = bores_obs_time_series[bores_obs_time_series['value'] > 5.]
        # Remove values that are false records, i.e. bores with a reading at 30/12/1899
        bores_obs_time_series = bores_obs_time_series[bores_obs_time_series['datetime']
                                                      != datetime.datetime(1899, 12, 30)]
        # Only consider head observations after 2010 to limit the size of the jacobian in PEST
        bores_obs_time_series = bores_obs_time_series[bores_obs_time_series['datetime']
                                                      > datetime.datetime(1980, 12, 30)]

        # Create outlier identifier, i.e multiplier for standard deviations from the mean
        # Using 4 here which is a rather large bound
        bores_obs_outliers_multiplier = 4.0

        bores_obs_group = bores_obs_time_series.groupby('name')
        bores_obs_means = bores_obs_group['value'].mean()
        bores_obs_stds = bores_obs_group['value'].std()

        bores_obs_cull = []
        for bore in bores_obs_time_series['name'].unique():
            bore_std = bores_obs_stds.loc[bore]
            mod_target = bores_obs_means.loc[bore]
            modifier = bores_obs_outliers_multiplier * bore_std
            upper = mod_target + modifier
            lower = mod_target - modifier
            for row in bores_obs_time_series[bores_obs_time_series['name'] == bore].iterrows():
                reading_value = row[1]['value']
                if reading_value < lower or reading_value > upper:
                    bores_obs_cull += [int(row[0])]

        # Now to remove those bores obs that have a reading greater or less than 5 * times sigma either side of the mean
        bores_obs_time_series = bores_obs_time_series[~bores_obs_time_series.index.isin(bores_obs_cull)]

        # Now to resample the time series to get the mean heads within a time interval
        bores_obs_time_series = resample_obs_time_series_to_model_data_index(
            bores_obs_time_series,
            date_index,
            frequencies,
            date_group,
            start,
            end,
            fill='none')

        # For the weights of observations we need to specify them as 1/sigma, where sigma is the standard deviation of measurement error
        tr_model.observations.set_as_observations(name,
                                                  bores_obs_time_series,
                                                  bore_points_3d,
                                                  domain='porous',
                                                  obs_type='head',
                                                  units='mAHD',
                                                  weights=1.0 / 0.2,
                                                  by_zone=True)

        return bores_obs_time_series


    bores_obs_time_series = bore_data_levels[bore_data_levels["HydroCode"].isin(
        final_bores["HydroCode"])]
    bores_obs_time_series = generate_bore_observations_for_model(bores_obs_time_series, 'farm_bores')

    bores_obs_time_series_policy = bore_data_levels[bore_data_levels["HydroCode"].isin(
        policy_bores)]
    bores_obs_time_series_policy = generate_bore_observations_for_model(bores_obs_time_series_policy, 'policy_bores')

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    # SETUP INITIAL CONDITIONS

    bores_obs_time_series_combined = pd.concat([bores_obs_time_series, bores_obs_time_series_policy])

    if not forecast_run:
        initial_heads_tr = generate_initial_head_from_bores(tr_model, bores_obs_time_series_combined,
                                                            final_bores,
                                                            time_max='1981',
                                                            interp_method='linear',
                                                            plot=False)
    else:
        initial_heads_tr = generate_initial_head_from_bores(tr_model, bores_obs_time_series_combined,
                                                            final_bores,
                                                            time_min='2014',
                                                            interp_method='linear',
                                                            plot=False)

    bores_in_layers = tr_model.map_points_to_raster_layers(
        bore_points, final_bores["depth"].tolist(), hu_raster_files_reproj)

    tr_model.initial_conditions.set_as_initial_condition("Head", initial_heads_tr)  # interp_heads[hu_raster_files[0]])

    if VERBOSE:
        print("************************************************************************")
        print(" Mapping pumping wells to grid ")

    wel = prepare_pumping_data_for_model(tr_model,
                                         pumps_points,
                                         start,
                                         start,
                                         end,
                                         date_index,
                                         pumping_data,
                                         frequencies,
                                         date_group
                                         )

    if VERBOSE:
        print("************************************************************************")
        print(" Creating pumping boundary ")

    tr_model.boundaries.create_model_boundary_condition('licenced_wells', 'wells', bc_static=True)
    tr_model.boundaries.assign_boundary_array('licenced_wells', wel)

    if VERBOSE:
        print("************************************************************************")
        print(" Mapping Campaspe river to grid")

    num_reaches = 20
    river_seg, reach_df, reach_data, known_points = \
        rivers.prepare_river_data_for_campaspe(tr_model,
                                               surface_raster_high_res,
                                               river_gauges,
                                               campaspe_river_poly_file,
                                               campaspe,
                                               campaspe_field_elevations,
                                               num_reaches=num_reaches,
                                               plot=True)

    campaspe_info = campaspe
    campaspe_info.index = campaspe_info['Site Id']

    camp_riv_cells = [x for x in zip(river_seg['i'], river_seg['j'])]

    criv, river_seg = rivers.create_riv_data_transient(tr_model,
                                                       river_seg,
                                                       river_stage_data,
                                                       campaspe_info,
                                                       date_index,
                                                       frequencies,
                                                       date_group,
                                                       start,
                                                       end)

    tr_model.river_mapping['Campaspe'] = river_seg

    if VERBOSE:
        print("************************************************************************")
        print(" Creating Campaspe river boundary")

    tr_model.boundaries.create_model_boundary_condition('Campaspe River',
                                                        'river',
                                                        bc_static=False)
    tr_model.boundaries.assign_boundary_array('Campaspe River', criv)

    if VERBOSE:
        print("************************************************************************")
        print(" Mapping Murray River to grid")

    riv, mriver_seg_ghb = \
        rivers.prepare_river_data_for_murray(tr_model,
                                             surface_raster_high_res_GSA,
                                             murray_river_poly_file,
                                             campaspe_relevant,
                                             river_stage_data,
                                             river_seg,
                                             plot=True)

    if VERBOSE:
        print("************************************************************************")
        print(" Creating Murray River boundary")

    tr_model.boundaries.create_model_boundary_condition('Murray River', 'river', bc_static=True)
    tr_model.boundaries.assign_boundary_array('Murray River', riv)

    if VERBOSE:
        print("************************************************************************")
        print(" Setting up Murray River GHB boundary")

    mriver_seg_ghb.loc[:, 'init_head'] = mriver_seg_ghb.apply(
        lambda x: initial_heads_tr[0][int(x['i'])][int(x['j'])], axis=1)

    ghb = prepare_ghb_boundary_from_murray_data(tr_model, mriver_seg_ghb)

    if VERBOSE:
        print("************************************************************************")
        print(" Creating GHB boundary")

    tr_model.boundaries.create_model_boundary_condition('GHB', 'general head', bc_static=True)
    tr_model.boundaries.assign_boundary_array('GHB', ghb)

    if VERBOSE:
        print("************************************************************************")
        print(" Mapping Drains to grid")

    drain = prepare_drain_data_for_model(tr_model,
                                         camp_riv_cells,
                                         start,
                                         date_index)

    if VERBOSE:
        print("************************************************************************")
        print(" Creating Drains boundary")

    tr_model.boundaries.create_model_boundary_condition('Drain', 'drain', bc_static=True)
    tr_model.boundaries.assign_boundary_array('Drain', drain)

    if VERBOSE:
        print("************************************************************************")
        print(" Collate observations")
    #
    tr_model.map_obs_loc2mesh3D(method='nearest', ignore=[-1, 7])
    tr_model.map_obs2model_times()
    tr_model.observations.collate_observations()

    obs_active_bores = bores_obs_time_series[bores_obs_time_series['zone'] != 'null']['name']
    obs_active_bores = obs_active_bores[obs_active_bores.isin(bores_in_top_layer)].tolist()
    obs_filter_bores = bore_points_3d[bore_points_3d.index.isin(obs_active_bores)]
    obs_bores_list = list(zip(obs_filter_bores['Easting'], obs_filter_bores['Northing']))

    stream_active = campaspe[campaspe.index.isin([int(x) for x in stream_gauges])]
    stream_gauges_list = list(zip(stream_active['Easting'], stream_active['Northing']))

    closest_bores_active = tr_model.find_closest_points_between_two_lists(obs_bores_list, stream_gauges_list)

    ecol_bores = []
    for ind in closest_bores_active:
        ecol_bores += [obs_filter_bores.index.tolist()[ind]]

    ecol_bores_df = obs_filter_bores[obs_filter_bores.index.isin(ecol_bores)]
    policy_found = [x for x in final_bores["HydroCode"] if x in policy_bores]

    # read in model link
    link_file = os.path.join(data_folder, "model_linking.csv")
    with open(link_file, 'r') as model_link:
        tmp = model_link.readlines()
        for i, line in enumerate(tmp):
            if "Ecology" in line:
                tmp[i] = "Ecology: {}\n".format(", ".join(ecol_bores))
            # End if

            if "Policy" in line:
                tmp[i] = "Policy: {}\n".format(", ".join(policy_found))
        # End for
    # End with

    # Write out ecology bore hydrocodes
    with open(os.path.join(data_folder, "model_linking.csv"), 'w') as model_link:
        model_link.writelines(tmp)
    # End with

    # Visuals checks on getting nearest mapped bore from top layer for the ecology part:
    if VERBOSE:
        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.scatter(obs_filter_bores['Easting'], obs_filter_bores['Northing'], label='bores')
        ax.scatter(stream_active['Easting'], stream_active[
                   'Northing'], color='r', label='stream gauges')
        for idx in stream_active.index:
            # Annotate the closest gauge
            x_y = stream_active.loc[stream_active.index == idx, ["Easting", "Northing"]]
            ax.annotate(idx, xy=x_y.values[0].tolist(), xycoords='data', xytext=(1.0, 0.8), textcoords='offset points')

        ax.scatter(ecol_bores_df['Easting'], ecol_bores_df['Northing'],
                   marker='+', color='orange', label='closest bores')
        for idx in ecol_bores_df.index:
            # Annotate closest bore
            x_y = ecol_bores_df.loc[ecol_bores_df.index == idx, ["Easting", "Northing"]]
            ax.annotate(idx, xy=x_y.values[0].tolist(), xycoords='data', xytext=(-120.0, 50.0),
                        textcoords='offset points',
                        arrowprops=dict(facecolor='orange', shrink=0.05))

        # plt.legend()
        # fig.suptitle('Finding nearest bores to stream gauges')
        # plt.xlabel('Easting')
        # plt.ylabel('Northing')
        # plt.axis('equal')
        #
        # plt.show()

        # set_ylabel('Northing')

        # fig = plt.figure()
        # ax = fig.add_subplot(111)
        # c = 'rgb'
        # for index, bore in enumerate(ecol_bores):
        #     bore_data_levels[bore_data_levels['HydroCode'] == bore][
        #         ['bore_date', 'result']].plot(x='bore_date', ax=ax, color=c[index])
        # # NOTE that the bores data is not going all the way to 2015, although bore filtering could include only those bores
        # # which have data that is recent .... can do this later!
        # # It is interesting to note that the distance can be quite far from gauge to bore
        # # Perhaps the restriction to top layer bores could be relaxed somewhat.
        #
        # plt.show()

    if VERBOSE:
        print("************************************************************************")
        print(" Mapping farm areas to grid")

    tr_model.map_polygon_to_grid(farms_poly, feature_name="ZoneID")

    if VERBOSE:
        print("************************************************************************")
        print(" Package up groundwater model builder object")

    tr_model.package_model()

    print(("Packaged into {}".format(tr_model.out_data_folder_grid)))
# End main()

if __name__ == '__main__':
    main()
