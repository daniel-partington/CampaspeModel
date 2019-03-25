import os

import logging
import numpy as np
import pandas as pd

from CampaspeModel.build_common import (get_bore_data, get_gw_licence_info,
                                        process_river_diversions,
                                        process_river_stations,
                                        process_weather_stations,
                                        read_hydrogeological_properties)


# Define the units for the project for consistency and to allow converions on input data
# ModelBuilderObject.length = 'm'
# ModelBuilderObject.time = 'd'


def process_custom_scripts_and_spatial_data(model_builder_object,
                                            campaspe_data_folder,
                                            verbose=True,
                                            model_boundary_file="GW_model_area2.shp",
                                            GW_link_Integrated=False):

    custom_data = {}
    mbo = model_builder_object
    # Set the model boundary using a polygon shapefile:

    if verbose:
        print "************************************************************************"
        print " Setting model boundary "

    mbo.set_model_boundary_from_polygon_shapefile(model_boundary_file,
                                                  shapefile_path=mbo.data_folder)

    # Set data boundary for model data
    if verbose:
        print "************************************************************************"
        print " Setting spatial data boundary "

    mbo.set_data_boundary_from_polygon_shapefile(mbo.boundary_poly_file,
                                                 shapefile_path=mbo.out_data_folder,
                                                 buffer_dist=20000)

    # Setup recharge:
    # ... read in climate data using Custom_Scripts
    weather_stations = ['Kyneton', 'Eppalock', 'Elmore', 'Rochester', 'Echuca']

    if verbose:
        print "************************************************************************"
        print " Executing custom script: processWeatherStations "

    custom_data['rain_gauges'] = mbo.read_points_data(
        os.path.join(campaspe_data_folder, "Climate", "Rain_gauges.shp"))

    rain_info_file = "rain_processed"
    # Check if this data has been processed and if not process it
    if os.path.exists(os.path.join(mbo.out_data_folder, rain_info_file + '.h5')):
        custom_data['long_term_historic_rainfall'] = \
            mbo.load_dataframe(os.path.join(mbo.out_data_folder, rain_info_file + '.h5'))
    else:
        custom_data['long_term_historic_rainfall'] = \
            process_weather_stations.process_weather_stations(weather_stations,
                                                              path=os.path.join(campaspe_data_folder, "Climate") + os.path.sep)
        mbo.save_dataframe(os.path.join(mbo.out_data_folder, rain_info_file),
                           custom_data['long_term_historic_rainfall'])

    rain_info_file2 = "rain_processed_transient"
    if os.path.exists(os.path.join(mbo.out_data_folder, rain_info_file2 + '.h5')):
        custom_data['long_term_historic_weather'] = mbo.load_dataframe(
            os.path.join(mbo.out_data_folder, rain_info_file2 + '.h5'))
    else:
        custom_data['long_term_historic_weather'] = \
            process_weather_stations.process_weather_stations(weather_stations,
                                                              path=os.path.join(campaspe_data_folder,
                                                                                "Climate") + os.path.sep,
                                                              frequency='M')
        mbo.save_dataframe(os.path.join(mbo.out_data_folder, rain_info_file2),
                           custom_data['long_term_historic_weather'])

    # Read in bore data:
    if verbose:
        print "************************************************************************"
        print " Executing custom script: getBoreData "

    bore_levels_file = "bore_levels"
    bore_salinity_file = "bore_salinity"
    bore_info_file = "bore_info"
    if os.path.exists(os.path.join(mbo.out_data_folder, bore_levels_file + ".h5")) & \
       os.path.exists(os.path.join(mbo.out_data_folder, bore_info_file + ".h5")) & \
       os.path.exists(os.path.join(mbo.out_data_folder, bore_salinity_file + ".h5")):
        custom_data['bore_data_levels'] = mbo.load_dataframe(
            os.path.join(mbo.out_data_folder, bore_levels_file + ".h5"))
        custom_data['bore_data_info'] = mbo.load_dataframe(os.path.join(mbo.out_data_folder, bore_info_file + ".h5"))
        custom_data['bore_data_salinity'] = mbo.load_dataframe(
            os.path.join(mbo.out_data_folder, bore_salinity_file + ".h5"))
    else:
        if not GW_link_Integrated:
            custom_data['bore_data_levels'], custom_data['bore_data_info'], custom_data['bore_data_salinity'] = \
                get_bore_data.get_bore_data(path=os.path.join(campaspe_data_folder, "ngis_shp_VIC"))
        else:
            custom_data['bore_data_levels'], custom_data['bore_data_info'], custom_data['bore_data_salinity'] = \
                get_bore_data.get_bore_data(path=os.path.join(
                    campaspe_data_folder, "ngis_shp_VIC"), construct_record_number_max=1000)

        mbo.save_dataframe(os.path.join(mbo.out_data_folder, bore_levels_file), custom_data['bore_data_levels'])
        mbo.save_dataframe(os.path.join(mbo.out_data_folder, bore_info_file), custom_data['bore_data_info'])
        mbo.save_dataframe(os.path.join(mbo.out_data_folder, bore_salinity_file), custom_data['bore_data_salinity'])
    # end if

    # getBoreDepth ... assuming that midpoint of screen interval is representative location and assign to layer accordingly
    custom_data['bore_data_info']['depth'] = \
        (custom_data['bore_data_info']['TopElev'] +
         custom_data['bore_data_info']['BottomElev']) / 2.0

    custom_data['bore_data_info']["HydroCode"] = custom_data['bore_data_info'].index

    if verbose:
        print "************************************************************************"
        print " Read in and filtering bore spatial data "

    bores_shpfile = mbo.read_points_data(os.path.join(
        campaspe_data_folder, "ngis_shp_VIC", "ngis_shp_VIC", "NGIS_Bores.shp"))
    custom_data['bores_shpfile'] = bores_shpfile
    bores_filtered_from_shpfile = mbo.points_shapefile_obj2dataframe(bores_shpfile, feature_id="HydroCode")

    # Get the intersection of bores_filtered_from_shpfile with bores_data_info

    final_bores = pd.merge(custom_data['bore_data_info'], bores_filtered_from_shpfile, how='inner', on="HydroCode")

    # Only consider bores whereby the measured values are above the bottom of screen
    final_bores = final_bores[final_bores['mean level'] > final_bores['BottomElev']]

    custom_data['final_bores'] = final_bores
    print 'Final number of bores within the data boundary that have level data and screen info: ', final_bores.shape[0]

    # final_bores.plot(kind='scatter', x="Easting", y="Northing", c="mean level", cmap="Spectral") # , edgecolor='None'

    # Load in the pumping wells data
    filename = "Groundwater licence information for Dan Partington bc301115.xlsx"
    # TODO: Check if these separators are necessary
    path = os.path.join(campaspe_data_folder, "GW", "Bore data") + os.path.sep
    out_path = os.path.join(campaspe_data_folder, "GW", "Bore data") + os.path.sep
    out_file = "pumping wells.shp"

    if verbose:
        print "************************************************************************"
        print " Executing custom script: get_gw_licence_info "

    custom_data['pumping_data'] = get_gw_licence_info.get_gw_licence_info(
        filename, path=path, out_file=out_file, out_path=out_path)
    custom_data['pumps_points'] = mbo.read_points_data(os.path.join(
        campaspe_data_folder, "GW", "Bore data", "pumping wells.shp"))

    if verbose:
        print "************************************************************************"
        print "Get the C14 data"

    c14_points = mbo.read_points_data(os.path.join(campaspe_data_folder, "Chemistry", "C14.shp"))
    c14data = os.path.join(campaspe_data_folder, "Chemistry", "C14_locs.xlsx")
    df_c14 = pd.read_excel(c14data)
    df_c14.drop_duplicates(subset=["Bore_id"], inplace=True)
    df_c14.dropna(inplace=True)
    c14_half_life = 5730.
    c14_lambda = np.log(2) / c14_half_life
    c13_carbonate_signature = -2  # Could be from 0 to -2
    c13_recharge_signature = -18.5  # Cartright 2010, -15 to -25
    df_c14.loc[:, 'a14C_corrected'] = df_c14['a14C(pMC)'] / (
        (df_c14['d13C'] - c13_carbonate_signature) /
        (c13_recharge_signature - c13_carbonate_signature))
    df_c14.loc[df_c14['a14C_corrected'] > 100., 'a14C_corrected'] = 100.

    custom_data['df_C14'] = df_c14
    custom_data['C14_points'] = c14_points

    if verbose:
        print "************************************************************************"
        print " Executing custom script: readHydrogeologicalProperties "

    hgu_file_location = os.path.join(campaspe_data_folder, "GW", "Aquifer properties", "Hydrogeologic_variables.xlsx")
    custom_data['HGU_props'] = read_hydrogeological_properties.get_hgu_properties(hgu_file_location)

    if verbose:
        print "************************************************************************"
        print " Executing custom script: processRiverStations "

    river_flow_file = "river_flow_processed"
    river_stage_file = "river_stage_processed"
    river_ec_file = "river_ec_processed"
    river_data_folder = os.path.join(campaspe_data_folder, "SW",
                                     "All_streamflow_Campaspe_catchment", "Updated", "June2017", "MOST RECENT")
    river_flow_file = os.path.join(mbo.out_data_folder, river_flow_file)
    river_stage_file = os.path.join(mbo.out_data_folder, river_stage_file)
    river_ec_file = os.path.join(mbo.out_data_folder, river_ec_file)

    site_details_file = "Site Details.csv"
    site_details = pd.read_csv(os.path.join(river_data_folder, site_details_file))
    # As all of the stream data for the whole of the Camaspe catchment is in the folder
    # to be processed, we can prefilter sites to examine by specifying sites.
    river_file = os.path.join(campaspe_data_folder, "Field_sampling", "BedElevation_GPS.csv")
    custom_data['Campaspe_field_elevations'] = \
        pd.read_csv(river_file, usecols=['Northing',
                                         'Easting',
                                         'Elevation',
                                         'Site',
                                         'Chainage (m)'])\
        .groupby('Site').mean().sort_values('Chainage (m)')

    custom_data['Campaspe_relevant'] = \
        site_details[
        site_details['Site Name'].str.contains("CAMPASPE RIVER") |
        site_details['Site Name'].str.contains("MURRAY RIVER") |
        site_details['Site Name'].str.contains("AXE CREEK") |
        site_details['Site Name'].str.contains("MOUNT PLEASANT")]

    campaspe = site_details[site_details['Site Name'].str.contains("CAMPASPE RIVER")]
    custom_data['Campaspe'] = campaspe[campaspe['Northing'] >=
                                       campaspe.loc[6]['Northing']]

    from geopandas import GeoDataFrame
    from shapely.geometry import Point

    geometry = [Point(xy) for xy in zip(campaspe.Easting, campaspe.Northing)]
    campaspe_geo = pd.DataFrame.copy(campaspe)
    campaspe_geo.drop(['Northing', 'Easting'], axis=1)
    crs = {'init': mbo.GISInterface.pcs_EPSG}
    Geo_df = GeoDataFrame(campaspe_geo, crs=crs, geometry=geometry)
    site_shapefile = os.path.join(river_data_folder, 'processed_river_sites_stage.shp')
    Geo_df.to_file(site_shapefile)
    sites = custom_data['Campaspe_relevant']['Site Id'].tolist()

    custom_data['Campaspe_gauge_zero'] = campaspe[(campaspe['Gauge Zero (Ahd)'] > 0.0) |
                                                  (campaspe['Cease to flow level'] > 10.0) |
                                                  (campaspe['Min value'] > 10.0)]

    # Check if this data has been processed and if not process it
    if os.path.exists(river_flow_file + '.pkl'):
        custom_data['river_flow_data'] = mbo.load_obj(river_flow_file + '.pkl')
    else:
        custom_data['river_flow_data'] = process_river_stations.get_flow(path=river_data_folder, sites=sites)
        mbo.save_obj(custom_data['river_flow_data'], river_flow_file)

    # Check if this data has been processed and if not process it
    if os.path.exists(river_stage_file + '.pkl'):
        custom_data['river_stage_data'] = mbo.load_obj(river_stage_file + '.pkl')
    else:
        custom_data['river_stage_data'] = process_river_stations.get_stage(path=river_data_folder, sites=sites)
        mbo.save_obj(custom_data['river_stage_data'], river_stage_file)

    # Check if this data has been processed and if not process it
    if os.path.exists(river_ec_file + '.pkl'):
        custom_data['river_ec_data'] = mbo.load_obj(river_ec_file + '.pkl')
    else:
        custom_data['river_ec_data'] = process_river_stations.get_EC(path=river_data_folder, sites=sites)
        mbo.save_obj(custom_data['river_ec_data'], river_ec_file)

    custom_data['river_gauges'] = mbo.read_points_data(site_shapefile)

    river_diversions_file = "river_diversions_processed"
    river_diversions_file = os.path.join(mbo.out_data_folder, river_diversions_file)

    diversions_data_file = os.path.join(campaspe_data_folder,
                                        "SW", "Campaspe_System_Data.xlsx")
    if os.path.exists(river_diversions_file + '.pkl'):
        custom_data['river_diversions_data'] = mbo.load_obj(river_diversions_file + '.pkl')
    else:
        custom_data['river_diversions_data'] = \
            process_river_diversions.get_diversions(diversions_data_file)
        mbo.save_obj(custom_data['river_diversions_data'], river_diversions_file)

    if verbose:
        print "************************************************************************"
        print "Load in the Campaspe river field sampling data"

    field_data_folder = os.path.join(campaspe_data_folder, "Chemistry", "Radon_interpreter")
    field_data_file = "data_Campaspe.csv"
    custom_data['FieldData'] = pd.read_csv(os.path.join(field_data_folder, field_data_file), skiprows=[
                                           1], index_col='Date', parse_dates=True, dayfirst=True, usecols=range(0, 15))

    if verbose:
        print "************************************************************************"
        print "Load in the river shapefiles"

    waterway_dir = os.path.join(mbo.data_folder, "Waterways")
    custom_data['Campaspe_river_poly_file'] = os.path.join(
        waterway_dir, "Campaspe_Riv.shp")
    custom_data['Campaspe_river_poly'] = mbo.read_poly(
        "Campaspe_Riv.shp", path=waterway_dir)
    custom_data['Murray_river_poly_file'] = os.path.join(
        waterway_dir, "River_Murray.shp")
    custom_data['Murray_river_poly'] = mbo.read_poly(
        "River_Murray.shp", path=waterway_dir)

    # custom_data['Axe_creek_poly_file'] = os.path.join(waterway_dir, "Axe_creek.shp")
    # custom_data['Axe_creek_poly'] = mbo.read_poly(
    #     "Axe_creek.shp", path=waterway_dir)
    # custom_data['MtPleasant_creek_poly_file'] = os.path.join(waterway_dir, "MtPleasant_creek.shp")
    # custom_data['MtPleasant_creek_poly'] = mbo.read_poly(
    #     "MtPleasant_creek.shp", path=waterway_dir)

    # This is disgusting, but works for now ... needs cleaning up through first testing
    # if raster is in the right projection and if not returning the name of the new
    # reprojected file

    custom_data['surface_raster_high_res_GSA'] = os.path.join(
        campaspe_data_folder, os.path.join("Surface_DEM_Geoscience_Australia", "Camp_1sHE_2026754", "Camp_1sHE_reproj.tif"))
    custom_data['surface_raster_high_res'] = os.path.join(campaspe_data_folder, "ESRI_GRID_raw", "ESRI_GRID", "sur_1t")

    if verbose:
        print "************************************************************************"
        print "Load in the recharge zones as determined by Yueqing Xie"

    custom_data['recharge_zones'] = os.path.join(campaspe_data_folder, "Linking_recharge", "Zones_24.tif")
    custom_data['recharge_zone_info'] = pd.read_csv(
        os.path.join(campaspe_data_folder, "Linking_recharge", "Curve_fitting_norm.dat"),
        delim_whitespace=True)
    custom_data['recharge_zone_info_detailed'] = pd.read_excel(os.path.join(
        campaspe_data_folder, r"Linking_recharge",
        r"zone_percentage_recharge_statistics_method3.xlsx"), sheet='Sheet1')

    if verbose:
        print "************************************************************************"
        print "Load in the shapefiles defining groundwater boundaries"

    custom_data['WGWbound_poly'] = mbo.read_poly("western_head.shp", path=mbo.data_folder)
    custom_data['EGWbound_poly'] = mbo.read_poly("eastern_head.shp", path=mbo.data_folder)

    if verbose:
        print "************************************************************************"
        print "Load in the farms shapefile"

    farms_path = os.path.join(campaspe_data_folder, "SW", "Farm")
    custom_data['farms_poly'] = mbo.read_poly("farm_v1_prj.shp", path=farms_path, poly_type='polygon')

    return custom_data


if __name__ == '__main__':
    pass
