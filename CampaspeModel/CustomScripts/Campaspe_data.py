import os
import pandas as pd
import numpy as np

from CampaspeModel.CustomScripts import processWeatherStations, getBoreData, get_GW_licence_info, processRiverStations, readHydrogeologicalProperties

# Define the units for the project for consistency and to allow converions on input data
# ModelBuilderObject.length = 'm'
# ModelBuilderObject.time = 'd'


def process_custom_scripts_and_spatial_data(ModelBuilderObject, 
                                            Campaspe_data_folder,
                                            verbose=True):
    
    custom_data = {}
    MBO = ModelBuilderObject
    # Set the model boundary using a polygon shapefile:
    if verbose:
        print "************************************************************************"
        print " Setting model boundary "
    
    MBO.set_model_boundary_from_polygon_shapefile("GW_model_area.shp", 
                                                       shapefile_path=MBO.data_folder)
    
    # Set data boundary for model data
    if verbose:
        print "************************************************************************"
        print " Setting spatial data boundary "

    MBO.set_data_boundary_from_polygon_shapefile(MBO.boundary_poly_file, 
                                                      shapefile_path=MBO.out_data_folder,
                                                      buffer_dist=20000)
    
    # Setup recharge:
    # ... read in climate data using Custom_Scripts
    weather_stations = ['Kyneton', 'Eppalock', 'Elmore', 'Rochester', 'Echuca']

    if verbose:
        print "************************************************************************"
        print " Executing custom script: processWeatherStations "

    custom_data['rain_gauges'] = MBO.read_points_data(os.path.join(Campaspe_data_folder, 
                                                    r"Climate\Rain_gauges.shp"))
    
    rain_info_file = "rain_processed"
    # Check if this data has been processed and if not process it
    if os.path.exists(MBO.out_data_folder + rain_info_file + '.h5'):
        custom_data['long_term_historic_rainfall'] = \
            MBO.load_dataframe(MBO.out_data_folder + rain_info_file + '.h5')
    else:
        custom_data['long_term_historic_rainfall'] = \
            processWeatherStations.processWeatherStations(weather_stations, 
                                                          path=os.path.join(Campaspe_data_folder,r"Climate\\"))
        MBO.save_dataframe(MBO.out_data_folder + rain_info_file, 
                           custom_data['long_term_historic_rainfall'])

    rain_info_file2 = "rain_processed_transient"
    if os.path.exists(MBO.out_data_folder + rain_info_file2 + '.h5'):
        custom_data['long_term_historic_weather'] = MBO.load_dataframe(MBO.out_data_folder + rain_info_file2 + '.h5')
    else:
        custom_data['long_term_historic_weather'] = \
            processWeatherStations.processWeatherStations(weather_stations, 
                                                          path=os.path.join(Campaspe_data_folder,r"Climate\\"), 
                                                          frequency='M')
        MBO.save_dataframe(MBO.out_data_folder + rain_info_file2, custom_data['long_term_historic_weather'])
        
    # Read in bore data:
    if verbose:
        print "************************************************************************"
        print " Executing custom script: getBoreData "
    
    bore_levels_file = "bore_levels"
    bore_salinity_file = "bore_salinity"
    bore_info_file = "bore_info"
    if os.path.exists(MBO.out_data_folder + bore_levels_file + ".h5") & \
       os.path.exists(MBO.out_data_folder + bore_info_file + ".h5") & \
       os.path.exists(MBO.out_data_folder + bore_salinity_file + ".h5"):
        custom_data['bore_data_levels'] = MBO.load_dataframe(MBO.out_data_folder + bore_levels_file + ".h5")
        custom_data['bore_data_info'] = MBO.load_dataframe(MBO.out_data_folder + bore_info_file + ".h5")
        custom_data['bore_data_salinity'] = MBO.load_dataframe(MBO.out_data_folder + bore_salinity_file + ".h5")
    else:
        custom_data['bore_data_levels'], custom_data['bore_data_info'], custom_data['bore_data_salinity'] = \
            getBoreData.getBoreData(path=os.path.join(Campaspe_data_folder, "ngis_shp_VIC"))
        MBO.save_dataframe(MBO.out_data_folder + bore_levels_file, custom_data['bore_data_levels'])
        MBO.save_dataframe(MBO.out_data_folder + bore_info_file, custom_data['bore_data_info'])
        MBO.save_dataframe(MBO.out_data_folder + bore_salinity_file, custom_data['bore_data_salinity'])
    # end if
    
    # getBoreDepth ... assuming that midpoint of screen interval is representative location and assign to layer accordingly
    custom_data['bore_data_info']['depth'] = \
        (custom_data['bore_data_info']['TopElev'] + \
         custom_data['bore_data_info']['BottomElev']) / 2.0
    
    custom_data['bore_data_info']["HydroCode"] = custom_data['bore_data_info'].index
    
    # For steady state model, only use bore details containing average level, not 
    #observation_bores = MBO.read_points_data(r"C:\Workspace\part0075\MDB modelling\ngis_shp_VIC\ngis_shp_VIC\NGIS_Bores.shp")
    
    if verbose:
        print "************************************************************************"
        print " Read in and filtering bore spatial data "
    
    bores_shpfile = MBO.read_points_data(os.path.join(Campaspe_data_folder,r"ngis_shp_VIC\ngis_shp_VIC\NGIS_Bores.shp"))
    custom_data['bores_shpfile'] = bores_shpfile
    bores_filtered_from_shpfile = MBO.points_shapefile_obj2dataframe(bores_shpfile, feature_id="HydroCode")
    
    # Get the intersection of bores_filtered_from_shpfile with bores_data_info
    
    final_bores = pd.merge(custom_data['bore_data_info'], bores_filtered_from_shpfile, how='inner', on="HydroCode")
    
    # Only consider bores whereby the measured values are above the bottom of screen
    final_bores = final_bores[final_bores['mean level'] > final_bores['BottomElev']]
    
    custom_data['final_bores'] = final_bores
    print 'Final number of bores within the data boundary that have level data and screen info: ', final_bores.shape[0]
    
    #final_bores.plot(kind='scatter', x="Easting", y="Northing", c="mean level", cmap="Spectral") # , edgecolor='None'
    
    
    # Load in the pumping wells data
    filename = "Groundwater licence information for Dan Partington bc301115.xlsx"
    path = os.path.join(Campaspe_data_folder, r"GW\Bore data\\")    
    out_path = os.path.join(Campaspe_data_folder, r"GW\Bore data\\")
    out_file = "pumping wells.shp"
    
    if verbose:
        print "************************************************************************"
        print " Executing custom script: get_GW_licence_info "
    
    custom_data['pumping_data'] = get_GW_licence_info.get_GW_licence_info(filename, path=path, out_file=out_file, out_path=out_path)
    custom_data['pumps_points'] = MBO.read_points_data(os.path.join(Campaspe_data_folder, r"GW\Bore data\pumping wells.shp"))        
    
    if verbose:
        print "************************************************************************"
        print "Get the C14 data"
    
    C14_points = MBO.read_points_data(os.path.join(Campaspe_data_folder, r"Chemistry\C14.shp"))    
    
    #C14_wells_info_file = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Chemistry\C14_bore_depth.csv"
    #df_C14_info = pd.read_csv(C14_wells_info_file)    
    #df_C14_info = df_C14_info.dropna()
    #df_C14_info = df_C14_info.set_index('Bore_id')    
    
    C14data = os.path.join(Campaspe_data_folder, r"Chemistry\C14_locs.xlsx")
    df_C14 = pd.read_excel(C14data)
    df_C14.drop_duplicates(subset=["Bore_id"], inplace=True)
    df_C14.dropna(inplace=True)
    
    C14_half_life = 5730
    C14_lambda = np.log(2) / 5730.
    C13_carbonate_signature = -2 # Could be from 0 to -2
    C13_recharge_signature = -18.5 # Cartright 2010, -15 to -25
    
    df_C14.loc[:, 'a14C_corrected'] = df_C14['a14C(pMC)'] / ( 
                               (df_C14['d13C'] - C13_carbonate_signature) / 
                               (C13_recharge_signature - C13_carbonate_signature))
    df_C14.loc[df_C14['a14C_corrected'] > 100., 'a14C_corrected'] = 100.
    
    custom_data['df_C14'] = df_C14    
    custom_data['C14_points'] = C14_points

    if verbose:
        print "************************************************************************"
        print " Executing custom script: readHydrogeologicalProperties "
    
    HGU_file_location = os.path.join(Campaspe_data_folder, r"GW\Aquifer properties\Hydrogeologic_variables.xlsx")
    custom_data['HGU_props'] = readHydrogeologicalProperties.getHGUproperties(HGU_file_location)
    
    if verbose:
        print "************************************************************************"
        print " Executing custom script: processRiverStations "
    
    river_flow_file = "river_flow_processed"
    river_stage_file = "river_stage_processed"
    river_ec_file = "river_ec_processed"
    river_data_folder = os.path.join(Campaspe_data_folder, r"SW\All_streamflow_Campaspe_catchment\Updated")
    river_flow_file = os.path.join(MBO.out_data_folder, river_flow_file)
    river_stage_file = os.path.join(MBO.out_data_folder, river_stage_file)
    river_ec_file = os.path.join(MBO.out_data_folder, river_ec_file)
    
    site_details_file = "Site Details.csv"
    site_details = pd.read_csv(os.path.join(river_data_folder, site_details_file))
    # As all of the stream data for the whole of the Camaspe catchment is in the folder
    # to be processed, we can prefilter sites to examine by specifying sites.
    custom_data['Campaspe_relevant'] = \
        site_details[
          site_details['Site Name'].str.contains("CAMPASPE RIVER") | \
          site_details['Site Name'].str.contains("MURRAY RIVER") | \
          site_details['Site Name'].str.contains("AXE CREEK") | \
          site_details['Site Name'].str.contains("MOUNT PLEASANT")]
    
    Campaspe = site_details[site_details['Site Name'].str.contains("CAMPASPE RIVER")]
    custom_data['Campaspe'] = Campaspe[Campaspe['Northing'] >= \
                        Campaspe.loc[6]['Northing']]
    
    from geopandas import GeoDataFrame
    from shapely.geometry import Point
    
    geometry = [Point(xy) for xy in zip(Campaspe.Easting, Campaspe.Northing)]
    Campaspe_geo = pd.DataFrame.copy(Campaspe)
    Campaspe_geo.drop(['Northing', 'Easting'], axis=1)
    crs = {'init': MBO.GISInterface.pcs_EPSG}
    Geo_df = GeoDataFrame(Campaspe_geo, crs=crs, geometry=geometry)
    site_shapefile = os.path.join(river_data_folder, 'processed_river_sites_stage.shp') 
    Geo_df.to_file(site_shapefile)
    sites = custom_data['Campaspe_relevant']['Site Id'].tolist()
    
    custom_data['Campaspe_gauge_zero'] = Campaspe[(Campaspe['Gauge Zero (Ahd)'] > 0.0) | \
                                   (Campaspe['Cease to flow level'] > 10.0) | \
                                   (Campaspe['Min value'] > 10.0)]
    
    # Check if this data has been processed and if not process it
    if os.path.exists(river_flow_file + '.pkl'):
        custom_data['river_flow_data'] = MBO.load_obj(river_flow_file + '.pkl')
    else:
        custom_data['river_flow_data'] = processRiverStations.getFlow(path=river_data_folder, sites=sites)
        MBO.save_obj(custom_data['river_flow_data'], river_flow_file)
    
    # Check if this data has been processed and if not process it
    if os.path.exists(river_stage_file + '.pkl'):
        custom_data['river_stage_data'] = MBO.load_obj(river_stage_file + '.pkl')
    else:
        custom_data['river_stage_data'] = processRiverStations.getStage(path=river_data_folder, sites=sites)
        MBO.save_obj(custom_data['river_stage_data'], river_stage_file)
    
    # Check if this data has been processed and if not process it
    if os.path.exists(river_ec_file + '.pkl'):
        custom_data['river_ec_data'] = MBO.load_obj(river_ec_file + '.pkl')
    else:
        custom_data['river_ec_data'] = processRiverStations.getEC(path=river_data_folder, sites=sites)
        MBO.save_obj(custom_data['river_ec_data'], river_ec_file)
    
    custom_data['river_gauges'] = MBO.read_points_data(site_shapefile)

    if verbose:
        print "************************************************************************"
        print "Load in the Campaspe river field sampling data"
    
    field_data_folder = os.path.join(Campaspe_data_folder, r"Chemistry\Radon_interpreter")
    field_data_file = "data_Campaspe.csv"
    custom_data['FieldData'] = pd.read_csv(os.path.join(field_data_folder, field_data_file), skiprows=[1], index_col='Date', parse_dates=True, dayfirst=True, usecols=range(0,15))
    
    if verbose:
        print "************************************************************************"
        print "Load in the river shapefiles"
    
    custom_data['Campaspe_river_poly_file'] = os.path.join(MBO.data_folder, r"Waterways\Campaspe_Riv.shp")
    custom_data['Campaspe_river_poly'] = MBO.read_poly("Campaspe_Riv.shp", path=os.path.join(MBO.data_folder, r"Waterways")) 
    custom_data['Murray_river_poly_file'] = os.path.join(MBO.data_folder, r"Waterways\River_Murray.shp") 
    custom_data['Murray_river_poly'] = MBO.read_poly("River_Murray.shp", path=os.path.join(MBO.data_folder, r"Waterways")) 
    
    # This is disgusting, but works for now ... needs cleaning up through first testing
    # if raster is in the right projection and if not returning the name of the new
    # reprojected file
    custom_data['surface_raster_high_res_GSA'] = os.path.join(Campaspe_data_folder, r"Surface_DEM_Geoscience_Australia\Camp_1sHE_2026754\Camp_1sHE_reproj.tif")
    custom_data['surface_raster_high_res'] = r"C:\Workspace\part0075\MDB modelling\ESRI_GRID_raw\ESRI_GRID\sur_1t"
    
    if verbose:
        print "************************************************************************"
        print "Load in the recharge zones as determined by Yueqing Xie"
    
    custom_data['recharge_zones'] = os.path.join(Campaspe_data_folder, r"Linking_recharge\Zones_24.tif")
    
    if verbose:
        print "************************************************************************"
        print "Load in the shapefiles defining groundwater boundaries"
    
    custom_data['WGWbound_poly'] = MBO.read_poly("western_head.shp", path=MBO.data_folder)
    custom_data['EGWbound_poly'] = MBO.read_poly("eastern_head.shp", path=MBO.data_folder)


    return custom_data
    
if __name__ == '__main__':
    pass