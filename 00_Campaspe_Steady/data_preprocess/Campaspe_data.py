import os

import pandas as pd

from CampaspeModel.CustomScripts import processWeatherStations, getBoreData, get_GW_licence_info, processRiverStations, readHydrogeologicalProperties

# Define the units for the project for consistency and to allow converions on input data
# ModelBuilderObject.length = 'm'
# ModelBuilderObject.time = 'd'

def custom_scripts_and_spatial_data(ModelBuilderObject):
    # Set the model boundary using a polygon shapefile:
    print "************************************************************************"
    print " Setting model boundary "
    
    ModelBuilderObject.set_model_boundary_from_polygon_shapefile("GW_model_area.shp", 
                                                          shapefile_path=ModelBuilderObject.data_folder)
    
    # Set data boundary for model data
    print "************************************************************************"
    print " Setting spatial data boundary "
    ModelBuilderObject.set_data_boundary_from_polygon_shapefile(ModelBuilderObject.boundary_poly_file, 
                                                        buffer_dist=20000)
    
    # Setup recharge:
    # ... read in climate data using Custom_Scripts
    weather_stations = ['Kyneton',  'Elmore', 'Rochester', 'Echuca']
    print "************************************************************************"
    print " Executing custom script: processWeatherStations "
    
    rain_info_file = "rain_processed"
    # Check if this data has been processed and if not process it
    if os.path.exists(ModelBuilderObject.out_data_folder + rain_info_file + '.h5'):
        long_term_historic_rainfall = ModelBuilderObject.load_dataframe(ModelBuilderObject.out_data_folder + rain_info_file + '.h5')
    else:
        long_term_historic_rainfall = processWeatherStations.processWeatherStations(weather_stations, path=r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Climate\\")
        ModelBuilderObject.save_dataframe(ModelBuilderObject.out_data_folder + rain_info_file, long_term_historic_rainfall)
    
    rain_gauges = ModelBuilderObject.read_points_data(r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Climate\Rain_gauges.shp")
    
    # $%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%
    # INCLUDE NSW bores in this next part too for better head representation at the border, i.e. Murray River
    
    # Read in bore data:
    print "************************************************************************"
    print " Executing custom script: getBoreData "
    
    bore_levels_file = "bore_levels"
    bore_info_file = "bore_info"
    if os.path.exists(ModelBuilderObject.out_data_folder + bore_levels_file + ".h5") & os.path.exists(ModelBuilderObject.out_data_folder + bore_info_file + ".h5"):
        bore_data_levels = ModelBuilderObject.load_dataframe(ModelBuilderObject.out_data_folder + bore_levels_file + ".h5")
        bore_data_info = ModelBuilderObject.load_dataframe(ModelBuilderObject.out_data_folder + bore_info_file + ".h5")
    else:
        bore_data_levels, bore_data_info = getBoreData.getBoreData()
        ModelBuilderObject.save_dataframe(ModelBuilderObject.out_data_folder + bore_levels_file, bore_data_levels)
        ModelBuilderObject.save_dataframe(ModelBuilderObject.out_data_folder + bore_info_file, bore_data_info)
    # end if
    
    # getBoreDepth ... assuming that midpoint of screen interval is representative location and assign to layer accordingly
    bore_data_info['depth'] = (bore_data_info['TopElev'] + bore_data_info['BottomElev'])/2.0
    
    bore_data_info["HydroCode"] = bore_data_info.index
    
    # For steady state model, only use bore details containing average level, not 
    #observation_bores = ModelBuilderObject.read_points_data(r"C:\Workspace\part0075\MDB modelling\ngis_shp_VIC\ngis_shp_VIC\NGIS_Bores.shp")
    
    print "************************************************************************"
    print " Read in and filtering bore spatial data "
    
    bores_shpfile = ModelBuilderObject.read_points_data(r"C:\Workspace\part0075\MDB modelling\ngis_shp_VIC\ngis_shp_VIC\NGIS_Bores.shp")
    
    bores_filtered_from_shpfile = ModelBuilderObject.points_shapefile_obj2dataframe(bores_shpfile, feature_id="HydroCode")
    
    # Get the intersection of bores_filtered_from_shpfile with bores_data_info
    
    final_bores = pd.merge(bore_data_info, bores_filtered_from_shpfile, how='inner', on="HydroCode")
    
    # Only consider bores whereby the measured values are above the bottom of screen
    final_bores = final_bores[final_bores['mean level'] > final_bores['BottomElev']]
    
    print 'Final number of bores within the data boundary that have level data and screen info: ', final_bores.shape[0]
    
    #final_bores.plot(kind='scatter', x="Easting", y="Northing", c="mean level", cmap="Spectral") # , edgecolor='None'
    
    
    # Load in the pumping wells data
    filename = "Groundwater licence information for Dan Partington bc301115.xlsx"
    path = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\GW\Bore data\\"    
    out_path = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\GW\Bore data\\"
    out_file = "pumping wells.shp"
    
    print "************************************************************************"
    print " Executing custom script: get_GW_licence_info "
    
    pumping_data = get_GW_licence_info.get_GW_licence_info(filename, path=path, out_file=out_file, out_path=out_path)
    pumps_points = ModelBuilderObject.read_points_data(r"C:\Workspace\part0075\MDB modelling\Campaspe_data\GW\Bore data\pumping wells.shp")
    
    print "************************************************************************"
    print " Executing custom script: readHydrogeologicalProperties "
    
    file_location = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\GW\Aquifer properties\Hydrogeologic_variables.xlsx"
    HGU_props = readHydrogeologicalProperties.getHGUproperties(file_location)
    
    print "************************************************************************"
    print "Get the C14 data"
    
    C14_points = ModelBuilderObject.read_points_data(r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Chemistry\C14.shp")    
    
    #C14_wells_info_file = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Chemistry\C14_bore_depth.csv"
    #df_C14_info = pd.read_csv(C14_wells_info_file)    
    #df_C14_info = df_C14_info.dropna()
    #df_C14_info = df_C14_info.set_index('Bore_id')    
    
    C14data = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Chemistry\C14_locs.xlsx"
    df_C14 = pd.read_excel(C14data)
    df_C14.drop_duplicates(subset=["Bore_id"], inplace=True)
    df_C14.dropna(inplace=True)
    
    print "************************************************************************"
    print " Executing custom script: processRiverStations "
    
    river_flow_file = "river_flow_processed"
    # Check if this data has been processed and if not process it
    if os.path.exists(ModelBuilderObject.out_data_folder + river_flow_file + '.h5'):
        river_flow_data = ModelBuilderObject.load_dataframe(ModelBuilderObject.out_data_folder + river_flow_file + '.h5')
    else:
        river_flow_data = processRiverStations.getFlow(path=r"C:\Workspace\part0075\MDB modelling\Campaspe_data\SW\All_streamflow_Campaspe_catchment\\")
        ModelBuilderObject.save_dataframe(ModelBuilderObject.out_data_folder + river_flow_file, river_flow_data)
    
    river_stage_file = "river_stage_processed"
    # Check if this data has been processed and if not process it
    if os.path.exists(ModelBuilderObject.out_data_folder + river_stage_file + '.h5'):
        river_stage_data = ModelBuilderObject.load_dataframe(ModelBuilderObject.out_data_folder + river_stage_file + '.h5')
    else:
        river_stage_data = processRiverStations.getStage(path=r"C:\Workspace\part0075\MDB modelling\Campaspe_data\SW\All_streamflow_Campaspe_catchment\\")
        ModelBuilderObject.save_dataframe(ModelBuilderObject.out_data_folder + river_stage_file, river_stage_data)
    
    #river_gauges = ModelBuilderObject.read_points_data(r"C:\Workspace\part0075\MDB modelling\Campaspe_data\20_gauges\Site_info.shp")
    river_gauges = ModelBuilderObject.read_points_data(r"C:\Workspace\part0075\MDB modelling\Campaspe_data\SW\All_streamflow_Campaspe_catchment\processed_river_sites_stage.shp")
    
    print "************************************************************************"
    print "Load in the river shapefiles"
    
    Campaspe_river_poly = ModelBuilderObject.read_polyline("Campaspe_Riv.shp", path=r"C:\Workspace\part0075\MDB modelling\Campaspe_model\GIS\GIS_preprocessed\Surface_Water\Streams\\") 
    Murray_river_poly = ModelBuilderObject.read_polyline("River_Murray.shp", path=r"C:\Workspace\part0075\MDB modelling\Campaspe_model\GIS\GIS_preprocessed\Surface_Water\Streams\\") 
    
    print "************************************************************************"
    print "Load in the shapefiles defining groundwater boundaries"
    
    WGWbound_poly = ModelBuilderObject.read_polyline("western_head.shp", path=r"C:\Workspace\part0075\MDB modelling\testbox\input_data\\") 
    EGWbound_poly = ModelBuilderObject.read_polyline("eastern_head.shp", path=r"C:\Workspace\part0075\MDB modelling\testbox\input_data\\") 

