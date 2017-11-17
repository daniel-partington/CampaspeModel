from osgeo import osr
import os
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt


from HydroModelBuilder.GWModelBuilder import GWModelBuilder
from HydroModelBuilder.GISInterface.GDALInterface.GDALInterface import GDALInterface

from CampaspeModel.CustomScripts import Campaspe_data
from CampaspeModel.build_common import Campaspe_mesh
from CampaspeModel.build_common import rivers
from CampaspeModel.build_common import groundwater_boundary

# Define basic model parameters:
Proj_CS = osr.SpatialReference()
Proj_CS.ImportFromEPSG(28355) # 28355 is the code for gda94 mga zone 55; http://spatialreference.org/ref/epsg/gda94-mga-zone-55/

Interface = GDALInterface()
Interface.projected_coordinate_system = Proj_CS 
Interface.pcs_EPSG = "EPSG:28355"

SS_model = GWModelBuilder(name="01_steady_state", 
                          data_folder=r"C:\Workspace\part0075\MDB modelling\testbox\input_data\\",
                          model_data_folder=r"C:\Workspace\part0075\MDB modelling\testbox\00_Campaspe_Cascade\01_steady_state\\",
                          out_data_folder=r"C:\Workspace\part0075\MDB modelling\testbox\data_build\\",
                          GISInterface=Interface,
                          model_type='Modflow',
                          mesh_type='structured')

# Cleanup
#SS_model.flush()

Campaspe_data_folder = r"C:\Workspace\part0075\MDB modelling\Campaspe_data"
hu_raster_path = r"C:\Workspace\part0075\MDB modelling\ESRI_GRID_raw\ESRI_GRID"

custom_data = \
    Campaspe_data.process_custom_scripts_and_spatial_data(SS_model, 
                                                          Campaspe_data_folder,
                                                          verbose=True)
HGU_props = custom_data['HGU_props']
rain_gauges = custom_data['rain_gauges']
long_term_historic_rainfall = custom_data['long_term_historic_rainfall'] 
recharge_zones = custom_data['recharge_zones']
surface_raster_high_res = custom_data['surface_raster_high_res'] 
surface_raster_high_res_GSA = custom_data['surface_raster_high_res_GSA'] 
river_gauges = custom_data['river_gauges']
Campaspe_river_poly_file = custom_data['Campaspe_river_poly_file']
Murray_river_poly_file = custom_data['Murray_river_poly_file']
river_flow_data = custom_data['river_flow_data']
river_stage_data = custom_data['river_stage_data']
Campaspe = custom_data['Campaspe']
Campaspe_relevant = custom_data['Campaspe_relevant']

#******************************************************************************
#******************************************************************************
#******************************************************************************
#******************************************************************************

print '########################################################################'
print '########################################################################'
print '## Mesh specific model building '
print '########################################################################'
print '########################################################################'

HGU, hu_raster_files_reproj = Campaspe_mesh.build_mesh_and_set_properties(SS_model,
                                                  hu_raster_path,
                                                  HGU_props,
                                                  resolution=1000,
                                                  create_basement=True)

SS_model.map_rasters_to_grid(os.path.basename(surface_raster_high_res), os.path.dirname(surface_raster_high_res))
surface_raster_high_res = os.path.join(SS_model.out_data_folder, os.path.basename(surface_raster_high_res) + '_clipped.tif')

print "************************************************************************"
print " Interpolating rainfall data to grid "

interp_rain = SS_model.interpolate_points2mesh(rain_gauges, long_term_historic_rainfall, feature_id='Name', method='linear')
# Adjust rainfall to m from mm and from year to day
interp_rain = interp_rain / 1000.0 / 365.0

SS_model.boundaries.create_model_boundary_condition('Rainfall', 'rainfall', bc_static=True)
SS_model.boundaries.assign_boundary_array('Rainfall', interp_rain)

# Replace interp_rain with a copy to prevent alteration of assigned boundary array
interp_rain = np.copy(interp_rain)

recharge_zone_array = SS_model.map_raster_to_regular_grid_return_array(recharge_zones)

rch_zone_dict = {i:x for i, x in enumerate(np.unique(recharge_zone_array))}
rch_zones = len(rch_zone_dict.keys())

SS_model.parameters.create_model_parameter_set('ssrch', 
                                               value=0.01,
                                               num_parameters=rch_zones)
SS_model.parameters.parameter_options_set('ssrch', 
                                      PARTRANS='fixed', 
                                      PARCHGLIM='factor', 
                                      PARLBND=1.0E-3, 
                                      PARUBND=0.5, 
                                      PARGP='ssrch', 
                                      SCALE=1, 
                                      OFFSET=0)

for i in range(rch_zones - 1):
    interp_rain[recharge_zone_array == rch_zone_dict[i + 1]] = \
        interp_rain[recharge_zone_array == rch_zone_dict[i + 1]] * \
        SS_model.parameters.param['ssrch{}'.format(i)]['PARVAL1']

interp_rain[recharge_zone_array==rch_zone_dict[0]] = interp_rain[recharge_zone_array == rch_zone_dict[0]] * 0.0
interp_rain[SS_model.model_mesh3D[1][0] == -1] = 0.
    
rch = {}
rch[0] = interp_rain

print "************************************************************************"
print " Creating recharge boundary "

SS_model.boundaries.create_model_boundary_condition('Rain_reduced', 'recharge', bc_static=True)
SS_model.boundaries.associate_zonal_array_and_dict('Rain_reduced', recharge_zone_array, rch_zone_dict)
SS_model.boundaries.assign_boundary_array('Rain_reduced', rch)

# Initial heads using uniform head over entire model making whole domain saturated:
initial_heads_SS = np.full(SS_model.model_mesh3D[1].shape, 500.)

# Initial heads using ground surface elevation:
#for i in range(len(hu_raster_files_reproj)/2):
#    initial_heads_SS[i] = SS_model.model_mesh3D[1][0] #(interp_heads[hu_raster_files[2]]) # 2*i

SS_model.initial_conditions.set_as_initial_condition("Head", initial_heads_SS)

print "************************************************************************"
print " Mapping Campaspe river to grid"

num_reaches = 20    
river_seg, reach_df, reach_data, known_points = \
    rivers.prepare_river_data_for_Campaspe(SS_model, 
                                    surface_raster_high_res,
                                    river_gauges,
                                    Campaspe_river_poly_file,
                                    Campaspe,
                                    num_reaches=num_reaches)

segment_data, seg_dict = \
    rivers.create_segment_data(SS_model, river_seg, river_flow_data)

SS_model.save_MODFLOW_SFR_dataframes('Campaspe', reach_df, segment_data)
SS_model.river_mapping['Campaspe'] = river_seg

print "************************************************************************"
print " Creating Campaspe river boundary"

SS_model.boundaries.create_model_boundary_condition('Campaspe River', 'river_flow', bc_static=True)
SS_model.boundaries.assign_boundary_array('Campaspe River', [reach_data, seg_dict])


print "************************************************************************"
print " Mapping Murray River to grid"

riv, mriver_seg_ghb = \
    rivers.prepare_river_data_for_Murray(SS_model, surface_raster_high_res_GSA,
                                         r"C:\Workspace\part0075\MDB modelling\test_model.shp", #Murray_river_poly_file,
                                         Campaspe_relevant,
                                         river_stage_data,
                                         river_seg,
                                         plot=True) 

print "************************************************************************"
print " Creating Murray River boundary"

SS_model.boundaries.create_model_boundary_condition('Murray River', 'river', bc_static=True)
SS_model.boundaries.assign_boundary_array('Murray River', riv)

print "************************************************************************"
print " Setting up Murray River GHB boundary"
  
ghb = groundwater_boundary.prepare_ghb_boundary_from_Murray_data(SS_model,
                                                                 mriver_seg_ghb)   

print "************************************************************************"
print " Creating GHB boundary"

SS_model.boundaries.create_model_boundary_condition('GHB', 'general head', bc_static=True)
SS_model.boundaries.assign_boundary_array('GHB', ghb)

print "************************************************************************"
print " Creating parameters for transport "

# General parameters for transport
for unit in HGU:
    SS_model.parameters.create_model_parameter('por_{}'.format(unit), value=0.25)
    SS_model.parameters.parameter_options('por_{}'.format(unit), 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=0.05, 
                                          PARUBND=0.4, 
                                          PARGP='transport', 
                                          SCALE=1, 
                                          OFFSET=0)

SS_model.parameters.create_model_parameter('disp', value=0.01)
SS_model.parameters.parameter_options('disp', 
                                      PARTRANS='log', 
                                      PARCHGLIM='factor', 
                                      PARLBND=1E-4, 
                                      PARUBND=100., 
                                      PARGP='transport', 
                                      SCALE=1, 
                                      OFFSET=0)

print "************************************************************************"
print " Package up groundwater model builder object"

SS_model.package_model()
