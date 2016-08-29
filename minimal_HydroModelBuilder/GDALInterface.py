import subprocess
from osgeo import ogr
import reproject
import pandas as pd

import utils.cache

import fishnet

def set_model_boundary_from_polygon_shapefile(shapefile_name, boundary_poly_file, projected_coordinate_system):
    """
    Function to set the model boundary using a shapefile, assuming polygon type:
    1. For structured mesh this will create a rectangle that creates an envelope on the polygon        
    
    """        
    ds = read(shapefile_name)
        
    layer = ds.GetLayer()
    # TODO use WKT rather than classes
    spatialRef = layer.GetSpatialRef()
    
    if projected_coordinate_system == spatialRef:
        print 'No transform required ... continuing'

    else:
        # TODO use ogr2ogr?
        reproject.reproject_layer(ds, 
                                  src_cs=spatialRef, 
                                  dst_cs=projected_coordinate_system, 
                                  geom_type=ogr.wkbPolygon,
                                  copy_dest=boundary_poly_file)
                                          
        ds = read(boundary_poly_file)
        
        layer = ds.GetLayer()

    if not layer.GetFeatureCount() == 1:
        raise Exception('Expected just on layer, got %d', layer.GetFeatureCount())

    # TODO could support unstructured meshes here
    feature = layer.GetFeature(0)       
    geom = feature.GetGeometryRef()
    xmin, xmax, ymin, ymax = geom.GetEnvelope()  
    # TODO use nice dictionary    
    boundary = [xmin, xmax, ymin, ymax, projected_coordinate_system]

    return boundary


def getXYpairs(points_obj, feature_id):

    Layer = points_obj.GetLayer()
    Layer.ResetReading()
    #points_list = []
    point_ids = {}
    for feat in Layer:
        geom = feat.GetGeometryRef()
        geom_str = geom.ExportToWkt()
        #points_list += geom_str.split('(')[1].split(')')[0]
        coords = geom_str.split('(')[1].split(')')[0]
        coords = coords.split(' ')
        coords = (float(coords[0]), float(coords[1]))
        point_ids[str(feat.GetField(feature_id))] = coords
    Layer = None
    return point_ids


def points_shapefile_obj2dataframe(points_obj, feature_id=None):
    """
    Converts GIS object to dataframe, which is designed for point style objects only
    """
    points_dict = getXYpairs(points_obj, feature_id=feature_id)
    df_points = pd.DataFrame.from_dict(points_dict, orient='index')
    df_points.columns = ['Easting', 'Northing']
    df_points[feature_id] = df_points.index
    return df_points


# TODO cache!
"""
see www.gdal.org/ogr2ogr.html
"""
@utils.cache.money()
def ogr2ogr(src_datasource_name, dst_datasource_name, t_srs, clipdst):
    ret = subprocess.check_output([
        "ogr2ogr",
        "-skipfailures",
        "-overwrite",
        "-t_srs", t_srs,
        "-clipdst", clipdst,
        # "-clipdstlayer", clipdstlayer,
        dst_datasource_name, 
        src_datasource_name, 
        # layer,
    ])

    return ret

def read(dst_datasource_name):
    ds = ogr.Open(dst_datasource_name)
    if ds is None:
        raise Exception('Could not open ' + dst_datasource_name)
    driver = ds.GetDriver()
    ds = driver.Open(dst_datasource_name, 0)
    return ds

def reproject_clip_read(src_datasource_name, dst_datasource_name, t_srs, clipdst):
    ret = ogr2ogr(src_datasource_name, dst_datasource_name, t_srs, clipdst)
    return read(dst_datasource_name)



"""
GWModelBuilder
"""


def set_temporal_components(steady_state=True, start_time=None, end_time=None, time_step=None, date_index=None):
    obj = {}
    obj['start_time'] = start_time
    obj['end_time'] = end_time
    obj['time_step'] = time_step
    obj['duration'] = obj['end_time'] - obj['start_time']
    
    date_index = pd.date_range(start=start_time, end=end_time, freq=time_step)
    obj['dateindex'] = date_index

    intervals = []
    old_time = date_index[0]
    for index, time in enumerate(date_index):
        if index > 0:
            intervals += [time - old_time]
            old_time = time
            
    obj['intervals'] = intervals
    obj['steps'] = len(obj['intervals'])
    return obj



def define_structured_mesh(model_boundary, gridHeight, gridWidth):
    out_data_folder = 'structured_model_grid_%im' % int(gridHeight)

    model_mesh = fishnet.create_fishnet({
                                                "xmin":model_boundary[0], 
                                                "xmax":model_boundary[1], 
                                                "ymin":model_boundary[2], 
                                                "ymax":model_boundary[3], 
                                                "gridHeight":gridHeight,
                                                "gridWidth":gridWidth
                                            }, 
                                            spatialRef=model_boundary[4], 
                                            copy_dest=out_data_folder)
    # End if
    # return self.model_mesh
    model_mesh_centroids = build_centroids_array(gridHeight)


    return out_data_folder_grid, model_mesh_centroids, model_mesh

