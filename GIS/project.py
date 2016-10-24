import json
import os
import subprocess
import numpy as np
from osgeo import ogr
from osgeo import osr
from utils import cache

"""
see www.gdal.org/ogr2ogr.html
"""
@cache.money()
def ogr2ogr(src_datasource_name, dst_datasource_name, layer_name_list=[], **kwargs):
    if os.path.exists(dst_datasource_name):
        os.remove(dst_datasource_name)

    extras = []
    for key in kwargs.keys():
        extras.append('-' + key)
        if isinstance(kwargs[key], (list, tuple)):
            for v in kwargs[key]:
                extras.append(str(v))
        else:
            extras.append(str(kwargs[key]))

    ret = subprocess.check_output([
        "ogr2ogr",
        "-skipfailures",
        "-overwrite",
    ] + extras + [
        dst_datasource_name, 
        src_datasource_name
    ] + layer_name_list)

    return ret


"""
    dest_text needs to work with ogr:SetFromUserInput 
    see http://www.gdal.org/classOGRSpatialReference.html#aec3c6a49533fe457ddc763d699ff8796
"""

def geojson_transformer(src_text, dst_text):
    dst_srs = osr.SpatialReference()
    dst_srs.SetFromUserInput(dst_text)

    src_srs = osr.SpatialReference()
    src_srs.SetFromUserInput(src_text)

    transform = osr.CoordinateTransformation(src_srs, dst_srs)

    def transform_geojson(geojson):
        geom = ogr.CreateGeometryFromJson(geojson)

        geom.Transform(transform)

        return geom.ExportToJson()

    return transform_geojson


def project_points(xs, ys, projector):

    reprojected = json.loads(projector(json.dumps({"type": "MultiPoint", "coordinates": 
        np.array([xs, ys]).T.tolist()
    })))
    return np.array(reprojected['coordinates']).T


"""
gets projection of layer_name of file_name
"""
def get_source_projection(file_name, layer_name):
    dataset = ogr.Open(file_name)
    layer = dataset.GetLayer(layer_name)
    projection_wkt = layer.GetSpatialRef().ExportToWkt()
    return projection_wkt

"""
reprojects a bounding box based 
"""

def get_source_bounding_box(bbox, bbox_projection, projection_wkt):

    min_x, max_x, min_y, max_y = bbox

    model_to_river = geojson_transformer(bbox_projection, projection_wkt)
    
    min_x, min_y = json.loads(model_to_river(json.dumps(
        {"type": "Point", "coordinates": [min_x, min_y]}
    )))['coordinates']
    max_x, max_y = json.loads(model_to_river(json.dumps(
        {"type": "Point", "coordinates": [max_x, max_y]}
    )))['coordinates']

    return [min_x, min_y, max_x, max_y]

