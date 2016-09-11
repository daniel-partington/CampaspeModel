from osgeo import ogr
from osgeo import osr
import subprocess

from utils import cache

"""
see www.gdal.org/ogr2ogr.html
"""
@cache.money()
def ogr2ogr(src_datasource_name, dst_datasource_name, layer_name_list=[], **kwargs):
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


"""
reprojects a bounding box based on the projection system of layer_name of file_name
"""

def get_source_bounding_box(bbox, bbox_projection, file_name, layer_name):
    import json
    from osgeo import ogr

    dataset = ogr.Open(file_name)
    layer = dataset.GetLayer(layer_name)
    river_projection = layer.GetSpatialRef().ExportToWkt()

    min_x, max_x, min_y, max_y = bbox

    model_to_river = geojson_transformer(bbox_projection, river_projection)
    
    min_x, min_y = json.loads(model_to_river(json.dumps(
        {"type": "Point", "coordinates": [min_x, min_y]}
    )))['coordinates']
    max_x, max_y = json.loads(model_to_river(json.dumps(
        {"type": "Point", "coordinates": [max_x, max_y]}
    )))['coordinates']

    return [min_x, min_y, max_x, max_y]

