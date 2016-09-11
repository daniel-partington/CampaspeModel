from osgeo import ogr
from osgeo import osr
import subprocess

"""
see www.gdal.org/ogr2ogr.html
"""
def ogr2ogr(src_datasource_name, dst_datasource_name, **kwargs):
    extras = []
    for key, val in kwargs.iteritems():
        extras.append('-' + key)
        extras.append(val)

    ret = subprocess.check_output([
        "ogr2ogr",
        "-skipfailures",
        "-overwrite",
    ] + extras + [
        dst_datasource_name, 
        src_datasource_name
    ])

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
