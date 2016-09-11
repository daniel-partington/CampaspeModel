from osgeo import gdal
from osgeo import ogr
import numpy as np

import io


def rasterize(geom_file, shape, geotransform, all_touched=False):
    # Thanks https://github.com/mapbox/rasterio/issues/45
    assert len(shape) == 2
    nrow, ncol = shape

    geom = io.load_geom(geom_file)
    sr = geom.GetSpatialReference()
    target_ds = gdal.GetDriverByName("MEM").Create("", ncol, nrow, gdal.GDT_Byte)
    target_ds.SetProjection(sr.ExportToWkt())
    target_ds.SetGeoTransform(geotransform)
    temp_features = ogr.GetDriverByName("Memory").CreateDataSource("")
    lyr = temp_features.CreateLayer("poly", srs=sr)
    feature = ogr.Feature(lyr.GetLayerDefn())
    feature.SetGeometryDirectly(geom)
    lyr.CreateFeature(feature)
    kwargs = {}
    if all_touched:
        kwargs['options'] = ["ALL_TOUCHED=TRUE"]
    gdal.RasterizeLayer(target_ds, [1], lyr, burn_values=[1], **kwargs)
    arr = target_ds.GetRasterBand(1).ReadAsArray(0, 0, ncol, nrow).astype(np.bool)
    del target_ds

    return arr