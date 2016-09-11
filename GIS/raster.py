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


"""
Converts source layer to raster (array) given a custom function 
    ```rasterize(layer, clipping poly)``` 
which takes the layer and polygon for once cell and returns cell's value.
More generic, but much slower than the above.
"""

def ogr_to_raster(source_file, shape, bbox, rasterize, filter=True):
    nrow, ncol = shape

    # open source (file)
    source_ds = ogr.Open(source_file)
    source_layer = source_ds.GetLayer()

    # to ensure cut geometries (like ibound) remain consistent
    min_x, max_x, min_y, max_y = bbox

    pixelWidth = (max_x - min_x) / ncol
    pixelHeight = (max_y - min_y) / nrow

    array = np.empty((nrow, ncol)).tolist()
    for row in range(nrow):
        for col in range(ncol):

            # Create clipping polygon
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(min_x + (col) * pixelWidth, min_y + (row) * pixelHeight)
            ring.AddPoint(min_x + (col + 1) * pixelWidth, min_y + (row) * pixelHeight)
            ring.AddPoint(min_x + (col + 1) * pixelWidth, min_y + (row + 1) * pixelHeight)
            ring.AddPoint(min_x + (col) * pixelWidth, min_y + (row + 1) * pixelHeight)
            ring.AddPoint(min_x + (col) * pixelWidth, min_y + (row) * pixelHeight)
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)
            
            # speeds things up but breaks them too:  source_layer.SetSpatialFilterRect(x_min+col*pixelWidth, x_min+(col+1)*pixelWidth, y_min+row*pixelHeight, y_min+(row+1)*pixelHeight)
            # speeds them up, doesn't seem to break them
            if filter:
                source_layer.SetSpatialFilter(poly)

            source_layer.ResetReading()

            array[row][col] = rasterize(source_layer, poly, row, col)

    return array
