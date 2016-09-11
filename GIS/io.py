from osgeo import ogr

def load_geom(file_name, driver_type='ESRI Shapefile'):
    driver = ogr.GetDriverByName(driver_type)
    dataSource = driver.Open(file_name, 0) # 0 means read-only. 1 means writeable.
    layer = dataSource.GetLayer()
    assert layer.GetFeatureCount() == 1
    feature = layer.GetFeature(0)
    geom = feature.GetGeometryRef()
    
    return geom.Clone()
