# Code derived from:
# https://jgomezdans.github.io/gdal_notes/reprojection.html
from osgeo import osr, ogr, gdal, gdalconst
import os

"""
This contains functions for reprojecting:
1. Rasters from one projected coordinate system to another, and allowing options:
    a. create_copy
    b. path_dest
    c. pixel_size
    d. new bounding box
    
2. Layer files, allowing the following options:
    a. create_copy
    b. path_dest    
"""

def coordinate_transform(pointX, pointY, coordTransform):
    """
    A function to transform a point on the horizontal plane
    using CoordinateTransformation object from GDAL osr.
    """
    # create a geometry from coordinates
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(pointX, pointY)
    # transform point
    point.Transform(coordTransform)

    return point.GetX(), point.GetY()

def reproject_raster(ds, # raster dataset
                     pixel_spacing=None, 
                     src_cs=None, 
                     dst_cs=None, 
                     reproj_method=gdal.GRA_Bilinear,
                     create_copy=False,
                     copy_dest=None,
                     raster_driver="GTiff",
                     set_bounds=None):
    
    """
    A function to reproject and resample a GDAL dataset from within 
    Python. The idea here is to reproject from one system to another, as well
    as to change the pixel size. The procedure is slightly long-winded, but
    goes like this:
    
    1. Set up the two Spatial Reference systems.
    2. Open the original dataset, and get the geotransform
    3. Calculate bounds of new geotransform by projecting the UL corners 
    4. Calculate the number of pixels with the new projection & spacing
    5. Create an in-memory raster dataset
    6. Perform the projection
    
    Creates a copy of the projected raster to file if create_copy == True

    Returns projected gdal raster object at the end    
    
    """
    
    # Define the projection for "from", and get from dataset if not provided
    if src_cs == None:    
        # Check raster GCS
        prj = ds.GetProjection()
        src_cs=osr.SpatialReference(wkt=prj)
        
    tx = osr.CoordinateTransformation(src_cs, dst_cs)
    
    # Get the Geotransform vector
    geo_t = ds.GetGeoTransform()
    x_size = ds.RasterXSize # Raster xsize
    y_size = ds.RasterYSize # Raster ysize
    NO_DATA = ds.GetRasterBand(1).GetNoDataValue()

    # Work out the boundaries of the new dataset in the target projection
    (ulx, uly) = coordinate_transform(geo_t[0], geo_t[3], tx)
    (lrx, lry) = coordinate_transform(geo_t[0] + geo_t[1]*x_size, \
                                          geo_t[3] + geo_t[5]*y_size, tx)
    # Now, we create an in-memory raster
    mem_drv = gdal.GetDriverByName( 'MEM' )
    # The size of the raster is given the new projection and pixel spacing
    # Using the values we calculated above. Also, setting it to store one band
    # and to use Float32 data type.

    if pixel_spacing == None: # use default spacing from input raster
        pixel_spacing = geo_t[1]
    
    if set_bounds != None:    
        # Define extent from set_bounds if it exists        
        (ulx, lrx, lry, uly) = set_bounds
        
    dest = mem_drv.Create('', int((lrx - ulx)/pixel_spacing), \
            int((uly - lry)/pixel_spacing), 1, gdal.GDT_Float32)
    # Calculate the new geotransform
    new_geo = (ulx, pixel_spacing, geo_t[2], \
               uly, geo_t[4], -pixel_spacing)
                           
    # Set the geotransform
    dest.SetGeoTransform(new_geo)
    dest.SetProjection(dst_cs.ExportToWkt())
    dest.GetRasterBand(1).SetNoDataValue(NO_DATA)
    # Perform the projection/resampling 
    gdal.ReprojectImage(ds, dest, src_cs.ExportToWkt(), dst_cs.ExportToWkt(), reproj_method)

    #reproj_ds = dest.ReadAsArray()
    # Let's save it as a GeoTIFF.
    if create_copy:     
        driver = gdal.GetDriverByName(raster_driver)
        metadata = driver.GetMetadata()
        # Create a copy of the transformed raster for external use        
        if metadata.has_key(gdal.DCAP_CREATECOPY) and metadata[gdal.DCAP_CREATECOPY] == 'YES':
            if os.path.exists(copy_dest):
                #driver.DeleteDataSource(copy_dest)            
                os.remove(copy_dest)
            driver.CreateCopy(copy_dest, dest, 0) #, [ 'TILED=YES', 'COMPRESS=PACKBITS' ] )
            #dst_ds = None
        else:
            print 'Driver %s does not support CreateCopy() method.' % raster_driver    

    return dest


                  
# Code derived from: 
# https://pcjericks.github.io/py-gdalogr-cookbook/projection.html

def reproject_layer(lyr_src, 
                    src_cs=None, 
                    dst_cs=None, 
                    create_copy=False,
                    copy_dest=None,
                    geom_type=None):
    
    """
    Reprojects a shapefile from one PCS to another. 
    Only works for single types of geomoetries at present
    
    :param geom_type: Type of geometry/geometries in the layer e.g. ogr.wkbMultiPolygon
    """    
    
    driver = ogr.GetDriverByName('ESRI Shapefile')
    
    # input SpatialReference
    #inSpatialRef = osr.SpatialReference()
    #inSpatialRef.ImportFromEPSG(epsg_from)
    if src_cs == None:
        inSpatialRef = lyr_src.GetSpatialRef()        
    else:
        inSpatialRef = src_cs    
    
    # output SpatialReference
    #outSpatialRef = osr.SpatialReference()
    #outSpatialRef.ImportFromEPSG(epsg_to)
    outSpatialRef = dst_cs    
    
    # create the CoordinateTransformation
    coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
    
    # get the input layer
    #inDataSet = lyr_src # driver.Open(r'c:\data\spatial\basemap.shp')
    inLayer = lyr_src.GetLayer()

    # create the output layer
    #outputShapefile = r'c:\data\spatial\basemap_4326.shp'
    if os.path.exists(copy_dest):
        os.remove(copy_dest)
    outDataSet = driver.CreateDataSource(copy_dest)
    outLayer = outDataSet.CreateLayer("", geom_type=geom_type)
    
    # add fields
    inLayerDefn = inLayer.GetLayerDefn()

    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)
    
    # get the output layer's feature definition
    outLayerDefn = outLayer.GetLayerDefn()
    
    # loop through the input features
    inFeature = inLayer.GetNextFeature()
    while inFeature:
        # get the input geometry
        geom = inFeature.GetGeometryRef()
        # reproject the geometry
        geom.Transform(coordTrans)
        # create a new feature
        outFeature = ogr.Feature(outLayerDefn)
        # set the geometry and attribute
        outFeature.SetGeometry(geom)
        for i in range(0, outLayerDefn.GetFieldCount()):
            outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
        # add the feature to the shapefile
        outLayer.CreateFeature(outFeature)
        # destroy the features and get the next input feature
        outFeature.Destroy()
        inFeature.Destroy()
        inFeature = inLayer.GetNextFeature()
    
    # close the shapefiles
    #outDataSet = None #.Destroy() 

    #create the prj projection file   
    outSpatialRef.MorphToESRI()   
    prj_file = open(copy_dest[:-4] + '.prj', 'w')   
    prj_file.write(outSpatialRef.ExportToWkt())   
    prj_file.close() 

    return outDataSet    

if __name__ == "__main__":
    # Test for raster reprojection
    dataset = r"C:\Workspace\part0075\MDB modelling\Campaspe_model\GIS\GIS_preprocessed\Hydrogeological_Unit_Layers\qa_1t_bb"
    pixel_spacing=10000 
    epsg_from=None 
    epsg_to=28355 
    reproj_method=gdal.GRA_Bilinear
    create_copy=True
    copy_dest= r"C:\Workspace\part0075\MDB modelling\integrated\Modules\Groundwater\GISInterface\GDALInterface\test.tif"
    raster_driver="EHdr" # "GTiff" 
    
    ds = gdal.Open(dataset, gdalconst.GA_ReadOnly)    

    Proj_CS = osr.SpatialReference()
    Proj_CS.ImportFromEPSG(28355) # 28355 is the code for gda94 mga zone 55; http://spatialreference.org/ref/epsg/gda94-mga-zone-55/

    #set_corner=[ulx, uly, 0, lrx, lry, 0]
    reproject_raster( ds, # raster dataset
                      pixel_spacing=pixel_spacing, 
                      src_cs=None, 
                      dst_cs=Proj_CS, 
                      reproj_method=gdal.GRA_NearestNeighbour,
                      create_copy=create_copy,
                      copy_dest=copy_dest,
                      raster_driver="GTiff",
                      set_bounds=None)

    # Test for layer reprojection
    driver = ogr.GetDriverByName("ESRI Shapefile")        
    ds = driver.Open(r"C:\Workspace\part0075\MDB modelling\Campaspe_model\GIS\GIS_preprocessed\Surface_Water\Streams\Campaspe_Riv.shp", 0)
    poly_obj = ds.GetLayer()    
    if poly_obj == None:
        print 'Could not open '
   
    srs = poly_obj.GetSpatialRef()
    #print srs.ExportToWkt()
    copy_dest= r"C:\Workspace\part0075\MDB modelling\integrated\Modules\Groundwater\GISInterface\GDALInterface\test_model.shp"

    ds = reproject_layer(ds, 
                         src_cs=srs, 
                         dst_cs=Proj_CS, 
                         create_copy=True,
                         copy_dest=copy_dest,    
                         geom_type=ogr.wkbMultiLineString)
