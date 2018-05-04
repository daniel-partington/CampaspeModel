import os

from GIS import project
from GIS import raster
from utils import cache


"""
BOM geofabric
    * Download from the [BOM FTP server](ftp://ftp.bom.gov.au/anon/home/geofabric/).
    * It can be explored using [BOM WMS connection](http://geofabric.bom.gov.au/documentation/) or the [National Map](http://nationalmap.nicta.com.au/).


ftp://ftp.bom.gov.au/anon/home/geofabric/SH_Cartography_GDB_V2_1_1.zip
AHGFCanalLine

ftp://ftp.bom.gov.au/anon/home/geofabric/SH_Network_GDB_V2_1_1.zip
AHGFNetworkStream

"""

def geofabric(river_file, river_layer, clipped_and_projected_folder, model_projection, bbox):
    """
    project and clip BOM geofabric data
    """
    projected_clipped_river_file = os.path.join(clipped_and_projected_folder, "projected_clipped_%s_file.shp" % river_layer)

    # clip first so we don't reproject the entire geodatabase, faster than clipdst
    river_projection = project.get_source_projection(file_name=river_file, layer_name=river_layer)
    river_bounding_box = project.get_source_bounding_box(bbox=bbox, bbox_projection=model_projection, 
                                                    projection_wkt=river_projection)
        
    # clipdst=projected_boundary_file removes some parts of the Murray, so better use the bbox
    project.ogr2ogr(src_datasource_name=river_file, 
        dst_datasource_name=projected_clipped_river_file,
        layer_name_list=[river_layer], 
        clipsrc=river_bounding_box,
        t_srs=model_projection,
    )

    # just for documentation
    project.ogr2ogr(src_datasource_name=river_file, 
        dst_datasource_name=os.path.join(clipped_and_projected_folder, "UNPROJECTED_%s.json" % river_layer),
        layer_name_list=[river_layer], 
        clipsrc=river_bounding_box,
        t_srs="EPSG:4326",
        f='GeoJSON'
    )

    return projected_clipped_river_file


@cache.money()
def rasterize_reaches(projected_clipped_river_file, bbox, shape, names):
    """
    cut reaches by mesh
    """
    reaches = []
    def rasterize_riv(source_layer, clipping_poly, row, col):
        properties = []
        for feature in source_layer:
            geom = feature.GetGeometryRef()
            if (feature.GetField('Name') in names)\
                and clipping_poly.Intersects(geom):

                clipped = clipping_poly.Clone().Intersection(geom)
                length = clipped.Length()
                properties.append(length)

                reach = {
                    "length": length,
                    "row": row,
                    "col": col,
                    "clipped": clipped.ExportToJson()
                }
                
                # NextDownID an HydroID can be used for routing
                reach.update(feature.items())
                reaches.append(reach)


        return sum(properties)

    nlay, nrow, ncol = shape

    riv_raster =  raster.ogr_to_raster(
        projected_clipped_river_file,
        (nrow, ncol),
        bbox,
        rasterize_riv) 

    return riv_raster, reaches


if __name__ == '__main__':
    import boundary
    import config
    import matplotlib.pylab as plt

    shape = (1, 30, 20)  # (nlay, nrow, ncol)
    boundary_params = boundary.boundary(boundary_file=config.boundary_file, 
                                    clipped_and_projected_folder=config.clipped_and_projected_folder,
                                    shape=shape)

    bbox = boundary_params['bbox']
    model_projection = boundary_params['model_projection']
    # seems to be higher resolution but similar to AHGFNetworkStream from SH_Network_GDB
    # TODO is network the same from SH_Network_GDB/AHGFNetworkStream and SH_Cartography_GDB/AHGFMappedStream
    projected_clipped_river_file = geofabric(river_file=os.path.join(config.data_folder, "SH_Cartography_GDB", "SH_Cartography.gdb"), 
                                            river_layer="AHGFMappedStream", 
                                            clipped_and_projected_folder=config.clipped_and_projected_folder, 
                                            model_projection=model_projection, bbox=bbox)

    river_length_raster, reaches = rasterize_reaches(projected_clipped_river_file=projected_clipped_river_file,
                                                    bbox=bbox, shape=shape,
                                                    names=['MURRAY RIVER', 'CAMPASPE RIVER'])


    projected_clipped_canal_file = geofabric(river_file=os.path.join(config.data_folder, "SH_Cartography_GDB", "SH_Cartography.gdb"), 
                                            river_layer="AHGFCanalLine", 
                                            clipped_and_projected_folder=config.clipped_and_projected_folder, 
                                            model_projection=model_projection, bbox=bbox)

    canal_length_raster, canal_reaches = rasterize_reaches(projected_clipped_river_file=projected_clipped_canal_file,
                                                                bbox=bbox, shape=shape,
                                                                names=["WARANGA WESTERN CHANNEL"])

    print reaches[:2]
    print canal_reaches[:2]

    plt.imshow(river_length_raster)
    plt.show()
    plt.imshow(canal_length_raster)
    plt.show()

