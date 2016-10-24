import os
import numpy as np
from GIS import io
from GIS import raster
from GIS import project
from GIS import interpolate

# Define basic model parameters
# EPSG:28355 is the code for gda94 mga zone 55; http://spatialreference.org/ref/epsg/gda94-mga-zone-55/
model_projection = "EPSG:28355"


def boundary(data_folder, clipped_and_projected_folder, shape):
    boundary_file = os.path.join(data_folder, "GW_model_area", "GW_model_area.shp")

    # this is just for documentation
    project.ogr2ogr(src_datasource_name=boundary_file,
        dst_datasource_name=os.path.join(clipped_and_projected_folder, "UNPROJECTED_GW_model_area.json"),
        t_srs="EPSG:4326",
        f='GeoJSON'
    )

    """
    project boundary file
    """
    projected_boundary_file = os.path.join(clipped_and_projected_folder, 'projected_' + os.path.basename(boundary_file))

    project.ogr2ogr(src_datasource_name=boundary_file, 
        dst_datasource_name=projected_boundary_file, 
        t_srs=model_projection)

    boundary_geom = io.load_geom(projected_boundary_file)

    """
    rasterize boundary file
    """

    min_x, max_x, min_y, max_y = boundary_geom.GetEnvelope()

    # buffer and round bounding box
    buffer_m = 5000.
    round_m = 1000.
    min_x = round_m * np.floor((min_x - buffer_m) / round_m)
    max_x = round_m * np.ceil((max_x + buffer_m) / round_m)
    min_y = round_m * np.floor((min_y - buffer_m) / round_m)
    max_y = round_m * np.ceil((max_y + buffer_m) / round_m)

    Lx = max_x - min_x
    Ly = max_y - min_y

    nlay, nrow, ncol = shape
    delr = Lx / ncol
    delc = Ly / nrow

    mesh_centroids = interpolate.make_mesh_centroids(bbox=(min_x, max_x, min_y, max_y), 
                                                    dels=(delr, delc), 
                                                    shape=(nlay, nrow, ncol))

    return {
        "bbox": (min_x, max_x, min_y, max_y),
        "dels": (delr, delc),
        "Ls": (Lx, Ly),
        "mesh_centroids": mesh_centroids,
        "model_projection": model_projection,
        "projected_boundary_file": projected_boundary_file,
    }


def load(projected_boundary_file, shape, bbox, dels):
    min_x, max_x, min_y, max_y = bbox
    delr, delc = dels
    nlay, nrow, ncol = shape
    ibound = raster.rasterize(geom_file=projected_boundary_file, 
                            shape=(nrow, ncol), 
                            geotransform=(min_x, delr, 0, min_y, 0, delc))

    return {
        "ibound": ibound,
        "strt": 10. * np.ones(shape, dtype=np.float32)
    }
