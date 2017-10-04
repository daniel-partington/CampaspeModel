"""
create east and west boundary line strings and rasters from boundary polgon
"""

import json
import numpy
import os

from GIS import raster
from GIS import project


# TODO this is boilerplate from example
def load(shape, hk, dels, zbot):
    nlay, nrow, ncol = shape
    delr, delc = dels
    # Make list for stress period 1
    stageleft = 10.
    stageright = 10.
    bound_sp1 = []
    for il in range(nlay):
        for ir in range(nrow):
            condleft = hk[il, ir, 0] * (stageleft - zbot) * delc
            condright = hk[il, ir,  ncol - 1] * (stageright - zbot) * delc
            bound_sp1.append([il, ir, 0, stageleft, condleft])
            bound_sp1.append([il, ir, ncol - 1, stageright, condright])
    print('Adding ', len(bound_sp1), 'GHBs for stress period 1.')

    # Make list for stress period 2
    stageleft = 10.
    stageright = 0.

    bound_sp2 = []
    for il in range(nlay):
        for ir in range(nrow):
            condleft = hk[il, ir, 0] * (stageleft - zbot) * delc
            condright = hk[il, ir, ncol - 1] * (stageright - zbot) * delc
            bound_sp2.append([il, ir, 0, stageleft, condleft])
            bound_sp2.append([il, ir, ncol - 1, stageright, condright])
    print('Adding ', len(bound_sp2), 'GHBs for stress period 2.')

    # We do not need to add a dictionary entry for stress period 3.
    # Flopy will automatically take the list from stress period 2 and apply it
    # to the end of the simulation
    stress_period_data = {0: bound_sp1, 1: bound_sp2}

    return {"stress_period_data": stress_period_data}
    


""" break points """
west_murray = [144.3136837132862, -35.75643626236702]
east_murray = [144.8098159275796, -36.118135894677835]
eppalock = [144.5294271300142, -36.8483087795658]


def cut_boundary(unprojected_boundary_geojson, model_projection, unprojected_from, unprojected_to, name):
    clipped_and_projected_folder = os.path.dirname(unprojected_boundary_geojson)

    """ cut boundary """
    with open(unprojected_boundary_geojson) as f:
        boundary_coordinates = json.loads(f.read())['features'][0]['geometry']['coordinates'][0]

    unprojected_from_i = numpy.argmin([(d[0] - unprojected_from[0]) ** 2 + (d[1] - unprojected_from[1]) ** 2 for d in boundary_coordinates])
    unprojected_to_i = numpy.argmin([(d[0] - unprojected_to[0]) ** 2 + (d[1] - unprojected_to[1]) ** 2 for d in boundary_coordinates])

    cut_boundary_collection = json.dumps({
        "type": "FeatureCollection",
        "features": [{
            "type": "Feature",
            "properties": {},
            "geometry": {
                "type": "LineString",
                "coordinates": boundary_coordinates[unprojected_from_i: unprojected_to_i] if unprojected_from_i < unprojected_to_i else boundary_coordinates[unprojected_from_i:] + boundary_coordinates[:unprojected_to_i]
            }
        }]
    })
    
    unprojected_cut_boundary_geojson = os.path.join(clipped_and_projected_folder, 'UNPROJECTED_%s_boundary.json' % name)
    
    """ for documentation """
    with open(unprojected_cut_boundary_geojson, 'w') as f:
        f.write(cut_boundary_collection)
    
    """ project """

    cut_boundary_file = os.path.join(clipped_and_projected_folder, '%s_boundary.shp' % name)
    
    project.ogr2ogr(src_datasource_name=unprojected_cut_boundary_geojson, 
        dst_datasource_name=cut_boundary_file,
        t_srs=model_projection,
        f='ESRI Shapefile'
    )

    return cut_boundary_file


def make_east_west_raster(bbox, dels, shape, model_projection, clipped_and_projected_folder):
    unprojected_boundary_geojson = os.path.join(clipped_and_projected_folder, 'UNPROJECTED_GW_model_area.json')
    
    west_file = cut_boundary(unprojected_boundary_geojson=unprojected_boundary_geojson,
                                model_projection=model_projection,
                                unprojected_from=eppalock,
                                unprojected_to=west_murray,
                                name="west")
    
    east_file = cut_boundary(unprojected_boundary_geojson=unprojected_boundary_geojson,
                                model_projection=model_projection,
                                unprojected_from=east_murray,
                                unprojected_to=eppalock,
                                name="east")


    min_x, max_x, min_y, max_y = bbox
    delr, delc = dels
    nlay, nrow, ncol = shape
    
    west_raster = raster.rasterize(geom_file=west_file, 
                            shape=(nrow, ncol), 
                            geotransform=(min_x, delr, 0, min_y, 0, delc))


    east_raster = raster.rasterize(geom_file=east_file, 
                            shape=(nrow, ncol), 
                            geotransform=(min_x, delr, 0, min_y, 0, delc))



    return west_raster, east_raster


if __name__ == '__main__':
    import matplotlib.pylab as plt
    import config
    import boundary

    shape = (1, 30, 20)  # (nlay, nrow, ncol)
    boundary_params = boundary.boundary(boundary_file=config.boundary_file,
                                    clipped_and_projected_folder=config.clipped_and_projected_folder,
                                    shape=shape)

    west_raster, east_raster = make_east_west_raster(bbox=boundary_params['bbox'], 
                            dels=boundary_params['dels'], 
                            shape=shape, 
                            model_projection=boundary_params['model_projection'], 
                            clipped_and_projected_folder=config.clipped_and_projected_folder)
    
    rows, cols = numpy.where(west_raster)    

    # stage = elev[lay, row, col]  + stage_parameter
    # conductance = gridXdim * (elev[lay, row, col] - elev[lay + 1, row, col]) * conductivity_param
    # [[lay, rows[i], cols[i], stage, conductance] for i in range(len(rows))]

    plt.imshow(west_raster)
    plt.show()
    plt.imshow(east_raster)
    plt.show()
