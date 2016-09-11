from dateutil.relativedelta import relativedelta
import datetime
import os
import numpy as np

from GIS import project
from GIS import raster

"""
The river model, given 
    * levels for each reach
    * inflows
    * outflow
will run and produce
    * new levels for each reach
    * predicted outflow

GW model will take 
    * new levels
and give 
    * net inflows and outflows
"""

def load(data_folder, clipped_and_projected_folder, 
        model_projection, bbox, 
        shape,
        hk):
    
    length_raster, reaches = rasterize_reaches(data_folder=data_folder, 
            clipped_and_projected_folder=clipped_and_projected_folder, 
            model_projection=model_projection, 
            bbox=bbox,
            shape=shape)

    return by_time(reaches, hk)


start_date = datetime.datetime(2012, 1, 1)

def by_time(reaches, hk):

    """
    set levels by stress period
    """

    riv = {}

    for month in range(12):
        first_of_month = start_date + relativedelta(months=month)
        
        day_of_year = first_of_month.timetuple().tm_yday
        
        time_reaches = []
        for reach in reaches:
            lay, row, col =(0, reach['row'], reach['col'])
            
            width = 10.
            # TODO conductance?
            cond = width * reach["length"] * hk[lay, row, col]  # the riverbed hydraulic conductance
            # TODO need elevation and surface water model inputs
            stage = np.abs(6 - month)  # the head in the river
            rbot = - 5.  # is the elevation of the bottom of the riverbed
            # TODO one based?
            time_reaches.append([lay, row, col, stage, cond, rbot])
        
        riv[day_of_year] = time_reaches

    return {"stress_period_data": riv}



"""
BOM geofabric
    * Download from the [BOM FTP server](ftp://ftp.bom.gov.au/anon/home/geofabric/).
    * It can be explored using [BOM WMS connection](http://geofabric.bom.gov.au/documentation/) or the [National Map](http://nationalmap.nicta.com.au/).

"""


def rasterize_reaches(data_folder, clipped_and_projected_folder, model_projection, bbox, shape):

    """
    project and clip BOM geofabric data
    """
    river_file = os.path.join(data_folder, "SH_Network_GDB", "SH_Network.gdb")
    river_layer = "AHGFNetworkStream"
    projected_clipped_river_file = os.path.join(clipped_and_projected_folder, "projected_clipped_river_file.shp")

    # clip first so we don't reproject the entire geodatabase, faster than clipdst
    river_bounding_box = project.get_source_bounding_box(bbox=bbox, bbox_projection=model_projection, 
                        file_name=river_file, layer_name=river_layer)

    # clipdst=projected_boundary_file removes some parts of the Murray, so better use the bbox
    project.ogr2ogr(src_datasource_name=river_file, 
        dst_datasource_name=projected_clipped_river_file,
        layer_name_list=[river_layer], 
        clipsrc=river_bounding_box,
        t_srs=model_projection,
    )

    """
    cut reaches by mesh
    """
    reaches = []
    def rasterize_riv(source_layer, clipping_poly, row, col):
        properties = []
        for feature in source_layer:
            geom = feature.GetGeometryRef()
            if (feature.GetField('Name') == 'MURRAY RIVER' or feature.GetField('Name') == 'CAMPASPE RIVER')\
                and clipping_poly.Intersects(geom):
                
                clipped = clipping_poly.Clone().Intersection(geom)
                length = clipped.Length()

                reaches.append({
                    "length": length,
                    "row": row,
                    "col": col,
                    'HydroID': feature.GetField('HydroID'),
                    'NextDownID': feature.GetField('NextDownID'),
                    'Name': feature.GetField('Name'),
                    })

                properties.append(length)

        return sum(properties)

    nlay, nrow, ncol = shape

    length_raster =  raster.ogr_to_raster(
        projected_clipped_river_file,
        (nrow, ncol),
        bbox,
        rasterize_riv) 

    return length_raster, reaches



if __name__ == '__main__':
    import boundary
    import config

    shape = (1, 30, 20)  # (nlay, nrow, ncol)
    boundary_params, boundary_inputs = boundary.load(boundary_file=os.path.join(config.data_folder, "GW_model_area", "GW_model_area.shp"), 
                                    clipped_and_projected_folder=config.clipped_and_projected_folder,
                                    shape=shape)

    length_raster, reaches = rasterize_reaches(data_folder=config.data_folder, 
            clipped_and_projected_folder=config.clipped_and_projected_folder, 
            model_projection=boundary.model_projection, 
            bbox=boundary_params['bbox'],
            shape=shape)

    hk = 1 * np.ones(shape)
    reaches_by_stress = by_time(reaches, hk)

    import matplotlib.pylab as plt
    plt.imshow(length_raster)
    plt.show()


