# TODO
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))


from dateutil.relativedelta import relativedelta
import datetime
import os
import json
import numpy as np
from osgeo import ogr

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
    
    # length_raster, reaches = rasterize_reaches(data_folder=data_folder, 
    #         clipped_and_projected_folder=clipped_and_projected_folder, 
    #         model_projection=model_projection, 
    #         bbox=bbox,
    #         shape=shape)
    start_date = datetime.datetime(2012, 1, 1)

    # return by_time(reaches, hk, start_date)
    return []




def by_time(reaches, hk, start_date):

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
            
            # TODO conductance?
            width = 10.
            cond = width * reach["length"] * hk[lay, row, col]  # the riverbed hydraulic conductance
            # TODO need elevation and surface water model inputs
            stage = np.abs(6 - month)  # the head in the river
            rbot = - 5.  # is the elevation of the bottom of the riverbed
            # TODO one based?
            time_reaches.append([lay, row, col, stage, cond, rbot])
        
        riv[day_of_year] = time_reaches

    return {"stress_period_data": riv}


def closest_hydroID_to_site(site, projected_clipped_river_file):
    source_ds = ogr.Open(projected_clipped_river_file)
    source_layer = source_ds.GetLayer()
    
    site_geom = ogr.CreateGeometryFromJson(json.dumps({"type": "Point", "coordinates": [site["x"], site["y"]]}))
    distances = [feature.GetGeometryRef().Distance(site_geom) for feature in source_layer]
    min_i = np.argmin(distances)
    
    return {
        "closest_HydroID": source_layer.GetFeature(min_i).GetField('HydroID'),
        "closest_HydroID_distance": distances[min_i]
    }


# def downstream(HydroID, reaches):
#     # note reaches have been cut by model cells so are not unique on HydroId

#     # network = collections.defaultdict(set)
#     network = {}
#     for reach in reaches:
#         # network[reach['HydroID']].add(reach['NextDownID'])
#         network[reach['HydroID']] = reach['NextDownID']

#     print len(network)

#     start = HydroID
#     print start, network
#     new_network = []
#     while start in network:
#         new_network.append(start)
#         start = network[start]

#     print len(new_network)


if __name__ == '__main__':
    import boundary
    import config
    
    import data_water_vic_gov_au
    import bom_geofabric

    shape = (1, 30, 20)  # (nlay, nrow, ncol)
    boundary_params = boundary.boundary(boundary_file=config.boundary_file,
                                    clipped_and_projected_folder=config.clipped_and_projected_folder,
                                    shape=shape)

    model_start = datetime.date(2015, 1, 1)
    model_end = datetime.date(2015, 12, 31)
    
    # data_folder = config.data_folder

    flow_and_level_sites = data_water_vic_gov_au.read_relevant_sites(sites_folder=os.path.join(config.data_folder, "data.water.vic", "download.20160923190321"),
                            clipped_and_projected_folder=config.clipped_and_projected_folder,
                            model_projection=boundary_params["model_projection"], 
                            projected_boundary_file=boundary_params['projected_boundary_file'],
                            model_start=model_start, model_end=model_end)

    bbox = boundary_params['bbox']
    model_projection = boundary_params['model_projection']

    # seems to be higher resolution but similar to AHGFNetworkStream from SH_Network_GDB
    projected_clipped_river_file = bom_geofabric.geofabric(river_file=os.path.join(config.data_folder, "SH_Cartography_GDB", "SH_Cartography.gdb"), 
                                            river_layer="AHGFMappedStream", 
                                            clipped_and_projected_folder=config.clipped_and_projected_folder, 
                                            model_projection=model_projection, bbox=bbox)

    murray_raster, murray_reaches = bom_geofabric.rasterize_reaches(projected_clipped_river_file=projected_clipped_river_file,
                                                    bbox=bbox, shape=shape, 
                                                    names=['MURRAY RIVER'])
    
    campaspe_raster, campaspe_reaches = bom_geofabric.rasterize_reaches(projected_clipped_river_file=projected_clipped_river_file,
                                                    bbox=bbox, shape=shape, 
                                                    names=['CAMPASPE RIVER'])

    projected_clipped_canal_file = bom_geofabric.geofabric(river_file=os.path.join(config.data_folder, "SH_Cartography_GDB", "SH_Cartography.gdb"), 
                                            river_layer="AHGFCanalLine", 
                                            clipped_and_projected_folder=config.clipped_and_projected_folder, 
                                            model_projection=model_projection, bbox=bbox)

    waranga_raster, waranga_reaches = bom_geofabric.rasterize_reaches(projected_clipped_river_file=projected_clipped_canal_file,
                                                                bbox=bbox, shape=shape,
                                                                names=["WARANGA WESTERN CHANNEL"])

    print closest_hydroID_to_site(
        site=[site for site in flow_and_level_sites if site['Site Name'] == "CAMPASPE RIVER @ EPPALOCK"][0], 
        projected_clipped_river_file=projected_clipped_river_file
    )
    print [site for site in flow_and_level_sites if site['Site Name'] == "CAMPASPE RIVER @ EPPALOCK"][0]

    print closest_hydroID_to_site(
        site=[site for site in flow_and_level_sites if site['Site Name'] == "MURRAY RIVER @ ECHUCA"][0], 
        projected_clipped_river_file=projected_clipped_river_file
    )
    print [site for site in flow_and_level_sites if site['Site Name'] == "MURRAY RIVER @ ECHUCA"][0]

    # TODO use historical average of "WARANGA WESTERN CHANNEL @ 0 MILE 40 CHAIN STRUCTURE"
    # TODO stage and level?
    print closest_hydroID_to_site(
        site=[site for site in flow_and_level_sites if site['Site Name'] == "DEAKIN MAIN DRAIN @ OUTFALL (ECHUCA)"][0], 
        projected_clipped_river_file=projected_clipped_canal_file
    )
    print [site for site in flow_and_level_sites if site['Site Name'] == "DEAKIN MAIN DRAIN @ OUTFALL (ECHUCA)"][0]


    # downstream_stages = downstream(campaspe_river[0]['closest_HydroID'], campaspe_reaches)
    
    # "CAMPASPE RIVER @ EPPALOCK"

    import matplotlib.pylab as plt
    plt.imshow(campaspe_raster)
    plt.show()
    
    plt.imshow(murray_raster)
    plt.show()

    plt.imshow(waranga_raster)
    plt.show()

    # for site in flow_and_level_sites:
        # print site['closest_HydroID_distance'], site['Site Name'], site.get('Mean stage (m)')

    # hk = 1 * np.ones(shape)
    # reaches_by_stress = by_time(reaches, hk, model_start)




