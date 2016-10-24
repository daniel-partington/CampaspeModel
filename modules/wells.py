def load_simple(shape):
    nlay, nrow, ncol = shape
    # Remember to use zero-based layer, row, column indices!
    pumping_rate = -500.
    wel_sp1 = [[0, nrow/2 - 1, ncol/2 - 1, 0.]]
    wel_sp2 = [[0, nrow/2 - 1, ncol/2 - 1, 0.]]
    wel_sp3 = [[0, nrow/2 - 1, ncol/2 - 1, pumping_rate]]
    stress_period_data = {0: wel_sp1, 1: wel_sp2, 200: wel_sp3}

    return {"stress_period_data": stress_period_data}



def load(shape):
    nlay, nrow, ncol = shape
    # Remember to use zero-based layer, row, column indices!
    pumping_rate = -500.
    wel_sp1 = [[0, nrow/2 - 1, ncol/2 - 1, 0.]]
    wel_sp2 = [[0, nrow/2 - 1, ncol/2 - 1, 0.]]
    wel_sp3 = [[0, nrow/2 - 1, ncol/2 - 1, pumping_rate]]
    stress_period_data = {0: wel_sp1, 1: wel_sp2, 200: wel_sp3}

    return {"stress_period_data": stress_period_data}


# TODO
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))


import pandas
import os 
import json

from utils import cache
from GIS import project


def read_licences(clipped_and_projected_folder, model_projection, licenses_file):

    cols = ["Entitlement Volume",
            "Daily Volume",
            "Annual Volume",
            "Zone",
            "Easting",
            "Northing",
            "Bore depth (m)",
            "Top screen depth (m)",
            "Bottom screen depth (m)",
            "Use 2005/06",
            "Use 2006/07",
            "Use 2007/08",
            "Use 2008/09",
            "Use 2009/10",
            "Use 2010/11",
            "Use 2011/12",
            "Use 2012/13",
            "Use 2013/14",
            "Use 2014/15"]

    licenses = pandas.read_excel(licenses_file, skiprows=10, usecols=cols)
    licenses = licenses.fillna(value=0.0)    


    licenses['x'], licenses['y'] = project.project_points(xs=licenses['Easting'].values, 
                                                    ys=licenses['Northing'].values,
                                                    projector=project.geojson_transformer(model_projection, "EPSG:4326"))
    
    with open(os.path.join(clipped_and_projected_folder, 'UNPROJECTED_bores.json'), 'w') as f:
        f.write(json.dumps({
            "features": [{
                "geometry": {
                            # "coordinates": [site["Easting"], site["Northing"]],
                            "coordinates": [site["x"], site["y"]],
                            "type": "Point"
                },
                "properties": site,
                "type": "Feature"
            } for site in licenses.T.to_dict().values()],
            "type": "FeatureCollection"
        }))


    # with open('UNPROJECTED_bores.json', 'w') as f:
    #     f.write(json.dumps({
    #         "features": [{
    #             "geometry": {
    #                         # "coordinates": [site["Easting"], site["Northing"]],
    #                         "coordinates": [site["x"], site["y"]],
    #                         "type": "Point"
    #             },
    #             "properties": site,
    #             "type": "Feature"
    #         } for site in licenses.T.to_dict().values()],
    #         "type": "FeatureCollection"
    #     }))



if __name__ == '__main__':

    import config
    import boundary


    read_licences(
        clipped_and_projected_folder=config.clipped_and_projected_folder,
        model_projection=boundary.model_projection,
        licenses_file=os.path.join(config.data_folder, "GW", "Bore data", "Groundwater licence information for Dan Partington bc301115.xlsx")
    )

