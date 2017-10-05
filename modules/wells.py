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


import pandas
import os
import collections
import json
import datetime
import numpy as np

from GIS import project
from GIS import interpolate


def day_range(start_date, end_date):
    for n in range((end_date - start_date).days + 1):
        yield start_date + datetime.timedelta(n)


def month_range(start_date, end_date):
    return [day for i, day in enumerate(day_range(start_date, end_date)) if i == 0 or day.day == 1]


def combined_licences(period_start, licenses):
    year = period_start.year
    if period_start.month < 8:
        year = year - 1
    # 'start_date': datetime.date(year, 8, 1),
    #  'end_date': datetime.date(year, 7, 30),
    cells = collections.defaultdict(float)
    for l in licenses:
        key = 'Use %4d/%02d' % (year, (year + 1) % 1000)
        # ML / year to m3 / day
        val = -1. * 1000. / 365. * l[key]
        cells[(l['layer'], l['row'], l['col'])] += val

    return list([
        list(k + (v,)) for k, v in cells.items()
    ])


def avg_licences(licenses_rows):
    cells = collections.defaultdict(float)
    years = list(range(2005, 2014))
    for l in licenses_rows:
        for year in years:
            key = 'Use %4d/%02d' % (year, (year + 1) % 1000)
            # ML / year to m3 / day
            val = -1. * 1000. / 365. * l[key] / len(years)
            cells[(l['layer'], l['row'], l['col'])] += val

    return list([
        list(k + (v,)) for k, v in cells.items()
    ])


def remove_outside(licenses, bbox):
    min_x, max_x, min_y, max_y = bbox
    idxs =  (min_x <= licenses['Easting'])  &  (licenses['Easting'] <= max_x) & (min_y <= licenses['Northing'])  &  (licenses['Northing'] <= max_y)
    return licenses.drop(licenses.index[~idxs])


def load(clipped_and_projected_folder, licenses_file, model_projection, bbox, dels, start_date, end_date):
    licenses = read_licences(clipped_and_projected_folder, model_projection, licenses_file)

    licenses = remove_outside(licenses, bbox)

    min_x, max_x, min_y, max_y = bbox
    delr, delc = dels

    licenses['row'] = (licenses['Northing'] - min_y) / delr
    licenses['row'] = licenses['row'].astype(np.int)
    licenses['col'] = (licenses['Easting'] - min_x) / delc
    licenses['col'] = licenses['col'].astype(np.int)
    licenses['layer'] = 0 

    licenses_rows = licenses.T.to_dict().values()

    stress_period_data = {
        0: avg_licences(licenses_rows)
    }
    for period_start in month_range(start_date, end_date):
        days_since_start_date = (period_start - start_date).days + 1
        # print("period_start", period_start, days_since_start_date)
        stress_period_data[days_since_start_date] = combined_licences(period_start, licenses_rows)
        print(period_start, sum([s[3] for s in stress_period_data[days_since_start_date]]))

    return {"stress_period_data": stress_period_data}


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

    licenses['x'], licenses['y'] = project.project_points(
        xs=licenses['Easting'].values, 
        ys=licenses['Northing'].values,
        projector=project.geojson_transformer(model_projection, "EPSG:4326")
    )

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

    return licenses


if __name__ == '__main__':

    import config
    import boundary

    # licences = read_licences(
    #     clipped_and_projected_folder=config.clipped_and_projected_folder,
    #     model_projection=boundary.model_projection,
    #     licenses_file=config.licenses_file
    # )

    start_date = datetime.datetime(2014, 1, 1)
    end_date = datetime.datetime(2014, 12, 31)

    shape = (1, 30, 20)  # (nlay, nrow, ncol)
    boundary_params = boundary.boundary(boundary_file=config.boundary_file, 
                                clipped_and_projected_folder=config.clipped_and_projected_folder,
                                shape=shape)

    sp = load(clipped_and_projected_folder=config.clipped_and_projected_folder, 
        licenses_file=config.licenses_file, 
        model_projection=boundary.model_projection, 
        bbox=boundary_params['bbox'], dels=boundary_params['dels'], 
        start_date=start_date, end_date=end_date)
    print(sp)