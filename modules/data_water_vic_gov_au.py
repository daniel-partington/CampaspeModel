import os
import json
import datetime
from osgeo import ogr
import pandas as pd

from GIS import io
from GIS import project
import utils.cache


"""
data.water.vic
    * river gauges
    * download from http://data.water.vic.gov.au/monitoring.htm
"""

def save_doc(sites, file_name):
    with open(file_name, 'w') as f:
        f.write(utils.cache.json.dumps({
            "type": "FeatureCollection",
            "features": [
                {
                    "type": "Feature",
                    "properties": site,
                    "geometry": {"type": "Point", "coordinates": [site['Longitude'], site['Latitude']]}
                }
                for site in sites
            ]
        }, cls=utils.cache.DateNumpyEncoder, ignore_nan=True))


def read_site_file(site_file_name, site_height, model_start, model_end):

    if not os.path.exists(site_file_name):
        print "doesn't exist", site_file_name
        return None

    site_df = pd.read_csv(site_file_name, skiprows=2, index_col='Date', parse_dates=True,
                                    date_parser=lambda d: datetime.datetime.strptime(d, '00:00:00 %d/%m/%Y'))

    site_df = site_df[model_start: model_end]
    site_df = site_df.fillna(0)        

    if site_df['Mean'].max() == 0.0 or pd.isnull(site_df['Mean'].max()):
        print "No data", site_file_name
        return None

    # more QC for MeanWaterLevel files only

    if 'Level' in site_file_name and set(site_df['Qual'].tolist()).intersection({8, 9}):
        print 'Quality code indicates poor or incorrect reading for: ', site_file_name        
        return None

    # hack if elevation is not AHD
    if 'Level' in site_file_name and site_df['Mean'].mean() < 10.0:
        if float(site_height) > 0:
            site_df['Mean'] += float(site_height)

    return site_df


def aggregate(ts, start_date):
    return sorted([{
        "date": date.to_datetime(), 
        "days_since_start_date": (date.to_datetime().date() - start_date).days + 1 + 1,
        "value": value
    } for date, value in 
    ts.resample("M", label='left').mean().iteritems()], key=lambda d: d['date'])


def read_relevant_sites(sites_folder, clipped_and_projected_folder, model_projection, projected_boundary_file, model_start, model_end):

    sites = pd.read_csv(os.path.join(sites_folder, "Site Details.csv"), 
        parse_dates=["Site Start Date", "Site Cease Date"])

    # project coordinates
    lng_lat_to_x_y = project.geojson_transformer("EPSG:4326", model_projection)

    sites['x'], sites['y'] = project.project_points(xs=sites['Longitude'].values, ys=sites['Latitude'].values,
                                                    projector=lng_lat_to_x_y)

    # filter by boundary
    model_boundary = io.load_geom(file_name=projected_boundary_file)
    for index, site in sites.iterrows():
        # within 5 km
        sites.loc[index, 'in_model_boundary'] = 5e3 > model_boundary.Distance(ogr.CreateGeometryFromJson(json.dumps({"coordinates": [site['x'], site['y']], "type": "Point"})))
    sites = sites[sites['in_model_boundary']]

    # filter out sites without incomplete data
    # not trustworthy, so just filter in read_site_file()
    # sites = sites[(sites["Site Start Date"] <= model_start.strftime('%Y-%m-%d')) & (sites["Site Cease Date"].isnull() | (model_end.strftime('%Y-%m-%d') <= sites["Site Cease Date"]))]

    # pandas -> list of dicts (more flexible)
    sites["Site Cease Date"] = sites["Site Cease Date"].astype(str)
    sites["Site Start Date"] = sites["Site Start Date"].astype(str)
    sites = sites.T.to_dict().values()

    save_doc(sites, file_name=os.path.join(clipped_and_projected_folder, 'UNPROJECTED_unfiltered_data.water.vic.gov.au.json'))

    for index, site in enumerate(sites):

        # load in MeanWaterFlow and MeanWaterLevel
        flow_file_name = os.path.join(sites_folder, "%d.%s.csv" % (site["Site Id"], "MeanWaterFlow"))
        flow_df = read_site_file(site_file_name=flow_file_name, site_height=site['Gauge Zero (Ahd)'], model_start=model_start, model_end=model_end)
        if flow_df is not None:
            flow_df['Mean flow (m3/s)'] = flow_df['Mean'] * 1000. / 86400  # 24 * 60 * 60
            
            site['flow_site'] = True
            site['flow file name'] = flow_file_name
            site['Mean flow (m3/s)'] = aggregate(flow_df['Mean flow (m3/s)'], start_date=model_start)

        level_file_name = os.path.join(sites_folder, "%d.%s.csv" % (site["Site Id"], "MeanWaterLevel"))
        level_df = read_site_file(site_file_name=level_file_name, site_height=site['Gauge Zero (Ahd)'], model_start=model_start, model_end=model_end)
        if level_df is not None:
            level_df['Mean stage (m)'] = level_df['Mean']
            
            site['level_site'] = True
            site['level file name'] = level_file_name
            site['Mean stage (m)'] = aggregate(level_df['Mean stage (m)'], start_date=model_start)

    save_doc(filter(lambda s: s.get('flow_site'), sites), file_name=os.path.join(clipped_and_projected_folder, 'UNPROJECTED_flow_data.water.vic.gov.au.json'))
    save_doc(filter(lambda s: s.get('level_site'), sites), file_name=os.path.join(clipped_and_projected_folder, 'UNPROJECTED_level_data.water.vic.gov.au.json'))

    return filter(lambda s: s.get('flow_site') or s.get('level_site'), sites)


if __name__ == '__main__':
    import config
    import boundary

    shape = (1, 30, 20)  # (nlay, nrow, ncol)
    boundary_params = boundary.boundary(data_folder=config.data_folder, 
                                    clipped_and_projected_folder=config.clipped_and_projected_folder,
                                    shape=shape)

    model_start = datetime.date(2015, 1, 1)
    model_end = datetime.date(2015, 12, 31)

    sites = read_relevant_sites(sites_folder=os.path.join(config.data_folder, "data.water.vic", "download.20160923190321"), 
                        clipped_and_projected_folder=config.clipped_and_projected_folder,
                        model_projection=boundary_params["model_projection"], 
                        projected_boundary_file=boundary_params['projected_boundary_file'],
                        model_start=model_start, model_end=model_end)

    print flow_sites, level_sites
