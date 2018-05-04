import os
import re
import json
import numpy as np
import pandas as pd

from GIS import project
from GIS import interpolate


def mm_day_to_m_day(rainfall):
    return rainfall / 1000.0


def load(silo_stations, bom_stations, model_projection, mesh_centroids, start_date, end_date):
    # load in stations
    stations = [read_bom_weather(**s) for s in bom_stations]\
                 + [read_silo_weather(**s) for s in silo_stations]

    # project coordinates
    lng_lat_to_x_y = project.geojson_transformer("EPSG:4326", model_projection)

    for station in stations:
        lon = station['lon']
        lat = station['lat']

        station['x'], station['y'] = json.loads(lng_lat_to_x_y(json.dumps(
            {"type": "Point", "coordinates": [lon, lat]}
        )))['coordinates']
        
    points = np.array([[station['x'], station['y']] for station in stations])

    # avg for first steady state
    period_mean = np.array([station['timeseries']['Rain'].mean() for station in stations])
    interp_rain = interpolate.interpolate(points, period_mean, mesh_centroids)
    rech = {0: mm_day_to_m_day(interp_rain)}

    # change recharge at the start of each month to that month's mean
    
    # resample each stations data
    station_means = [station['timeseries']['Rain'][start_date: end_date].resample("M", label='left').mean() for station in stations] 
    for period_start in station_means[0].index:
        # TODO zero indexing?
        days_since_start_date = (period_start - start_date).days + 1 + 1
        # print("period_start", period_start, days_since_start_date)
        period_mean = np.array([station_mean[period_start] for station_mean in station_means])
        interp_rain = interpolate.interpolate(points, period_mean, mesh_centroids)

        # Adjust rainfall from mm/year to m/day
        rech[days_since_start_date] = mm_day_to_m_day(interp_rain)
        print(period_start, np.nanmean(rech[days_since_start_date]))

    return {
        "rech": rech
    }


"""
BOM 
    * stations.txt contains a list of BOM observation sites
        - search for sites [manually]( http://www.bom.gov.au/climate/data ) or
        - (download list)[ftp://ftp.bom.gov.au/anon2/home/ncc/metadata/sitelists/stations.zip] -> stations.txt
    * climate observations at a site eg. http://www.bom.gov.au/jsp/ncc/cdio/weatherData/av?p_nccObsCode=136&p_display_type=dailyDataFile&p_stn_num=080015&p_startYear=
    * climate averages eg. http://www.bom.gov.au/clim_data/cdio/tables/text/IDCJCM0037_066196.csv
"""


def read_bom_weather(station_name, meta_file, data_file):
    with open(meta_file, 'r') as f:
        text = f.read()

    station_number = re.search('Bureau of Meteorology station number: (\S+)', text).group(1)

    station_lat = re.search('Latitude \(decimal degrees, south negative\): (\S+)', text).group(1)
    station_long = re.search('Longitude \(decimal degrees, east positive\): (\S+)', text).group(1)
    
    station_elev = re.search('Height of station above mean sea level \(metres\): (\S+)', text).group(1)

    # Read in time series data:
    df = pd.read_csv(data_file, 
                   index_col=0,
                   parse_dates=[[2,3,4]],
                   infer_datetime_format=True,
                   comment='"', 
                   skipinitialspace=True)

    df['Rain'] = df['Rainfall amount (millimetres)']

    return {
        "name": station_name,
        "timeseries": df,
        "number": station_number,
        "lat": float(station_lat),
        "lon": float(station_long),
        "elev": float(station_elev)
    }


"""
SILO
* SILO contains Australian climate data from 1889 (current to yesterday)
* https://www.longpaddock.qld.gov.au/silo/
"""

def read_silo_weather(station_name, station_file):
    with open(station_file, 'r') as f:
        text = f.read()

    station_number = re.search('Patched Point data for station: (\S+)', text).group(1)

    station_lat, station_long = re.search('Lat: ([-+]?\d*\.\d+|\d+)\s+Long: ([-+]?\d*\.\d+|\d+)', text).groups()
    station_elev = re.search('Elevation:\s+([-+]?\d*\.\d+|\d+)m', text).group(1)
        
    df = pd.read_csv(station_file, 
                   index_col=0, 
                   skiprows=[41], 
                   parse_dates=True, 
                   infer_datetime_format=True, 
                   delim_whitespace=True, 
                   comment='"', 
                   skipinitialspace=True, 
                   usecols=[0,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])

    return {
        "name": station_name,
        "timeseries": df,
        "number": station_number,
        "lat": float(station_lat),
        "lon": float(station_long),
        "elev": float(station_elev)
    }


if __name__ == '__main__':
    import datetime
    import config
    import boundary
    import matplotlib.pylab as plt

    # print read_silo_weather(config.data_folder, silo_stations[0])['timeseries']['Rain']
    # print read_bom_weather(config.data_folder, bom_stations[0])['timeseries']['Rain']

    shape = (1, 30, 20)  # (nlay, nrow, ncol)
    start_date = datetime.datetime(2014, 1, 1)
    end_date = datetime.datetime(2014, 12, 31)
    boundary_params = boundary.boundary(boundary_file=config.boundary_file,
                                    clipped_and_projected_folder=config.clipped_and_projected_folder,
                                    shape=shape)

    rech = load(bom_stations=config.bom_stations,
                            silo_stations=config.silo_stations, 
                            model_projection=boundary.model_projection, 
                            mesh_centroids=boundary_params['mesh_centroids'],
                            start_date=start_date, end_date=end_date)
    for rain in rech['rech'].values()[:3]:
        plt.imshow(rain)
        plt.colorbar()
        plt.savefig('/shared/plot/rech.png')
        plt.close()


