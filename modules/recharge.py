import os
import re
import json
import datetime
import numpy as np
from dateutil.relativedelta import relativedelta
import pandas as pd

from utils import cache

from GIS import project
from GIS import interpolate

silo_stations = ['Kyneton',  'Elmore', 'Rochester', 'Echuca']

bom_stations = [{
        "name": "Rainfall at Eppalock",
        "meta_file": 'IDCJAC0009_081083_1800_Note.txt',
        "data_file": 'IDCJAC0009_081083_1800_Data.csv',
}]

start_date = datetime.datetime(2012, 1, 1)
end_date = datetime.datetime(2012, 12, 31)


@cache.money()
def load(data_folder, model_projection, mesh_centroids):
    # load in stations
    stations = [read_bom_weather(data_folder, s) for s in bom_stations]\
                 + [read_silo_weather(data_folder, s) for s in silo_stations]

    # project coordinates
    lng_lat_to_x_y = project.geojson_transformer("EPSG:4326", model_projection)

    for station in stations:
        lon = station['lon']
        lat = station['lat']

        station['x'], station['y'] = json.loads(lng_lat_to_x_y(json.dumps(
            {"type": "Point", "coordinates": [lon, lat]}
        )))['coordinates']
        
    points = np.array([[station['x'], station['y']] for station in stations])

    # change recharge at the start of each month to that month's mean
    rech = {}
    
    for month in range(12):
        first_of_month = start_date + relativedelta(months=month)
        last_of_month = start_date + relativedelta(months=month + 1, days=-1)
        
        monthly_mean = np.array([station['timeseries']['Rain'][first_of_month: last_of_month].mean()
                                for station in stations])
        
        day_of_year = first_of_month.timetuple().tm_yday
        
        interp_rain = interpolate.interpolate(points, monthly_mean, mesh_centroids)
        
        # Adjust rainfall to m from mm and from year to days
        rech[day_of_year] = interp_rain / 1000.0 / 365.0

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


def read_bom_weather(data_folder, bom_station):
    with open(os.path.join(data_folder, 'Climate', bom_station['name'], bom_station['meta_file']), 'r') as f:
        text = f.read()

    station_number = re.search('Bureau of Meteorology station number: (\S+)', text).group(1)

    station_lat = re.search('Latitude \(decimal degrees, south negative\): (\S+)', text).group(1)
    station_long = re.search('Longitude \(decimal degrees, east positive\): (\S+)', text).group(1)
    
    station_elev = re.search('Height of station above mean sea level \(metres\): (\S+)', text).group(1)

    # Read in time series data:
    df = pd.read_csv(os.path.join(data_folder, 'Climate', bom_station['name'], bom_station['data_file']), 
                   index_col=0,
                   parse_dates=[[2,3,4]],
                   infer_datetime_format=True,
                   comment='"', 
                   skipinitialspace=True)

    df['Rain'] = df['Rainfall amount (millimetres)']

    return {
        "name": bom_station['name'],
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

def read_silo_weather(data_folder, station):
    station_file = os.path.join(data_folder, 'Climate', station + '.txt')
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
        "name": station,
        "timeseries": df,
        "number": station_number,
        "lat": float(station_lat),
        "lon": float(station_long),
        "elev": float(station_elev)
    }


if __name__ == '__main__':
    
    import config
    import boundary
    import matplotlib.pylab as plt

    # print read_silo_weather(config.data_folder, silo_stations[0])['timeseries']['Rain']
    # print read_bom_weather(config.data_folder, bom_stations[0])['timeseries']['Rain']

    shape = (1, 30, 20)  # (nlay, nrow, ncol)
    boundary_params, boundary_inputs = boundary.load(boundary_file=os.path.join(config.data_folder, "GW_model_area", "GW_model_area.shp"), 
                                    clipped_and_projected_folder=config.clipped_and_projected_folder,
                                    shape=shape)

    rech = load(config.data_folder, boundary.model_projection, boundary_params['mesh_centroids'])
    for rain in rech['rech'].values()[:3]:
        plt.imshow(rain)
        plt.colorbar()
        plt.show()


