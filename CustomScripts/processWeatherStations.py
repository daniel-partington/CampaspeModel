import pandas as pd
import matplotlib.pyplot as plt
#import matplotlib.colors as colors
import re
import os

import utils.cache

def add_alternative_format(path, other_stations, weather_dfs, weather_station_details):
    #  Lake Eppalock is in a different format so need something else to process
    for other_station in other_stations:
        other_station_name = other_station['name']
        with open(os.path.join(path, other_station['name'], other_station['meta_file']), 'r') as f:
            text = f.read()
            
        station_number = re.search('Bureau of Meteorology station number: (\S+)', text).group(1)

        # Get station Lat Long which corresponds to GDA94:

        station_lat = re.search('Latitude \(decimal degrees, south negative\): (\S+)', text).group(1)
        station_long = re.search('Longitude \(decimal degrees, east positive\): (\S+)', text).group(1)
        station_latlong = station_lat + " " + station_long

        # Get elevation of station:
        station_elev = re.search('Height of station above mean sea level \(metres\): (\S+)', text).group(1)

        weather_station_details[other_station_name] = [station_number, station_latlong , station_elev]

        #Read in time series data:
        weather_dfs[other_station_name] = pd.read_csv(os.path.join(path, other_station['name'], other_station['data_file']), 
                                               index_col=0, 
                                               skiprows=[41], 
                                               parse_dates=[[2,3,4]], #True, 
                                               #index_col = [[2,3,4]],
                                               infer_datetime_format=True, 
                                               comment='"', 
                                               skipinitialspace=True 
                                               )

        weather_dfs[other_station_name]['Rain'] = weather_dfs[other_station_name]['Rainfall amount (millimetres)']

    return 0


@utils.cache.money()
def processWeatherStations(weather_stations, other_stations, path='', frequency='A'):
    
    weather_station_details = {}
    weather_dfs = {}
    
    for station in weather_stations:
        station_file = os.path.join(path, station + '.txt')
        with open(station_file, 'r') as f:
            text = f.read()         
    
            # Get station number:
            station_number = re.search('Patched Point data for station: (\S+)', text).group(1)
    
            # Get station Lat Long which corresponds to GDA94:
            station_latlong = re.search('Lat: (\S+) Long: (\S+)', text).group().strip('"')
    
            # Get elevation of station:
            station_elev = re.search('Elevation:\s+(\w+)', text).group()
            
        weather_station_details[station] = [station_number, station_latlong , station_elev]
    
        #Read in time series data:
        weather_dfs[station] = pd.read_csv(station_file, 
                                           index_col=0, 
                                           skiprows=[41], 
                                           parse_dates=True, 
                                           infer_datetime_format=True, 
                                           delim_whitespace=True, 
                                           comment='"', 
                                           skipinitialspace=True, 
                                           usecols=[0,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])
 
    add_alternative_format(path, other_stations, weather_dfs, weather_station_details)

    if frequency == 'M':
        monthly_rain_df = pd.DataFrame() 
        for station in weather_stations + [d['name'] for d in other_stations]:
            monthly_rain_df[station]= weather_dfs[station]['Rain'].resample("M", how='mean')

        # monthly_rain_df[other_station]= weather_dfs[other_station]['Rainfall amount (millimetres)'].resample("M", how='mean')    
    
        return monthly_rain_df
    elif frequency == 'A':
        annual_rain_df = pd.DataFrame() 
        for station in weather_stations + [d['name'] for d in other_stations]:
            annual_rain_df[station]= weather_dfs[station]['Rain'].resample("A", how='sum')

        # annual_rain_df[other_station]= weather_dfs[other_station]['Rainfall amount (millimetres)'].resample("A", how='sum')    
        return annual_rain_df.mean()
    else:
        raise Exception("Invalid frequency %s. Should be M or A" % frequency)

    # Look at long term mean values for each
    
    # def cm2inch(*tupl):
    #     inch = 2.54
    #     if isinstance(tupl[0], tuple):
    #         return tuple(i/inch for i in tupl[0])
    #     else:
    #         return tuple(i/inch for i in tupl)
    
    #plt.figure(figsize=cm2inch(18,8))
    #plt.ylabel("Annual Rainfall [mm]")
    
    #for station in weather_stations:
        #weather_dfs[station]['Rain'].plot()        
        #weather_dfs[station]['Rain'].resample("M", how='sum').plot()    
    #    weather_dfs[station]['Rain'].resample("A", how='sum'). \
    #    plot(legend=True, 
    #         label=station + ', ' 
    #         + weather_station_details[station][0] + ', '  
    #         + weather_station_details[station][2] + ', Average: '  
    #         + str(weather_dfs[station]['Rain'].resample("A", how='sum').mean())[:5] + 'mm')
        
    #plt.xlabel("Year")
    #plt.legend(bbox_to_anchor=(0, 1), loc='upper left', ncol=1)
    

    
    #annual_rain_df.plot(kind='box')
    #plt.ylabel("Annual Rainfall [mm]")
    
    #Months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    #month_avg = pd.groupby(monthly_rain_df,by=[monthly_rain_df.index.month]).mean()
    #month_avg['Months'] = Months
    
    #month_avg.plot(kind='bar',x='Months',y=weather_stations)    
    
    #plt.ylabel('Average Monthly Rainfall [mm]')
    #plt.xlabel("")
    #plt.tight_layout()
    #plt.legend(bbox_to_anchor=(0, 1), loc='upper left', ncol=1)



    
if __name__ == "__main__":
    
    weather_stations = ['Kyneton',  'Elmore', 'Rochester', 'Echuca']
    weather = processWeatherStations(weather_stations, path=r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Climate\\", frequency='M')