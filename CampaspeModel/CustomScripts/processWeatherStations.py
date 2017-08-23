<<<<<<< HEAD
import pandas as pd
import matplotlib.pyplot as plt
#import matplotlib.colors as colors
import re

def processWeatherStations(weather_stations, path='', frequency='A', \
                           plot_monthly_pattern=False, \
                           plot_yearly_rainfall=False):
    
    weather_station_details = {}
    weather_dfs = {}
    
    for station in weather_stations:
        with open(path + station + '.txt', 'r') as f:
            text = f.read()         
    
            # Get station number:
            station_number = re.search('Patched Point data for station: (\S+)', text).group(1)
    
            # Get station Lat Long which corresponds to GDA94:
            station_latlong = re.search('Lat: (\S+) Long: (\S+)', text).group().strip('"')
    
            # Get elevation of station:
            station_elev = re.search('Elevation:\s+(\w+)', text).group()
            
        weather_station_details[station] = [station_number, station_latlong , station_elev]
    
        #Read in time series data:
        weather_dfs[station] = pd.read_csv(path + station + '.txt', 
                                           index_col=0, 
                                           skiprows=[41], 
                                           parse_dates=True, 
                                           infer_datetime_format=True, 
                                           delim_whitespace=True, 
                                           comment='"', 
                                           skipinitialspace=True, 
                                           usecols=[0,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@ This little chuck below was used for processing BOM data for Eppalock  @@
# @@@ but it has been updated with SILO data ... going forward would be good @@
# @@@ to make general SILO and BOM functions for processing and import them  @@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
#    #  Lake Eppalock is in a different format so need something else to process
#    other_station = 'Eppalock Reservoir'
#    other_path = path + r"Rainfall at Eppalock\\"
#    with open(other_path + 'IDCJAC0009_081083_1800_Note.txt', 'r') as f:
#        text = f.read()
#        
#        station_number = re.search('Bureau of Meteorology station number: (\S+)', text).group(1)
#    
#        # Get station Lat Long which corresponds to GDA94:
#     
#        station_lat = re.search('Latitude \(decimal degrees, south negative\): (\S+)', text).group(1)
#        station_long = re.search('Longitude \(decimal degrees, east positive\): (\S+)', text).group(1)
#        station_latlong = station_lat + " " + station_long
#    
#        # Get elevation of station:
#        station_elev = re.search('Height of station above mean sea level \(metres\): (\S+)', text).group(1)
#    
#        weather_station_details[other_station] = [station_number, station_latlong , station_elev]
#        
#        #Read in time series data:
#        weather_dfs[other_station] = pd.read_csv(other_path + 'IDCJAC0009_081083_1800_Data.csv', 
#                                               index_col=0, 
#                                               skiprows=[41], 
#                                               parse_dates=[[2,3,4]], #True, 
#                                               #index_col = [[2,3,4]],
#                                               infer_datetime_format=True, 
#                                               comment='"', 
#                                               skipinitialspace=True 
#                                               )

    # Get info:
    
    # Look at long term mean values for each
    
    def cm2inch(*tupl):
        inch = 2.54
        if isinstance(tupl[0], tuple):
            return tuple(i/inch for i in tupl[0])
        else:
            return tuple(i/inch for i in tupl)
    
   
    annual_weather_df = pd.DataFrame()
    for station in weather_stations:
        annual_weather_df[station]= weather_dfs[station]['Rain'].resample("A").sum()
        annual_weather_df[station + '_ET']= weather_dfs[station]['Evap'].resample("A").sum()

#    annual_weather_df[other_station]= weather_dfs[other_station]['Rainfall amount (millimetres)'].resample("A").sum()    
    
    monthly_weather_df = pd.DataFrame() 
    for station in weather_stations:
        monthly_weather_df[station]= weather_dfs[station]['Rain'].resample("M").mean()
        monthly_weather_df[station + '_ET']= weather_dfs[station]['Evap'].resample("M").mean()

#    monthly_weather_df[other_station]= weather_dfs[other_station]['Rainfall amount (millimetres)'].resample("M").mean()    
    
    if plot_yearly_rainfall:
        plt.figure(figsize=cm2inch(18,8))
        plt.ylabel("Annual Rainfall [mm]")
        
        for station in weather_stations:
            weather_dfs[station]['Rain'].plot()        
            weather_dfs[station]['Rain'].resample("M", how='sum').plot()    
            weather_dfs[station]['Rain'].resample("A", how='sum'). \
            plot(legend=True, 
                 label=station + ', ' 
                 + weather_station_details[station][0] + ', '  
                 + weather_station_details[station][2] + ', Average: '  
                 + str(weather_dfs[station]['Rain'].resample("A", how='sum').mean())[:5] + 'mm')
            
        plt.xlabel("Year")
        plt.legend(bbox_to_anchor=(0, 1), loc='upper left', ncol=1)

        annual_weather_df.plot(kind='box')
        plt.ylabel("Annual Rainfall [mm]")

    if plot_monthly_pattern:
        Months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        month_avg = pd.groupby(monthly_weather_df,by=[monthly_weather_df.index.month]).mean()
        month_avg['Months'] = Months
        
        month_avg.plot(kind='bar',x='Months',y=weather_stations)    
        
        plt.ylabel('Average Monthly Rainfall [mm]')
        plt.xlabel("")
        plt.tight_layout()
        plt.legend(bbox_to_anchor=(0, 1), loc='upper left', ncol=1)

    if frequency == 'A':
        return annual_weather_df.mean()
    if frequency == 'M':
        return monthly_weather_df
    
if __name__ == "__main__":
    
    weather_stations = ['Kyneton', 'Eppalock', 'Elmore', 'Rochester', 'Echuca']
    weather = processWeatherStations(weather_stations, 
                                     path=r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Climate\\", 
                                     frequency='M',
                                     plot_monthly_pattern=True)
=======
#import matplotlib.colors as colors
import os
import re

import matplotlib.pyplot as plt
import pandas as pd


def processWeatherStations(weather_stations, path='', frequency='A'):

    weather_station_details = {}
    weather_dfs = {}

    p_j = os.path.join
    for station in weather_stations:
        with open(p_j(path, station + ".txt"), 'r') as f:
            text = f.read()

            # Get station number:
            station_number = re.search('Patched Point data for station: (\S+)', text).group(1)

            # Get station Lat Long which corresponds to GDA94:
            station_latlong = re.search('Lat: (\S+) Long: (\S+)', text).group().strip('"')

            # Get elevation of station:
            station_elev = re.search('Elevation:\s+(\w+)', text).group()

        weather_station_details[station] = [station_number, station_latlong, station_elev]

        # Read in time series data:
        weather_dfs[station] = pd.read_csv(os.path.join(path, station + ".txt"),
                                           index_col=0,
                                           skiprows=[41],
                                           parse_dates=True,
                                           infer_datetime_format=True,
                                           delim_whitespace=True,
                                           comment='"',
                                           skipinitialspace=True,
                                           usecols=[0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17])

    #  Lake Eppalock is in a different format so need something else to process
    other_station = 'Eppalock Reservoir'
    other_path = p_j(path, r"Rainfall at Eppalock\\")

    with open(p_j(other_path, "IDCJAC0009_081083_1800_Note.txt"), 'r') as f:
        text = f.read()

        station_number = re.search('Bureau of Meteorology station number: (\S+)', text).group(1)

        # Get station Lat Long which corresponds to GDA94:

        station_lat = re.search(
            'Latitude \(decimal degrees, south negative\): (\S+)', text).group(1)
        station_long = re.search(
            'Longitude \(decimal degrees, east positive\): (\S+)', text).group(1)
        station_latlong = station_lat + " " + station_long

        # Get elevation of station:
        station_elev = re.search(
            'Height of station above mean sea level \(metres\): (\S+)', text).group(1)

        weather_station_details[other_station] = [station_number, station_latlong, station_elev]

        # Read in time series data:
        weather_dfs[other_station] = pd.read_csv(other_path + 'IDCJAC0009_081083_1800_Data.csv',
                                                 index_col=0,
                                                 skiprows=[41],
                                                 parse_dates=[[2, 3, 4]],  # True,
                                                 #index_col = [[2,3,4]],
                                                 infer_datetime_format=True,
                                                 comment='"',
                                                 skipinitialspace=True
                                                 )

    # Get info:

    # Look at long term mean values for each

    def cm2inch(*tupl):
        inch = 2.54
        if isinstance(tupl[0], tuple):
            return tuple(i / inch for i in tupl[0])
        else:
            return tuple(i / inch for i in tupl)

    # plt.figure(figsize=cm2inch(18,8))
    #plt.ylabel("Annual Rainfall [mm]")

    # for station in weather_stations:
        # weather_dfs[station]['Rain'].plot()
        #weather_dfs[station]['Rain'].resample("M", how='sum').plot()
    #    weather_dfs[station]['Rain'].resample("A", how='sum'). \
    #    plot(legend=True,
    #         label=station + ', '
    #         + weather_station_details[station][0] + ', '
    #         + weather_station_details[station][2] + ', Average: '
    #         + str(weather_dfs[station]['Rain'].resample("A", how='sum').mean())[:5] + 'mm')

    # plt.xlabel("Year")
    #plt.legend(bbox_to_anchor=(0, 1), loc='upper left', ncol=1)

    annual_rain_df = pd.DataFrame()
    for station in weather_stations:
        annual_rain_df[station] = weather_dfs[station]['Rain'].resample("A").sum()

    annual_rain_df[other_station] = weather_dfs[other_station][
        'Rainfall amount (millimetres)'].resample("A").sum()

    # annual_rain_df.plot(kind='box')
    #plt.ylabel("Annual Rainfall [mm]")

    monthly_rain_df = pd.DataFrame()
    for station in weather_stations:
        monthly_rain_df[station] = weather_dfs[station]['Rain'].resample("M").mean()

    monthly_rain_df[other_station] = weather_dfs[other_station][
        'Rainfall amount (millimetres)'].resample("M").mean()

    #Months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    #month_avg = pd.groupby(monthly_rain_df,by=[monthly_rain_df.index.month]).mean()
    #month_avg['Months'] = Months

    # month_avg.plot(kind='bar',x='Months',y=weather_stations)

    #plt.ylabel('Average Monthly Rainfall [mm]')
    # plt.xlabel("")
    # plt.tight_layout()
    #plt.legend(bbox_to_anchor=(0, 1), loc='upper left', ncol=1)

    if frequency == 'A':
        return annual_rain_df.mean()
    if frequency == 'M':
        return monthly_rain_df

if __name__ == "__main__":

    weather_stations = ['Kyneton',  'Elmore', 'Rochester', 'Echuca']
    weather = processWeatherStations(
        weather_stations, path=r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Climate\\", frequency='M')
>>>>>>> GW_int_dev
