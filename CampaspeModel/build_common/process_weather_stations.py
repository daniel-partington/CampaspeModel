import pandas as pd
import matplotlib.pyplot as plt
import os
#import matplotlib.colors as colors
import re

def process_weather_stations(weather_stations, path='', frequency='A', \
                           plot_monthly_pattern=False, \
                           plot_yearly_rainfall=False):
    
    assert frequency in ['A', 'M'], "Frequency must be either: 'A' or 'M'"
    
    weather_station_details = {}
    weather_dfs = {}
    
    for station in weather_stations:
        with open(os.path.join(path, station + '.txt'), 'r') as f:
            text = f.read()         
    
            # Get station number:
            station_number = re.search('Patched Point data for station: (\S+)', text).group(1)
    
            # Get station Lat Long which corresponds to GDA94:
            station_latlong = re.search('Lat: (\S+) Long: (\S+)', text).group().strip('"')
    
            # Get elevation of station:
            station_elev = re.search('Elevation:\s+(\w+)', text).group()
            
        weather_station_details[station] = [station_number, station_latlong , station_elev]
    
        #Read in time series data:
        weather_dfs[station] = pd.read_csv(os.path.join(path, station + '.txt'), 
                                           index_col=0, 
                                           skiprows=[41], 
                                           parse_dates=True, 
                                           infer_datetime_format=True, 
                                           delim_whitespace=True, 
                                           comment='"', 
                                           skipinitialspace=True, 
                                           usecols=[0,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])

    def cm2inch(*tupl):
        inch = 2.54
        if isinstance(tupl[0], tuple):
            return tuple(i/inch for i in tupl[0])
        else:
            return tuple(i/inch for i in tupl)
    

    def get_rain_and_ET_from_df(df, stations, freq, how='sum'):
        new_df = pd.DataFrame()
        for station in stations:
            if how == 'mean':
                new_df.loc[:, station] = df[station]['Rain'].resample(freq).mean()
                new_df.loc[:, station + '_ET'] = df[station]['Evap'].resample(freq).mean()
            elif how == 'sum':
                new_df.loc[:, station] = df[station]['Rain'].resample(freq).sum()
                new_df.loc[:, station + '_ET'] = df[station]['Evap'].resample(freq).sum()
            # end if
        #end for
        return new_df
        
    annual_weather_df = get_rain_and_ET_from_df(weather_dfs, weather_stations,
                                                'A', how='sum')
    monthly_weather_df = get_rain_and_ET_from_df(weather_dfs, weather_stations,
                                                'M', how='mean') 

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
        # Keeping this as is for now but should not calculate mean here 
        return annual_weather_df.mean()
    if frequency == 'M':
        return monthly_weather_df
    
if __name__ == "__main__":
    
    weather_stations = ['Kyneton', 'Eppalock', 'Elmore', 'Rochester', 'Echuca']
    weather = process_weather_stations(weather_stations, 
                                     path=r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Climate\\", 
                                     frequency='M',
                                     plot_monthly_pattern=True)

