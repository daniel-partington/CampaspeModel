import os
import pandas as pd
import re

def processRiverGauges(cwd=None):
    if cwd == None:    
        cwd = os.getcwd() + '\\'
    
    listdir = os.listdir(cwd)
    
    dfs = {}
    
    # Get site info and data and place in dictionary with site no and name as key
    for item in listdir:
        print item        
        if not os.path.isdir(cwd+item):
            continue
        path = cwd + item + '\\'
        files = os.listdir(path)
        for all_file in files:
            if '.csv' in all_file:
                with open(path + all_file, 'r') as csvfile:            
                    skip_lines = 0
                    lines = csvfile.readlines()
                    for line in lines:
                        if 'Datetime' in line:
                            break
                        skip_lines += 1
    
                with open(path + all_file, 'r') as csvfile:
                    text = csvfile.read()
                    #quality_codes = re.search('', raw_data)                        
                    site_info = re.search('Site.*Elev:\S+', text).group()
                    site_info = {'site_no':re.search('Site (\S+)', text).group(1),                
                                 'site_lat':re.search('Lat:(\S+)', text).group(1),
                                 'site_long':re.search('Long:(\S+)', text).group(1),
                                 'site_elev':re.search('Elev:(\d+)', text).group(1)}                 
                    site_info['site_name'] = re.search(site_info['site_no']+'(.*)Lat:', text).group(1),
                    #  site_data = io.StringIO(re.search('"Datetime.*', raw_data).group())
                            #site_data = site_data.write(re.search('"Datetime.*', raw_data).group())  
                            #site_data = '"Datetime' + raw_data.split('Datetime')[1]                        
                dfs[item] = [site_info, pd.io.parsers.read_csv(open(path + all_file, 'r'), index_col=0, infer_datetime_format=True, parse_dates=True, skiprows=skip_lines, dayfirst=True)]                        
    
    max_elevation = 0
    min_elevation = 150
    for item in dfs:
        if dfs[item][0]['site_elev'] > max_elevation:
            max_elevation = dfs[item][0]['site_elev']
        if dfs[item][0]['site_elev'] < min_elevation:
            min_elevation = dfs[item][0]['site_elev']
        
    return dfs        
        
if __name__ == "__main__":
    
    river_gauge_folder = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\20_gauges\\"
    riv_data = processRiverGauges(cwd=river_gauge_folder)