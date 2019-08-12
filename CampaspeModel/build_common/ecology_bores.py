
import pandas as pd
import xlrd
import matplotlib.pyplot as plt
import os, sys
import re
from osgeo import ogr, osr  # to get this install the GDAL python library

###############################################################################
#
# This is a script to: 
#   1. read in stream gauge files and extract their details
#   2. read in bore data and find closest bores to the stream gauges
#   3. filter bores by those having salinity and level data
#   4. resmaple bore data at daily time step and output to csv file
#
###############################################################################


# Location for output files:
filepath = None  # Change if you want to specify where to place the csv files

# Set stream gauges of interest:
stream_gauges_of_interest = ['406201', '406202', '406207', '406218', '406265']

# Get info for stream gauges:
stream_gauge_info_dir = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\20_gauges\\"

listdir = os.listdir(stream_gauge_info_dir)

dfs = {}

# Get site info and data and place in dictionary with site no and name as key
for item in listdir:
    if not os.path.isdir(stream_gauge_info_dir+item):
        continue
    path = stream_gauge_info_dir + item + '\\'
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


# We need this function to get from spherical to cartesion coordinates to aid
# the calculation of distance between the stream gauge and bore locations

def coordinate_transform(pointX, pointY, coordTransform):
    """
    A function to transform a point on the horizontal plane
    using CoordinateTransformation object from GDAL osr.
    """
    if type(pointX) != float:
        (pointX, pointY) = (float(pointX), float(pointY))
    # create a geometry from coordinates
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(pointX, pointY)
    # transform point
    point.Transform(coordTransform)

    return point.GetX(), point.GetY()

# EPSG is a code system for coordinate systems 

epsg_from = 4283 # http://spatialreference.org/ref/epsg/4283/
epsg_to = 28355 # http://spatialreference.org/ref/epsg/gda94-mga-zone-55/

# input SpatialReference
inSpatialRef = osr.SpatialReference()
inSpatialRef.ImportFromEPSG(epsg_from)

# output SpatialReference
outSpatialRef = osr.SpatialReference()
outSpatialRef.ImportFromEPSG(epsg_to)

# create the CoordinateTransformation
coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)

relevant_stream_gauges = {}

# Do the transformations to get easting and northings
for key in dfs:
    if any(x in key for x in stream_gauges_of_interest):     
        x, y = coordinate_transform(dfs[key][0]['site_long'], dfs[key][0]['site_lat'], coordTrans)
        relevant_stream_gauges[dfs[key][0]['site_no']] = (x, y)

#colors = list("rgbcmyk")
#for gauge in relevant_stream_gauges.values():
#    x, y = gauge[0], gauge[1]
#    plt.scatter(x,y,color=colors.pop())
#plt.legend(relevant_stream_gauges.keys())

# Now load in the bore data to pandas dataframes

fname = r"Shallow monitoring bores bc301115.xlsx"
fname2 = r"State Observation Bores bc271115.xlsx"

with xlrd.open_workbook(fname, on_demand=True) as xls:
    sheets = xls.sheet_names()

with xlrd.open_workbook(fname2, on_demand=True) as xls:
    sheets2 = xls.sheet_names()

df_set = {}
df_set2 = {}

for sheet in sheets:
    df_set[sheet] = pd.read_excel(fname, sheetname=sheet, parse_dates=True)

for sheet in sheets2:
    df_set2[sheet] = pd.read_excel(fname2, sheetname=sheet, parse_dates=True)

for key in df_set:
    df_set[key] = pd.concat([df_set[key], df_set2[key]])


# Filter construction details to get screens:
# Index by bore ID for the filter construction drails
Screen_info = df_set['Bore Construction'].loc[df_set['Bore Construction']['Component'] == 'Screen']

# Filter lab chem sheet to get all of the salinities at bores
WaterLevel = df_set['Water Levels']
water_level_bores = pd.unique(WaterLevel['Bore ID'])

# Filter lab chem sheet to get all of the salinities at bores
SalinityEC = df_set['Lab Chem'].loc[df_set['Lab Chem']['Parameter name'] == 'EC (field) @ S/T, microS/cm']
SalinityEC2 = df_set['Field Chem'].loc[df_set['Field Chem']['EC (uS/cm)'] != 0]

# Get unique values of bores with salinity data
salinity_bores = pd.unique(SalinityEC['Bore ID'])
salinity_bores2 = pd.unique(SalinityEC2['Bore ID'])
salinity_bores = list(set(salinity_bores) | set(salinity_bores2))

# Define converter for EC to TDS ppm
def EC2ppm(ECvalue, conversion_factor = 0.7):
    # Assuming EC is in microsiemens per cm
    return ECvalue*conversion_factor

water_level_and_salinity = list(set(water_level_bores) & set(salinity_bores))

# Get the distance between bores and stream gauges
#for bore in water_level_and_salinity:

relevant_bores = df_set['Site Details'].loc[df_set['Site Details']['Bore ID'].isin(water_level_and_salinity)]

def points_dist(x1,y1,x2,y2):
    return ((x2-x1)**2+(y2-y1)**2)**0.5
    
gauges = {}    
gauge_x = []
gauge_y = []
bore_x = []
bore_y = []
bore_linking = {}
stream_linking = {}

for gauge in relevant_stream_gauges:
    gauge_loc = relevant_stream_gauges[gauge]    
    gauge_x.append(gauge_loc[0])
    gauge_y.append(gauge_loc[1])
    gauges[gauge] = [] 
    min_dist = 99999    
    for row in relevant_bores.iterrows():
        dist = points_dist(row[1]['Easting'], row[1]['Northing'], gauge_loc[0], gauge_loc[1])        
        if dist < min_dist:
            min_dist = dist
            min_bore = row[1]['Bore ID']
        gauges[gauge].append(dist)
        bore_x.append(row[1]['Easting'])
        bore_y.append(row[1]['Northing'])
    if min_dist == 99999:
        print('Default minimum distance used, looks to be problems!')
        sys.exit()
    bore_linking[gauge] = [min_dist, min_bore]
    stream_linking[min_bore] = gauge
    
    # Couldn't get this to work ... would be much more efficient ... might revisit:    
    #relevant_bores[gauge] = relevant_bores.apply(lambda row: points_dist(row['Easting'], row['Northing'], gauge_loc[0], gauge_loc[1]), axis=1)        

plt.scatter(gauge_x, gauge_y,color='r', label='Stream gauges')
plt.scatter(bore_x, bore_y,color='g', label='Bores with EC and DTW')
plt.legend()

# Now get salinity and level data from nearest bores:

relevant_water_levels = {}
relevant_water_salinity = {}

for values in list(bore_linking.values()):
    # get bore ID    
    bore_near = values[1]
    # Filter water levels by bore ID    
    water_levels = df_set['Water Levels'].loc[df_set['Water Levels']['Bore ID'] == bore_near]
    # Create new dataframe for nearest bore but only include date and depth to 
    # water, set the index as the date, resample it to daily data, and infill 
    # missing values with linear interpolation    
    relevant_water_levels[bore_near] = water_levels[['Reading date','Depth to water (m)']].set_index(['Reading date']).resample('1D').interpolate(method='linear')
    # Write data to csv file with name using stream gauge number first 
    # and bore ID second
    relevant_water_levels[bore_near].to_csv(str(stream_linking[bore_near])+'_'+str(bore_near)+'_level.csv')
    water_salinity = df_set['Field Chem'].loc[df_set['Field Chem']['Bore ID'] == bore_near]
    relevant_water_salinity[bore_near] = water_salinity[['Date','EC (uS/cm)']].set_index(['Date']).resample('1D').interpolate(method='linear')
    # Convert salinity from EC to TDS ppm    
    relevant_water_salinity[bore_near] = relevant_water_salinity[bore_near].apply(EC2ppm)
    # Change column name
    relevant_water_salinity[bore_near].columns = ['TDS (ppm)']
    relevant_water_salinity[bore_near].to_csv(str(stream_linking[bore_near])+'_'+str(bore_near)+'_salinity.csv')        

# plot up available data
plt.figure()      
for key in relevant_water_salinity:
    plt.plot(relevant_water_salinity[key], '-', label=key)
plt.xlabel('Date')
plt.ylabel(relevant_water_salinity[key].columns[0])
plt.legend()

plt.figure()      
for key in relevant_water_levels:
    plt.plot(relevant_water_levels[key], '-', label=key)      
plt.xlabel('Date')
plt.ylabel('Depth to water (m)')
plt.legend()

      