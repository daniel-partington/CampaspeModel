import os
import datetime
import pandas as pd

# If files are in another directory, change the following line to the file path 
# of you files, e.g. r"C:\Home\my_files"
working_directory = os.getcwd()

## Get the flow between the 31/3/2016 and 5/04/2016
#
#start_date = datetime.date(2016,03,31)
#end_date = datetime.date(2016,04,02)


def getFlow(path=None, start_date=None, end_date=None, summary=False):
    if path == None:
        path = working_directory
    #end if
    site_details_file = "Site Details.csv"
    site_details = pd.read_csv(os.path.join(path, site_details_file))
    flow_files = [x for x in os.listdir(path) if "MeanWaterFlow" in x]
    
    flow_stations = [int(x.split('.')[0]) for x in flow_files]
    
    relevant_data = {'Site ID':[], 'Site Name':[], 'Easting':[], 'Northing':[], 'Mean flow (m3/s)':[], 'Max flow (m3/s)':[], '5th percentile flow (m3/s)':[], 'Min flow (m3/s)':[]}
    
    for index, flow_file in enumerate(flow_files):
        flow_file_df = pd.read_csv(os.path.join(path, flow_file), skiprows=2, index_col='Date', parse_dates=True, dayfirst=True)
        if start_date != None:
            flow_file_df = flow_file_df.ix[start_date:end_date]
        #end if
        if flow_file_df.empty:
            continue
        #end if
        flow_file_df = flow_file_df.fillna(0)        
        if flow_file_df['Mean'].max() == 0:
            continue
        #end if
        relevant_site_details = site_details[site_details['Site Id']==[flow_stations[index]]]

        relevant_data['Site ID'] += [relevant_site_details['Site Id'].values[0]]
        relevant_data['Site Name'] += [relevant_site_details['Site Name'].values[0]]
        relevant_data['Easting'] += [float(relevant_site_details['Easting'])]
        relevant_data['Northing'] += [float(relevant_site_details['Northing'])]
        relevant_data['Mean flow (m3/s)'] += [flow_file_df['Mean'].mean() * 1000. / 86400]
        relevant_data['Max flow (m3/s)'] += [flow_file_df['Mean'].max() * 1000. / 86400]
        relevant_data['Min flow (m3/s)'] += [flow_file_df['Mean'].min() * 1000. / 86400]
        relevant_data['5th percentile flow (m3/s)'] += [flow_file_df['Mean'].quantile(q=0.05) * 1000. / 86400]
    # end for    
    
    processed_river_sites = pd.DataFrame(relevant_data)
    return processed_river_sites    
    #processed_river_sites.to_csv(os.path.join(working_directory,'processed_river_sites_flow.csv'))
    
###############################################################################

def getStage(path=None, start_date=None, end_date=None, summary=False):
    if path == None:
        path = working_directory
    #end if
    site_details_file = "Site Details.csv"
    site_details = pd.read_csv(os.path.join(path, site_details_file))

    stage_files = [x for x in os.listdir(path) if "MeanWaterLevel" in x]
    #stage_files = [r"406265.MeanWaterLevel.csv"]
                   
    stage_stations = [int(x.split('.')[0]) for x in stage_files]
    
    relevant_data = {'Site ID':[], 'Site Name':[], 'Easting':[], 'Northing':[], 
                     'Mean stage (m)':[], 'High stage (m)':[], '5th percentile stage (m)':[], 
                     'Low stage (m)':[], 'Gauge Zero (Ahd)':[]}
    
    for index, stage_file in enumerate(stage_files):
        stage_file_df = pd.read_csv(os.path.join(path, stage_file), skiprows=2, index_col='Date', parse_dates=True, dayfirst=True)
        if start_date != None:
            stage_file_df = stage_file_df.ix[start_date:end_date]    
        #end if
        if stage_file_df.empty:
            print 'No data at provided dates for: ', stage_file        
            continue
        #end if 
        stage_file_df = stage_file_df.fillna(0)        
        if stage_file_df['Mean'].max() == 0:
            print 'Level reading is 0, so ignoring for: ', stage_file        
            continue
        
        qual_codes_to_ignore = [8, 9, 21, 100, 101, 120, 149, 150, 151, 152, 
                                153, 154, 155, 156, 160, 161, 165, 180, 190, 
                                200, 201, 237, 250, 254, 255]
        #qual_codes_to_include = [1, 2, 82]
        stage_file_df.drop(stage_file_df[stage_file_df['Qual'].isin(qual_codes_to_ignore)].index, inplace=True)
        #stage_file_df = stage_file_df[stage_file_df['Qual'].isin(qual_codes_to_include)]
        
        relevant_site_details = site_details[site_details['Site Id']==[stage_stations[index]]]

        #if float(relevant_site_details['Gauge Zero (Ahd)']) == float(0):
        #    print 'Gauge zero reading is 0, so ignoring for: ', stage_file
        #    continue
    
        mean_stage = stage_file_df['Mean'].mean()
        high_stage = stage_file_df['Mean'].max()
        low_fifth_stage = stage_file_df['Mean'].quantile(q=0.05)
        low_stage = stage_file_df['Mean'].min()
        if mean_stage < 10.0:
            if float(relevant_site_details['Gauge Zero (Ahd)']) > 0:
                mean_stage += float(relevant_site_details['Gauge Zero (Ahd)'])
            else:
                print 'Mean value of stage less than 10m and Gauge Zero not known, so ignoring for: ', stage_file
                continue
                
        relevant_data['Site ID'] += [relevant_site_details['Site Id'].values[0]]
        relevant_data['Site Name'] += [relevant_site_details['Site Name'].values[0]]
        relevant_data['Easting'] += [float(relevant_site_details['Easting'])]
        relevant_data['Northing'] += [float(relevant_site_details['Northing'])]
        relevant_data['Mean stage (m)'] += [mean_stage]
        relevant_data['High stage (m)'] += [high_stage]
        relevant_data['5th percentile stage (m)'] += [low_fifth_stage]
        relevant_data['Low stage (m)'] += [low_stage]
        relevant_data['Gauge Zero (Ahd)'] += [float(relevant_site_details['Gauge Zero (Ahd)'])]
    # end for

    processed_river_sites_stage = pd.DataFrame(relevant_data)
    processed_river_sites_stage.to_csv(os.path.join(working_directory,'processed_river_sites_stage.csv'))

    return processed_river_sites_stage

if __name__ == "__main__":
    run = getStage(path=r"C:\Workspace\part0075\MDB modelling\Campaspe_data\SW\All_streamflow_Campaspe_catchment\\")
    