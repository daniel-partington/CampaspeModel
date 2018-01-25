import os
import pandas as pd

# If files are in another directory, change the following line to the file path
# of you files, e.g. r"C:\Home\my_files"
working_directory = os.getcwd()


def getFlow(path=None, start_date=None, end_date=None, summary=False, sites=None):

    if path is None:
        path = working_directory
    # End if
    site_details_file = "Site Details.csv"
    site_details = pd.read_csv(os.path.join(path, site_details_file))
    flow_files = [x for x in os.listdir(path) if "MeanWaterFlow" in x]

    if sites is not None:
        flow_files = [x for x in flow_files if int(x.split('.')[0]) in sites]
    # End if

    flow_stations = [int(x.split('.')[0]) for x in flow_files]

    relevant_data = {'Site ID': [], 'Site Name': [], 'Easting': [], 'Northing': [],
                     'Mean flow (m3/s)': [], 'Max flow (m3/s)': [], '5th percentile flow (m3/s)': [],
                     'Min flow (m3/s)': []}

    processed_river_sites_ts = {}

    for index, flow_file in enumerate(flow_files):
        flow_file_df = pd.read_csv(os.path.join(path, flow_file), skiprows=2,
                                   index_col='Date', parse_dates=True, dayfirst=True)
        if start_date is not None:
            flow_file_df = flow_file_df.ix[start_date:end_date]
        # end if

        qual_codes_to_ignore = [8, 9, 21, 100, 101, 120, 149, 150, 151, 152,
                                153, 154, 155, 156, 160, 161, 165, 180, 190,
                                200, 201, 237, 250, 254, 255]

        flow_file_df.drop(flow_file_df[flow_file_df['Qual'].isin(qual_codes_to_ignore)].index, inplace=True)

        if flow_file_df.empty or flow_file_df['Mean'].max() == 0:
            continue
        # End if

        relevant_site_details = site_details[site_details['Site Id'] == [flow_stations[index]]]

        relevant_data['Site ID'] += [relevant_site_details['Site Id'].values[0]]
        relevant_data['Site Name'] += [relevant_site_details['Site Name'].values[0]]
        relevant_data['Easting'] += [float(relevant_site_details['Easting'])]
        relevant_data['Northing'] += [float(relevant_site_details['Northing'])]
        relevant_data['Mean flow (m3/s)'] += [flow_file_df['Mean'].mean() * 1000. / 86400]
        relevant_data['Max flow (m3/s)'] += [flow_file_df['Mean'].max() * 1000. / 86400]
        relevant_data['Min flow (m3/s)'] += [flow_file_df['Mean'].min() * 1000. / 86400]
        relevant_data['5th percentile flow (m3/s)'] += [flow_file_df['Mean'].quantile(q=0.05) * 1000. / 86400]

        processed_river_sites_ts[flow_stations[index]] = flow_file_df
    # end for

    processed_river_sites = pd.DataFrame(relevant_data)
    if summary:
        return processed_river_sites, processed_river_sites_ts
    else:
        return processed_river_sites_ts
    # end if

###############################################################################


def getStage(path=None, start_date=None, end_date=None, summary=True, sites=None):
    '''
    Function to process flow data from the "Water Measurement Information System"
    located at http://data.water.vic.gov.au/monitoring.htm

    '''

    if path is None:
        path = working_directory
    # end if
    site_details_file = "Site Details.csv"
    site_details = pd.read_csv(os.path.join(path, site_details_file))

    stage_files = [x for x in os.listdir(path) if "MeanWaterLevel" in x]

    if sites is not None:
        stage_files = [x for x in stage_files if int(x.split('.')[0]) in sites]

    relevant_data = {'Site ID': [], 'Site Name': [], 'Easting': [], 'Northing': [],
                     'Mean stage (m)': [], 'High stage (m)': [], '5th percentile stage (m)': [],
                     'Low stage (m)': [], 'Gauge Zero (Ahd)': []}

    processed_river_sites_stage_ts = {}
    for index, stage_file in enumerate(stage_files):
        stage_file_df = pd.read_csv(os.path.join(path, stage_file), skiprows=2, index_col='Date',
                                    parse_dates=True, dayfirst=True)
        if start_date is not None:
            stage_file_df = stage_file_df.ix[start_date:end_date]
        # end if
        if stage_file_df.empty:
            print 'No data at provided dates for: ', stage_file
            continue
        # end if
        stage_file_df = stage_file_df.fillna(0)
        if stage_file_df['Mean'].max() == 0:
            print 'Level reading is 0, so ignoring for: ', stage_file
            continue

        qual_codes_to_ignore = [8, 9, 21, 100, 101, 120, 149, 150, 151, 152,
                                153, 154, 155, 156, 160, 161, 165, 180, 190,
                                200, 201, 237, 250, 254, 255]

        stage_file_df.drop(stage_file_df[stage_file_df['Qual'].isin(qual_codes_to_ignore)].index, inplace=True)

        stage_id = int(stage_file.split('.')[0])

        relevant_site_details = site_details[site_details['Site Id'] == stage_id]

        if summary:
            mean_stage = stage_file_df['Mean'].mean()
            high_stage = stage_file_df['Mean'].max()
            low_fifth_stage = stage_file_df['Mean'].quantile(q=0.05)
            low_stage = stage_file_df['Mean'].min()
            if mean_stage < 10.0:
                if float(relevant_site_details['Gauge Zero (Ahd)']) > 0.:
                    mean_stage += float(relevant_site_details['Gauge Zero (Ahd)'])
                    stage_file_df['Mean'] = stage_file_df['Mean'] + float(relevant_site_details['Gauge Zero (Ahd)'])
                elif float(relevant_site_details['Cease to flow level']) > 10.:
                    mean_stage += float(relevant_site_details['Cease to flow level'])
                    stage_file_df['Mean'] = stage_file_df['Mean'] + float(relevant_site_details['Cease to flow level'])
                elif float(relevant_site_details['Min value']) > 10.:
                    mean_stage += float(relevant_site_details['Min value'])
                    stage_file_df['Mean'] = stage_file_df['Mean'] + float(relevant_site_details['Min value'])
                else:
                    print 'Mean value of stage less than 10m and Gauge Zero not known, so ignoring for: ', stage_file
                    continue
                # End if
            # End if

            mean_stage = stage_file_df['Mean'].mean()
            high_stage = stage_file_df['Mean'].max()
            low_fifth_stage = stage_file_df['Mean'].quantile(q=0.05)
            low_stage = stage_file_df['Mean'].min()

            relevant_data['Site ID'] += [relevant_site_details['Site Id'].values[0]]
            relevant_data['Site Name'] += [relevant_site_details['Site Name'].values[0]]
            relevant_data['Easting'] += [float(relevant_site_details['Easting'])]
            relevant_data['Northing'] += [float(relevant_site_details['Northing'])]
            relevant_data['Mean stage (m)'] += [mean_stage]
            relevant_data['High stage (m)'] += [high_stage]
            relevant_data['5th percentile stage (m)'] += [low_fifth_stage]
            relevant_data['Low stage (m)'] += [low_stage]
            relevant_data['Gauge Zero (Ahd)'] += [float(relevant_site_details['Gauge Zero (Ahd)'])]
        # end if
        processed_river_sites_stage_ts[stage_id] = stage_file_df
    # end for

    processed_river_sites_stage = pd.DataFrame(relevant_data)
    processed_river_sites_stage.to_csv(os.path.join(working_directory, 'processed_river_sites_stage.csv'))

    if summary:
        return processed_river_sites_stage_ts, processed_river_sites_stage
    else:
        return processed_river_sites_stage_ts
    # End if

###############################################################################


def getEC(path=None, start_date=None, end_date=None, summary=False, sites=None):

    if path is None:
        path = working_directory
    # end if
    site_details_file = "Site Details.csv"
    site_details = pd.read_csv(os.path.join(path, site_details_file))
    ec_files = [x for x in os.listdir(path) if "Conductivity" in x]

    if sites is not None:
        ec_files = [x for x in ec_files if int(x.split('.')[0]) in sites]

    ec_stations = [int(x.split('.')[0]) for x in ec_files]

    relevant_data = {'Site ID': [], 'Site Name': [], 'Easting': [], 'Northing': [], 'Mean EC (??)': [],
                     'Max EC (??)': [], '5th percentile EC (??)': [], 'Min EC (??)': []}

    processed_ec_sites_ts = {}

    for index, ec_file in enumerate(ec_files):
        ec_file_df = pd.read_csv(os.path.join(path, ec_file), skiprows=2,
                                 index_col='Date', parse_dates=True, dayfirst=True)
        if start_date is not None:
            ec_file_df = ec_file_df.ix[start_date:end_date]
        # End if
        if ec_file_df.empty or ec_file_df['Mean'].max() == 0:
            continue
        # End if

        relevant_site_details = site_details[site_details['Site Id'] == [ec_stations[index]]]

        relevant_data['Site ID'] += [relevant_site_details['Site Id'].values[0]]
        relevant_data['Site Name'] += [relevant_site_details['Site Name'].values[0]]
        relevant_data['Easting'] += [float(relevant_site_details['Easting'])]
        relevant_data['Northing'] += [float(relevant_site_details['Northing'])]
        relevant_data['Mean EC (??)'] += [ec_file_df['Mean'].mean()]
        relevant_data['Max EC (??)'] += [ec_file_df['Mean'].max()]
        relevant_data['Min EC (??)'] += [ec_file_df['Mean'].min()]
        relevant_data['5th percentile EC (??)'] += [ec_file_df['Mean'].quantile(q=0.05)]

        processed_ec_sites_ts[ec_stations[index]] = ec_file_df
    # End for

    processed_ec_sites = pd.DataFrame(relevant_data)
    if summary:
        return processed_ec_sites, processed_ec_sites_ts
    else:
        return processed_ec_sites_ts
    # end if


if __name__ == "__main__":
    run = getStage(path=r"C:\Workspace\part0075\MDB modelling\Campaspe_data\SW\All_streamflow_Campaspe_catchment\Updated\June2017\MOST RECENT")
    run2 = getFlow(path=r"C:\Workspace\part0075\MDB modelling\Campaspe_data\SW\All_streamflow_Campaspe_catchment\Updated\June2017\MOST RECENT")
    run3 = getEC(path=r"C:\Workspace\part0075\MDB modelling\Campaspe_data\SW\All_streamflow_Campaspe_catchment\Updated\June2017\MOST RECENT")
