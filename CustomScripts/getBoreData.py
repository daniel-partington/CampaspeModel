import os
import pandas as pd
import dbf2df
import xlrd
"""
NOTES ON NGIS DATA
The quality codes for water levels are:

Quality code	Description
quality-A	Best available given the technologies, techniques and monitoring objectives at the time of classification
quality-B	Compromised in its ability to truly represent the parameter
quality-C	Estimate
quality-E	Ability to truly represent the monitored parameter is not known
quality-F	Not of release quality or contains missing data
Y	Deemed anomalous by the data provider
N	Deemed not anomalous by the data provider
1	Unedited data; no data in the period have been adjusted from sensing system measurement in any way. Maximum discrepancy between the sensor reading and the reference point is +/- 5 mm.
43	State Observation Bore Network verified data
44	Data from water corporations, other government agencies or consultants
45	Data from universities or catchment management authorities
47	Historic data from the Victorian Groundwater Management System
100	Data for which quality, accuracy, and/or derivation is not known, or data suspect or data supplied by other authorities
157	Data unreliable; issue with bore
165	Suspect or bad data supplied by other authority
201	Rating not available for this stage
255	No data exists
~	Value is approximate
>	Measured value is approximate; true value is greater than measured
<	Measured value is approximate; true value is less than measured
 	No code; assumed ok
 

The quality codes for salinity are:

Quality code	Description
quality-E	Ability to truly represent the monitored parameter is not known
Y	Deemed anomalous by the data provider
N	Deemed not anomalous by the data provider
32	Good result - imported from Water Quality Branch spreadsheets
43	State Observation Bore Network verified data
44	Data from water corporations, other government agencies or consultants
45	Data from universities or catchment management authorities
47	Historic data from the Victorian Groundwater Management System
81	Satisfactory result - uncontrolled sampling, lab results manually entered
82	Satisfactory result - uncontrolled sampling, lab analysis
83	Satisfactory result - uncontrolled sampling, lab results rechecked against lab sheets
85	Satisfactory result - uncontrolled sampling, lab results transferred from Laboratory Information Management System electronically
100	Discrete data - recorded value represents the true value
185	Transferred from Laboratory Information Management System electronically
 	No code; assumed ok
"""

def getBoreDataGMW(path=""):
    '''
    Function to process bore data received from Goulburn Murray Water
    
    '''
    # Now load in the bore data to pandas dataframes
    
    fname = os.path.join(path, r"Shallow monitoring bores bc301115.xlsx")
    fname2 = os.path.join(path, r"State Observation Bores bc271115.xlsx")
    
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

def getBoreData(get='transient', path=""): 
    '''
    Function to process National Groudnwater Information System (NGIS) data
    to extract bores with level readings and that have clear info on 
    the construction, i.e. top and bottom of screen.
    '''

    VIC_level_data = os.path.join(path, "level_VIC.csv")
    VIC_salinity_data = os.path.join(path, "salinity_VIC.csv")
    #NSW_level_data = r"C:\Workspace\part0075\MDB modelling\ngis_shp_NSW\level_NSW.csv"
    #NSW_salinity_data = r"C:\Workspace\part0075\MDB modelling\ngis_shp_NSW\salinity_NSW.csv"
    
    fields_level = ['bore_id', 'bore_date', 'obs_point_datum', 'result', 'quality_flag', 'hydroid']
    fields_salinity = ['bore_id', 'bore_date', 'uom', 'result']
    
    dfVIC_level = pd.read_csv(VIC_level_data, sep=r',', usecols=fields_level, dtype={fields_level[0]:str, fields_level[4]:str})
    #dfNSW_level = pd.read_csv(NSW_level_data, sep=r',', usecols=fields_level, dtype={fields_level[0]:str})
    dfVIC_salinity = pd.read_csv(VIC_salinity_data, sep=r',', usecols=fields_salinity, dtype={fields_salinity[0]:str})
    #dfNSW_salinity = pd.read_csv(NSW_salinity_data, sep=r',', usecols=fields_salinity, dtype={fields_salinity[0]:str})
    
    #df_level = pd.concat([dfVIC_level, dfNSW_level])
    #df_salinity = pd.concat([dfVIC_salinity,dfNSW_salinity])
    
    #del dfVIC_level
    #del dfNSW_level
    #del dfVIC_salinity
    #del dfNSW_salinity

    df_ConstructionLog_VIC = dbf2df.dbf2df(os.path.join(path, r"ngis_shp_VIC\NGIS_ConstructionLog.dbf"), cols=["BoreID", "HydroCode", "TopElev", "BottomElev", "Constructi"])
    df_HydrogeologicUnit_VIC = dbf2df.dbf2df(os.path.join(path, r"ngis_shp_VIC\NGIS_HydrogeologicUnit.dbf"), cols=["HGUNumber", "HGCCode"])
    df_BoreholeLog_VIC = dbf2df.dbf2df(os.path.join(path, r"ngis_shp_VIC\NGIS_BoreholeLog.dbf"), cols=["HydroCode", "HGUNumber"])

    df_ConstructionLog_VIC["BoreID"] = df_ConstructionLog_VIC["BoreID"].astype(str)
    df_BoreholeLog_VIC["HydroCode"] = df_BoreholeLog_VIC["HydroCode"].astype(str) 
    
    
    #df_ConstructionLog_NSW = dbf2df.dbf2df(r"C:\Workspace\part0075\MDB modelling\ngis_shp_NSW\ngis_shp_NSW\NGIS_ConstructionLog.dbf", cols=["BoreID","TopElev", "BottomElev", "Constructi"])
    #df_HydrogeologicUnit_NSW = dbf2df.dbf2df(r"C:\Workspace\part0075\MDB modelling\ngis_shp_NSW\ngis_shp_NSW\NGIS_HydrogeologicUnit.dbf", cols=["HGUNumber", "HGCCode"])
    #df_BoreholeLog_NSW = dbf2df.dbf2df(r"C:\Workspace\part0075\MDB modelling\ngis_shp_NSW\ngis_shp_NSW\NGIS_BoreholeLog.dbf", cols=["BoreID", "HGUNumber"])

    #df_ConstructionLog = pd.concat([df_ConstructionLog_VIC, df_ConstructionLog_NSW])
    #df_HydrogeologicUnit = pd.concat([df_HydrogeologicUnit_VIC, df_HydrogeologicUnit_NSW])
    #df_BoreholeLog = pd.concat([df_BoreholeLog_VIC, df_BoreholeLog_NSW])

    #del df_ConstructionLog_VIC
    #del df_HydrogeologicUnit_VIC
    #del df_BoreholeLog_VIC
    #del df_ConstructionLog_NSW
    #del df_HydrogeologicUnit_NSW
    #del df_BoreholeLog_NSW
    
    # Only use reading in AHD ... would be nice to later convert the other ones
    print 'Total level records: ', dfVIC_level.shape[0]

    dfVIC_level = dfVIC_level[dfVIC_level['obs_point_datum'] == "RSWL (mAHD)"]
    df_ConstructionLog_VIC = df_ConstructionLog_VIC[df_ConstructionLog_VIC['Constructi'] == "INLT"]

    # Get rid of unnecessary columns:
    dfVIC_level = dfVIC_level.drop(dfVIC_level[['obs_point_datum', 'hydroid']], axis=1)

    # Only use data from the state observation bores:
    dfVIC_level = dfVIC_level[dfVIC_level['quality_flag'].isin(['43', '47'])]

    # For the salinity data only three quality codes appear which are all acceptable data,
    # and the rest of the data is missing quality codes.
    #dfVIC_salinity = dfVIC_salinity[dfVIC_salinity['quality_flag'].isin(['43', '47'])]

    # Group bores by ID and get the mean of the heads
    dfVIC_level_summary = dfVIC_level.groupby('bore_id').count()
    dfVIC_level_summary['mean level'] = dfVIC_level.groupby('bore_id').mean()
   
    print 'Total number of unique bores with level readings: ', dfVIC_level_summary.shape[0]
    # Filter out bores with less than obs_num_min records
    obs_num_min = 1
    dfVIC_level_summary = dfVIC_level_summary[dfVIC_level_summary['result'] > obs_num_min]

    print 'Total number of unique bores with at least %i readings: ' %(obs_num_min + 1), dfVIC_level_summary.shape[0]
    # Get column with index
    dfVIC_level_summary['HydroCode'] = dfVIC_level_summary.index

    # Filter original dataset
    dfVIC_level = dfVIC_level[dfVIC_level['bore_id'].isin(dfVIC_level_summary.index)]
    
    # Rename column id of 'bore_id' to bring inline with dbf files 'HydroCode'
    dfVIC_level.rename(columns={'bore_id':'HydroCode'}, inplace=True) 
    dfVIC_salinity.rename(columns={'bore_id':'HydroCode'}, inplace=True) 

    # Get bore construction info
    df_bore_construction_info = pd.merge(dfVIC_level_summary, df_ConstructionLog_VIC, how='inner', on=['HydroCode'])
    
    # For bores with multiple entries, they are ambiguous, so remove
    df_bores_clear = df_bore_construction_info.groupby('HydroCode').count()
 
    print 'Total number of bores with levels and screen info: ', df_bores_clear.shape[0] 
 
    # Filter bores by those with only one construction record as multiscreened wells are ambiguous with respect to observations in NGIS database   
    df_bores_clear = df_bores_clear[df_bores_clear['result'] < 3]

    print 'Total number of bores with levels and screen info non-ambiguous: ', df_bores_clear.shape[0] 
    
    # Assume bottom is the screened part and that well is not multi-screened    
    df_bores_clear['mean level'] = df_bore_construction_info.groupby('HydroCode').min()['mean level']
    df_bores_clear['BottomElev'] = df_bore_construction_info.groupby('HydroCode').min()['BottomElev']
    df_bores_clear['TopElev'] = df_bore_construction_info.groupby('HydroCode').min()['TopElev']

    # There is probably a cleaner way to do this ... but ...
    # Remove unnecessary columns

    df_bores_clear = df_bores_clear[['mean level', 'BottomElev', 'TopElev']]
    
    #print 'Total number of bores with levels and screen info that is non-ambiguous: ', df_bores_clear.shape[0]
 
    df_level_ordered = dfVIC_level.sort_values(['HydroCode', 'bore_date']) 
    df_salinity_ordered = dfVIC_salinity.sort_values(['HydroCode', 'bore_date'])
    # Need to kill bore obs for which the obervation is below the bottom of the screen.
    
    
    
    #df_bores['HGUNumber'] = df_bores.lookup(df_BoreholeLog_VIC["BoreID"], df_BoreholeLog_VIC['HGUNumber'])

    #if get == 'transient':
    #    return dfVIC_level[df_level['HydroCode'].isin()]
    #else:
    #    return df_level_ordered, df_bores_clear #, df_ConstructionLog_VIC #, df_HydrogeologicUnit, df_level, df_salinity    

    return df_level_ordered, df_bores_clear, df_salinity_ordered #, df_ConstructionLog_VIC #, df_HydrogeologicUnit, df_level, df_salinity    

    #filter_wells = getWells(df_level, wellslist)
    #date_sorted_wells = filter_wells.sort_values('bore_date')
    #further_sorting = date_sorted_wells.loc[date_sorted_wells['obs_point_datum']=='RSWL (mAHD)']
    #ax = further_sorting.plot(marker ='o', x='bore_date', y='result')
    #ax.set_ylabel('Head (mAHD)')
    #ax.plot(title="Bore @" + wellslist[0])

if __name__ == "__main__":
    path = r"C:\Workspace\part0075\MDB modelling\ngis_shp_VIC_2016"

    df_level, df_bores, df_salinity = getBoreData(get='transient', path=path)    
    #getBoreDataGMW()
