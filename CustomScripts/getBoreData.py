import pandas as pd
import dbf2df
import utils.cache

@utils.cache.money()
def getBoreData(level_file, construction_file): 
    # TODO ngis_shp_NSW\salinity_NSW.csv, fields_salinity = ['bore_id', 'bore_date', 'uom', 'result']
    # TODO NGIS_HydrogeologicUnit.dbf, cols=["HGUNumber", "HGCCode"]
    # TODO NGIS_BoreholeLog.dbf", cols=["HydroCode", "HGUNumber"]

    fields_level = ['bore_id', 'bore_date', 'obs_point_datum', 'result', 'hydroid']

    dfVIC_level = pd.read_csv(level_file, sep=r',', usecols=fields_level, dtype={fields_level[0]:str, fields_level[4]:str})
    
    df_ConstructionLog_VIC = dbf2df.dbf2df(construction_file, cols=["BoreID", "HydroCode", "TopElev", "BottomElev", "Constructi"])

    df_ConstructionLog_VIC["BoreID"] = df_ConstructionLog_VIC["BoreID"].astype(str)
    
    # Only use reading in AHD ... would be nice to later convert the other ones
    print 'Total level records: ', dfVIC_level.shape[0]

    dfVIC_level = dfVIC_level[dfVIC_level['obs_point_datum'] == "RSWL (mAHD)"]
    df_ConstructionLog_VIC = df_ConstructionLog_VIC[df_ConstructionLog_VIC['Constructi'] == "INLT"]

    # Get rid of unnecessary columns:
    dfVIC_level = dfVIC_level.drop(dfVIC_level[['obs_point_datum', 'hydroid']], axis=1)

    # Group bores by ID and get the mean of the heads
    dfVIC_level_summary = dfVIC_level.groupby('bore_id').count()
    dfVIC_level_summary['mean level'] = dfVIC_level.groupby('bore_id').mean()
   
    print 'Total number of unique bores with level readings: ', dfVIC_level_summary.shape[0]
    # Filter out bores with less than obs_num_min records
    obs_num_min = 1
    dfVIC_level_summary = dfVIC_level_summary[dfVIC_level_summary['result'] > obs_num_min]

    print 'Total number of unique bores with at least %i readings: ' %(obs_num_min), dfVIC_level_summary.shape[0]
    # Get column with index
    dfVIC_level_summary['HydroCode'] = dfVIC_level_summary.index

    # Filter original dataset
    dfVIC_level = dfVIC_level[dfVIC_level['bore_id'].isin(dfVIC_level_summary.index)]
    
    # Rename column id of 'bore_id' to bring inline with dbf files 'HydroCode'
    dfVIC_level.rename(columns={'bore_id':'HydroCode'}, inplace=True) 

    # Get bore construction info
    df_bore_construction_info = pd.merge(dfVIC_level_summary, df_ConstructionLog_VIC, how='inner', on=['HydroCode'])
    
    # For bores with multiple entries, they are ambiguous, so remove
    df_bores_clear = df_bore_construction_info.groupby('HydroCode').count()
 
    print 'Total number of bores with levels and screen info: ', df_bores_clear.shape[0] 
 
    # Filter bores by those with only one construction record as multiscreened wells are ambiguous with respect to observations in NGIS database   
    # df_bores_clear = df_bores_clear[df_bores_clear['result'] < 3]
    
    # Assume bottom is the screened part and that well is not multi-screened    
    df_bores_clear['mean level'] = df_bore_construction_info.groupby('HydroCode').min()['mean level']
    df_bores_clear['BottomElev'] = df_bore_construction_info.groupby('HydroCode').min()['BottomElev']
    df_bores_clear['TopElev'] = df_bore_construction_info.groupby('HydroCode').min()['TopElev']

    # There is probably a cleaner way to do this ... but ...
    # Remove unnecessary columns

    df_bores_clear = df_bores_clear[['mean level', 'BottomElev', 'TopElev']]
    
    #print 'Total number of bores with levels and screen info that is non-ambiguous: ', df_bores_clear.shape[0]
 
    df_level_ordered = dfVIC_level.sort_values(['HydroCode', 'bore_date']) 

    return df_level_ordered, df_bores_clear #, df_ConstructionLog_VIC #, df_HydrogeologicUnit, df_level, df_salinity    
