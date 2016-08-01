import pandas as pd
import dbf2df

def getBoreData(): 
    VIC_level_data = r"C:\Workspace\part0075\MDB modelling\ngis_shp_VIC_2016\level_VIC.csv"
    VIC_salinity_data = r"C:\Workspace\part0075\MDB modelling\ngis_shp_VIC_2016\salinity_VIC.csv"
    #NSW_level_data = r"C:\Workspace\part0075\MDB modelling\ngis_shp_NSW\level_NSW.csv"
    #NSW_salinity_data = r"C:\Workspace\part0075\MDB modelling\ngis_shp_NSW\salinity_NSW.csv"
    
    fields_level = ['bore_id', 'bore_date', 'obs_point_datum', 'result', 'hydroid']
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

    df_ConstructionLog_VIC = dbf2df.dbf2df(r"C:\Workspace\part0075\MDB modelling\ngis_shp_VIC_2016\ngis_shp_VIC\NGIS_ConstructionLog.dbf", cols=["BoreID", "HydroCode", "TopElev", "BottomElev", "Constructi"])
    df_HydrogeologicUnit_VIC = dbf2df.dbf2df(r"C:\Workspace\part0075\MDB modelling\ngis_shp_VIC_2016\ngis_shp_VIC\NGIS_HydrogeologicUnit.dbf", cols=["HGUNumber", "HGCCode"])
    df_BoreholeLog_VIC = dbf2df.dbf2df(r"C:\Workspace\part0075\MDB modelling\ngis_shp_VIC_2016\ngis_shp_VIC\NGIS_BoreholeLog.dbf", cols=["HydroCode", "HGUNumber"])

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

    # Group bores by ID and get the mean of the heads
    dfVIC_level_summary = dfVIC_level.groupby('bore_id').count()
    dfVIC_level_summary['mean level'] = dfVIC_level.groupby('bore_id').mean()
   
    print 'Total number of unique bores with level readings: ', dfVIC_level_summary.shape[0]
    # Filter out bores with less than 10 records
    obs_num_min = 10
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

    
    #df_bores['HGUNumber'] = df_bores.lookup(df_BoreholeLog_VIC["BoreID"], df_BoreholeLog_VIC['HGUNumber'])

        
    return df_level_ordered, df_bores_clear #, df_ConstructionLog_VIC #, df_HydrogeologicUnit, df_level, df_salinity    

    #wellslist = ['82999']
    
    #def getWells(df, wellslist):
    #    wells = df.loc[df['bore_id'].isin(wellslist)]
    #    return wells    
        
    #def plotWells(well):    
    #    return well.result.plot(marker ='o', linestyle="None")

    #filter_wells = getWells(df_level, wellslist)
    #date_sorted_wells = filter_wells.sort_values('bore_date')
    #further_sorting = date_sorted_wells.loc[date_sorted_wells['obs_point_datum']=='RSWL (mAHD)']
    #ax = further_sorting.plot(marker ='o', x='bore_date', y='result')
    #ax.set_ylabel('Head (mAHD)')
    #ax.plot(title="Bore @" + wellslist[0])

if __name__ == "__main__":
    df_level, df_bores = getBoreData()    
    
