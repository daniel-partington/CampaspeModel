import pysal as ps
import pandas as pd

def dbf2df(dbf_path, index=None, cols=False, incl_index=False):
    '''
    Read a dbf file as a pandas.DataFrame, optionally selecting the index
    variable and which columns are to be loaded.
    __author__  = "Dani Arribas-Bel <darribas@asu.edu> "
    ...
    Arguments
    ---------
    dbf_path    : str
                  Path to the DBF file to be read
    index       : str
                  Name of the column to be used as the index of the DataFrame
    cols        : list
                  List with the names of the columns to be read into the
                  DataFrame. Defaults to False, which reads the whole dbf
    incl_index  : Boolean
                  If True index is included in the DataFrame as a
                  column too. Defaults to False
    Returns
    -------
    df          : DataFrame
                  pandas.DataFrame object created
    '''
    db = ps.open(dbf_path)
    if cols:
        if incl_index:
            cols.append(index)
        vars_to_read = cols
    else:
        vars_to_read = db.header
    data = dict([(var, db.by_col(var)) for var in vars_to_read])
    if index:
        index = db.by_col(index)
        db.close()
        return pd.DataFrame(data, index=index)
    else:
        db.close()
        return pd.DataFrame(data)
        
if __name__ == "__main__":
    
    df_ConstructionLog = dbf2df(r"C:\Workspace\part0075\MDB modelling\ngis_shp_VIC\ngis_shp_VIC\NGIS_ConstructionLog.dbf", cols=["BoreID","TopElev","BottomElev","Constructi"])
    df_HydrogeologicUnit = dbf2df(r"C:\Workspace\part0075\MDB modelling\ngis_shp_VIC\ngis_shp_VIC\NGIS_HydrogeologicUnit.dbf")
    