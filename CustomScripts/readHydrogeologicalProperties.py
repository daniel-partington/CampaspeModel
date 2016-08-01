"""
Script to load up the hydrogeological properties database for all of the layers
where data is available. 
"""

import pandas as pd

def getHGUproperties(fname):
    
    HGUdata = pd.read_excel(fname, sheetname='Hydrogeologic properties summar', skiprows=1, index_col=1)

    def two_col_avg(row, a ,b):
        if row[a] == '-' and row[b] == '-':
            new = '-'
            unknown = {'T':, 'Kh':, 'Kz':, 'Sy':, 'SS':}
            return new
        elif row[a] == '-' and row[b] != '-':   
            new = row[b]
            return new
        elif row[a] != '-' and row[b] == '-':    
            new = row[a]
            return new
        else:
            new = (row[a] + row[b])/2.0
            return new
        
    HGUdata['T mean'] = HGUdata.apply(lambda row: two_col_avg(row, 'T Lower', 'T Upper'), axis=1)
    HGUdata['Kh mean'] = HGUdata.apply(lambda row: two_col_avg(row, 'Kh Lower', 'Kh Upper'), axis=1)
    HGUdata['Kz mean'] = HGUdata.apply(lambda row: two_col_avg(row, 'Kz Lower', 'Kz Upper'), axis=1)
    HGUdata['Sy mean'] = HGUdata.apply(lambda row: two_col_avg(row, 'Sy Lower', 'Sy Upper'), axis=1)
    HGUdata['SS mean'] = HGUdata.apply(lambda row: two_col_avg(row, 'SS Lower', 'SS Upper'), axis=1)

    return HGUdata

if __name__ == "__main__":

    file_location = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\GW\Aquifer properties\Hydrogeologic_variables.xlsx"
    HGUdata = getHGUproperties(file_location)