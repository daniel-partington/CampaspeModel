"""
Script to load up the hydrogeological properties database for all of the layers
where data is available. 
"""

import pandas as pd

def get_hgu_properties(fname):
    
    hgu_data = pd.read_excel(fname, sheetname='Hydrogeologic properties summar', skiprows=1, index_col=1)

    def two_col_avg(row, a ,b):
        if row[a] == '-' and row[b] == '-':
            unknown = {'T ':100, 'Kh':1, 'Kz':0.1, 'Sy':0.25, 'SS':0.00001}
            new = unknown[a[0:2]] #'-'
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
        
    hgu_data['T mean'] = hgu_data.apply(lambda row: two_col_avg(row, 'T Lower', 'T Upper'), axis=1)
    hgu_data['Kh mean'] = hgu_data.apply(lambda row: two_col_avg(row, 'Kh Lower', 'Kh Upper'), axis=1)
    hgu_data['Kz mean'] = hgu_data.apply(lambda row: two_col_avg(row, 'Kz Lower', 'Kz Upper'), axis=1)
    hgu_data['Sy mean'] = hgu_data.apply(lambda row: two_col_avg(row, 'Sy Lower', 'Sy Upper'), axis=1)
    hgu_data['SS mean'] = hgu_data.apply(lambda row: two_col_avg(row, 'SS Lower', 'SS Upper'), axis=1)

    return hgu_data

if __name__ == "__main__":

    file_location = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\GW\Aquifer properties\Hydrogeologic_variables.xlsx"
    hgu_data = get_hgu_properties(file_location)