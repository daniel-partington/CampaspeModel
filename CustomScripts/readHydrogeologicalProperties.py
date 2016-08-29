"""
Script to load up the hydrogeological properties database for all of the layers
where data is available. 
"""

import pandas as pd
import numpy as np


def getHGUproperties(fname):
    
    HGUdata = pd.read_excel(fname, sheetname='Hydrogeologic properties summar', skiprows=1, index_col=1)

    DEFAULT = {'T': 100, 'Kh': 1, 'Kz': 0.1, 'Sy': 0.25, 'SS': 0.00001}
        
    for prop in DEFAULT.keys():
        HGUdata['%s mean' % prop] = np.nanmean([HGUdata['%s Lower' % prop], HGUdata['%s Upper' % prop]], axis=0)
        HGUdata['%s mean' % prop][np.isnan(HGUdata['%s mean' % prop])] = DEFAULT[prop]
    
    return HGUdata

if __name__ == "__main__":

    # file_location = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\GW\Aquifer properties\Hydrogeologic_variables.xlsx"
    # HGUdata = getHGUproperties(file_location)
    import os
    data_folder = "/Users/masher/Desktop/phd/data_campaspe"
    HGU_file = os.path.join(data_folder, "GW", "Aquifer properties", "Hydrogeologic_variables.xlsx")
    HGUdata = getHGUproperties(HGU_file)
    print HGUdata