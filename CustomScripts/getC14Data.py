import pandas as pd

#C14data = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Chemistry\bores_with_Carbon-14_&_new_field_data_added.xls"
C14data = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\Chemistry\C14_locs.xlsx"
df_C14 = pd.read_excel(C14data)
df_C14.drop_duplicates(subset=["Bore_id"], inplace=True)

# Fix ages that are negative to be zero
#df_C14[df_C14['Age (years)'] <= 0] = 0.

# Make sure easting and northing values are floats:
#df_C14[['zone55_easting', 'zone55_northing']] = df_C14[['zone55_easting', 'zone55_northing']]

# Get rid of duplicates


#import matplotlib.pyplot as plt
#
#ax = plt.subplot(1,1,1, aspect='equal')
#
#vmin = 0
#vmax = 5000
#
#scat = ax.scatter(df_C14['zone55_easting'], df_C14['zone55_northing'], s=50, c=df_C14['Age (years)'], vmin=vmin, vmax=vmax, alpha=0.8)
#
#print len(df_C14['zone55_easting'])
#
#plt.colorbar(scat)