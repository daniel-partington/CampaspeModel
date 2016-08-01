import os

import pandas
from osgeo import gdal, ogr, osr

def get_GW_licence_info(filename, path=None, out_file=None, out_path=None):
    
    if path == None:
        filename = filename
    else:
        filename = path + filename
    
    df = pandas.read_excel(filename, skiprows=10)
    
    cols = ['RECORD ID', 'Water System Source', 'Trading Zone', 'Use of Water', 'Entitlement Volume']
    
    df[cols] = df[cols].fillna(method='ffill')

    df = df.fillna(value=0.0)    
    
    #for trade_zone in pandas.Series.unique(df['Trading Zone']):
    #    print 'Total volume in %s: ' % (str(trade_zone))
    #    print df.loc[df['Trading Zone'] == trade_zone]['Use 2014/15'].sum()
    #print df.groupby('Trading Zone').loc['Bamawn - Lower Campaspe Valley'] #.sum() #.isin(trade_zone)]   

    # set up the shapefile driver
    driver = ogr.GetDriverByName("ESRI Shapefile")
    
    # create the data source
    if out_path == None:
        out_file = out_file
    else:
        out_file = out_path + out_file

    if os.path.exists(out_file):
        os.remove(out_file)

    data_source = driver.CreateDataSource(out_file)
    
    # create the spatial reference, WGS84
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(28355)

    # create the layer
    layer = data_source.CreateLayer("pumping wells", srs, ogr.wkbPoint)
    
    # Add the fields we're interested in
    field_name = ogr.FieldDefn("OLD ID", ogr.OFTString)
    field_name.SetWidth(24)
    layer.CreateField(field_name)
    field_name2 = ogr.FieldDefn("Works ID", ogr.OFTString)
    field_name2.SetWidth(24)
    layer.CreateField(field_name2)
    field_region = ogr.FieldDefn("TradeZone", ogr.OFTString)
    field_region.SetWidth(24)
    #field_region.SetName("Trading Zone")
    layer.CreateField(field_region)
    layer.CreateField(ogr.FieldDefn("Northing", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("Easting", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("TopScreen", ogr.OFTReal))

    bot_screen = ogr.FieldDefn("BotScreen", ogr.OFTReal)
    #bot_screen.SetName("Bottom Screen")
    layer.CreateField(bot_screen)
    
    # Process the text file and add the attributes and features to the shapefile
    for row in df.iterrows():
        row = row[1]
        # create the feature
        feature = ogr.Feature(layer.GetLayerDefn())
        # Set the attributes using the values from the delimited text file
        feature.SetField("OLD ID", str(row['OLD ID']))
        feature.SetField("Works ID", str(row['Works ID']))
        feature.SetField("TradeZone", str(row['Trading Zone']))
        feature.SetField("Northing", float(row['Northing']))
        feature.SetField("Easting", float(row['Easting']))
        feature.SetField("TopScreen", float(row['Top screen depth (m)']))
        feature.SetField("BotScreen", float(row['Bottom screen depth (m)']))
        
        # create the WKT for the feature using Python string formatting
        wkt = "POINT(%f %f)" %  (float(row['Easting']) , float(row['Northing']))
    
        # Create the point from the Well Known Txt
        point = ogr.CreateGeometryFromWkt(wkt)
    
        # Set the feature geometry using the point
        feature.SetGeometry(point)
        # Create the feature in the layer (shapefile)
        layer.CreateFeature(feature)
        # Destroy the feature to free resources
        feature = None #.Destroy()
    
    # Destroy the data source to free resources
    data_source = None

    df = df[['Annual Volume', 'OLD ID', 'Works ID','Top screen depth (m)']]
    df['OLD ID'] = df['OLD ID'].astype(str)    
    df['Works ID'] = df['Works ID'].astype(str)    
    df = df.set_index('OLD ID')
    #df2['Annual Volume'] = df2['Annual Volume']/365.
    return df #['OLD ID', 'Annual Volume']
    

if __name__ == "__main__":

    filename = "Groundwater licence information for Dan Partington bc301115.xlsx"
    path = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\GW\Bore data\\"    
    out_path = r"C:\Workspace\part0075\MDB modelling\Campaspe_data\GW\Bore data\\"
    out_file = "pumping wells.shp"
    Original, GW_pumping_bores = get_GW_licence_info(filename, path=path, out_file=out_file, out_path=out_path)    