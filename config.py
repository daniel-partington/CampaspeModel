"""
Sets up useful parameters for running debug scripts in each module
"""

import os
import multiprocessing

# TODO shouldn't do this here
# os.environ['DATA_CAMPASPE'] = '/Users/masher/Downloads/data_campaspe/'
# os.environ['DATA_CAMPASPE'] = '/Volumes/KINGSTON/downloaded_data/'
os.environ['DATA_CAMPASPE'] = '/shared/campaspeim_data/Groundwater/2017-08-21/GW_data'
data_folder = os.environ['DATA_CAMPASPE']

boundary_file = os.path.join(data_folder, 'input_data', 'GW_model_area.shp')
silo_stations = [{
    'station_name': station,
    'station_file': os.path.join(data_folder, 'Campaspe_data', 'Climate', station + '.txt'),
} for station in ['Kyneton', 'Elmore', 'Rochester', 'Echuca']]
bom_stations = [{
    "station_name": "Rainfall at Eppalock",
    "meta_file": os.path.join(data_folder, 'Campaspe_data', 'Climate', 'Rainfall at Eppalock', 'IDCJAC0009_081083_1800_Note.txt'),
    "data_file": os.path.join(data_folder, 'Campaspe_data', 'Climate', 'Rainfall at Eppalock', 'IDCJAC0009_081083_1800_Data.csv'),
}]
licenses_file=os.path.join(data_folder, "Campaspe_data", "GW", "Bore data", "Groundwater licence information for Dan Partington bc301115.xlsx")




modelname = 'campaspe'
model_dir = os.path.join(os.path.dirname(__file__), 'too_big_to_git', modelname, str(multiprocessing.current_process().name))
clipped_and_projected_folder = os.path.join(os.path.dirname(__file__), 'too_big_to_git', 'clipped_and_projected_folder')

if not os.path.exists(clipped_and_projected_folder):
	os.makedirs(clipped_and_projected_folder)