"""
Sets up useful parameters for running debug scripts in each module
"""

import os
import multiprocessing

# TODO shouldn't do this here
# os.environ['DATA_CAMPASPE'] = '/Users/masher/Downloads/data_campaspe/'
os.environ['DATA_CAMPASPE'] = '/Volumes/KINGSTON/downloaded_data/'
data_folder = os.environ['DATA_CAMPASPE']

modelname = 'campaspe'
model_dir = os.path.join(os.path.dirname(__file__), 'too_big_to_git', modelname, str(multiprocessing.current_process().name))
clipped_and_projected_folder = os.path.join(os.path.dirname(__file__), 'too_big_to_git', 'clipped_and_projected_folder')

if not os.path.exists(clipped_and_projected_folder):
	os.makedirs(clipped_and_projected_folder)