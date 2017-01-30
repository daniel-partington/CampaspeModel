import os
import sys
import shutil
import glob

with open("packaged_folder.txt", 'r') as f:
    path = f.read().strip().strip('"')

# Check that directory can be found
if not os.path.exists(path):
    sys.exit("{} does not exist".format(path))

# Now go through all items in the directory and copy those with package.    
files = [x for x in [y for y in glob.glob(path + '/*/**/***')] if x.endswith('packaged.pkl')]       
#[shutil.copy(z, os.getcwd()) for z in files]

cwd = os.getcwd()  
for fil in files:
    folder = fil.split('\\')[-2]
    to_dir = os.path.join(cwd, folder, "")
    if not os.path.exists(to_dir):
        os.mkdir(to_dir)
    shutil.copy(fil, to_dir)