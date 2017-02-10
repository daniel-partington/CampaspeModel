import os
import sys
import shutil
import glob

with open("packaged_folder.txt", 'r') as f:
    path = f.read().strip().strip('"')

# Check that directory can be found
if not os.path.exists(path):
    sys.exit("{} does not exist".format(path))

# Now go through all items in the directory and copy those with "package.pkl"    
files = [x for x in [y for y in glob.glob(path + '/*/**/***')] if x.endswith('packaged.pkl')]       
#[shutil.copy(z, os.getcwd()) for z in files]

cwd = os.getcwd()  
for fil in files:
    folder = fil.split('\\')[-2]
    to_dir = os.path.join(cwd, folder, "")
    if not os.path.exists(to_dir):
        os.mkdir(to_dir)
    shutil.copy(fil, to_dir)
    
# Now go through all items in the directory and copy folder pilot_points "package.pkl"    

pp_folds = [x for x in [y for y in glob.glob(path + '/*/**/***')] if 'pilot_points' in x]       
    
for folder in pp_folds:
    folder2 = folder.split('\\')[-2]
    to_dir = os.path.join(cwd, folder2, "pilot_points")
    if os.path.exists(to_dir):
        shutil.rmtree(to_dir)
        
    ignore_patterns = ('*.in','grid.spc','reg*.dat','*.ref','*.fig', '*.ref', '*.inf', 'struct.dat', 'debug1.dat')
    shutil.copytree(folder, to_dir, ignore=shutil.ignore_patterns(*ignore_patterns))                        