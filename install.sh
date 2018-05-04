# compile modflow
mkdir too_big_to_git
cd too_big_to_git
git clone git@github.com:mjasher/pymake.git
# on windows this may work better
# git clone git@github.com:modflowpy/pymake.git
mkdir temp
cd temp
curl -O http://water.usgs.gov/ogw/modflow/MODFLOW-2005_v1.11.00/mf2005v1_11_00_unix.zip
unzip mf2005v1_11_00_unix.zip
python ../pymake/makebin.py -i Unix/src/ -o mf2005mja
cd ../
sudo rm /usr/bin/mf2005
sudo ln -s $(pwd)/temp/mf2005mja /usr/bin/mf2005

# install normal python dependencies
pip install -r ../requirements.txt

# install flopy
git clone git@github.com:modflowpy/flopy.git
export PYTHONPATH=$PYTHONPATH:$(pwd)/flopy

# install depency HydroModelBuilder
git clone git@github.com:daniel-partington/HydroModelBuilder.git
export PYTHONPATH=$PYTHONPATH:$(pwd)/HydroModelBuilder

# as a test, run flopy example, you should see some plots
python flopy/examples/Tutorials/Tutorial02/tutorial02.py

cd ..