# CampaspeModel
Model building scripts for the Lower Campaspe Valley Surface Water-Groundwater model.


Installation
------------------
In addition to more common Python libraries listed in [requirements.txt](requirements.txt), CampaspeModel depends on 

* [HydroModelBuilder](https://github.com/daniel-partington/HydroModelBuilder)
* [MODFLOW](http://water.usgs.gov/ogw/modflow/MODFLOW.html)
* [Flopy](https://github.com/modflowpy/flopy)

These are installed into a ```too_big_to_git``` folder by running
```
source install.sh
```


Data
------------------
Email maintainers for a link to download ```data_campaspe.zip```, which contains the required data.
The model will create intemediary data files```too_big_to_git/clipped_and_projected_folder```.


Running the model
------------------
Set the environment variable ```DATA_CAMPASPE``` to the path of the unzipped data directory. Check that ```config.py``` runs.
Then
```
python run.py
```


Quantities of Interest
------------------

1. average depth to GW in the farming area
2. depth to GW at key ecological locations
3. GW level at key bores with trigger levels (relating to policy)
4. SW-GW interaction along significant reaches of the river


Structure
------------------
* daily timesteps
* meter horizontal resolution