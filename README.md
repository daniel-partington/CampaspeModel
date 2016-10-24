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
The model will create intemediary data files in ```too_big_to_git/clipped_and_projected_folder```.


Running the model
------------------
Set the environment variable ```DATA_CAMPASPE``` to the path of the unzipped data directory. 
Check that ```config.py``` runs.
Then
```
python run.py
```
The first run will take a long time. Slow functions are cached, so subsequent runs should be quicker.


Quantities of Interest
------------------

1. average depth to GW in the farming area
2. depth to GW at key ecological locations
3. GW level at key bores with trigger levels (relating to policy)
4. SW-GW interaction along significant reaches of the river




Structure
------------------
* model boundary (active cells) [Campaspe river catchment](docs/GW_model_area.json)
* all geometries projected to EPSG:28355 [gda94 mga zone 55](http://spatialreference.org/ref/epsg/gda94-mga-zone-55)
* ? years, daily timesteps
* ? by ? horizontal grid
* recharge
	- rainfall timeseries from [SILO](https://www.longpaddock.qld.gov.au/silo/) - 'Kyneton',  'Elmore', 'Rochester', 'Echuca'
	- rainfall timeseries from [BOM](http://www.bom.gov.au/climate/data/) - "Rainfall at Eppalock"
	- monthly mean
	- interpolated using cubic spline
* river
	- Shapes
		- "AHGFMappedStream" layer of "SH_Cartography.gdb" - "MURRAY RIVER", "CAMPASPE RIVER"
		- "AHGFCanalLine" layer of "SH_Cartography.gdb" - "WARANGA WESTERN CHANNEL"
	- guages
		- WARANGA WESTERN CHANNEL
			- discontinued 2007 - 405705 WARANGA WESTERN CHANNEL @ 0 MILE 40 CHAIN STRUCTURE
			- best proxy ? - D406704 DEAKIN MAIN DRAIN @ OUTFALL (ECHUCA)
		- MURRAY RIVER
			- 409200 MURRAY RIVER @ ECHUCA
		- CAMPASPE RIVER
			- 406207 CAMPASPE RIVER @ EPPALOCK
			- 406201 CAMPASPE RIVER @ BARNADOWN
	
* wells
 	- Bore data\pumping wells.shp
* irrigation recharge
* hydraulic properties
* edges (general head boundary)
	- Western GW
	- Eastern GW

