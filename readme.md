Follow these steps to build and run the Campaspe model.

### Download data
```
# Data will be publicly hosted at some point, but for now...
wget <SECRET URL ASK TAKUYA>
tar -xvf campaspeim_data.tar
```

### Download code
```
git clone git@github.com:daniel-partington/HydroModelBuilder.git
cd HydroModelBuilder
git checkout dev-working
# git reset --hard 32f8bb6bc8885b3ac4685f4459d772220c242b30
cd ..

git clone git@github.com:daniel-partington/CampaspeModel.git
cd CampaspeModel
git checkout Amalgamation
cd ..
```

### Build dockerfile
```
docker build -t campaspe CampaspeModel/docker_build/.
```

### Prepare folder for output
```
mkdir data_build
```

### Build model

Note that the build process requires a minimum of 8GB of memory.

#### Linux
```
docker run -it \
-e PYTHONPATH=/shared/CampaspeModel:/shared/HydroModelBuilder \
-v $PWD:/shared \
campaspe \
python /shared/CampaspeModel/CampaspeModel/GW_link_Integrated/GW_link_Integrated_build.py
```

#### Windows

On Windows systems, you may encounter an issue with drives/directories not being shareable.
Error messages relate to `firewall blocking access`, or `docker error occurred`.

The suggested workaround is to create a new user and use this newly created user's credentials to share the drive in question via `Settings`->`Shared Drives`.

```
docker run -it ^
-e PYTHONPATH=/shared/CampaspeModel:/shared/HydroModelBuilder ^
-v %cd%:/shared ^
campaspe ^
python /shared/CampaspeModel/CampaspeModel/GW_link_Integrated/GW_link_Integrated_build.py
```

#### Powershell
```
docker run -it `
-e PYTHONPATH=/shared/CampaspeModel:/shared/HydroModelBuilder `
-v ${pwd}:/shared `
campaspe `
python /shared/CampaspeModel/CampaspeModel/GW_link_Integrated/GW_link_Integrated_build.py
```

### Message on successful build

`Packaged into /shared/campaspeim_data/Groundwater/2017-08-21/GW_data/structured_model_grid_5000m`

### Example model run

#### Linux
```
docker run -it \
-e PYTHONPATH=/shared/CampaspeModel:/shared/HydroModelBuilder \
-v $PWD:/shared \
campaspe \
python /shared/CampaspeModel/CampaspeModel/GW_link_Integrated/GW_link_Integrated.py
```

#### Windows
```
docker run -it ^
-e PYTHONPATH=/shared/CampaspeModel:/shared/HydroModelBuilder ^
-v %cd%:/shared ^
campaspe ^
python /shared/CampaspeModel/CampaspeModel/GW_link_Integrated/GW_link_Integrated.py
```

#### Powershell
```
docker run -it `
-e PYTHONPATH=/shared/CampaspeModel:/shared/HydroModelBuilder `
-v ${pwd}:/shared `
campaspe `
python /shared/CampaspeModel/CampaspeModel/GW_link_Integrated/GW_link_Integrated.py
```

### Output on successful run
```
swgw_exchanges [(-1014.78955078125, 12442.4189453125, 10813.69921875, 2225.982177734375)]
avg_depth_to_gw [ (9.734519958496094, 11.03411865234375, 13.479754130045572, 9.599544525146484, 9.240921020507812, 6.667
884826660156, 18.53559398651123, 13.302096198586856, 13.344627380371094, 7.946220397949219, 7.979389391447368, 12.024791
717529297)]
ecol_depth_to_gw [(13.04412841796875, -1.1079940795898438, 12.796051025390625)]
```
