### download data
```
wget <SECRET URL ASK TAKUYA>
tar -xvf campaspeim_data.tar 
```

### download code
```
git clone git@github.com:daniel-partington/HydroModelBuilder.git
cd HydroModelBuilder
git checkout dev-working
cd ..

git clone git@github.com:daniel-partington/CampaspeModel.git
cd CampaspeModel
git checkout Amalgamation
cd ..
```

### build dockerfile
```
docker build -t campaspe CampaspeModel/docker_build/.
```

### prepare folder for output
```
mkdir data_build
```

### build model
```
docker run -it \
-e PYTHONPATH=/shared/CampaspeModel:/shared/HydroModelBuilder \
-v $PWD:/shared \
campaspe \
python /shared/CampaspeModel/CampaspeModel/GW_link_Integrated/GW_link_Integrated_build.py
```

### TODO run model
```
# just edit paths in config?
mkdir /app/CampaspeModel/CampaspeModel/testbox/integrated/data/structured_model_grid_5000m
cp /shared/campaspeim_data/Groundwater/2017-08-21/structured_model_grid_5000m/* /app/CampaspeModel/CampaspeModel/testbox/integrated/data/structured_model_grid_5000m/


docker run -it \
-e PYTHONPATH=/shared/CampaspeModel:/shared/HydroModelBuilder \
-v $PWD:/shared \
campaspe \
python /shared/CampaspeModel/CampaspeModel/GW_link_Integrated/GW_link_Integrated.py
```

