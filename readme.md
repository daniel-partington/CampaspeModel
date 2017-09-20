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

### run
```
docker run -it \
-e PYTHONPATH=/shared/CampaspeModel \
-v $PWD:/shared \
campaspe \
python /shared/CampaspeModel/CampaspeModel/GW_link_Integrated/GW_link_Integrated_build.py
```


### TODO
cp /shared/campaspeim_data/Groundwater/2017-08-21/GW_data/model_*.csv /app/CampaspeModel/CampaspeModel/testbox/integrated/data/

# edit /app/CampaspeModel/CampaspeModel/config/model_config.json

# BUILD
python /app/CampaspeModel/CampaspeModel/GW_link_Integrated/GW_link_Integrated_build.py

# RUN
python /app/CampaspeModel/CampaspeModel/GW_link_Integetrated/GW_link_Integrated.py

mkdir /app/CampaspeModel/CampaspeModel/testbox/integrated/data/structured_model_grid_5000m
cp /shared/campaspeim_data/Groundwater/2017-08-21/structured_model_grid_5000m/* /app/CampaspeModel/CampaspeModel/testbox/integrated/data/structured_model_grid_5000m/
