Follow these steps to build and run the Campaspe model.

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
# git reset --hard 32f8bb6bc8885b3ac4685f4459d772220c242b30
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

### run model
```
docker run -it \
-e PYTHONPATH=/shared/CampaspeModel:/shared/HydroModelBuilder \
-v $PWD:/shared \
campaspe \
python /shared/CampaspeModel/CampaspeModel/GW_link_Integrated/GW_link_Integrated.py
```

### TODO view output
