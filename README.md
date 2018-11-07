# CampaspeModel

CampaspeModel comprises a series of MODFLOW-based surface water-groundwater [models](https://github.com/daniel-partington/CampaspeModel/tree/master/CampaspeModel/models) generated with the help of [HydroModelBuilder](https://github.com/daniel-partington/HydroModelBuilder).
It was developed for the Lower Campaspe catchment in Victoria, Australia.

## Development setup
We use the Anaconda Python Distribution to ensure consistent runtime environments.

Quick setup can be achieved by downloading or cloning this repo and running the following:

`conda env create --file dev_environment.yml`

This will create a new environment called 'gwdev', which can be activated by running

* `activate gwdev` on Windows
* `source activate gwdev` on *nix-like systems

An alternative environment name may be specified using the `-n` switch, e.g.

`conda env create -n my_env_name --file dev_environment.yml`

The environment may be removed with:

`conda env remove -n gwdev`

After environment setup run `conda develop .` from within the repository folder to allow use of CampaspeModel as a Python package (yes, the dot matters).

### Flow and transport modelling

To model flow and transport, CampaspeModel uses [MODFLOW-NWT](https://water.usgs.gov/ogw/modflow-nwt/) and [MT3D-USGS](https://water.usgs.gov/ogw/mt3d-usgs/) which may be obtained (compiled) using [pymake](https://github.com/modflowpy/pymake)


## Other information
CampaspeModel requires additional data files for its initial build which, unfortunately, consists of some proprietary data. Example data suitable for demonstrating the build process is currently not available, but is intended to be included in the future.
