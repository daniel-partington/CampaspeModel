import config
import datetime

import flopy

from modules import boundary
from modules import discrete
from modules import acquifer
from modules import edges
from modules import wells
from modules import recharge
from modules import river
from modules import output

"""
parameters should go here
"""
shape = (1, 30, 20)  # (nlay, nrow, ncol)
zbot = -50.

# inclusive
start_date = datetime.datetime(2014, 1, 1)
end_date = datetime.datetime(2014, 12, 31)

boundary_params = boundary.boundary(data_folder=config.data_folder, 
								clipped_and_projected_folder=config.clipped_and_projected_folder,
								shape=shape)

"""
MODFLOW/flopy *_inputs objects
all data inter-dependencies should be visible here (ie. modules shouldn't import one another)
"""

boundary_inputs = boundary.load(projected_boundary_file=boundary_params['projected_boundary_file'],
								shape=shape,
								bbox=boundary_params['bbox'],
								dels=boundary_params['dels'])

discrete_inputs = discrete.load(shape=shape, 
							bbox=boundary_params["bbox"], 
							dels=boundary_params["dels"],
							zbot=zbot)

acquifer_inputs = acquifer.load(shape=shape)

edge_inputs = edges.load(shape=shape, 
						hk=acquifer_inputs['hk'], 
						dels=boundary_params['dels'], 
						zbot=zbot)

well_inputs = wells.load(shape=shape)

recharge_inputs = recharge.load(data_folder=config.data_folder, 
							model_projection=boundary.model_projection, 
							mesh_centroids=boundary_params['mesh_centroids'],
							start_date=start_date, end_date=end_date)

river_inputs = river.load(data_folder=config.data_folder, 
				clipped_and_projected_folder=config.clipped_and_projected_folder, 
				model_projection=boundary.model_projection,
				bbox=boundary_params['bbox'],
				shape=shape,
				hk=acquifer_inputs['hk'])

output_control_inputs = output.load()

"""
flopy admin
"""

# Flopy objects
mf = flopy.modflow.Modflow(config.modelname, exe_name='mf2005', model_ws=config.model_dir)
dis = flopy.modflow.ModflowDis(mf, **discrete_inputs)
bas = flopy.modflow.ModflowBas(mf, **boundary_inputs)
lpf = flopy.modflow.ModflowLpf(mf, **acquifer_inputs)
pcg = flopy.modflow.ModflowPcg(mf)
ghb = flopy.modflow.ModflowGhb(mf, **edge_inputs)
well = flopy.modflow.ModflowWel(mf, **well_inputs)
rch = flopy.modflow.ModflowRch(mf, **recharge_inputs)
# riv = flopy.modflow.ModflowRiv(mf, **river_inputs)
oc = flopy.modflow.ModflowOc(mf, **output_control_inputs)

# Write the model input files
mf.write_input()

# Run the model
success, mfoutput = mf.run_model(silent=False, pause=False)
if not success:
    raise Exception('MODFLOW did not terminate normally.')

"""
Read output and calculate QOIs
"""

output.read_and_plot(mf, shape, boundary_params["dels"], boundary_params["Ls"])
