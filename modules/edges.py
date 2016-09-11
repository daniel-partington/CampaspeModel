def load(shape, hk, dels, zbot):
	nlay, nrow, ncol = shape
	delr, delc = dels
	# Make list for stress period 1
	stageleft = 10.
	stageright = 10.
	bound_sp1 = []
	for il in range(nlay):
	    condleft = hk * (stageleft - zbot) * delc
	    condright = hk * (stageright - zbot) * delc
	    for ir in range(nrow):
	        bound_sp1.append([il, ir, 0, stageleft, condleft])
	        bound_sp1.append([il, ir, ncol - 1, stageright, condright])
	print('Adding ', len(bound_sp1), 'GHBs for stress period 1.')

	# Make list for stress period 2
	stageleft = 10.
	stageright = 0.
	condleft = hk * (stageleft - zbot) * delc
	condright = hk * (stageright - zbot) * delc
	bound_sp2 = []
	for il in range(nlay):
	    for ir in range(nrow):
	        bound_sp2.append([il, ir, 0, stageleft, condleft])
	        bound_sp2.append([il, ir, ncol - 1, stageright, condright])
	print('Adding ', len(bound_sp2), 'GHBs for stress period 2.')

	# We do not need to add a dictionary entry for stress period 3.
	# Flopy will automatically take the list from stress period 2 and apply it
	# to the end of the simulation
	stress_period_data = {0: bound_sp1, 1: bound_sp2}

	return {"stress_period_data": stress_period_data}
	