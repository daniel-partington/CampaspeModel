import numpy as np

def load(shape, bbox, dels, zbot):
	nlay, nrow, ncol = shape
	min_x, max_x, min_y, max_y = bbox
	delr, delc = dels
	
	# Time step parameters
	days = 365
	perlen = [1] + [1 for i in range(days)]
	nper = len(perlen)
	nstp = [1] + [1 for i in range(days)]
	steady = [True] + [False for i in range(days)]
	
	ztop = 0.

	botm = np.linspace(ztop, zbot, nlay + 1)

	return {
		"nlay": nlay,
		"nrow": nrow,
		"ncol": ncol,
		"delr": delr,
		"delc": delc,
		"top": ztop,
		"botm": botm[1:],
		"nper": nper,
		"perlen": perlen,
		"nstp": nstp,
		"steady": steady
	}