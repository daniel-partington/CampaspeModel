def load(shape):
	nlay, nrow, ncol = shape
	# Remember to use zero-based layer, row, column indices!
	pumping_rate = -500.
	wel_sp1 = [[0, nrow/2 - 1, ncol/2 - 1, 0.]]
	wel_sp2 = [[0, nrow/2 - 1, ncol/2 - 1, 0.]]
	wel_sp3 = [[0, nrow/2 - 1, ncol/2 - 1, pumping_rate]]
	stress_period_data = {0: wel_sp1, 1: wel_sp2, 200: wel_sp3}

	return {"stress_period_data": stress_period_data}