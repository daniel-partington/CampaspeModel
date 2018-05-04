import numpy as np

def load(shape):
	hk = 1. * np.ones(shape)
	vka = 1.
	sy = 0.1
	ss = 1.e-4
	laytyp = 1

	return {
		"hk": hk,
		"vka": vka,
		"sy": sy,
		"ss": ss,
		"laytyp": laytyp
	}
