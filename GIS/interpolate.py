import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate import Rbf


def make_mesh_centroids(bbox, dels, shape):
	min_x, max_x, min_y, max_y = bbox
	nlay, nrow, ncol = shape
	delr, delc = dels

	x = np.linspace(min_x + delr / 2.0, max_x - delr / 2.0, ncol)
	y = np.linspace(min_y + delc / 2.0, max_y - delc / 2.0, nrow)

	X, Y = np.meshgrid(x, y)
	mesh_centroids = (X, Y)

	return mesh_centroids


def interpolate(points, values, mesh_centroids):
	"""
	Options for griddata method are ['nearest', 'linear', 'cubic']
	This has been modded to use nearest method outside of interpolation points
	convex hull.    
	"""

	z0 = griddata(points, values, mesh_centroids, method='nearest')
	z1 = griddata(points, values, mesh_centroids, method='cubic')
	z1[np.isnan(z1)] = z0[np.isnan(z1)]

	return z1


def rbf_interpolate(points, values, mesh_centroids):
    """
    Available types of interpolation are:
    ['multiquadric', 'inverse', 'gaussian', 'linear', 'cubic', 'quintic', 'thin_plate']
    """
    rbf = Rbf(np.array(points)[:,0], np.array(points)[:,1], np.array(values), 
    		epsilon=2, function='multiquadric')
    ZI = rbf(mesh_centroids[0], mesh_centroids[1])
    return ZI