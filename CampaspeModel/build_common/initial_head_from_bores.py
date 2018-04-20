import numpy as np
from HydroModelBuilder.Utilities import interpolation
import matplotlib.pyplot as plt

def generate_initial_head_from_bores(model_builder_object, df_ts, df_info, 
                                     time_min=None, time_max=None,
                                     interp_method='nearest', plot=False):
    '''
    Generate an initial head array for model mesh using time series data of 
    bore heads.
    
    :param df_ts: Pandas dataframe, time series of bore levels at all bores of interest
    :param df_info: Pandas dataframe, contains easting and northing for each of the bores
    :param time_max: str, maximum allowable datetime to keep in the time series df_ts
    :param time_min: str, minimum allowable datetime to keep in the time series df_ts
    :param interp_method: str, interpolation method to use in generating the heads
    :param plot: boolean, define whether to output plots 
    
    returns 2D numpy array of initial head in shape
    '''
    
    mbo = model_builder_object
    # Get policy bores with data before time_max and after time_min:
    if time_max:
        df_ts_filter = df_ts[df_ts['datetime'] < time_max]
    if time_min:
        df_ts_filter = df_ts[df_ts['datetime'] >= time_min]

    df_ts_filter = df_ts_filter[df_ts_filter['active'] == True]
    df_ts_filter.loc[:, 'HydroCode'] = df_ts_filter['name']
    df_ts_filter_loc_merge = df_ts_filter.merge(df_info, on='HydroCode')

    points = zip(df_ts_filter_loc_merge['Easting'].tolist(), df_ts_filter_loc_merge['Northing'].tolist())
    values = df_ts_filter_loc_merge['value'].tolist()
    xi = mbo.model_mesh_centroids

    grid_z0 = interpolation.Interpolator('structured', np.array(points), np.array(values), xi, method=interp_method)

    if plot:
        extent = (np.min(xi[0]), np.max(xi[0]), np.min(xi[1]), np.max(xi[1]))
        plt.subplot(111)
        plt.imshow(grid_z0.T, extent=extent, origin='lower', interpolation='none')
        plt.title("Heads interpolated using {}".format(interp_method))
        plt.gcf().set_size_inches(6, 6)
        plt.show()

    initial_heads = np.full(mbo.model_mesh3D[1].shape, 0.)
    for i in range(mbo.model_mesh3D[1].shape[0]):
        initial_heads[i] = grid_z0  

    return initial_heads

def generate_initial_head_from_surface_elevation(model_builder_object, plot=False):        
    mbo = model_builder_object
    initial_heads = np.full(mbo.model_mesh3D[1].shape, 0.)
    for i in range(mbo.model_mesh3D[1].shape[0]):
        initial_heads[i] = (mbo.model_mesh3D[0][i] + mbo.model_mesh3D[0][i + 1]) / 2.

    return initial_heads