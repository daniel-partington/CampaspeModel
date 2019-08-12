import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import colors as mpl_colors
import numpy as np

def obs_raw_resample_plots(df_raw, df_resample, ob_type, units, save_location='', inset=False, locations=None, model_object=None):
    if inset:
        zone2D_info = _zone_array2layers(model_object.model_mesh3D[1])
    
    for name in df_raw['name'].unique():
        fig = plt.figure()
        ax = fig.add_subplot(111)
        df_raw_name = df_raw[df_raw['name'] == name]
        ax.plot(df_raw_name['datetime'], df_raw_name['value'], label='Raw', 
                marker='o', linewidth=0, markersize=5, markeredgecolor='none', 
                alpha=0.5)
        df_reasmple_name = df_resample[df_resample['name'] == name]
        ax.plot(df_reasmple_name['datetime'], df_reasmple_name['value'], 
                label='Resampled', drawstyle='steps-post', color='red', 
                linewidth=2)
        ax.legend()
        ax.set_xlabel('Date')
        ax.set_ylabel('{} {}'.format(ob_type, units))
        ax.set_title(name)
        if inset:
            left, bottom, width, height = [0.0, 0.7, 0.4, 0.2]
            ax2 = fig.add_axes([left, bottom, width, height])
            plot_stream_reaches_basic(model_object, ax2, zone2D_info)
            loc = locations.loc[name, ['Easting', 'Northing']]
            ax2.scatter(x=loc['Easting'], y=loc['Northing'], color='orange', s=30)
            ax2.set_aspect('equal')
            ax2.set_xticklabels('')
            ax2.set_xlabel('')
            ax2.set_ylabel('')
            ax2.get_xaxis().set_visible(False)
            ax2.get_yaxis().set_visible(False)
            ax2.axis('off')
            ax2.patch.set_alpha(0.1)
            ax2.set_yticklabels('')
        plt.savefig(os.path.join(save_location, 
            '{}_data_resampling_{}.png').format(ob_type, name), dpi=300)
        plt.close()

def _zone_array2layers(zone_array, plots=False):
    '''
    Function to generate 2D masked layer arrays for each zone
    '''
    zones = [x + 1 for x in range(len(np.unique(zone_array)) - 1)]
    layers = zone_array.shape[0]
    zone_mask2D = {}
    zone_top2D = {}

    for index, zone in enumerate(zones):
        zone_mask = np.ma.masked_array(zone_array, 
                                       zone_array == zone).mask
    
        zone_mask2D[index] = np.full_like(zone_mask[0], False, dtype=bool)
        zone_top2D[index] = np.full_like(zone_mask[0], 0, dtype=int)
        for layer in range(layers):
            zone_mask2D[index] |= zone_mask[layer] 
    
        for layer in range(layers, 0, -1):
            zone_top2D[index][zone_mask[layer - 1]] = layer - 1

    if plots:
        import matplotlib.pyplot as plt
        for index, zone in enumerate(zones):
            plt.figure()
            plt.imshow(np.ma.masked_array(zone_top2D[index], ~zone_mask2D[index]), interpolation='none')

    return zone_mask2D, zone_top2D            
    
def plot_stream_reaches_basic(m, ax, zone2D_info, bounds=None, zones=[1, 2, 3, 4, 5]):
    model_mesh3D = m.model_mesh3D
    model_boundary = m.model_boundary
    river_seg = m.river_mapping['Campaspe']
    
    river_segs_reach = [river_seg['iseg'][river_seg['reach'] == x].tolist() for x in river_seg['reach'].unique()]
    nrow = model_mesh3D[0].shape[1]
    ncol = model_mesh3D[0].shape[2]
    delr = m.gridHeight
    delc = m.gridWidth
    xul = model_boundary[0]
    yul = model_boundary[3]
    
    x = np.linspace(xul, xul + ncol * delc, ncol)
    y = np.linspace(yul - nrow * delr, yul, nrow)
    X, Y = np.meshgrid(x, y)

    flatten = lambda l: [item for sublist in l for item in sublist]
    cmap_grey_white = mpl_colors.ListedColormap(['white', 'lightgrey'])
    
    if not zones:
        ax.pcolormesh(X, Y, np.flipud(zone2D_info[0][list(zone2D_info[0].keys())[-1]]), cmap=cmap_grey_white)
    else:
        for i, zone in enumerate(zones):
            if i == 0:
                zone_temp = zone2D_info[0][zone]
            else:
                zone_temp = zone_temp | zone2D_info[0][zone]
        ax.pcolormesh(X, Y, np.flipud(zone_temp), cmap=cmap_grey_white)
        
    for index, reach in enumerate(river_segs_reach[:]):
        reach_river = river_seg[river_seg['iseg'].isin(reach)]
        points = [i for i in reach_river['amalg_riv_points_collection']]
        points = flatten(points)
        reach_dist = 0
        for ind, point in enumerate(points[:-1]):
            reach_dist += np.linalg.norm(np.array(point) - np.array(points[ind + 1]))
        x_points = [i[0] for i in points]
        y_points = [j[1] for j in points]    
        ax.plot(x_points, y_points, c='blue', alpha=0.5)
    #Campaspe_info.plot(kind='scatter', x='Easting', y='Northing', ax=ax)#, label='Site Id')    
    start_ax, end_ax = ax.get_xlim()
    start_ax = start_ax // 1000 * 1000 + 1000
    end_ax = end_ax // 1000 * 1000 - 1000
    ax.xaxis.set_ticks(np.arange(start_ax, end_ax, 20000.))
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
    if bounds:
        ax.set_xlim(bounds[0], bounds[1])
        ax.set_ylim(bounds[2], bounds[3])    
    ax.set_xlabel('Easting', fontsize=10)
    ax.set_ylabel('Northing', fontsize=10)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)        