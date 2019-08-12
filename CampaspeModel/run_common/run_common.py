import numpy as np
import time

def update_HGU_parameters(m, data_folder, verbose=True):
    if verbose:
        print("************************************************************************")
        print(" Updating HGU parameters ")
    
        # This needs to be automatically generated with the map_raster2mesh routine ...
    zone_map = {1: 'qa', 2: 'utb', 3: 'utqa', 4: 'utam', 5: 'utaf', 6: 'lta', 7: 'bse'}

    default_array = m.model_mesh3D[1].astype(float)
    Zone = np.copy(default_array)
    Kh = np.copy(default_array)
    Kv = np.copy(default_array)
    Sy = np.copy(default_array)
    SS = np.copy(default_array)
    
    def create_pp_points_dict(zone_map, Zone, prop_array, prop_name, m):
        points_values_dict = {}
        for index, key in enumerate(zone_map.keys()):
            for index2, param in enumerate(m.parameters.param_set[prop_name + zone_map[key]]):
                if index2 == 0:
                    points_values_dict[index] = [m.parameters.param[param]['PARVAL1']]
                else: 
                    points_values_dict[index] += [m.parameters.param[param]['PARVAL1']]
        return points_values_dict    
        
    def update_pilot_points(zone_map, Zone, prop_array, par_name, prop_name, prop_folder, m, prop_array_fname):
        points_values_dict = create_pp_points_dict(zone_map, Zone, prop_array, prop_name, m)
        p = m.pilot_points[par_name]
        zones = len(list(zone_map.keys()))
        p.output_directory = os.path.join(data_folder, prop_folder)
        p.update_pilot_points_files_by_zones(zones, points_values_dict)
        time.sleep(3)
        p.run_pyfac2real_by_zones(zones) 
        #p.save_mesh3D_array(filename=os.path.join(data_folder, prop_array_fname))
        return p.val_array

    Kh = update_pilot_points(zone_map, Zone, Kh, 'hk', 'kh', 'hk_pilot_points',
                             m, 'hk_val_array')  
    m.save_array(os.path.join(data_folder, 'Kh'), Kh)

    print(("Erroneous K pilot cells: {}".format(len(Kh[Kh > 10000.]))))
    Kh[Kh > 10000.] = 25.
    Kv = Kh * 0.1
    m.save_array(os.path.join(data_folder, 'Kv'), Kv)

    Sy = update_pilot_points(zone_map, Zone, Sy, 'sy', 'sy', 'sy_pilot_points',
                             m, 'sy_val_array')
    print(("Erroneous Sy pilot cells: {}".format(len(Sy[Sy > 0.5]))))
    Sy[Sy > 0.5] = 0.25
    m.save_array(os.path.join(data_folder, 'Sy'), Sy)
    
    SS = update_pilot_points(zone_map, Zone, SS, 'ss', 'ss', 'ss_pilot_points',
                             m, 'ss_val_array')
    print(("Erroneous Ss pilot cells: {}".format(len(SS[SS > 0.01]))))
    SS[SS > 0.01] = 1E-5
    m.save_array(os.path.join(data_folder, 'SS'), SS)
    
    m.properties.update_model_properties('Kh', Kh)
    m.properties.update_model_properties('Kv', Kv)
    m.properties.update_model_properties('Sy', Sy)
    m.properties.update_model_properties('SS', SS)
    
def update_SFR_river_parameters(m, data_folder, verbose=True):
    if verbose:
        print("************************************************************************")
        print(" Updating river parameters ")
    
    reach_df = m.mf_sfr_df['Campaspe']['reach_df'] 
    segment_data = m.mf_sfr_df['Campaspe']['seg_df']

    num_reaches = m.pilot_points['Campaspe'].num_points #4
    known_points = m.pilot_points['Campaspe'].points
    # Create reach data
    river_seg = m.river_mapping['Campaspe']
    
    strcond_val = [m.parameters.param['kv_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    river_seg.loc[river_seg['strhc1'] != 0.0, 'strhc1'] = np.interp(
            river_seg[river_seg['strhc1'] != 0.0]['Cumulative Length'].tolist(), 
            known_points, strcond_val)
    
    strthick_val = [m.parameters.param['bedthck{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    river_seg['strthick'] = np.interp(river_seg['Cumulative Length'].tolist(), known_points, strthick_val)
    
    reach_df = river_seg[['k','i','j','iseg','ireach','rchlen','strtop','slope','strthick','strhc1']]
    reach_data = reach_df.to_records(index=False)
    
    # Set the roughness for the channel
    roughch_val = [m.parameters.param['mn_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    roughch = np.interp(river_seg['Cumulative Length'].tolist(), 
                                    known_points, roughch_val)
    # Set the roughness for the banks
    roughbk_val = [m.parameters.param['mn_riv{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    roughbk = np.interp(river_seg['Cumulative Length'].tolist(), 
                                    known_points, roughbk_val)
    river_seg['roughch'] = roughch
    river_seg['roughbk'] = roughbk
    
    width1_val = [m.parameters.param['rivwdth{}'.format(x)]['PARVAL1'] for x in range(num_reaches)] 
    width1 = np.interp(river_seg['Cumulative Length'].tolist(), 
                                    known_points, width1_val)

    segment_data1 = {}
    for key in list(segment_data.keys()):
        segment_data[key]['width2'] = segment_data[key]['width1'] = width1
        segment_data1[key] = segment_data[key].to_records(index=False)
    
    seg_dict = segment_data1   
    
    if verbose:
        print("************************************************************************")
        print(" Updating Campaspe river boundary")
    
    m.boundaries.update_boundary_array('Campaspe River', [reach_data, seg_dict])
    
