import os


def process_line(line):
    return [x.strip() for x in line.split(':')[1].strip().split(',')]
# end process_line


def update_recharge(vals, interp_rain, rch_zones, recharge_zone_array, rch_zone_dict, mesh_1):
    for key in list(interp_rain.keys()):
        for i in range(rch_zones - 1):
            matching_zone = recharge_zone_array == rch_zone_dict[i + 1]
            interp_rain[key][matching_zone] = interp_rain[key][matching_zone] * vals[i]
        # End for

        matching_zone = recharge_zone_array == rch_zone_dict[0]
        interp_rain[key][matching_zone] = interp_rain[key][matching_zone] * 0.0
        interp_rain[key][mesh_1[0] == -1] = 0.0
    # End for

    return interp_rain
# End update_recharge()


def create_pp_points_dict(zone_map, zone, prop_array, prop_name, m, use_alt_vals=False):
    points_values_dict = {}
    if use_alt_vals is False:
        for index, key in enumerate(zone_map.keys()):
            for index2, param in enumerate(m.parameters.param_set[prop_name + zone_map[key]]):
                if index2 == 0:
                    points_values_dict[index] = [m.parameters.param[param]['PARVAL1']]
                else:
                    points_values_dict[index] += [m.parameters.param[param]['PARVAL1']]
                # End if
            # End for
        # End for
    else:
        # Hacky model param changes
        alt_k_vals = {'khutqa': 1.0,
                      'khutb': 0.1,  # 1.,
                      'khqa': 20.,
                      'khutam': 10.,  # 0.1,
                      'khutaf': 50.,  # 170
                      'khlta': 50.,  # 170
                      'khbse': 0.05}

        alt_ss_vals = {'ssutqa': 1E-5,
                       'ssutb': 1E-5,
                       'ssqa': 1E-5,
                       'ssutam': 1E-5,
                       'ssutaf': 1E-3,
                       'sslta': 1E-3,
                       'ssbse': 1E-5}

        alt_sy_vals = {'syutqa': 0.10,
                       'syutb': 0.08,
                       'syqa': 0.22,
                       'syutam': 0.25,
                       'syutaf': 0.25,
                       'sylta': 0.25,
                       'sybse': 0.0009}

        k_factor = 1.0  # 2.5
        sy_factor = 1  # 2.5
        for key in alt_k_vals:
            alt_k_vals[key] = alt_k_vals[key] * k_factor
        for key in alt_sy_vals:
            alt_sy_vals[key] = alt_sy_vals[key] * sy_factor

        for index, key in enumerate(zone_map.keys()):
            for index2, param in enumerate(m.parameters.param_set[prop_name + zone_map[key]]):
                if index2 == 0:
                    if prop_name == 'ss':
                        m.parameters.param[param]['PARVAL1'] = alt_ss_vals[prop_name + zone_map[key]]
                    elif prop_name == 'sy':
                        m.parameters.param[param]['PARVAL1'] = alt_sy_vals[prop_name + zone_map[key]]
                    elif prop_name == 'kh':
                        m.parameters.param[param]['PARVAL1'] = alt_k_vals[prop_name + zone_map[key]]
                    # End if
                    points_values_dict[index] = [m.parameters.param[param]['PARVAL1']]
                else:
                    if prop_name == 'ss':
                        m.parameters.param[param]['PARVAL1'] += alt_ss_vals[prop_name + zone_map[key]]
                    elif prop_name == 'sy':
                        m.parameters.param[param]['PARVAL1'] += alt_sy_vals[prop_name + zone_map[key]]
                    elif prop_name == 'kh':
                        m.parameters.param[param]['PARVAL1'] += alt_k_vals[prop_name + zone_map[key]]
                    # End if
                    points_values_dict[index] += [m.parameters.param[param]['PARVAL1']]
                # End if
            # End for
        # End for
    return points_values_dict
# End create_pp_points_dict()


def update_pilot_points(zone_map, zone, prop_array, par_name, prop_name, prop_folder, m,
                        prop_array_fname, model_folder, use_alt_vals):
    points_values_dict = create_pp_points_dict(zone_map, zone, prop_array, prop_name, m, use_alt_vals=use_alt_vals)
    p = m.pilot_points[par_name]
    zones = len(list(zone_map.keys()))
    p.output_directory = os.path.join(model_folder, prop_folder)
    p.update_pilot_points_files_by_zones(zones, points_values_dict)
    p.run_pyfac2real_by_zones(zones)
    return p.val_array
# End update_pilot_points()


def update_campaspe_pilot_points(model, model_folder, use_alt_vals=False):
    import numpy as np
    # This needs to be automatically generated with the map_raster2mesh routine ...
    zone_map = {1: 'qa', 2: 'utb', 3: 'utqa', 4: 'utam', 5: 'utaf', 6: 'lta', 7: 'bse'}

    default_array = model.model_mesh3D[1].astype(float)
    zone = np.copy(default_array)
    kh = np.copy(default_array)  # horizonal and vertical hydraulic conductivity
    kv = np.copy(default_array)
    sy = np.copy(default_array)  # Specific yield and storage
    ss = np.copy(default_array)

    kh = update_pilot_points(zone_map, zone, kh, 'hk', 'kh', 'hk_pilot_points',
                             model, 'hk_val_array', model_folder, use_alt_vals)
    kh[kh > 10000.0] = 25.0
    kv = kh * 0.1

    sy = update_pilot_points(zone_map, zone, sy, 'sy', 'sy', 'sy_pilot_points',
                             model, 'sy_val_array', model_folder, use_alt_vals)
    sy[sy > 0.5] = 0.5

    ss = update_pilot_points(zone_map, zone, ss, 'ss', 'ss', 'ss_pilot_points',
                             model, 'ss_val_array', model_folder, use_alt_vals)

    model.properties.update_model_properties('Kh', kh)
    model.properties.update_model_properties('Kv', kv)
    model.properties.update_model_properties('Sy', sy)
    model.properties.update_model_properties('SS', ss)

# End update_campaspe_pilot_points()


def load_obj(filename):
    try:
        import pickle as pickle
    except ImportError:
        import pickle

    if filename[-4:] == '.pkl':
        with open(filename, 'r') as f:
            return pickle.load(f)
    else:
        print('File type not recognised as "pkl"')
    # End if
# End load_obj()


def load_model_config():
    from HydroModelBuilder.Utilities.Config.ConfigLoader import ConfigLoader

    this_file_loc = os.path.realpath(__file__)
    pkg_path = this_file_loc[0:this_file_loc.index('models')]
    config_path = os.path.join(pkg_path, "config", "model_config.json")
    CONFIG = ConfigLoader(config_path).set_environment("GW_link_Integrated")

    return CONFIG
# End load_model_config()
