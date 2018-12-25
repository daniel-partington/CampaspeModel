import numpy as np


def build_mesh_and_set_properties(model_builder_object,
                                  hu_raster_path,
                                  HGU_props,
                                  resolution=5000,
                                  create_basement=True,
                                  pilot_points_YX=False,
                                  verbose=True,
                                  hu_raster_files=None,
                                  hu_to_process=['qa', 'utb', 'utqa', 'utam', 
                                                 'utaf', 'lta', 'bse']):

    mbo = model_builder_object
    # Define the grid width and grid height for the model mesh which is stored as a multipolygon shapefile GDAL object
    print "************************************************************************"
    print " Defining structured mesh"
    resolution = resolution
    mbo.define_structured_mesh(resolution, resolution)

    # Read in hydrostratigraphic raster info for layer elevations:

    # Build basement file ... only need to do this once as it is time consuming so commented out for future runs
    if create_basement:
        mbo.create_basement_bottom(hu_raster_path, "sur_1t", "bse_1t", "bse_2b", hu_raster_path)
    # end if

    if not hu_raster_files:
        hu_raster_files = []
        for hu in hu_to_process:
            if hu != 'bse':
                hu_raster_files += ["{}_1t".format(hu), "{}_2b".format(hu)]
            else:
                hu_raster_files += ["{}_1t".format(hu), "{}_2b.tif".format(hu)]
            # end if
        # end for
        
    # This loads in the raster files and transforms them into the correct coordinate
    # sytstem.
    mbo.read_rasters(hu_raster_files, path=hu_raster_path)

    hu_raster_files_reproj = [x + "_reproj.tif" for x in hu_raster_files]

    # Map HGU's to grid
    print "************************************************************************"
    print " Mapping hydrogeological unit (HGU) rasters to grid "

    # Build 3D grid
    model_grid_raster_files = [x + "_model_grid.tif" for x in hu_raster_files]
    mbo.map_rasters_to_grid(hu_raster_files, hu_raster_path)

    # First two arguments of next function are arbitrary and not used ... need to rework module
    print "************************************************************************"
    print " Building 3D mesh based on HGU rasters"
    mbo.build_3D_mesh_from_rasters(model_grid_raster_files,
                                   mbo.out_data_folder_grid, minimum_thickness=1.0, maximum_thickness=1000.0)

    # Cleanup any isolated cells:
    mbo.reclassIsolatedCells(assimilate=True)

    print "************************************************************************"
    print " Assign properties to mesh based on pilot points and zonal information"

    # create list of HGU's from hu_raster_files
    HGU = list(set([x.split('_')[0] for x in hu_raster_files]))

    # NOTE *** utam is mapping to Shepparton Sands but it represents the
    # Loxton-Parilla Sand ... the HGU database needs updating to include this.
    HGU_map = {'bse': 'Bedrock',
               'utb': 'Newer Volcanics Basalts',
               'utaf': 'Calivil',
               'utam': 'Shepparton Sands',
               'lta': 'Renmark',
               'qa': 'Coonambidgal Sands',
               'utqa': 'Shepparton Sands'}

    HGU_zone = {'qa': 0,
                'utb': 1,
                'utqa': 2,
                'utam': 3,
                'utaf': 4,
                'lta': 5,
                'bse': 6}

    zone_HGU = {HGU_zone[x]: x for x in HGU_zone.keys()}

    pilot_points = True

    if pilot_points:
        # Set up pilot points:
        pilot_point_groups = ['hk', 'ss', 'sy']
        pp_group_dict = {}
        # Create some references to data inside the model builder object
        mesh_array = mbo.model_mesh3D
        cell_centers = mbo.model_mesh_centroids
        model_boundary = mbo.model_boundary
        # Define list of active zones as denoted by unique positive integers in
        # model mesh zone array
        zones_active = np.unique(mbo.model_mesh3D[1]).astype(int).tolist()
        zones_active.remove(-1)
        
        for pilot_points_group in pilot_point_groups:
            mbo.create_pilot_points(pilot_points_group)
            pp_group_dict[pilot_points_group] = mbo.pilot_points[pilot_points_group]
            # create alias for brevity ...
            pp_grp = pp_group_dict[pilot_points_group]
            # Create dict of zones and properties
            zone_prop_dict = {zone: HGU_props['Kh mean'][HGU_map[HGU[index]]] for index, zone in enumerate(zones_active)}
            # Define some parameters for pilot point distribution
            def reduce_list(lst, indices):
                return[lst[x - 1] for x in indices]
            
            if resolution == 1000:
                skip = reduce_list([3, 0, 6, 0, 6, 6, 6], zones_active)
                skip_active = reduce_list([3, 20, 0, 34, 0, 0, 0], zones_active)
            elif resolution == 500:
                skip = reduce_list([0, 0, 12, 0, 12, 12, 12], zones_active)
                skip_active = reduce_list([100, 40, 0, 70, 0, 0, 0], zones_active)
            elif resolution == 100:
                skip = reduce_list([0, 0, 120, 0, 120, 120, 120], zones_active)
                skip_active = reduce_list([1000, 400, 0, 700, 0, 0, 0], zones_active)
            elif resolution == 20000:
                skip = reduce_list([0, 0, 0, 0, 0, 0, 0], zones_active)
                skip_active = reduce_list([0, 0, 0, 0, 0, 0, 0], zones_active)
            else:
                skip = reduce_list([0,  0, 3, 0, 2, 3, 3], zones_active)
                skip_active = reduce_list([3, 20, 0, 4, 0, 0, 0], zones_active)

            # Generate the pilot points
            points_dict, points_zone_dict, points_val_dict, zone_mask2D, zone_active = \
                pp_grp.generate_points_from_mesh(mesh_array, cell_centers,
                                                 skip=skip,
                                                 skip_active=skip_active,
                                                 zone_prop_dict=zone_prop_dict,
                                                 zones=zones_active)
            # Create some necessary files from pilot points utilities
            pp_grp.write_settings_fig()
            pp_grp.write_grid_spec(mesh_array, model_boundary, delc=resolution, delr=resolution)
            pp_grp.write_struct_file(mesh_array, nugget=0.0,
                                     transform='log', numvariogram=1, variogram=0.15,
                                     vartype=2, bearing=0.0, a=20000.0, anisotropy=1.0)

            # These search_radius values have been tested on the 1000m grid, would be good
            # to add in other resolution lists as they are developed
            if resolution == 1000:
                search_radius = reduce_list([40000, 20000, 15000, 20000, 20000, 20000, 20000], zones_active)
            else:
                search_radius = reduce_list([30000, 20000, 40000, 20000, 40000, 50000, 20000], zones_active)

            prefixes = ['{}_{}'.format(pilot_points_group, zone_HGU[x - 1]) for x in zones_active]
            pp_grp.setup_pilot_points_by_zones(mesh_array, zones_active, search_radius, prefixes=prefixes)

            pp_grp.generate_cov_mat_by_zones(zones_active)

            #print("Running pyfac2real")
            import time
            time.sleep(3)
            pp_grp.run_pyfac2real_by_zones(zones_active)

        mbo.save_pilot_points()
        hk = pp_group_dict['hk']
        ss = pp_group_dict['ss']
        sy = pp_group_dict['sy']

    for index, unit in enumerate(HGU):
        if HGU_zone[unit] + 1 in zones_active:
            if pilot_points and not pilot_points_YX:
                mbo.parameters.create_model_parameter_set('kh' + unit,
                                                          value=HGU_props['Kh mean'][HGU_map[unit]],
                                                          num_parameters=hk.num_ppoints_by_zone[HGU_zone[unit] + 1])
                mbo.parameters.parameter_options_set('kh' + unit,
                                                     PARTRANS='log',
                                                     PARCHGLIM='factor',
                                                     PARLBND=HGU_props['Kh mean'][HGU_map[unit]] / 10.,
                                                     PARUBND=HGU_props['Kh mean'][HGU_map[unit]] * 10.,
                                                     PARGP='k' + unit,
                                                     SCALE=1,
                                                     OFFSET=0)
                mbo.parameters.create_model_parameter_set('sy' + unit,
                                                          value=HGU_props['Sy mean'][HGU_map[unit]],
                                                          num_parameters=ss.num_ppoints_by_zone[HGU_zone[unit] + 1])
                mbo.parameters.parameter_options_set('sy' + unit,
                                                     PARTRANS='log',
                                                     PARCHGLIM='factor',
                                                     PARLBND=1.0E-2,
                                                     PARUBND=0.44,
                                                     PARGP='sy' + unit,
                                                     SCALE=1,
                                                     OFFSET=0)
                mbo.parameters.create_model_parameter_set('ss' + unit,
                                                          value=HGU_props['SS mean'][HGU_map[unit]],
                                                          num_parameters=sy.num_ppoints_by_zone[HGU_zone[unit] + 1])
                mbo.parameters.parameter_options_set('ss' + unit,
                                                     PARTRANS='log',
                                                     PARCHGLIM='factor',
                                                     PARLBND=1E-6, #HGU_props['SS mean'][HGU_map[unit]] / 10.,
                                                     PARUBND=1E-4, #HGU_props['SS mean'][HGU_map[unit]] * 10.,
                                                     PARGP='ss' + unit,
                                                     SCALE=1,
                                                     OFFSET=0)
            elif pilot_points_YX:
                mbo.parameters.create_model_parameter_set('kh' + unit,
                                                          value=HGU_props['Kh mean'][HGU_map[unit]],
                                                          num_parameters=hk.num_ppoints_by_zone[HGU_zone[unit] + 1])
                mbo.parameters.parameter_options_set('kh' + unit,
                                                     PARTRANS='log',
                                                     PARCHGLIM='factor',
                                                     PARLBND=HGU_props['Kh mean'][HGU_map[unit]] / 10.,
                                                     PARUBND=HGU_props['Kh mean'][HGU_map[unit]] * 10.,
                                                     PARGP='k' + unit,
                                                     SCALE=1,
                                                     OFFSET=0)
                mbo.parameters.create_model_parameter_set('sy' + unit,
                                                          value=HGU_props['Sy mean'][HGU_map[unit]],
                                                          num_parameters=ss.num_ppoints_by_zone[HGU_zone[unit] + 1])
                mbo.parameters.parameter_options_set('sy' + unit,
                                                     PARTRANS='fixed',
                                                     PARCHGLIM='factor',
                                                     PARLBND=1.0E-3,
                                                     PARUBND=0.8,
                                                     PARGP='sy' + unit,
                                                     SCALE=1,
                                                     OFFSET=0)
                mbo.parameters.create_model_parameter_set('ss' + unit,
                                                          value=HGU_props['SS mean'][HGU_map[unit]],
                                                          num_parameters=sy.num_ppoints_by_zone[HGU_zone[unit] + 1])
                mbo.parameters.parameter_options_set('ss' + unit,
                                                     PARTRANS='fixed',
                                                     PARCHGLIM='factor',
                                                     PARLBND=HGU_props['SS mean'][HGU_map[unit]] / 10.,
                                                     PARUBND=HGU_props['SS mean'][HGU_map[unit]] * 10.,
                                                     PARGP='ss' + unit,
                                                     SCALE=1,
                                                     OFFSET=0)
    
            else:
                mbo.parameters.create_model_parameter('kh' + unit,
                                                      value=HGU_props['Kh mean'][HGU_map[unit]])
                mbo.parameters.parameter_options('kh' + unit,
                                                 PARTRANS='log',
                                                 PARCHGLIM='factor',
                                                 PARLBND=HGU_props['Kh mean'][HGU_map[unit]] / 10.,
                                                 PARUBND=HGU_props['Kh mean'][HGU_map[unit]] * 10.,
                                                 PARGP='k' + unit,
                                                 SCALE=1,
                                                 OFFSET=0)
                mbo.parameters.create_model_parameter('sy' + unit, value=HGU_props['Sy mean'][HGU_map[unit]])
                mbo.parameters.parameter_options('sy' + unit,
                                                 PARTRANS='log',
                                                 PARCHGLIM='factor',
                                                 PARLBND=1.0E-3,
                                                 PARUBND=0.8,
                                                 PARGP='sy' + unit,
                                                 SCALE=1,
                                                 OFFSET=0)
                mbo.parameters.create_model_parameter('ss' + unit, value=HGU_props['SS mean'][HGU_map[unit]])
                mbo.parameters.parameter_options('ss' + unit,
                                                 PARTRANS='log',
                                                 PARCHGLIM='factor',
                                                 PARLBND=HGU_props['SS mean'][HGU_map[unit]] / 10.,
                                                 PARUBND=HGU_props['SS mean'][HGU_map[unit]] * 10.,
                                                 PARGP='ss' + unit,
                                                 SCALE=1,
                                                 OFFSET=0)
    
    #        mbo.parameters.create_model_parameter('kv' + unit, value=HGU_props['Kz mean'][HGU_map[unit]])
    #        mbo.parameters.parameter_options('kv' + unit,
    #                                         PARTRANS='log',
    #                                         PARCHGLIM='factor',
    #                                         PARLBND=HGU_props['Kz mean'][HGU_map[unit]] / 10.,
    #                                         PARUBND=HGU_props['Kz mean'][HGU_map[unit]] * 10.,
    #                                         PARGP='k' + unit,
    #                                         SCALE=1,
    #                                         OFFSET=0)
    
            mbo.parameters.create_model_parameter('kv' + unit, value=0.1)
            mbo.parameters.parameter_options('kv' + unit,
                                             PARTRANS='log',
                                             PARCHGLIM='factor',
                                             PARLBND=1E-2,
                                             PARUBND=1.0,
                                             PARGP='k' + unit,
                                             SCALE=1,
                                             OFFSET=0)
        
    # This needs to be automatically generated from with the map_raster2mesh routine ...
    zone_map = {1: 'qa', 2: 'utb', 3: 'utqa', 4: 'utam', 5: 'utaf', 6: 'lta', 7: 'bse'}

    Kh = mbo.model_mesh3D[1].astype(float)
    Kv = mbo.model_mesh3D[1].astype(float)
    Sy = mbo.model_mesh3D[1].astype(float)
    SS = mbo.model_mesh3D[1].astype(float)
    for key in zones_active:
        if key not in np.unique(mbo.model_mesh3D[1]):
            continue
        if not pilot_points:
            Kh[Kh == key] = mbo.parameters.param['kh' + zone_map[key]]['PARVAL1']
            Sy[Sy == key] = mbo.parameters.param['sy' + zone_map[key]]['PARVAL1']
            SS[SS == key] = mbo.parameters.param['ss' + zone_map[key]]['PARVAL1']
        Kv[Kv == key] = Kh[Kv == key] * mbo.parameters.param['kv' + zone_map[key]]['PARVAL1']

    if pilot_points:
        Kh = hk.val_array
        Kh[mbo.model_mesh3D[1] == -1] = 0.0
        for key in zones_active:
            Kv[Kv == key] = Kh[Kv == key] * mbo.parameters.param['kv' + zone_map[key]]['PARVAL1']
        Sy = sy.val_array
        SS = ss.val_array

    mbo.properties.assign_model_properties('Kh', Kh)
    mbo.properties.assign_model_properties('Kv', Kv)
    mbo.properties.assign_model_properties('Sy', Sy)
    mbo.properties.assign_model_properties('SS', SS)

    return HGU, hu_raster_files_reproj
