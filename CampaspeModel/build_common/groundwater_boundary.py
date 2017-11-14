def prepare_ghb_boundary_from_Murray_data(ModelBuilderObject,
                                          mriver_seg_ghb):
    MBO = ModelBuilderObject
    MBO.parameters.create_model_parameter('mghb_stage', value=0.01)
    MBO.parameters.parameter_options('mghb_stage', 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=-20.0, 
                                          PARUBND=50, 
                                          PARGP='ghb', 
                                          SCALE=1, 
                                          OFFSET=0)
    MBO.parameters.create_model_parameter('mghbk', value=10)
    MBO.parameters.parameter_options('mghbk', 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=1E-8, 
                                          PARUBND=50, 
                                          PARGP='ghb', 
                                          SCALE=1, 
                                          OFFSET=0)
    
    # First find which cells should make up the boundary based on the mapping 
    # from the Murray river polyline to the grid
    
    MurrayGHB = []
    Active_MurrayGHB_cells = []
    Murray_df_ind = []
    checked = []
    for mrow in mriver_seg_ghb.iterrows():
        ind = mrow[0]
        mrow = mrow[1]
        row = int(mrow['i'])
        col = int(mrow['j'])
        for lay in range(MBO.model_mesh3D[1].shape[0]):    
            if [lay, row, col] in checked:
                continue
            checked += [lay, row, col]
            if MBO.model_mesh3D[1][0][row][col] == -1:
                continue
            MurrayGHBstage = mrow['stage'] + MBO.parameters.param['mghb_stage']['PARVAL1']
            #if MurrayGHBstage < SS_model.model_mesh3D[0][lay][row][col]:
            #    continue
            if lay <= mrow['k']:
                continue
    
            Murray_df_ind += [ind]        
            Active_MurrayGHB_cells += [[lay, row, col]]
    
    # Now make sure that no cells are being caught surrounded by other GHB cells to prevent short circuiting
    def active_check(check_param, check_val, ref_cell, zone_cell, Active_MurrayGHB_cells):
        """
        Check target cell if it meets criteria to mark it as active
    
        :param check_param: int, target parameter to check
        :param check_val: int, value to check against, target parameter should not equal this
        :param ref_cell: list, reference to neighbouring cell ([layer, row, column])
        :param zone_cell: list, reference to HGU of cell with -1 indicating inactive cell 
        :param Active_MurrayGHB_cells: list, list of cell locations ([layer, row, column])
    
        :returns: bool, cell is active or not active
        """
        if (check_param != check_val) and \
                (zone_cell != -1) and \
                (ref_cell not in Active_MurrayGHB_cells):
            return True
    
        return False
    # End active_check()
    
    Murray_df_ind2 = []
    Final_MurrayGHB_cells = []
    zone = MBO.model_mesh3D[1]
    shape = zone.shape
    lay, row, col = 0, 1, 2
    for index, ac in enumerate(Active_MurrayGHB_cells):
        # check if active GHB cell has any active non GHB cells N,E,S,W, above or below
        active_non_GHB = False
        acl, acr, acc = int(ac[lay]), int(ac[row]), int(ac[col])
    
        # Check north:
        ref_cell = [acl, acr + 1, acc]
        active_non_GHB = active_check(acr, 0, ref_cell, zone[acl, acr + 1, acc], Active_MurrayGHB_cells)
    
        # Check east:
        if not active_non_GHB:
            ref_cell = [acl, acr, acc + 1]
            active_non_GHB = active_check(acc, shape[col] - 1, ref_cell, zone[acl, acr, acc + 1], Active_MurrayGHB_cells)
    
        # Check south:
        if not active_non_GHB:
            ref_cell = [acl, acr - 1, acc]
            active_non_GHB = active_check(acr, shape[row] - 1, ref_cell, zone[acl, acr - 1, acc], Active_MurrayGHB_cells)
    
        # Check west:
        if not active_non_GHB:
            ref_cell = [acl, acr, acc - 1]
            active_non_GHB = active_check(acc, 0, ref_cell, zone[acl, acr, acc - 1], Active_MurrayGHB_cells)
    
        if active_non_GHB:
            Final_MurrayGHB_cells += [ac]
            Murray_df_ind2 += [Murray_df_ind[index]]
    
    #for MurrayGHB_cell in SS_model.polyline_mapped['River_Murray_model.shp']:
    for index, MurrayGHB_cell in enumerate(Final_MurrayGHB_cells):
    
        lay = MurrayGHB_cell[0]
        row = MurrayGHB_cell[1]
        col = MurrayGHB_cell[2]
            
        MurrayGHBstage = mriver_seg_ghb['stage'].loc[Murray_df_ind2[index]] + MBO.parameters.param['mghb_stage']['PARVAL1']
        #if MurrayGHBstage < SS_model.model_mesh3D[0][0][row][col]:
        #    continue
        dx = MBO.gridHeight
        dz = MBO.model_mesh3D[0][lay][row][col] - MBO.model_mesh3D[0][lay + 1][row][col]
        MGHBconductance = dx * dz * MBO.parameters.param['mghbk']['PARVAL1']
        MurrayGHB += [[lay, row, col, MurrayGHBstage, MGHBconductance]]
    
    ghb = {}
    ghb[0] = MurrayGHB

    return ghb