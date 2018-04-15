import os
import numpy as np
import pandas as pd

def prepare_drain_data_for_model(ModelBuilderObject,
                                 Camp_riv_cells,
                                 start_irrigation,
                                 date_index,
                                 pilot_points_YX=False):
    MBO = ModelBuilderObject
    drain_poly = MBO.read_poly("Drain_Clip.shp", path=os.path.join(MBO.data_folder, r"SW\\")) 
    MBO.map_polyline_to_grid(drain_poly)
    
    MBO.parameters.create_model_parameter('drndrp', value=0.01)
    MBO.parameters.parameter_options('drndrp', 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=0.001, 
                                          PARUBND=0.1, 
                                          PARGP='drain', 
                                          SCALE=1, 
                                          OFFSET=0)
    MBO.parameters.create_model_parameter('kv_drn', value=5E-3)
    MBO.parameters.parameter_options('kv_drn', 
                                          PARTRANS='log', 
                                          PARCHGLIM='factor', 
                                          PARLBND=1E-8, 
                                          PARUBND=20, 
                                          PARGP='drain', 
                                          SCALE=1, 
                                          OFFSET=0)
    
    simple_drain = []
    drain_width_avg = 3.0 #m
    drain_bed_thickness = 0.10 #m
    for drain_cell in MBO.polyline_mapped['Drain_Clip_model.shp']:
        row = drain_cell[0][0]
        col = drain_cell[0][1]
        if MBO.model_mesh3D[1][0][row][col] == -1:
            continue
        if (row, col) in Camp_riv_cells:
            continue
        #print tr_model.model_mesh3D
        drain_bed = MBO.model_mesh3D[0][0][row][col] - MBO.parameters.param['drndrp']['PARVAL1']
        drain_cond = drain_cell[1] * drain_width_avg * MBO.parameters.param['kv_drn']['PARVAL1'] / drain_bed_thickness
        simple_drain += [[0, row, col, drain_bed, drain_cond]]
    
    def findInterval(row, times):
        key_time = pd.to_datetime(row)
        lower_time = times[0]
        for period, time in enumerate(times):
            if period > 0:
                if lower_time <= key_time < time:
                    return period - 1
            lower_time = time
        return 0
    
    drain_start = findInterval(start_irrigation, date_index)
    drain = {}
    if not pilot_points_YX and len(date_index) > 1:
        drain[drain_start + 1] = simple_drain
    elif not pilot_points_YX and len(date_index) == 1:
        drain[drain_start] = simple_drain
    elif pilot_points_YX:
        drain[0] = simple_drain

    return drain