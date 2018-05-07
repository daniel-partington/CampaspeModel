import os


def prepare_channel_data_for_model(ModelBuilderObject,
                                   start_irrigation,
                                   date_index,
                                   Camp_riv_cells):

    MBO = ModelBuilderObject
    channel_poly = MBO.read_poly("Channel_Clip.shp",
                                 path=os.path.join(MBO.data_folder, "SW") + os.path.sep)

    MBO.map_polyline_to_grid(channel_poly)

    MBO.parameters.create_model_parameter('chndrp', value=0.01)
    MBO.parameters.parameter_options('chndrp',
                                     PARTRANS='log',
                                     PARCHGLIM='factor',
                                     PARLBND=0.001,
                                     PARUBND=0.1,
                                     PARGP='channl',
                                     SCALE=1,
                                     OFFSET=0)
    MBO.parameters.create_model_parameter('kv_ch', value=5E-3)
    MBO.parameters.parameter_options('kv_ch',
                                     PARTRANS='log',
                                     PARCHGLIM='factor',
                                     PARLBND=1E-8,
                                     PARUBND=20,
                                     PARGP='channl',
                                     SCALE=1,
                                     OFFSET=0)

    simple_channel = []
    channel_width_avg = 5.0  # m
    channel_bed_thickness = 0.10  # m
    for channel_cell in MBO.polyline_mapped['Channel_Clip_model.shp']:
        row = channel_cell[0][0]
        col = channel_cell[0][1]
        if MBO.model_mesh3D[1][0][row][col] == -1:
            continue
        if (row, col) in Camp_riv_cells:
            continue
        channel_stage = MBO.model_mesh3D[0][0][row][col]
        channel_bed = MBO.model_mesh3D[0][0][row][col] - MBO.parameters.param['chndrp']['PARVAL1']
        channel_cond = channel_cell[1] * channel_width_avg * \
            MBO.parameters.param['kv_ch']['PARVAL1'] / channel_bed_thickness
        simple_channel += [[0, row, col, channel_stage, channel_cond, channel_bed]]

    channel_start = findInterval(start_irrigation, date_index)

    channel = {}
    channel[channel_start + 1] = simple_channel
