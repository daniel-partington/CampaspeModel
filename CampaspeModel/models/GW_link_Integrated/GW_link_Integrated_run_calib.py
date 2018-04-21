"""
Run Campaspe GW model using inputs from Farm model and SW model and
then return SW/GW exchanges, avg depth to GW, depth to GW at ecology sites and
head at trigger bores.
"""

import cPickle as pickle
import os
import sys
import time
import warnings

import flopy.utils.binaryfile as bf
import numpy as np

from common_run_funcs import *
from HydroModelBuilder.GWModelManager import GWModelManager
from HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface


# TODO: Set the stream gauges, ecology bores, policy bores at the start in some
# other class or in here but so they are available in the run function.


def process_line(line):
    return [x.strip() for x in line.split(':')[1].strip().split(',')]
# end process_line


def run(model_folder, data_folder, mf_exe_folder, farm_zones=None, param_file=None, riv_stages=None,
        rainfall_irrigation=None, pumping=None, verbose=True, MM=None, recharge_amt=0.03, is_steady=False, cache={}):
    """
    GW Model Runner.

    :param model_folder: str, path to model folder
    :param data_folder: str, path to data
    :param mf_exe_folder: str, path to MODFLOW executable
    :param farm_zones: list[str], farm zone IDs to map ground water level values to (in order)
    :param param_file: str, path to parameter file
    :param riv_stages: np.recarray, gauge numbers and stage
    :param rainfall_irrigation: np.ndarray, array representing rainfall and irrigation input.
                                Must match the model extent.
    :param pumping: float, daily pumping amount in m^3/day (ML/day to m^3/day => ML * 1000)

    :returns: tuple[np.recarray], four elements
              xchange: exchange for each gauge by gauge ID
              avg_gw_depth: average depth for each zone
              ecol_depth_to_gw: average gw depth at each cell that contains an ecological bore of interest.
              trigger_head: Policy trigger well heads
    """
    p_j = os.path.join

    # DEBUG: Setting constant high pumping rate to see if gw levels go down
    # pumping = 1.0

    if MM is None:
        MM = GWModelManager()
        MM.load_GW_model(p_j(model_folder, "GW_link_Integrated_packaged.pkl"))
    # End if

    # Complimentary models requirements, i.e. bore and gauge data that should be
    # referenceable to this model for parsing specific outputs and receiving inputs:
    if not hasattr(MM, 'ext_linkage_bores') or MM.ext_linkage_bores is None:

        # Read in bores that relate to external models
        model_linking = p_j(data_folder, "model_linking.csv")

        # Need to check given data folder and its parent directory
        # Dirty hack, I know :(
        if not os.path.exists(model_linking):
            model_linking = p_j(data_folder, "..", "model_linking.csv")
            if not os.path.exists(model_linking):
                raise IOError("Could not find bore linkages information (`model_linking.csv`)")
            # End if
        # End if

        with open(model_linking, 'r') as f:
            lines = f.readlines()

            for line in lines:
                if line.split(':')[0] == 'Ecology':
                    Ecology_bores = process_line(line)
                elif line.split(':')[0] == 'Policy':
                    Policy_bores = process_line(line)
                elif line.split(':')[0] == 'SW_stream_gauges':
                    Stream_gauges = process_line(line)
            # End for
        # End with

        MM.ext_linkage_bores = {
            "Ecology_bores": Ecology_bores,
            "Policy_bores": Policy_bores,
            "Stream_gauges": Stream_gauges
        }

        cache['Ecology_bores'] = Ecology_bores
        cache['Policy_bores'] = Policy_bores
        cache['Stream_gauges'] = Stream_gauges
    else:
        Ecology_bores = cache['Ecology_bores']
        Policy_bores = cache['Policy_bores']
        Stream_gauges = cache['Stream_gauges']
    # End if

    name = MM.GW_build.keys()[0]
    this_model = MM.GW_build[name]
    mesh = this_model.model_mesh3D
    mesh_0, mesh_1 = mesh[0], mesh[1]

    # Load in the new parameters based on parameters.txt or dictionary of new parameters
    if param_file:
        if 'parameters.txt' not in param_file:
            param_file = p_j(param_file, 'parameters.txt')
        # End if
        this_model.updateModelParameters(param_file, verbose=verbose)
    # End if

    model_params = this_model.parameters.param
    model_boundaries = this_model.boundaries
    model_boundaries_bc = model_boundaries.bc

    #+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_
    #_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_
    #+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_
    #_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_
    #+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_

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

    ghb_k_factor = 1.0  # 10.
    riv_factor = 0.03

#    red_red = 0.75
#    non_irrig_red = 0.2 * red_red
#    irrig_red = 0.50 * red_red

    if is_steady:
        red_red = 1.
        non_irrig_red = 1. * red_red
        irrig_red = 1. * red_red
    else:
        red_red = 1.
        non_irrig_red = 1. * red_red
        irrig_red = 1. * red_red

    this_model.parameters.param['mghbst']['PARVAL1'] = -10.

    #+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_
    #_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_
    #+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_
    #_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_
    #+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_

    if verbose:
        print "************************************************************************"
        print " Updating HGU parameters "

    kh = this_model.properties.properties['Kh']

    if verbose:
        print "************************************************************************"
        print " Updating Campaspe River boundary"

    river_seg = this_model.river_mapping['Campaspe']
    num_reaches = this_model.pilot_points['Campaspe'].num_points
    known_points = this_model.pilot_points['Campaspe'].points

    for x in xrange(num_reaches):
        model_params['kv_riv{}'.format(x)]['PARVAL1'] = model_params['kv_riv{}'.format(x)]['PARVAL1'] * riv_factor

    strcond_val = [model_params['kv_riv{}'.format(x)]['PARVAL1'] for x in xrange(num_reaches)]
    river_seg['strhc1'] = np.interp(river_seg['Cumulative Length'].values.tolist(), known_points, strcond_val)
    river_seg['multi'] = river_seg['strhc1'] * river_seg['rchlen'] * river_seg['width1']
    river_bc = this_model.boundaries.bc['Campaspe River']['bc_array']

    if not is_steady:
        for key in river_bc:
            riv_stress_period = river_bc[key]
            for idx, rc in enumerate(riv_stress_period):
                river_bc[key][idx] = [rc[0], rc[1], rc[2], rc[3], river_seg['multi'].tolist()[idx], rc[5]]
    else:
        mean_stage = {}
        for key in river_bc:
            riv_stress_period = river_bc[key]
            for idx, rc in enumerate(riv_stress_period):
                if key == 0:
                    mean_stage[idx] = rc[3]
                else:
                    mean_stage[idx] = mean_stage[idx] * (1. - 1. / (key + 1.)) + rc[3] / (key + 1.)

        for key in river_bc:
            riv_stress_period = river_bc[key]
            for idx, rc in enumerate(riv_stress_period):
                river_bc[key][idx] = [rc[0], rc[1], rc[2], mean_stage[idx], river_seg['multi'].tolist()[idx], rc[5]]

     #simple_river = river_seg[['k', 'i', 'j', 'stage', 'multi', 'strtop']].values.tolist()

    model_boundaries.update_boundary_array('Campaspe River', river_bc)

    if verbose:
        print "************************************************************************"
        print " Updating Murray River boundary"

    # Updating Murray River
    mriver_seg = this_model.river_mapping['Murray']
    mriver_seg['strhc1'] = model_params['kv_rm']['PARVAL1']
    mriver_seg['multi'] = mriver_seg['strhc1'] * mriver_seg['rchlen'] * mriver_seg['width1']

    mriver_bc = this_model.boundaries.bc['Murray River']['bc_array']
    for key in mriver_bc:
        mriv_stress_period = mriver_bc[key]
        for idx, rc in enumerate(mriv_stress_period):
            mriver_bc[key][idx] = [rc[0], rc[1], rc[2], rc[3], mriver_seg['multi'].tolist()[idx], rc[5]]
    #msimple_river = mriver_seg[['k', 'i', 'j', 'stage', 'multi', 'strtop']].values.tolist()

    model_boundaries.update_boundary_array('Murray River', mriver_bc)

    if verbose:
        print "************************************************************************"
        print " Updating recharge boundary "

    # Adjust rainfall to recharge using zoned rainfall reduction parameters
    # Need to make copies of all rainfall arrays
    interp_rain = model_boundaries_bc['Rainfall']['bc_array']
    for key in interp_rain.keys():
        interp_rain[key] = np.copy(interp_rain[key])

    if is_steady:
        interp_rain[0] = np.mean([interp_rain[x] for x in interp_rain], axis=0)

    recharge_zone_array = model_boundaries_bc['Rain_reduced']['zonal_array']
    rch_zone_dict = model_boundaries_bc['Rain_reduced']['zonal_dict']

    rch_zones = len(rch_zone_dict.keys())

    par_rech_vals = [model_params['rchred{}'.format(i)]['PARVAL1']
                     for i in xrange(rch_zones - 1)]
    par_rech_vals = [model_params['rchred{}'.format(i)]['PARLBND']
                     for i in xrange(rch_zones - 1)]

    def update_recharge(vals):
        for key in interp_rain.keys():
            for i in xrange(rch_zones - 1):
                #                interp_rain[key][recharge_zone_array == rch_zone_dict[i + 1]] = \
                #                    interp_rain[key][recharge_zone_array == rch_zone_dict[i + 1]] * \
                #                    vals[i]
                if i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]:
                    interp_rain[key][recharge_zone_array == rch_zone_dict[i + 1]] = \
                        interp_rain[key][recharge_zone_array == rch_zone_dict[i + 1]] * \
                        vals[i] * non_irrig_red
#                elif i in [15]:
#                    interp_rain[key][recharge_zone_array == rch_zone_dict[i + 1]] = \
#                        interp_rain[key][recharge_zone_array == rch_zone_dict[i + 1]] * \
#                        vals[i] * 0.
                else:
                    interp_rain[key][recharge_zone_array == rch_zone_dict[i + 1]] = \
                        interp_rain[key][recharge_zone_array == rch_zone_dict[i + 1]] * \
                        vals[i] * irrig_red

            interp_rain[key][recharge_zone_array == rch_zone_dict[0]] = \
                interp_rain[key][recharge_zone_array == rch_zone_dict[0]] * 0.0
            # Killing recharge on inactive cells
            interp_rain[key][this_model.model_mesh3D[1][0] == -1] = 0.
            # Killing recharge across the outcropping bedrock
            interp_rain[key][this_model.model_mesh3D[1][0] == 7] = 0.
        return interp_rain

    interp_rain = update_recharge(par_rech_vals)
    rch = interp_rain

    model_boundaries.update_boundary_array('Rain_reduced', rch)
    #model_boundaries.assign_boundary_array('Rain_reduced', {0: interp_rain})

    if verbose:
        print "************************************************************************"
        print " Updating pumping boundary "

    pumpy = model_boundaries_bc['licenced_wells']['bc_array']
    wel = {key: [[b[0], b[1], b[2], b[3] * pumping] for b in a] for key, a in pumpy.iteritems()}
    if is_steady:
        wel[0] = wel[wel.keys()[-1]]

    model_boundaries.assign_boundary_array('licenced_wells', wel)

    if verbose:
        print "************************************************************************"
        print " Updating GHB boundary "

    MurrayGHB = []

    MurrayGHB_cells = [[x[0], x[1], x[2], x[3]] for x in this_model.boundaries.bc['GHB']['bc_array'][0]]
    for MurrayGHB_cell in MurrayGHB_cells:
        lay, row, col = MurrayGHB_cell[:3]
        MurrayGHBstage = MurrayGHB_cell[3] + this_model.parameters.param['mghbst']['PARVAL1']
        if MurrayGHBstage < this_model.model_mesh3D[0][lay + 1][row][col]:
            continue
        # end except
        dx = this_model.gridHeight
        dz = this_model.model_mesh3D[0][lay][row][col] - \
            this_model.model_mesh3D[0][lay + 1][row][col]
        MGHBconductance = dx * dz * kh[lay][row][col] * ghb_k_factor  # m.parameters.param['mghbk']['PARVAL1']
        MurrayGHB += [[lay, row, col, MurrayGHBstage, MGHBconductance]]

    ghb = {}
    ghb[0] = MurrayGHB

    model_boundaries.assign_boundary_array('GHB', {0: MurrayGHB})

    if verbose:
        print "************************************************************************"
        print " Initialise head start values "

    if not is_steady:
        fname = 'model_{}'.format(name)
        try:
            headobj = bf.HeadFile(p_j(data_folder, fname, name + '.hds'))
            times = headobj.get_times()
            head = headobj.get_data(totim=times[-1])
            this_model.initial_conditions.set_as_initial_condition("Head", head)
        except IndexError:
            raise IndexError(
                "Corrupted MODFLOW hds file - check, replace, or clear {}".format(
                    p_j(data_folder, fname, name + '.hds')))
    #    except IOError:
    #        warnings.warn("MODFLOW hds file not found. Recreating head state from existing values.")
    #        head = np.stack([this_model.model_mesh3D[0][0] for i in range(7)], axis=0)
    #        this_model.initial_conditions.set_as_initial_condition("Head", head)
        # End try

    # Currently using flopyInterface directly rather than running from the ModelManager ...
    modflow_model = flopyInterface.ModflowModel(this_model, data_folder=p_j(data_folder, "model_{}".format(name)))

    # Override temporal aspects of model build:
    modflow_model.steady = is_steady  # This is to tell FloPy that is a transient model
    modflow_model.executable = mf_exe_folder

    modflow_model.buildMODFLOW()
    modflow_model.runMODFLOW()

    if not is_steady:
        #modflow_model.viewHeadsByZone2(1, head_name='policy_bores')
        modflow_model.compareAllObs('policy_bores')
        modflow_model.waterBalanceTS()
        modflow_model.viewHeads2()

        # TODO: Check that all of the wells listed were mapped to the model mesh and
    # are available for inspection
    # return swgw_exchanges, avg_depth_to_gw, ecol_depth_to_gw, trigger_heads
    return modflow_model
# End run()


def main():
    print("Running from: " + os.getcwd())

    CONFIG = load_model_config()

    # Example river level data (to be inputted from SW Model)
    fname = "initial_river_levels.pkl"
    riv_stages = load_obj(os.path.join(CONFIG.settings['data_folder'], fname))

    args = sys.argv
    if len(args) > 1:
        model_folder = sys.argv[1]
        data_folder = sys.argv[2]
        mf_exe_folder = sys.argv[3]

        if len(args) > 4:
            param_file = sys.argv[4]
    else:
        model_config = CONFIG.model_config
        model_folder = model_config['model_folder'] + model_config['grid_resolution'] + os.path.sep
        data_folder = os.path.join(model_config['data_folder'], 'hindcast')
        mf_exe_folder = model_config['mf_exe_folder']
        param_file = model_config['param_file']
    # End if

    model_folder = model_folder.replace("structured_model_grid_5000m", "hindcast/structured_model_grid_5000m")
    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder, r"GW_link_Integrated_packaged.pkl"))

    this_model = MM.GW_build[MM.name]
    update_campaspe_pilot_points(this_model, model_folder, use_alt_vals=True)

    # SS first
    run_params = {
        "model_folder": model_folder,
        "data_folder": data_folder,
        "mf_exe_folder": mf_exe_folder,
        "farm_zones": ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'],
        "param_file": param_file if param_file else None,
        "riv_stages": riv_stages,
        "rainfall_irrigation": None,
        "pumping": 0.0,  # m^3/day
        "MM": MM,
        "verbose": False,
        "is_steady": True
    }

    model_results_ss = run(**run_params)
    model_results_ss.model_data.write_parameters_to_file(os.path.join(
        data_folder.replace('hindcast', 'forecast'), 'ss_parameters.txt'))

    model_name = this_model.name
    model_dir = 'model_{}'.format(model_name)
    heads_file_loc = os.path.join(data_folder, model_dir, '{}.hds'.format(model_name))
    # dst_loc = os.path.join(data_folder.replace('hindcast', 'forecast'), model_dir, '{}.hds'.format(model_name))
    dst_loc = os.path.join(data_folder.replace('hindcast', ''), 'hindcast_initial.hds')

    from shutil import copyfile
    print('Copying {} to {}'.format(heads_file_loc, dst_loc))
    copyfile(heads_file_loc, dst_loc)

    # TR second
    MM = GWModelManager()
    MM.load_GW_model(os.path.join(model_folder, r"GW_link_Integrated_packaged.pkl"))
    run_params = {
        "model_folder": model_folder,
        "data_folder": data_folder,
        "mf_exe_folder": mf_exe_folder,
        "farm_zones": ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'],
        "param_file": param_file if param_file else None,
        "riv_stages": riv_stages,
        "rainfall_irrigation": None,
        "pumping": 1.5,  # m^3/day
        "MM": MM,
        "verbose": False,
        "is_steady": False
    }

    model_results = run(**run_params)
    model_results.model_data.write_parameters_to_file(os.path.join(
        data_folder.replace('hindcast', 'forecast'), 'tr_parameters.txt'))

    dst_loc = os.path.join(data_folder.replace('hindcast', ''), 'forecast_initial.hds')
    print('Copying {} to {}'.format(heads_file_loc, dst_loc))
    copyfile(heads_file_loc, dst_loc)

#    swgw_exchanges, avg_depth_to_gw, ecol_depth_to_gw, trigger_heads = run(**run_params)
#    print("swgw_exchanges", swgw_exchanges)
#    print("avg_depth_to_gw", avg_depth_to_gw)
#    print("ecol_depth_to_gw", ecol_depth_to_gw)
#    print("trigger_heads", trigger_heads)
    return model_results_ss, model_results


if __name__ == "__main__":
    model_results_ss, model_results = main()
