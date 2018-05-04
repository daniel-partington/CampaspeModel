"""
docker run -it \
     -v /Users/masher/Documents:/shared \
     -e PYTHONPATH=/shared/benjamin-holt:/shared/masher \
     -e CACHE_DIR=/shared/masher/too_big_to_git/cache\
     mjasher/heat \
    python /shared/masher/simple_campaspe/run.py

docker run -it \
-v /Users/masher/Documents:/shared \
mjasher/heat \
python /shared/modflow_heat_docker/too_big_to_git/heat/pyheat/models/heat_equation/kle.py


TODO
- kriging
- paramaterize all
- discontinuity


"""


"""
rainfall + irrigation R (36, 23)
sw stages (8)
gw pumping (12)
initial heads (36, 23, 7)
--->
sw-gw exchange (4)
policy bore level(2), ecology bore level (3)
average gw level in farm zone (12)


point -> grid (kriging)
grid -> point (take value at points)


1. parameterize each input
2. write a function for each output


"""



"""
Estimating groundwater recharge using an integrated surface and groundwater modelling approach
https://ac.els-cdn.com/002216949290216I/1-s2.0-002216949290216I-main.pdf?_tid=0aa9b6bc-a991-11e7-8446-00000aab0f02&acdnat=1507182787_22b46417dbd258ccb67d94b57489fd29
"""


import numpy as np
import os
import multiprocessing
import flopy


def extent():
    import tempfile
    import ogr
    import json
    from GIS import project
    # this is just for documentation
    with tempfile.NamedTemporaryFile() as tempfile:
        project.ogr2ogr(src_datasource_name='/shared/CampaspeModel/docs/GW_model_area.json',
            dst_datasource_name=tempfile.name,
            t_srs='EPSG:28355',
            f='GeoJSON'
        )
        with open(tempfile.name) as f:
            geojson = json.dumps(json.loads(f.read())['features'][0]['geometry'])
        geom = ogr.CreateGeometryFromJson(geojson)
        (minX, maxX, minY, maxY) = geom.GetEnvelope()
        print(minX, maxX, minY, maxY)


def locations(minX, maxX, minY, maxY, Lx, Ly, **kwargs):
    np.random.seed(989765)
    data = np.random.rand(10, 3)
    data[:, 0] = minX + Lx * data[:, 0]
    data[:, 1] = minY + Ly * data[:, 1]
    return data.T


def interpolation(layer_avgs, plot, minX, maxX, minY, maxY, Lx, Ly, ncol, nrow, delr, delc, **kwargs):
    # pip install pykrige
    from pykrige.ok import OrdinaryKriging

    xs, ys, zs = locations(minX, maxX, minY, maxY, Lx, Ly)

    gridx = np.linspace(minX + 0.5 * delr, maxX - 0.5 * delr, ncol)
    gridy = np.flip(np.linspace(minY + 0.5 * delc, maxY - 0.5 * delc, nrow), axis=0)

    # data = np.array([[0.3, 1.2, 0.47],
    #                  [1.9, 0.6, 0.56],
    #                  [1.1, 3.2, 0.74],
    #                  [3.3, 4.4, 1.47],
    #                  [4.7, 3.8, 1.74]])

    # gridx = np.arange(0.0, 5.5, 0.5)
    # gridy = np.arange(0.0, 5.5, 0.5)

    # Create the ordinary kriging object. Required inputs are the X-coordinates of
    # the data points, the Y-coordinates of the data points, and the Z-values of the
    # data points. If no variogram model is specified, defaults to a linear variogram
    # model. If no variogram model parameters are specified, then the code automatically
    # calculates the parameters by fitting the variogram model to the binned
    # experimental semivariogram. The verbose kwarg controls code talk-back, and
    # the enable_plotting kwarg controls the display of the semivariogram.
    OK = OrdinaryKriging(xs, ys, zs, variogram_model='linear',
                         verbose=False, enable_plotting=False)

    # Creates the kriged grid and the variance grid. Allows for kriging on a rectangular
    # grid of points, on a masked rectangular grid of points, or with arbitrary points.
    # (See OrdinaryKriging.__doc__ for more information.)
    z, ss = OK.execute('grid', gridx, gridy)

    if plot:
        import matplotlib.pyplot as plt
        z_min, z_max = np.percentile(zs, 20), np.percentile(zs, 80)
        print(np.percentile(z, np.linspace(0, 100, 10)))
        plt.imshow(z, extent=[minX, maxX, minY, maxY], vmin=z_min, vmax=z_max)
        plt.scatter(xs, ys, 100, zs, vmin=z_min, vmax=z_max)
        # plt.plot(data[:, 0], data[:, 1], '.')
        plt.colorbar()
        plt.savefig('/shared/plot/base/kriging.png')
        plt.close()

        plt.imshow(ss, extent=[minX, maxX, minY, maxY])
        plt.colorbar()
        plt.savefig('/shared/plot/base/kriging_var.png')
        plt.close()

    return [l_avg * (z + 0.5) for l_avg in layer_avgs]


def get_riv(hk, botm, nrow, ncol, delr, delc, riv_level):
    # [lay, row, col, stage, cond, rbot]
    # Make list for stress period 1
    # TODO river and well locations based on coordinates
    stage_murray = 3.
    stage_campaspe = riv_level *  1.
    rbot = -5.
    cond_murray = hk[0] * (stage_murray - botm[0]) * delc
    cond_campaspe = hk[0] * (stage_campaspe - botm[0]) * delc
    bound_sp1 = []
    for ic in range(ncol):
        bound_sp1.append([0, 0, ic, stage_murray, cond_murray[0, ic], rbot])
    for ir in range(nrow):
        bound_sp1.append([0, ir, ncol/2, stage_campaspe, cond_campaspe[ir, ncol/2], rbot])
    
    print('Adding ', len(bound_sp1), 'RIVs for stress period 1.')
    # print('BOTM', botm[:, 0, ncol/2], 'to', botm[:, nrow - 1, ncol/2])

    # Make list for stress period 2
    stage_murray = 3.
    stage_campaspe = 0.
    cond_murray = hk[0] * (stage_murray - botm[0]) * delc
    cond_campaspe = hk[0] * (stage_campaspe - botm[0]) * delc
    bound_sp2 = []
    for ic in range(ncol):
        bound_sp2.append([0, 0, ic, stage_murray, cond_murray[0, ic], rbot])
    for ir in range(nrow):
        bound_sp1.append([0, ir, ncol/2, stage_campaspe, cond_campaspe[ir, ncol/2], rbot])
    print('Adding ', len(bound_sp2), 'RIVs for stress period 2.')

    # We do not need to add a dictionary entry for stress period 3.
    # Flopy will automatically take the list from stress period 2 and apply it
    # to the end of the simulation
    return {'stress_period_data': {0: bound_sp1, 1: bound_sp2}}


def get_well(nrow, ncol, well_pumping):
    # Remember to use zero-based layer, row, column indices!
    # TODO -8000, -4000 what the
    ml_year_to_m3_day = 1. * 1000 / 365
    well_l_r_c = [1, nrow / 2 - 1, ncol / 2 - 1]
    wel_sp1 = [well_l_r_c + [well_pumping * -80000 * ml_year_to_m3_day]]
    wel_sp2 = [well_l_r_c + [1.]]
    wel_sp3 = [well_l_r_c + [well_pumping * -40000 * ml_year_to_m3_day]]
    return {'stress_period_data': {0: wel_sp1, 1: wel_sp2, 2: wel_sp3}}


def get_chd(nrow, chd_boundary_head):
    peak = chd_boundary_head * 15.
    return {
        'stress_period_data': {
            0: [[0, ir, 0, 10., 10.] for ir in range(nrow)],
            1: [[0, ir, 0, 10., peak] for ir in range(nrow)],
            2: [[0, ir, 0, peak, 10.] for ir in range(nrow)],
        }
    }


def get_dis():
    # Model domain and grid definition
    minX, maxX, minY, maxY = (247220.78719887036, 323214.2661659557, 5905095.566330256, 6040228.733922092)
    Lx = maxX - minX
    Ly = maxY - minY
    print('Lx, Ly', Lx, Ly)
    (nlay, nrow, ncol) = (2, 18, 10)
    delr = Lx / ncol
    delc = Ly / nrow
    return {
        'minX': minX, 'maxX': maxX, 'minY': minY, 'maxY': maxY,
        'Lx': Lx, 'Ly': Ly,
        'delr': delr, 'delc': delc,
        'ncol': ncol, 'nrow': nrow,
    }


def run_modflow(well_pumping, riv_level, chd_boundary_head, plot=False):

    # Model domain and grid definition
    minX, maxX, minY, maxY = (247220.78719887036, 323214.2661659557, 5905095.566330256, 6040228.733922092)
    Lx = maxX - minX
    Ly = maxY - minY
    print('Lx, Ly', Lx, Ly)
    (nlay, nrow, ncol) = (2, 18, 10)
    delr = Lx / ncol
    delc = Ly / nrow
    # delv = (ztop - zbot) / nlay
    # botm = [0 * np.ones((nrow, ncol)), -35. * np.ones((nrow, ncol)), -55. * np.ones((nrow, ncol))]
    botm = np.array([
        np.meshgrid(np.ones(ncol), np.linspace(-40, 0, nrow))[1],
        np.meshgrid(np.ones(ncol), np.linspace(-100, -20, nrow))[1],
        np.meshgrid(np.ones(ncol), np.linspace(-140, -23, nrow))[1],
    ])
    
    dis_args = get_dis()

    hk = interpolation(layer_avgs=[4., 1000.], plot=plot, **dis_args)  # hydraulic conductivity in meters/day
    vka = 1.
    sy = [0.1, 1e-2] # specific yield
    ss = 1.e-4
    laytyp = 1

    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    strt = 10. * np.ones((nlay, nrow, ncol), dtype=np.float32)

    # Time step parameters
    nper = 3
    perlen = [1, 100, 100]
    nstp = perlen
    steady = [True, False, False]

    # Flopy objects
    model_name = 'base'
    MODFLOW_executable = '/mfnwt'
    model_dir = os.path.join(os.path.dirname(__file__), 'too_big_to_git', model_name, multiprocessing.current_process().name)
    mf = flopy.modflow.Modflow(model_name, exe_name=MODFLOW_executable, model_ws=model_dir)
    dis = flopy.modflow.ModflowDis(mf, nlay, nrow, ncol, delr=delr, delc=delc,
                                   top=botm[0], botm=botm[1:],
                                   nper=nper, perlen=perlen, nstp=nstp,
                                   steady=steady)
    y_g, x_g, z_g = dis.get_node_coordinates()
    assert np.allclose(x_g + minX, np.linspace(minX + 0.5 * delr, maxX - 0.5 * delr, ncol))
    assert np.allclose(np.flip(y_g + minY, axis=0), np.linspace(minY + 0.5 * delc, maxY - 0.5 * delc, nrow))

    bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt)
    lpf = flopy.modflow.ModflowLpf(mf, hk=hk, vka=vka, sy=sy, ss=ss, laytyp=laytyp, ipakcb=53)
    pcg = flopy.modflow.ModflowPcg(mf)

    chd_data = get_chd(nrow, chd_boundary_head)
    # Specified Head boundaries (which are examples of  Dirichlet or first-type boundary conditions)
    chd = flopy.modflow.ModflowChd(mf, **chd_data)

    mm_year_to_m_day = 1. / 1000 / 365
    rch = flopy.modflow.ModflowRch(mf, rech={0: 45 * mm_year_to_m_day, 1: 70 * mm_year_to_m_day, 2: 30 * mm_year_to_m_day})

    riv_data = get_riv(hk, botm, nrow, ncol, delr, delc, riv_level)

    # Create the flopy ghb object
    # Head-Dependent Flux boundaries (which are examples of Robin or mixed boundary conditions)
    riv = flopy.modflow.ModflowRiv(mf, **riv_data)

    well_data = get_well(nrow, ncol, well_pumping)
    # Create the well package
    # Specified Flux boundaries (which are examples of  Neumann or second-type boundary conditions)
    wel = flopy.modflow.ModflowWel(mf, **well_data)

    # Output control
    stress_period_data = {(0, 0): ['save head',
                                   'save drawdown',
                                   'save budget',
                                   'print head',
                                   'print budget']}
    save_head_every = 1
    oc = flopy.modflow.ModflowOc(mf, stress_period_data=stress_period_data, compact=True)

    # Write the model input files
    mf.write_input()
    # Run the model
    success, mfoutput = mf.run_model(silent=False, pause=False)
    if not success:
        raise Exception('MODFLOW did not terminate normally.')

    # Create the headfile and budget file objects
    headobj = flopy.utils.binaryfile.HeadFile(os.path.join(model_dir, model_name+'.hds'))
    cbb = flopy.utils.binaryfile.CellBudgetFile(os.path.join(model_dir, model_name+'.cbc'))
    # flows = cbb.get_data(text='FLOW RIGHT FACE', kstpkper=(99,2))[0] # zero based (kstp, kper)

    if plot:
        plot_outputs(mf, headobj, cbb)
    return mf, headobj, cbb



def plot_outputs(mf, headobj, cbb):
    times = headobj.get_times()
    
    dis = mf.get_package('DIS')
    nrow = dis.nrow
    ncol = dis.ncol
    nlay = dis.nlay
    delr = dis.delr.array
    delc = dis.delc.array

    assert np.min(delr) == np.max(delr)
    delr = delr[0]
    assert np.min(delc) == np.max(delc)
    delc = delc[0]

    Lx = delr * ncol
    Ly = delc * nrow
    
    wel = mf.get_package('WEL')
    w_l, w_r, w_c, w_f = wel.stress_period_data[0][0]
    # TODO why 1.5
    wpt = ((w_c + 0.5)*delr, (w_r + 1.5)*delc)

    import matplotlib.pyplot as plt

    # Setup contour parameters
    levels = np.linspace(0, 10, 11)
    print('Levels: ', levels)
    # extent = (delr/2., Lx - delr/2., delc/2., Ly - delc/2.)
    # print('Extent: ', extent)

    # Make the plots
    print('times', times)
    mytimes = [1.0, 101.0, 201.0]
    for iplot, time in enumerate(mytimes):
        print('*****Processing time: ', time)
        head = headobj.get_data(totim=time)
        #Print statistics
        print('Head statistics')
        print('  min: ', head.min())
        print('  max: ', head.max())
        print('  std: ', head.std())

        # Extract flow right face and flow front face
        frf = cbb.get_data(text='FLOW RIGHT FACE', totim=time)[0]
        fff = cbb.get_data(text='FLOW FRONT FACE', totim=time)[0]

        #Create the plot
        #plt.subplot(1, len(mytimes), iplot + 1, aspect='equal')
        plt.subplot(1, 1, 1, aspect='equal')
        plt.title('stress period ' + str(iplot + 1))

        modelmap = flopy.plot.ModelMap(model=mf, layer=0)
        qm = modelmap.plot_ibound()
        lc = modelmap.plot_grid()
        qm = modelmap.plot_bc('RIV', alpha=0.5)
        cs = modelmap.contour_array(head, levels=levels)
        plt.clabel(cs, inline=1, fontsize=10, fmt='%1.1f', zorder=11)
        quiver = modelmap.plot_discharge(frf, fff, head=head)


        mfc = 'None'
        if (iplot+1) == len(mytimes):
            mfc='black'
        plt.plot(wpt[0], wpt[1], lw=0, marker='o', markersize=8, 
                 markeredgewidth=0.5,
                 markeredgecolor='black', markerfacecolor=mfc, zorder=9)
        plt.text(wpt[0]+25, wpt[1]-25, 'well', size=12, zorder=12)
        plt.savefig('/shared/plot/base/%d_tut.png' % iplot)

    # Plot the head versus time
    idx = (0, nrow/2 - 1, ncol/2 - 1)
    ts = headobj.get_ts(idx)
    plt.subplot(1, 1, 1)
    ttl = 'Head at cell ({0},{1},{2})'.format(idx[0] + 1, idx[1] + 1, idx[2] + 1)
    plt.title(ttl)
    plt.xlabel('time')
    plt.ylabel('head')
    plt.plot(ts[:, 0], ts[:, 1])
    plt.show()
    plt.savefig('/shared/plot/base/ts.png')


def snap_to_grid(x, y, minX, minY, delr, delc, **kwargs):
    return (0, int((y - minY) / delc), int((x - minX) / delr))


def final_head_at_locations(headobj):
    times = headobj.get_times()
    head = headobj.get_data(totim=times[-1])
    dis = get_dis()
    xs, ys, zs = locations(**dis)
    return [head[snap_to_grid(x, y, **dis)] for x, y in zip(xs, ys)]


if __name__ == '__main__':
    # extent()
    x = [1, 1, 1]
    mf, headobj, cbb = run_modflow(well_pumping=x[0], riv_level=x[1], chd_boundary_head=x[2], plot=True)

    print(final_head_at_locations(headobj))

