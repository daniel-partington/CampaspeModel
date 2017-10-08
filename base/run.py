"""
docker run -it \
     -v /Users/masher/Documents:/shared \
     -e PYTHONPATH=/shared/benjamin-holt:/shared/CampaspeModel \
     -e CACHE_DIR=/shared/CampaspeModel/too_big_to_git/cache\
     mjasher/heat \
    python /shared/CampaspeModel/base/run.py
"""

"""
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


def run_modflow(plot=False):

    # Model domain and grid definition
    Lx = 75993.
    Ly = 135133.
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
    zbot = botm[0]

    hk = [4., 1000.]  # hydraulic conductivity in meters/day
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
                                   top=zbot, botm=botm[1:],
                                   nper=nper, perlen=perlen, nstp=nstp,
                                   steady=steady)
    bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt)
    lpf = flopy.modflow.ModflowLpf(mf, hk=hk, vka=vka, sy=sy, ss=ss, laytyp=laytyp, ipakcb=53)
    pcg = flopy.modflow.ModflowPcg(mf)
    mm_year_to_m_day = 1. / 1000 / 365
    rch = flopy.modflow.ModflowRch(mf, rech={0: 45 * mm_year_to_m_day, 1: 70 * mm_year_to_m_day, 2: 30 * mm_year_to_m_day})

    # [lay, row, col, stage, cond, rbot]

    # Make list for stress period 1
    stage_murray = 3.
    stage_campaspe = 1.
    rbot = -5.
    cond_murray = hk[0] * (stage_murray - zbot) * delc
    cond_campaspe = hk[0] * (stage_campaspe - zbot) * delc
    bound_sp1 = []
    for ic in range(ncol):
        bound_sp1.append([0, 0, ic, stage_murray, cond_murray[0, ic], rbot])
    for ir in range(nrow):
        bound_sp1.append([0, ir, ncol/2, stage_campaspe, cond_campaspe[ir, ncol/2], rbot])
    
    print(botm[:, 0, ncol/2], 'to', botm[:, nrow - 1, ncol/2])

    print('Adding ', len(bound_sp1), 'RIVs for stress period 1.')

    # Make list for stress period 2
    stage_murray = 3.
    stage_campaspe = 0.
    cond_murray = hk[0] * (stage_murray - zbot) * delc
    cond_campaspe = hk[0] * (stage_campaspe - zbot) * delc
    bound_sp2 = []
    for ic in range(ncol):
        bound_sp2.append([0, 0, ic, stage_murray, cond_murray[0, ic], rbot])
    for ir in range(nrow):
        bound_sp1.append([0, ir, ncol/2, stage_campaspe, cond_campaspe[ir, ncol/2], rbot])
    print('Adding ', len(bound_sp2), 'RIVs for stress period 2.')


    # We do not need to add a dictionary entry for stress period 3.
    # Flopy will automatically take the list from stress period 2 and apply it
    # to the end of the simulation
    riv_data = {0: bound_sp1, 1: bound_sp2}


    # Create the flopy ghb object
    riv = flopy.modflow.ModflowRiv(mf, stress_period_data=riv_data)


    # Create the well package
    # Remember to use zero-based layer, row, column indices!
    ml_year_to_m3_day = 1. * 1000 / 365
    wel_sp1 = [[1, nrow/2 - 1, ncol/2 - 1, -8000 * ml_year_to_m3_day]]
    wel_sp2 = [[1, nrow/2 - 1, ncol/2 - 1, 1.]]
    wel_sp3 = [[1, nrow/2 - 1, ncol/2 - 1, -4000 * ml_year_to_m3_day]]
    well_data = {0: wel_sp1, 1: wel_sp2, 2: wel_sp3}
    wel = flopy.modflow.ModflowWel(mf, stress_period_data=well_data)

    # Output control
    stress_period_data = {(0, 0): ['save head',
                                   'save drawdown',
                                   'save budget',
                                   'print head',
                                   'print budget']}
    save_head_every = 1
    oc = flopy.modflow.ModflowOc(mf, stress_period_data=stress_period_data,
                                 compact=True)

    # Write the model input files
    mf.write_input()

    # Run the model
    success, mfoutput = mf.run_model(silent=False, pause=False)
    if not success:
        raise Exception('MODFLOW did not terminate normally.')

    # Create the headfile and budget file objects
    headobj = flopy.utils.binaryfile.HeadFile(os.path.join(model_dir, model_name+'.hds'))
    times = headobj.get_times()
    cbb = flopy.utils.binaryfile.CellBudgetFile(os.path.join(model_dir, model_name+'.cbc'))

    if plot:
        import matplotlib.pyplot as plt

        # Setup contour parameters
        levels = np.linspace(0, 10, 11)
        extent = (delr/2., Lx - delr/2., delc/2., Ly - delc/2.)
        print('Levels: ', levels)
        print('Extent: ', extent)

        # Well point
        wpt = ((float(ncol/2)-0.5)*delr, (float(nrow/2-1)+0.5)*delc)

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

    return headobj.get_alldata(), cbb.get_data(text='FLOW RIGHT FACE', kstpkper=(99,2))[0]

if __name__ == '__main__':
    heads, flows = run_modflow(plot=True)
