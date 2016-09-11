import os
import numpy as np

import matplotlib.pyplot as plt
import flopy
import flopy.utils.binaryfile as bf

def control():
    stress_period_data = {(0, 0): ['save head',
                               'save drawdown',
                               'save budget',
                               'print head',
                               'print budget']}

    return {
        "stress_period_data":stress_period_data,
        "compact": True
     }



def read_and_plot(mf, shape, dels, Ls):
    model_dir = mf.model_ws
    modelname = mf.name

    # Create the headfile and budget file objects
    headobj = bf.HeadFile(os.path.join(model_dir, modelname+'.hds'))
    # times = headobj.get_times()
    cbb = bf.CellBudgetFile(os.path.join(model_dir, modelname+'.cbc'))

    nlay, nrow, ncol = shape
    delr, delc = dels
    Lx, Ly = Ls

    # Setup contour parameters
    levels = np.linspace(0, 10, 11)
    extent = (delr/2., Lx - delr/2., delc/2., Ly - delc/2.)
    print('Levels: ', levels)
    print('Extent: ', extent)

    # Well point
    wpt = ((float(ncol/2)-0.5)*delr, (float(nrow/2-1)+0.5)*delc)
    wpt = (450., 550.)

    # Make the plots
    mytimes = [1.0, 101.0, 365.0]
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
        qm = modelmap.plot_bc('GHB', alpha=0.5)
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
        plt.show()

    plt.show()

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
