"""

Das Program dient der Veranschaulichung der 5 Minutigen RX Radolan Daten!
Es werden durch den Zeitstempel Deutschlandweite Niederschlags und
Reflektivitaeten dargestellt!

"""



import numpy as np
import matplotlib.pyplot as plt
import wradlib
import h5py



############################################### Zeitstempel nach YYYYMMDDhhmmss

from pcc import zeitschleife as zt

zeit = zt(2014,10,07,02,30,00,
          2014,10,07,02,50,00,
          steps=5)


for ij in range(len(zeit)):

    print 'Zeitstempel ',ij,' von  ', len(zeit)

    ZP = zeit[ij]
    print 'Zeitpunkt: ',ZP

    year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
    ye = ZP[2:4]



    if year < '2015':
        print 'archive'
        sc = 'radar-archiv'
    else:
        sc = 'radar'
    ################################################### Read RADOLAN GK Koordinaten

    iii = 0
    pfad = ('/automount/'+sc+'/scans/'+str(year)+'/'+str(year)+'-'+str(m)+'/'+
            str(year)+'-'+str(m)+'-'+str(d)+'/n_rhi_lacros/'+str(year)+'-'+str(m)+'-'+
            str(d)+'--'+str(ht)+':'+str(mt)+':01,00.mvol')

    #Todo: die an der Uhrzeit am nachliegensten Uhrzeit nehmen

    boxpol_filename = wradlib.util.get_wradlib_data_file(pfad)



    ppi=h5py.File(boxpol_filename,'r')
    data, attrs = wradlib.io.read_GAMIC_hdf5(boxpol_filename)

    print data[u'SCAN0'].keys()
    print attrs['VOL'].keys()
    print attrs['SCAN0'].keys()


    zh = data['SCAN0'][u'ZH']['data']
    phidp = data['SCAN0'][u'PHIDP']['data']
    rhohv = data['SCAN0'][u'RHOHV']['data']
    zv = data['SCAN0'][u'ZV']['data']
    zdr = data['SCAN0'][u'ZDR']['data']
    kdp = data['SCAN0'][u'KDP']['data']

    r = attrs['SCAN0']['r']
    az = attrs['SCAN0']['az']
    th = attrs['SCAN0']['el']

    lon_ppi = attrs['VOL']['Longitude']
    lat_ppi = attrs['VOL']['Latitude']
    alt_ppi = attrs['VOL']['Height']

    from mpl_toolkits.axisartist.grid_finder import FixedLocator, DictFormatter
    # reading in GAMIC hdf5 file

    mask_ind = np.where(data <= np.nanmin(zh))
    zh[mask_ind] = np.nan
    ma = np.ma.array(zh, mask=np.isnan(zh))


    fft = 15
    import pcc as pcc
    cmap2 = pcc.get_miub_cmap()

    sb = [231,232,233,234,235,236]
    radpar = [zh, zv, zdr,phidp, rhohv, kdp]
    vmini = [0, 0, -1, -100, .7,-0.75]
    vmaxi = [50, 50, 3, 50, .99, 2]
    name = ['Zh (dBZ)','Zv (dBZ)','Zdr (db)','PhiDP','RhoHV (/)', r"kdp ( $^\circ$ $km^{-1}$)"]
    #cc = [cmap2, cmap2, cmap2, 'jet', 'jet', 'jet']
    fig = plt.figure(figsize=(16,12))

    for iii in range(len(radpar)):

        cgax, caax, paax, pm, xxx, yyy = pcc.pcc_plot_cg_rhi(radpar[iii], r=r,
                                                             th=th,
                                                             rf=1e3,
                                                             subplot=sb[iii],
                                                             autoext=True,
                                                             vmin=vmini[iii],
                                                             vmax= vmaxi[iii],
                                                             cmap=cmap2)
        cgax.set_ylim(0,7)
        cbar = plt.gcf().colorbar(pm, shrink=0.8, extend='both')#, pad=0.05)
        plt.gca().invert_xaxis()
        cbar.set_label(name[iii],fontsize=fft)
        caax.set_xlabel('x (km)',fontsize=fft)
        caax.set_ylabel('z (km)',fontsize=fft)

        cbar.ax.tick_params(labelsize=fft)

        caax.tick_params(labelsize=fft)

        plt.title(' BoXPol RHI: '+ ZP, fontsize=fft)


        gh = cgax.get_grid_helper()
        locs = [0.]
        gh.grid_finder.grid_locator1 = FixedLocator(locs)
        gh.grid_finder.tick_formatter1 = DictFormatter(dict([(i, r"${0:.0f}^\circ$".format(i)) for i in locs]))
        plt.yticks(fontsize=fft)
        plt.xticks(fontsize=fft)



    plt.tight_layout()
    plt.show()


'''
    #### PLOT
    fig = plt.figure(figsize=(16,12))
    fig.suptitle('Date: '+ ZP)
    fig.add_subplot(231, aspect='equal')
    ax1, pm1 = wradlib.vis.plot_ppi(zh,r,az, vmin=0, vmax=50)
    cbar = plt.colorbar(pm1, shrink=0.75)
    cbar.set_label("zh (dBz)")
    plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
    plt.xlim((-101000,101000))
    plt.ylim((-101000,101000))

    fig.add_subplot(232, aspect='equal')
    ax2, pm2 = wradlib.vis.plot_ppi(zv,r,az, vmin=0, vmax=50)
    cbar = plt.colorbar(pm2, shrink=0.75)
    cbar.set_label("zv (dBz)")
    plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
    plt.xlim((-101000,101000))
    plt.ylim((-101000,101000))

    fig.add_subplot(233, aspect='equal')
    ax3, pm3 = wradlib.vis.plot_ppi(zdr,r,az, vmin=-1, vmax=6)
    cbar = plt.colorbar(pm3, shrink=0.75)
    cbar.set_label("zdr (dB)")
    plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
    plt.xlim((-101000,101000))
    plt.ylim((-101000,101000))

    fig.add_subplot(234, aspect='equal')
    ax4, pm4 = wradlib.vis.plot_ppi(phidp,r,az)
    cbar = plt.colorbar(pm4, shrink=0.75)
    cbar.set_label("phidp")
    plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
    plt.xlim((-101000,101000))
    plt.ylim((-101000,101000))

    fig.add_subplot(235, aspect='equal')
    ax5, pm5 = wradlib.vis.plot_ppi(rhohv,r,az, vmin=.7, vmax=.99)
    cbar = plt.colorbar(pm5, shrink=0.75)
    cbar.set_label("rhohv ()")
    plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
    plt.xlim((-101000,101000))
    plt.ylim((-101000,101000))


    fig.add_subplot(236, aspect='equal')
    ax6, pm6 = wradlib.vis.plot_ppi(kdp,r,az, vmin=-0.75, vmax=2)
    cbar = plt.colorbar(pm6, shrink=0.75)
    cbar.set_label(r"kdp $^\circ$ $km^{-1}$")
    plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
    plt.xlim((-101000,101000))
    plt.ylim((-101000,101000))

    plt.tight_layout()
    plt.show()

'''