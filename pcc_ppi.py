"""

Das Program dient der Veranschaulichung der 5 Minutigen RX Radolan Daten!
Es werden durch den Zeitstempel Deutschlandweite Niederschlags und
Reflektivitaeten dargestellt!

"""



import numpy as np
import matplotlib.pyplot as plt
import wradlib
import h5py
import wradlib as wrl
import pandas as pd
import datetime as dt
import pcc as pcc

from pcc import boxpol_pos
from pcc import plot_radar

bonn_pos = boxpol_pos()
bx, by = bonn_pos['gkx_ppi'], bonn_pos['gky_ppi']
blat, blon = bonn_pos['lat_ppi'], bonn_pos['lon_ppi']


############################################### Zeitstempel nach YYYYMMDDhhmmss

from pcc import zeitschleife as zt

zeit = zt(2014,10,7,2,45,0,
          2014,10,7,2,50,0,
          steps=30)


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
            str(year)+'-'+str(m)+'-'+str(d)+'/ppi_1p5deg/'+str(year)+'-'+str(m)+'-'+
            str(d)+'--'+str(ht)+':'+str(mt)+':00,00.mvol')

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

    lon_ppi = attrs['VOL']['Longitude']
    lat_ppi = attrs['VOL']['Latitude']
    alt_ppi = attrs['VOL']['Height']


    #Z = wradlib.trafo.idecibel(ZH)
    #R = wradlib.zr.z2r(Z, a=200., b=1.6)

    #### PLOT
    fig = plt.figure(figsize=(16,12))

    fig.add_subplot(231, aspect='equal')
    ax1, pm1 = wradlib.vis.plot_ppi(zh,r,az)
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
    ax2, pm2 = wradlib.vis.plot_ppi(zv,r,az)
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
    ax3, pm3 = wradlib.vis.plot_ppi(zdr,r,az)
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
    ax5, pm5 = wradlib.vis.plot_ppi(rhohv,r,az)
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
    ax6, pm6 = wradlib.vis.plot_ppi(kdp,r,az)
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
    radolan_zeit = rwattrs['datetime'].strftime("%Y.%m.%d -- %H:%M:%S")
    radolan_zeit_sav = rwattrs['datetime'].strftime("%Y%m%d-%H%M%S")

    radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
    x = radolan_grid_xy[:,:,0]
    y = radolan_grid_xy[:,:,1]

    ZZ = rwdata
    Z = wradlib.trafo.idecibel(rwdata)

    # Marshall and Palmer 1948
    rwdata = wradlib.zr.z2r(Z, a=200., b=1.6)


    my_cmap = pcc.get_my_cmap()
    cmap2 = pcc.get_miub_cmap() #' bei Reflektivitat'
    ########################################################################## PLOT

    ff = 15
    fig = plt.figure(figsize=(14,10))

    ax1 = fig.add_subplot(121, aspect='equal')
    plt.pcolormesh(x, y, rwdata, cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
    #plt.scatter(x, y, rwdata, cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
    cb = plt.colorbar(shrink=0.5)
    cb.set_label("Rainrate (mm/h)",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    plot_borders(ax1)
    plt.title('RADOLAN Rainrate: \n'+ radolan_zeit + 'UTC',fontsize=ff)

    #plot_ocean(ax1)
    plt.xlabel("x [km] ",fontsize=ff)
    plt.ylabel("y [km]  ",fontsize=ff)
    #plt.xticks(fontsize=0)
    #plt.yticks(fontsize=0)
    plt.grid(color='r')
    plt.xlim(-420,390)
    plt.ylim(-4700, -3700)
    plt.tight_layout()

    plot_radar(blon, blat, ax1, reproject=True)

    ax2 = fig.add_subplot(122, aspect='equal')
    plt.pcolormesh(x, y, ZZ,vmin=-30,vmax=50, cmap=cmap2, zorder=2)
    #plt.scatter(x, y, rwdata, cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
    cb = plt.colorbar(shrink=0.5, extend='both')
    cb.set_label("Reflectivity (dBZ)",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    plot_borders(ax2)
    plt.title('RADOLAN Reflectivity: \n'+ radolan_zeit + 'UTC',fontsize=ff)

    #plot_ocean(ax1)
    plt.xlabel("x [km] ",fontsize=ff)
    plt.ylabel("y [km]  ",fontsize=ff)
    #plt.xticks(fontsize=0)
    #plt.yticks(fontsize=0)
    plt.grid(color='r')
    plt.xlim(-420,390)
    plt.ylim(-4700, -3700)
    plt.tight_layout()

    plot_radar(blon, blat, ax2, reproject=True)


    plt.savefig('/home/velibor/shkgpm/plot/radolan/rx_'+ radolan_zeit_sav+ '.png')
    plt.show()
    plt.close()
'''
