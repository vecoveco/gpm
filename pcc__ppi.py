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

cmap2 = pcc.get_miub_cmap()

from pcc import boxpol_pos
from pcc import plot_radar

bonn_pos = boxpol_pos()
bx, by = bonn_pos['gkx_ppi'], bonn_pos['gky_ppi']
blat, blon = bonn_pos['lat_ppi'], bonn_pos['lon_ppi']


############################################### Zeitstempel nach YYYYMMDDhhmmss

from pcc import zeitschleife as zt

zeit = zt(2014,10,07,02,30,0,
          2014,10,07,02,35,0,
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
            str(year)+'-'+str(m)+'-'+str(d)+'/ppi_1p5deg/'+str(year)+'-'+str(m)+'-'+
            str(d)+'--'+str(ht)+':'+str(mt)+':00,00.mvol')
    print sc
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
    fft = 15
    fig.suptitle('Date: '+ ZP)
    fig.add_subplot(231, aspect='equal')
    ax1, pm1 = wradlib.vis.plot_ppi(zh,r,az, vmin=0, vmax=50, cmap=cmap2)
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
    plt.title(' BoXPol PPI: '+ ZP ,fontsize=fft)
    plt.xlim((-101000,101000))
    plt.ylim((-101000,101000))

    fig.add_subplot(232, aspect='equal')
    ax2, pm2 = wradlib.vis.plot_ppi(zv,r,az, vmin=0, vmax=50, cmap=cmap2)
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
    plt.title(' BoXPol PPI: '+ ZP ,fontsize=fft)
    plt.xlim((-101000,101000))
    plt.ylim((-101000,101000))

    fig.add_subplot(233, aspect='equal')
    ax3, pm3 = wradlib.vis.plot_ppi(zdr,r,az, vmin=-1, vmax=3, cmap=cmap2)
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
    plt.title(' BoXPol PPI: '+ ZP ,fontsize=fft)
    plt.xlim((-101000,101000))
    plt.ylim((-101000,101000))

    fig.add_subplot(234, aspect='equal')
    ax4, pm4 = wradlib.vis.plot_ppi(phidp,r,az, cmap=cmap2, vmin=-100, vmax=0)
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
    plt.title(' BoXPol PPI: '+ ZP ,fontsize=fft)
    plt.xlim((-101000,101000))
    plt.ylim((-101000,101000))

    fig.add_subplot(235, aspect='equal')
    ax5, pm5 = wradlib.vis.plot_ppi(rhohv,r,az, vmin=.7, vmax=.99, cmap=cmap2)
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
    plt.title(' BoXPol PPI: '+ ZP ,fontsize=fft)
    plt.xlim((-101000,101000))
    plt.ylim((-101000,101000))


    fig.add_subplot(236, aspect='equal')
    ax6, pm6 = wradlib.vis.plot_ppi(kdp,r,az, vmin=-0.75, vmax=2, cmap=cmap2)
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
    plt.title(' BoXPol PPI: '+ ZP ,fontsize=fft)
    plt.xlim((-101000,101000))
    plt.ylim((-101000,101000))

    plt.tight_layout()
    plt.show()

