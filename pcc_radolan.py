"""

Das Program dient der Veranschaulichung der 5 Minutigen RX Radolan Daten!
Es werden durch den Zeitstempel Deutschlandweite Niederschlags und
Reflektivitaeten dargestellt!


"""



import numpy as np
import matplotlib.pyplot as plt
import wradlib  
import wradlib as wrl
import pandas as pd
import datetime as dt
import pcc as pcc

from pcc import boxpol_pos
from pcc import plot_radar

bonn_pos = boxpol_pos()
bx, by = bonn_pos['gkx_ppi'], bonn_pos['gky_ppi']
blat, blon = bonn_pos['lat_ppi'], bonn_pos['lon_ppi']
from pcc import plot_borders

############################################### Zeitstempel nach YYYYMMDDhhmmss

from pcc import zeitschleife as zt

zeit = zt(2017,7,1,00,00,0,
          2017,7,31,23,55,0,
          steps=5)
#zeit = zt(2014,10,07,02,35,0, 2014,10,07,02,40,0, steps=5)

for ij in range(len(zeit)):

    print 'Zeitstempel ',ij,' von  ', len(zeit)

    ZP = zeit[ij]
    print 'Zeitpunkt: ',ZP

    year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
    ye = ZP[2:4]



    ################################################### Read RADOLAN GK Koordinaten

    iii = 0
    r_pro = 'rx'
    pfad = ('/automount/radar/dwd/'+ r_pro +'/'+str(year)+'/'+str(year)+'-'+str(m)+'/'+
            str(year)+'-'+str(m)+'-'+str(d)+'/raa01-'+r_pro+'_10000-'+str(ye)+str(m)+
            str(d)+str(ht)+str(mt)+'-dwd---bin.gz')

    pfad_radolan = pfad[:-3]


    try:
        rw_filename = wradlib.util.get_wradlib_data_file(pfad)
    except EnvironmentError:
        rw_filename = wradlib.util.get_wradlib_data_file(pfad_radolan)


    rwdata, rwattrs = wradlib.io.read_RADOLAN_composite(rw_filename)

    rwdata = np.ma.masked_equal(rwdata, -9999) / 2 - 32.5

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
    fig = plt.figure(figsize=(16,10))

    ax1 = fig.add_subplot(121, aspect='equal')
    plt.pcolormesh(x, y, rwdata, cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
    #plt.scatter(x, y, rwdata, cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
    cb = plt.colorbar(shrink=0.5)
    cb.set_label("Rainrate (mm/h)",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    plot_borders(ax1)

    #plt.plot(gpm_x[:,0], gpm_y[:,0], color='black')
    #plt.plot(gpm_x[:,-1], gpm_y[:,-1], color='black')
    #plt.plot(gpm_x[:,23], gpm_y[:,23], color='black', ls='--')

    plt.title('RADOLAN RY Rainrate: \n'+ radolan_zeit + 'UTC',fontsize=20)
    plt.tick_params(
        axis='both',
        which='both',
        bottom='off',
        top='off',
        labelbottom='off',
        right='off',
        left='off',
        labelleft='off')
    #plot_ocean(ax1)
    #plt.xlabel("x [km] ",fontsize=ff)
    #plt.ylabel("y [km]  ",fontsize=ff)
    #plt.xticks(fontsize=0)
    #plt.yticks(fontsize=0)
    plt.grid(color='r')
    plt.xlim(-420,390)
    plt.ylim(-4700, -3700)

    plt.tight_layout()


    from pcc import get_radar_locations
    radar = get_radar_locations()

    #for i in range(len(radar.keys())):
    #    plot_radar(radar[radar.keys()[i]]['lon'],
    #               radar[radar.keys()[i]]['lat'],
    #               ax1, reproject=True, cband=True, col='red')
    plot_radar(blon, blat, ax1, reproject=True, cband=False,col='black')

    ax2 = fig.add_subplot(122, aspect='equal')
    plt.pcolormesh(x, y, ZZ,vmin=0,vmax=50, cmap=cmap2, zorder=2)
    #plt.scatter(x, y, rwdata, cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
    cb = plt.colorbar(shrink=0.5, extend='both')
    cb.set_label("Reflectivity (dBZ)",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)

    #plt.plot(gpm_x[:,0], gpm_y[:,0], color='black')
    #plt.plot(gpm_x[:,-1], gpm_y[:,-1], color='black')
    #plt.plot(gpm_x[:,23], gpm_y[:,23], color='black', ls='--')

    plot_borders(ax2)
    plt.title('RADOLAN RX Reflectivity: \n'+ radolan_zeit + 'UTC',fontsize=20)

    plt.tick_params(
        axis='both',
        which='both',
        bottom='off',
        top='off',
        labelbottom='off',
        right='off',
        left='off',
        labelleft='off')

    #plot_ocean(ax1)
    #plt.xlabel("x [km] ",fontsize=ff)
    #plt.ylabel("y [km]  ",fontsize=ff)
    #plt.xticks(fontsize=0)
    #plt.yticks(fontsize=0)
    plt.grid(color='r')
    plt.xlim(-420,390)
    plt.ylim(-4700, -3700)
    plt.tight_layout()

    plot_radar(blon, blat, ax2, reproject=True, cband=False,col='black')
    #for i in range(len(radar.keys())):
    #    plot_radar(radar[radar.keys()[i]]['lon'],
    #               radar[radar.keys()[i]]['lat'],
    #               ax2, reproject=True, cband=True, col='red')


    #plt.savefig('/automount/ags/velibor/plot/radolan/RealPEP/lauf2/r_'+ radolan_zeit_sav+ '.png')
    plt.savefig('/automount/ags/velibor/plot/radolan/1-31_07_2017/r_'+ radolan_zeit_sav+ '.png')

    #plt.show()
    plt.close()

