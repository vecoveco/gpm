"""

Einlesen und darstellen von GPM DPR und Radolan Dateien

Radolanpfad:

"""


import h5py
import numpy as np
import matplotlib.pyplot as plt
import wradlib
import glob
from scipy import stats, linspace
import wradlib as wrl
from osgeo import osr
from pcc import get_time_of_gpm
from pcc import cut_the_swath

## Landgrenzenfunktion
## -------------------
from pcc import boxpol_pos
bonn_pos = boxpol_pos()
bx, by = bonn_pos['gkx_ppi'], bonn_pos['gky_ppi']
bonnlat, bonnlon = bonn_pos['lat_ppi'], bonn_pos['lon_ppi']
from pcc import plot_borders
from pcc import plot_radar
from pcc import get_miub_cmap
my_cmap = get_miub_cmap()

from pcc import get_my_cmap
my_cmap2 = get_my_cmap()
from pcc import get_4_cmap


#zz = np.array([20140609, 20140610, 20140629, 20140826, 20140921, 20141007,
#               20141016, 20150128, 20150227, 20150402, 20150427, 20160405,
#               20160607, 20160805, 20160904, 20160917, 20161001, 20161024,
#               20170113, 20170203,])

zz = np.array(['20141007'])
for i in range(len(zz)):
    ZP = str(zz[i])
    #year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
    year, m, d = ZP[0:4], ZP[4:6], ZP[6:8]

    ye = ZP[2:4]

    ## Read GPM Data
    ## -------------

    pfad2 = ('/home/velibor/shkgpm/data/'+str(year)+str(m)+str(d)+'/dpr/*.HDF5')
    pfad_gpm = glob.glob(pfad2)
    pfad_gpm_g = pfad_gpm[0]

    gpmdpr = h5py.File(pfad_gpm_g, 'r')
    lat=np.array(gpmdpr['NS']['Latitude'])
    lon=np.array(gpmdpr['NS']['Longitude'])

    gprof_pp=np.array(gpmdpr['NS']['SLV']['zFactorCorrectedNearSurface'])

    gprof_pp[gprof_pp<=0]= np.nan


    gpm_time = gpmdpr['NS']['ScanTime']
    gpm_zeit = get_time_of_gpm(lon, lat, gpm_time)

    gprof_phase=np.array(gpmdpr['NS']['SLV']['phaseNearSurface'],dtype=int)
    #Phase 0-solid, 1-mixed, 2-liquid, 255-missing
    #gprof_phase[gprof_phase==255.]= np.nan

    #Precipitation water vertically integrate
    gprof_wi=np.array(gpmdpr['NS']['SLV']['precipRateNearSurface'])
    gprof_wi[gprof_wi==-9999]=np.nan


    # Precip Typ
    gprof_typ=np.array(gpmdpr['NS']['CSF']['typePrecip'])
    gprof_typ[gprof_typ==-9999] = 0
    gprof_typ[gprof_typ==-1111] = 0

    from pcc import get_3_cmap
    from pcc import plot_world
    x1,x2,y3,y4 = 5, 15, 45, 55
    ff = 15
    asp = 'equal'
    fig = plt.figure(figsize=(15,15))
    ax1 = fig.add_subplot(2,2,1, aspect=asp)
    plt.pcolormesh(lon, lat, np.ma.masked_invalid(gprof_pp), cmap=my_cmap)
    cb = plt.colorbar(shrink=0.8)
    cb.set_label("Ref [dbz]",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    plot_world(ax1,x1,x2,y3,y4)

    plt.grid()
    plt.xlim(x1,x2)
    plt.ylim(y3,y4)

    ax2 = fig.add_subplot(2,2,2, aspect=asp)
    plt.pcolormesh(lon, lat, np.ma.masked_invalid(gprof_phase/100),cmap=get_3_cmap())
    print np.unique(gprof_phase)
    cb = plt.colorbar(shrink=0.8)
    cb.set_label("Phase 0-solid, 1-mixed, 2-liquid",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    plot_world(ax2,x1,x2,y3,y4)

    plt.grid()
    plt.xlim(x1,x2)
    plt.ylim(y3,y4)

    ax3 = fig.add_subplot(2,2,3, aspect=asp)
    plt.pcolormesh(lon, lat, np.ma.masked_invalid(gprof_wi), vmin=0.1, vmax=10, cmap=get_my_cmap())
    cb = plt.colorbar(shrink=0.8)
    cb.set_label("RR (mm/h)",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    plot_world(ax3,x1,x2,y3,y4)

    plt.grid()
    plt.xlim(x1,x2)
    plt.ylim(y3,y4)

    ax4 = fig.add_subplot(2,2,4, aspect=asp)
    plt.pcolormesh(lon, lat, np.ma.masked_invalid(gprof_typ/10000000),cmap=get_4_cmap())
    cb = plt.colorbar(shrink=0.8)
    cb.set_label("Typ, (1=strat, 2=conv, 3=other)",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    plot_world(ax4,x1,x2,y3,y4)
    print np.unique(gprof_typ/10000000)
    plt.grid()
    plt.xlim(x1,x2)
    plt.ylim(y3,y4)

    plt.show()
    #plt.savefig('/home/velibor/shkgpm/plot/pro_dpr_'+ZP + '.png' )
    #fig.clf() # CLEAR FIGURE
    #plt.close()
