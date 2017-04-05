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

#zz = np.array(['20141007'])
pfad2 = ('/automount/ags/velibor/gpmdata/dpr/*.HDF5')
pfad_gpm = glob.glob(pfad2)
for i in range(len(pfad_gpm)):

    pfad_gpm_g = pfad_gpm[i]

    try:

        gpmdpr = h5py.File(pfad_gpm_g, 'r')
        lat=np.array(gpmdpr['NS']['Latitude'])
        lon=np.array(gpmdpr['NS']['Longitude'])

        gprof_pp=np.array(gpmdpr['NS']['SLV']['zFactorCorrectedNearSurface'])

        gprof_pp[gprof_pp<=0]= np.nan


        gpm_time = gpmdpr['NS']['ScanTime']

        try:
            gpm_zeit = get_time_of_gpm(lon, lat, gpm_time)

        except ValueError:

            gpm_zeit = ' '

        gprof_phase=np.array(gpmdpr['NS']['SLV']['phaseNearSurface'],dtype=int)
        gprof_wi=np.array(gpmdpr['NS']['CSF']['flagShallowRain'],dtype=float)

        gprof_bb=np.array(gpmdpr['NS']['CSF']['heightBB'])
        gprof_bb[gprof_bb==-1111.1] = np.nan
        gprof_wi[gprof_wi==-1111.] = np.nan
        #Phase 0-solid, 1-mixed, 2-liquid, 255-missing
        #gprof_phase[gprof_phase==255.]= np.nan

        #Precipitation water vertically integrate
        gprof_rr=np.array(gpmdpr['NS']['SLV']['precipRateNearSurface'])
        gprof_rr[gprof_rr==-9999]=np.nan


        # Precip Typ
        gprof_typ=np.array(gpmdpr['NS']['CSF']['typePrecip'])
        gprof_typ[gprof_typ==-9999] = 0
        gprof_typ[gprof_typ==-1111] = 0

        from pcc import get_3_cmap
        from pcc import plot_world
        x1,x2,y3,y4 = 5, 15, 45, 55
        ff = 15
        sh = 0.5
        asp = 'equal'#'auto'
        ori = 'vertical'

        ####################### PLOT

        fig = plt.figure(figsize=(20,15))

        ax1 = fig.add_subplot(2,3,1, aspect=asp)
        plt.pcolormesh(lon, lat, np.ma.masked_invalid(gprof_pp), cmap=my_cmap)
        plt.plot(lon[:,0],lat[:,0], color='black',lw=1)
        plt.plot(lon[:,-1],lat[:,-1], color='black',lw=1)
        cb = plt.colorbar(shrink=sh, orientation=ori)
        cb.set_label("Reflectivity (dbz)",fontsize=ff)
        plt.title('GPM DPR Reflectivity: ',fontsize=ff)
        cb.ax.tick_params(labelsize=ff)
        plot_world(ax1,x1,x2,y3,y4)

        plt.grid()
        plt.xlim(x1,x2)
        plt.ylim(y3,y4)

        ax2 = fig.add_subplot(2,3,3, aspect=asp)
        plt.pcolormesh(lon, lat, np.ma.masked_invalid(gprof_phase/100),cmap=get_3_cmap())
        plt.plot(lon[:,0],lat[:,0], color='black',lw=1)
        plt.plot(lon[:,-1],lat[:,-1], color='black',lw=1)
        cb = plt.colorbar(shrink=sh, orientation=ori)
        cb.set_label("Phase",fontsize=ff)
        cb.ax.tick_params(labelsize=ff)
        plt.title('GPM DPR Phase: \n Phase: 0=solid, 1=mixed, 2=liquid',fontsize=ff)
        plot_world(ax2,x1,x2,y3,y4)

        plt.grid()
        plt.xlim(x1,x2)
        plt.ylim(y3,y4)

        ax3 = fig.add_subplot(2,3,2, aspect=asp)
        plt.pcolormesh(lon, lat, np.ma.masked_invalid(gprof_wi), cmap='copper')
        plt.plot(lon[:,0],lat[:,0], color='black',lw=1)
        plt.plot(lon[:,-1],lat[:,-1], color='black',lw=1)
        cb = plt.colorbar(shrink=sh, orientation=ori)
        cb.set_label("Flag",fontsize=ff)
        cb.ax.tick_params(labelsize=ff)
        plt.title('GPM DPR flagShallowRain: ',fontsize=ff)
        plot_world(ax3,x1,x2,y3,y4)

        plt.grid()
        plt.xlim(x1,x2)
        plt.ylim(y3,y4)

        ax4 = fig.add_subplot(2,3,4, aspect=asp)
        plt.pcolormesh(lon, lat, np.ma.masked_invalid(gprof_typ/10000000),cmap=get_4_cmap())
        plt.plot(lon[:,0],lat[:,0], color='black',lw=1)
        plt.plot(lon[:,-1],lat[:,-1], color='black',lw=1)
        cb = plt.colorbar(shrink=sh, orientation=ori)
        cb.set_label("Type",fontsize=ff)
        cb.ax.tick_params(labelsize=ff)
        plot_world(ax4,x1,x2,y3,y4)
        plt.title('GPM DPR Raintype:  \n Type: 1=strat, 2=conv, 3=other',fontsize=ff)
        plt.grid()
        plt.xlim(x1,x2)
        plt.ylim(y3,y4)

        ax33 = fig.add_subplot(2,3,5, aspect=asp)
        plt.pcolormesh(lon, lat, np.ma.masked_invalid(gprof_bb), cmap='terrain')
        plt.plot(lon[:,0],lat[:,0], color='black',lw=1)
        plt.plot(lon[:,-1],lat[:,-1], color='black',lw=1)
        cb = plt.colorbar(shrink=sh, orientation=ori)
        cb.set_label("High in m",fontsize=ff)
        cb.ax.tick_params(labelsize=ff)
        plt.title('GPM DPR BrightBandHigh: ',fontsize=ff)
        plot_world(ax33,x1,x2,y3,y4)
        plt.grid()
        plt.xlim(x1,x2)
        plt.ylim(y3,y4)

        ax333 = fig.add_subplot(2,3,6, aspect=asp)
        plt.pcolormesh(lon, lat, np.ma.masked_invalid(gprof_rr), cmap=my_cmap2, vmin=0.1, vmax=10)
        plt.plot(lon[:,0],lat[:,0], color='black',lw=1)
        plt.plot(lon[:,-1],lat[:,-1], color='black',lw=1)
        cb = plt.colorbar(shrink=sh, orientation=ori)
        cb.set_label("Rain Rate in mm/h",fontsize=ff)
        cb.ax.tick_params(labelsize=ff)
        plt.title('GPM DPR Rain Rate : ',fontsize=ff)
        plot_world(ax333,x1,x2,y3,y4)
        plt.grid()
        plt.xlim(x1,x2)
        plt.ylim(y3,y4)

        plt.suptitle('Date: ' + str(gpm_zeit) + ' UTC', fontsize=ff)
        #plt.tight_layout(pad=0.9, w_pad=0.9, h_pad=0.9)
        #plt.tight_layout()
        #plt.show()
        plt.savefig('/automount/ags/velibor/plot/alledprskill/pro_dpr_'+str(gpm_zeit) + '.png')#, transparent=True)
        #fig.clf() # CLEAR FIGURE
        plt.close()
        #plt.show()

        del(pfad_gpm_g, gpmdpr, lat, lon, gprof_pp, gpm_time,
        gpm_zeit, gprof_phase, gprof_wi, gprof_typ, fig, ax1, ax2, ax3, ax4, cb, ax33, ax333)

    except:
        print '_________Fehler__________'

