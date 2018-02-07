"""

Das Program dient der Veranschaulichung der 5 Minutigen RX Radolan Daten!
Es werden durch den Zeitstempel Deutschlandweite Niederschlags und
Reflektivitaeten dargestellt!

Made by Velibor Pejcic


"""



import numpy as np
import h5py
import matplotlib.pyplot as plt
import wradlib
import pandas as pd
import pcc as pcc
from pcc import boxpol_pos
from pcc import plot_radar
from pcc import plot_borders
from time import *
my_cmap = pcc.get_my_cmap()
cmap2 = pcc.get_miub_cmap()
from pcc import get_radar_locations
radar = get_radar_locations()
from pcc import zeitschleife as zt
bonn_pos = boxpol_pos()
bx, by = bonn_pos['gkx_ppi'], bonn_pos['gky_ppi']
blat, blon = bonn_pos['lat_ppi'], bonn_pos['lon_ppi']
import os

t1 = clock()

zeit = zt(2016,01,01,0,00,0,
          2016,12,31,23,55,0,
          steps=5)

RADOLAN = []

att = ['RX','RY','zgrids_me', 'zgrids_st', 'zgrids_co',
                'rgrids_>0', 'rgrids_0-1',
                'rgrids_1-5', 'rgrids_5-10', 'rgrids_>10', 'rsum']

df0 = pd.DataFrame(columns=att)
csv_name = '/automount/ags/velibor/plot/radolan/test/rxy_'+ zeit[0]+'-'+zeit[-1]+ '.csv'
df0.to_csv(csv_name)


for ij in range(len(zeit)):

    print 'Zeitstempel ',ij,' von  ', len(zeit)

    ZP = zeit[ij]
    print 'Zeitpunkt: ',ZP

    # Zeitstempel definieren
    year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
    ye = ZP[2:4]

    ################################################### Read RADOLAN GK Koordinaten

    iii = 0
    #Try Read RX Data
    r_pro = 'rx'

    try:
        pfad = ('/automount/radar/dwd/'+ r_pro +'/'+str(year)+'/'+str(year)+'-'+str(m)+'/'+
                str(year)+'-'+str(m)+'-'+str(d)+'/raa01-'+r_pro+'_10000-'+str(ye)+str(m)+
                str(d)+str(ht)+str(mt)+'-dwd---bin.gz')
        pfad_radolan = pfad[:-3]



        try:
            rw_filename = wradlib.util.get_wradlib_data_file(pfad)
        except EnvironmentError:
            rw_filename = wradlib.util.get_wradlib_data_file(pfad_radolan)

        if os.path.getsize(rw_filename) != 0:

            rwdata, rwattrs = wradlib.io.read_RADOLAN_composite(rw_filename)

            rwdata = np.ma.masked_equal(rwdata, -9999) / 2 - 32.5
            radolan_zeit = rwattrs['datetime'].strftime("%Y.%m.%d -- %H:%M:%S")
            radolan_zeit_sav = rwattrs['datetime'].strftime("%Y%m%d-%H%M%S")

            radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
            x = radolan_grid_xy[:,:,0]
            y = radolan_grid_xy[:,:,1]
            ZZ = rwdata
            RX=1
            #print ('RX gefunden')
        else:
            RX=0
            print ('RX nicht gefunden')

    except:
        RX=0
        print ('RX nicht gefunden')

    #Try read RY Data
    r_pro = 'ry'
    try:

        pfad = ('/automount/radar/dwd/'+ r_pro +'/'+str(year)+'/'+str(year)+'-'+str(m)+'/'+
                str(year)+'-'+str(m)+'-'+str(d)+'/raa01-'+r_pro+'_10000-'+str(ye)+str(m)+
                str(d)+str(ht)+str(mt)+'-dwd---bin.gz')

        pfad_radolan = pfad[:-3]

        try:
            ry_filename = wradlib.util.get_wradlib_data_file(pfad)
        except EnvironmentError:
            ry_filename = wradlib.util.get_wradlib_data_file(pfad_radolan)

        if os.path.getsize(ry_filename) != 0:
            rydata, ryattrs = wradlib.io.read_RADOLAN_composite(ry_filename)

            rydata = np.ma.masked_equal(rydata, -9999)

            radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
            x = radolan_grid_xy[:,:,0]
            y = radolan_grid_xy[:,:,1]
            RR = rydata
            RY=1
            #print ('RY gefunden!')
        else:
            RY=0
            print ('RY nicht gefunden!')

    except:
        RY=0
        print ('RY nicht gefunden!')


    if ((RX==1) & (RY==1)):
        ploty = 1

    elif ((RX==1) & (RY==0)):
        # Wenn kein RY vorhanden MP
        # Marshall and Palmer 1948
        RR = wradlib.zr.z2r(wradlib.trafo.idecibel(rwdata), a=200., b=1.6)
        ploty = 1

    elif ((RX==0) & (RY==1)):
        # Wenn kein RY vorhanden MP
        # Marshall and Palmer 1948
        ZZ = wradlib.zr.r2z(wradlib.trafo.idecibel(rydata), a=200., b=1.6)
        ploty = 1

    else:# (RX==0 & RY==0):
        # Wenn RX und RY nicht vorhanden sind
        zgrids_me, zgrids_st, zgrids_co= -9999., -9999., -9999.
        rgrids_0,rgrids_01,rgrids_15 = -9999., -9999., -9999.
        rsum, rgrids_510, rgrids_10 = -9999.,-9999., -9999.
        ploty = 0
        #alles daten auf nan

    # PLOT an/aus
    '''
    ##################################################################### PLOT
    if ploty != 0:

        ff = 15
        fig = plt.figure(figsize=(14,8))

        ax1 = fig.add_subplot(121, aspect='equal')
        plt.pcolormesh(x, y, RR, cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
        cb = plt.colorbar(shrink=0.5, extend='both')
        cb.set_label("Rainrate (mm/h)",fontsize=ff)
        cb.ax.tick_params(labelsize=ff)
        plot_borders(ax1)

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

        plt.grid(color='r')
        plt.xlim(-420,390)
        plt.ylim(-4700, -3700)

        plt.tight_layout()




        plot_radar(blon, blat, ax1, reproject=True, cband=False,col='black')

        ax2 = fig.add_subplot(122, aspect='equal')
        plt.pcolormesh(x, y, ZZ,vmin=0,vmax=50, cmap=cmap2, zorder=2)
        cb = plt.colorbar(shrink=0.5, extend='both')
        cb.set_label("Reflectivity (dBZ)",fontsize=ff)
        cb.ax.tick_params(labelsize=ff)
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

        plt.grid(color='r')
        plt.xlim(-420,390)
        plt.ylim(-4700, -3700)
        plt.tight_layout()

        plot_radar(blon, blat, ax2, reproject=True, cband=False,col='black')

        plt.savefig('/automount/ags/velibor/plot/radolan/test/r_'+ radolan_zeit_sav+ '.png')
        plt.close()
        ##########################################################################





    else:
        print ('kein Plot!')
        pass
    '''

    # Auswahl bestimmter Grid Parameter
    zgrids_me = ZZ[np.where(ZZ>=-32.5)].shape[0]
    zgrids_st = ZZ[np.where((ZZ>=5)&(ZZ<=30))].shape[0]
    zgrids_co = ZZ[np.where(ZZ>30)].shape[0]

    rgrids_0 = RR[np.where(RR>0)].shape[0]
    rgrids_01 = RR[np.where((RR>0)&(RR<=1))].shape[0]
    rgrids_15 = RR[np.where((RR>=1)&(RR<=5))].shape[0]
    rgrids_510 = RR[np.where((RR>=5)&(RR<=10))].shape[0]
    rgrids_10 = RR[np.where(RR>10)].shape[0]

    # Summe mm/h ueber Regenflaeche
    rsum = np.sum(RR[np.where(RR>0)])

    #Time Index Setzen
    time_index = pd.to_datetime(([ZP]), format='%Y%m%d%H%M%S')

    # Benutzte Grid Parameter als Array speichern
    radolan_array = np.array([[RX,RY,zgrids_me, zgrids_st, zgrids_co,rgrids_0, rgrids_01,
                    rgrids_15, rgrids_510, rgrids_10,rsum]])

    # Pandas Dataframe erstellen
    df_new = pd.DataFrame(radolan_array, index=time_index, columns=att)

    # In die erstellte CSV schreiben
    df_new.to_csv(csv_name, mode='a', header=False)
'''
    del(ZZ,RR,RX,RY,zgrids_me, zgrids_st, zgrids_co,rgrids_0,
        rgrids_01,rgrids_15, rgrids_510, rgrids_10,rsum, radolan_array,
        time_index,ploty)

    try:
        del(rwdata,rydata)
    except NameError:
        pass

'''


# Berechnung der Laufzeit
t2 = clock()
print ('___________________________________')
print ('')
print ('Laufzeit: ', (t2 -t1), 'Sekunden')
print ('')
print ('___________________________________')



