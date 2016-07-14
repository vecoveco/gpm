
'''Dieses Program soll dazu dienen die
Radardaten von BoxPol mit den GPM Daten
hinsichtlich der Reflektivitat zu validieren.
Hier werden mehrere Ueberflug analysiert'''

#!/usr/bin/env python

import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import wradlib
import glob
import math
import pandas as pd
from scipy import stats
# ftp://ftp.meteo.uni-bonn.de/pub/pablosaa/gpmdata/


# Pfad mit String
# ---------------
hi = 1 # Hohe von DPR
TH = 0.01 #Threshold um Nullen fuer Niederschlag raus zu filtern
corra = []
error = []
corra2 = []
error2 = []
time = []
ipoli = [wradlib.ipol.Idw, wradlib.ipol.Linear, wradlib.ipol.Nearest, wradlib.ipol.OrdinaryKriging]
offset = 2
#LISTE der Ueberfluege des GPM mit Niederschlagevents
LISTE = ("20140629145000","20140629145925","20140921070500","20140921071058","20141007023744")#,"20141007023000")#bei 3 ohne "20150128171500", bei 2 ohne ,"20141016001500" ,schlecht:"20140826150322","20141016001500","20140826145000","20141016002458"
LISTE=sorted(LISTE)

for i in range(0,len(LISTE)):
    ZP=LISTE[i]
    print "initialize ::: " + ZP

    year = ZP[0:4]
    m = ZP[4:6]
    d = ZP[6:8]
    ht = ZP[8:10]
    mt = ZP[10:12]
    st = ZP[12:14]



    pfad = ('/home/velibor/shkgpm/data/' + year + m + d + '/radar/*.HDF5')
    pfad_radar = sorted(glob.glob(pfad))
    pfad_radar_Ku = pfad_radar[2]
    pfad_radar_DPR = pfad_radar[0]
    pfad_radar_DPRGMI = pfad_radar[3]
    pfad_radar_Ka = pfad_radar[1]

    pfad1 = ('/home/velibor/shkgpm/data/' + year + m + d + '/boxpol/gpm_rhi_01/*.mvol')
    pfad_boxpol = glob.glob(pfad1)
    pfad_boxpol_rhi01 = pfad_boxpol[0]

    pfad2 = ('/home/velibor/shkgpm/data/' + year + m + d + '/gprof/*.HDF5')
    pfad_gprof = glob.glob(pfad2)
    pfad_gprof_g = pfad_gprof[0]

    ppi_datapath=('/automount/radar-archiv/scans/' + year+ "/" + year +"-"+ m + "/" + year+ "-" + m +"-"+ d + "/ppi_1p5deg/"+ year + "-" + m +"-"+ d + "--" +ht +":"+mt+":"+st+",00.mvol")


    # Rhi BoxPol Daten einlesen
    # --------------------------

    gpmrhi01 = h5py.File(pfad_boxpol_rhi01, 'r')

    # PPI BoxPol Daten einlesen
    #---------------------------

    ppi=h5py.File(ppi_datapath,'r')
    data, attrs = wradlib.io.read_GAMIC_hdf5(ppi_datapath)

    ZH = data['SCAN0']['ZH']['data']
    PHIDP = data['SCAN0']['PHIDP']['data']
    r = attrs['SCAN0']['r']
    az = attrs['SCAN0']['az']
    lon_ppi = attrs['VOL']['Longitude']
    lat_ppi = attrs['VOL']['Latitude']
    alt_ppi = attrs['VOL']['Height']
    #Umwandeln von Z in RR Marshal-Palmer Z(R)
    ZH = ZH + offset
    Z = wradlib.trafo.idecibel(ZH)
    R = wradlib.zr.z2r(Z, a=200., b=1.6)
    R[151:165]=np.nan



    # DPR Einlesen
    # ------------
    gpmdpr = h5py.File(pfad_radar_DPR, 'r')
    gpmdpr_HS=gpmdpr['HS']['SLV']
    dpr_lat=np.array(gpmdpr['HS']['Latitude'])			#(7934, 24)
    dpr_lon=np.array(gpmdpr['HS']['Longitude'])			#(7934, 24)
    #dpr_pp=np.array(gpmdpr_HS['zFactorCorrected'])   	    #(7934, 24, 88)
    dpr_pp=np.array(gpmdpr_HS['precipRateNearSurface'])
    # Lon Lat Bestimmung
    # ------------------

    bonn_lat1 = 49.9400
    bonn_lat2 = 51.3500
    bonn_lon1 = 6.40000
    bonn_lon2 = 8.10000

    ilat= np.where((dpr_lat>49.9400) & (dpr_lat<51.3500))
    ilon= np.where((dpr_lon>6.40000) & (dpr_lon<8.10000))
    lonstart = ilon[0][0]
    lonend = ilon[0][-1]
    latstart = ilat[0][0]
    latend = ilat[0][-1]
    dpr_pp[dpr_pp==-9999] = np.nan
    #dpr_pp = dpr_pp[:,:,hi]  # Untersete Schicht

    radar_location = (lon_ppi, lat_ppi, alt_ppi)
    elevation = 1.5
    azimuths = az
    ranges = r
    polargrid = np.meshgrid(ranges, azimuths)
    lon, lat, alt = wradlib.georef.polar2lonlatalt_n(polargrid[0], polargrid[1], elevation, radar_location)

    gk3 = wradlib.georef.epsg_to_osr(31467)
    x, y = wradlib.georef.reproject(lon, lat, projection_target=gk3)
    xgrid, ygrid = wradlib.georef.reproject(dpr_lon[latstart:latend], dpr_lat[latstart:latend], projection_target=gk3)

    grid_xy = np.vstack((xgrid.ravel(), ygrid.ravel())).transpose()

    xy=np.concatenate([x.ravel()[:,None],y.ravel()[:,None]], axis=1)
    gridded = wradlib.comp.togrid(xy, grid_xy, ranges[-1], np.array([x.mean(), y.mean()]), R.ravel(), ipoli[0],nnearest=5,p=2)
    gridded = np.ma.masked_invalid(gridded).reshape(xgrid.shape)


    # Plot
    # ----
    fig = plt.figure(figsize=(13,10))
    maxv = np.max([np.max(np.ma.masked_invalid(dpr_pp)[latstart:latend]),np.nanmax(gridded)])

    plt.subplot(221)  # ==== Scatterplot GPM/boxpol ==== #
    A = gridded
    B = np.ma.masked_invalid(dpr_pp)[latstart:latend]
    A[A<TH]=np.nan
    B[B<TH]=np.nan
    #A[A>130]=np.nan
    #B[B>130]=np.nan

    mask = ~np.isnan(B) & ~np.isnan(A)
    slope, intercept, r_value, p_value, std_err = stats.linregress(B[mask], A[mask])
    line = slope*B+intercept
    plt.plot(B,line,'r-',B,A,'ob')
    maxAB = np.nanmax([np.nanmax(A),np.nanmax(B)])
    plt.xlim(0,maxAB + 1)
    plt.ylim(0,maxAB + 1)
    plt.xlabel("DPR RR [mm/h]")
    plt.ylabel("BoxPol RR [mm/h]")
    plt.grid(True)
    plt.title("Scatterplot (gprof/ppi), cor: " + str(r_value))

    plt.subplot(222)  # ==== RainRate boxpol ==== #
    ax1, pm2 = wradlib.vis.plot_ppi(R,r,az,vmin=0,vmax=maxv)
    cbar = plt.colorbar(pm2, shrink=0.75)
    cbar.set_label("RainRate Boxpol [mm/h]")
    plt.xlim((-101000,101000))
    plt.ylim((-101000,101000))
    plt.xticks(())
    plt.yticks(())
    plt.xlabel("X Range [km]")
    plt.ylabel("Y Range [km]")
    plt.title(ppi_datapath[-28:-8])

    plt.subplot(224)
    pm2 = plt.pcolormesh(dpr_lon[latstart:latend], dpr_lat[latstart:latend], gridded,vmin=0,vmax=maxv)
    plt.xlim((bonn_lon1,bonn_lon2))
    plt.ylim((bonn_lat1,bonn_lat2))
    plt.title(ppi_datapath[-28:-8])
    cbar = plt.colorbar(pm2, shrink=0.75)
    cbar.set_label("Boxpol  interpolated [mm/h]")
    plt.xlabel("Easting (m)")
    plt.ylabel("Northing (m)")

    plt.subplot(223)
    pm2 = plt.pcolormesh(dpr_lon[latstart:latend], dpr_lat[latstart:latend],np.ma.masked_invalid(dpr_pp[latstart:latend])
                   ,vmin=0,vmax=maxv)
    plt.xlim((bonn_lon1,bonn_lon2))
    plt.ylim((bonn_lat1,bonn_lat2))
    plt.title(pfad_boxpol_rhi01[-28:-6])
    cbar = plt.colorbar(pm2, shrink=.75)
    cbar.set_label("DPR  [mm/h]")
    plt.xlabel("Easting (m)")
    plt.ylabel("Northing (m)")

    plt.show()


