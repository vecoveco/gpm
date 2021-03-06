 
'''Dieses Program soll dazu dienen die 
Radardaten von BoxPol mit den GPM Daten
hinsichtlich der Regenraten zu validieren.
Hier werden mehrere Ueberflug analysiert'''

#!/usr/bin/env python

#--------------------------------------------------------------------------------------------------------
'''Einlesen von Modulen'''
#--------------------------------------------------------------------------------------------------------
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import wradlib
import glob
import wradlib as wrl

'''ftp://ftp.meteo.uni-bonn.de/pub/pablosaa/gpmdata/'''

#--------------------------------------------------------------------------------------------------------
### Pfad mit String ##
#--------------------------------------------------------------------------------------------------------
TH = 0.01 #Threshold um Nullen fuer Niederschlag raus zu filtern

ipoli = [wradlib.ipol.Idw, wradlib.ipol.Linear, wradlib.ipol.Nearest, wradlib.ipol.OrdinaryKriging]

#LISTE der Ueberfluege des GPM mit Niederschlagevents
#BITTE ZEITRAUM DES PPI EINTRAGEN!
LISTE = ("20140629145000","20140629145925","20140921070500","20140921071058","20141007023744","20141007023000")#bei 3 ohne "20150128171500", bei 2 ohne ,"20141016001500" ,schlecht:"20140826150322","20141016001500","20140826145000","20141016002458"
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

#ppi_1p5deg,ppi_2p4deg, ppi_3p4deg
#--------------------------------------------------------------------------------------------------------
### Rhi BoxPol Daten einlesen ##
#--------------------------------------------------------------------------------------------------------

    gpmrhi01 = h5py.File(pfad_boxpol_rhi01, 'r')

#--------------------------------------------------------------------------------------------------------
### PPI BoxPol Daten einlesen ##
#--------------------------------------------------------------------------------------------------------

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
    Z = wradlib.trafo.idecibel(ZH)
    R = wradlib.zr.z2r(Z, a=200., b=1.6) 
    # Todo: Verbesserung Regenraten Berechnung
    # Todo: Daempfung enthalten ?
    #PHIDP = np.deg2rad(PHIDP)
    #kdp = wradlib.dp.kdp_from_phidp_linregress(PHIDP)
    #R_kdp = wradlib.trafo.kdp2r(kdp,10, a=129.0, b=0.85)
    #R = R_kdp#TESTWEISE!
#--------------------------------------------------------------------------------------------------------
### gprof Daten einlesen ##
#--------------------------------------------------------------------------------------------------------

    gpmgmi = h5py.File(pfad_gprof_g, 'r')

#--------------------------------------------------------------------------------------------------------
###---- GPROF ---- Einlesen von Oberflachenniederschlag und Lat/Lon von Gprof ##
#--------------------------------------------------------------------------------------------------------

    gpmgmi.keys()
    gpmgmi_S1=gpmgmi['S1']				
    gprof_lat=gpmgmi_S1['Latitude']			#(2962, 221)
    gprof_lon=gpmgmi_S1['Longitude']			#(2962, 221)
    gprof_pp=gpmgmi_S1['surfacePrecipitation']   	#(2962, 221)

#--------------------------------------------------------------------------------------------------------
### ------- In Arrays umwandeln ----------- ##
#--------------------------------------------------------------------------------------------------------

    gprof_pp_a = np.array(gprof_pp)
    gprof_lon_a = np.array(gprof_lon)
    gprof_lat_a = np.array(gprof_lat)

#--------------------------------------------------------------------------------------------------------
### ------- Lon Lat Bestimmung ----------- ##
#--------------------------------------------------------------------------------------------------------
# Indexzugriff
#   gprof_lat_a[[Zeile][Zeilenelement]]
# Indexsuche
#   itemin= np.where((gprof_lat_a<58) & (gprof_lat_a>57))
    bonn_lat1 = 49.9400
    bonn_lat2 = 51.3500
    bonn_lon1 = 6.40000
    bonn_lon2 = 8.10000

    ilat= np.where((gprof_lat_a>49.9400) & (gprof_lat_a<51.3500))
    ilon= np.where((gprof_lon_a>6.40000) & (gprof_lon_a<8.10000))
    lonstart = ilon[0][0]	#erstes Element
    lonend = ilon[0][-1]	#letztes Element
    latstart = ilat[0][0]
    latend = ilat[0][-1]



    gprof_pp_a[gprof_pp_a==-9999] = np.nan #NAN!
    radar_location = (lon_ppi, lat_ppi, alt_ppi) # (lon, lat, alt) in decimal degree and meters
    elevation = 1.5 # in degree
    azimuths = az # in degrees
    ranges = r # in meters
    polargrid = np.meshgrid(ranges, azimuths)
    lon, lat, alt = wradlib.georef.polar2lonlatalt_n(polargrid[0], polargrid[1], elevation, radar_location)

    gk3 = wradlib.georef.epsg_to_osr(31467)
    x, y = wradlib.georef.reproject(lon, lat, projection_target=gk3)
    xgrid, ygrid = wradlib.georef.reproject(gprof_lon_a[latstart:latend], gprof_lat_a[latstart:latend], projection_target=gk3)

    grid_xy = np.vstack((xgrid.ravel(), ygrid.ravel())).transpose()

    xy=np.concatenate([x.ravel()[:,None],y.ravel()[:,None]], axis=1)
    gridded = wradlib.comp.togrid(xy, grid_xy, ranges[-1], np.array([x.mean(), y.mean()]), R.ravel(), ipoli[0],nnearest=500,p=2)
    gridded = np.ma.masked_invalid(gridded).reshape(xgrid.shape)
    # ========== PLOT ========== #
    fig = plt.figure(figsize=(13,10))

    maxvl = np.max([np.max(np.log10(gridded)),np.max(np.log10(np.ma.masked_invalid(gprof_pp_a)[latstart:latend]))])
    maxv = np.max([np.max(gridded),np.max(np.ma.masked_invalid(gprof_pp_a)[latstart:latend])])

    plt.subplot(221) # ==== Scatterplot GPM/boxpol ==== #
    A = gridded
    B = np.ma.masked_invalid(gprof_pp_a)[latstart:latend]
    A[A<TH]=np.nan
    B[B<TH]=np.nan

    from scipy import stats
    mask = ~np.isnan(B) & ~np.isnan(A)
    slope, intercept, r_value, p_value, std_err = stats.linregress(B[mask], A[mask])
    line = slope*B+intercept

    #Todo:RMSE und BIAS berechene
    #from sklearn.metrics import mean_squared_error
    #RMSE = np.sqrt(np.mean((B[mask], A[mask])**2))
    #print('RMSE:', RMSE)
    from scipy import signal
    corr = signal.correlate(B[mask], A[mask], mode='same')


    plt.plot(B[mask], label='Gprof RR [mm/h]')
    plt.plot(A[mask], label='BoxPol RR [mm/h]')

    legend = plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=2, fancybox=True, shadow=True,
                        fontsize='small', title="GPROF_vs_BoxPol r: " + str(round(r_value,3)) +
                                                r'$\pm$' + str(round(std_err,3)) + ', p:' +  str(round(p_value,3)))
    plt.xlabel("Pixel Position")
    plt.ylabel("RR [mm/h]")
    plt.grid(True)

    plt.subplot(222) # ==== RainRate boxpol ==== #

    plt.scatter(B,A, color='gray', label='RR [mm/h]')
    plt.plot(B,line,'r-')
    plt.plot(B,line,'r-')
    maxAB = np.nanmax([np.nanmax(A),np.nanmax(B)])
    plt.xlim(0,maxAB + 1)
    plt.ylim(0,maxAB + 1)
    legend = plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=2, fancybox=True, shadow=True,
                        fontsize='small', title="GPROF_vs_BoxPol  Slope: " + str(round(slope,3)) + ' Intercept: '
                                                +  str(round(intercept,3)))
    plt.xlabel("GPROF RR [mm/h]")
    plt.ylabel("BoxPol RR [mm/h]")
    plt.title(" .")
    plt.grid(True)

    plt.subplot(223) # ==== RainRate Gprof ==== #
    pm2 = plt.pcolormesh(gprof_lon_a[latstart:latend], gprof_lat_a[latstart:latend], np.ma.masked_invalid(
        gprof_pp_a)[latstart:latend],vmin=0,vmax=maxv)
    plt.xlim((bonn_lon1,bonn_lon2))
    plt.ylim((bonn_lat1,bonn_lat2))
    plt.title(pfad_boxpol_rhi01[-28:-6])
    cbar = plt.colorbar(pm2, shrink=.75)
    cbar.set_label("GPROF RainRate [mm/h]")
    plt.xlabel("Easting ")
    plt.ylabel("Northing ")

    plt.subplot(224) # ==== RainRate Boxpol interpolation in GPROF Grid  ==== #
    pm2 = plt.pcolormesh(gprof_lon_a[latstart:latend], gprof_lat_a[latstart:latend], gridded,vmin=0,vmax=maxv)
    plt.xlim((bonn_lon1,bonn_lon2))
    plt.ylim((bonn_lat1,bonn_lat2))
    plt.title(ppi_datapath[-28:-8])
    cbar = plt.colorbar(pm2, shrink=0.75)
    cbar.set_label("Boxpol RainRate interpolated [mm/h]")
    plt.xlabel("Easting ")
    plt.ylabel("Northing ")
    plt.tight_layout()
    plt.savefig('/user/velibor/SHKGPM/data/plot/' + ppi_datapath[-28:-8] + '_Vergleich.png')
    plt.close()







