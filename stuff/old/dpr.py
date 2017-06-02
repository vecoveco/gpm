
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

import matplotlib.cm as cm
my_cmap = cm.get_cmap('jet',40)
my_cmap.set_under('lightgrey')
my_cmap.set_over('darkred')

# Pfad mit String
# ---------------

# Hohe von DPR
TH = 0.1 #Threshold um Nullen fuer Niederschlag raus zu filtern

ipoli = [wradlib.ipol.Idw, wradlib.ipol.Linear, wradlib.ipol.Nearest, wradlib.ipol.OrdinaryKriging]
offset = 2


ZP = '20140921070500'


year = ZP[0:4]
m = ZP[4:6]
d = ZP[6:8]
ht = ZP[8:10]
mt = ZP[10:12]
st = ZP[12:14]



pfad = ('/home/velibor/shkgpm/data/' + year + m + d + '/radar/*.HDF5')
print pfad
pfad_radar = sorted(glob.glob(pfad))
print pfad_radar

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

#Clutter
#clutter = wradlib.clutter.filter_gabella(Z, tr1=12, n_p=6, tr2=1.1)
#Z_no_clutter = wradlib.ipol.interpolate_polar(Z, clutter)

R = wradlib.zr.z2r(Z, a=200., b=1.6)
R[151:165]=np.nan



# DPR Einlesen
# ------------
gpmdprs = h5py.File(pfad_radar_DPR, 'r')
gpmdprs_HS=gpmdprs['HS']['SLV']
dprs_lat=np.array(gpmdprs['HS']['Latitude'])			#(7934, 24)
dprs_lon=np.array(gpmdprs['HS']['Longitude'])			#(7934, 24)
dprs_pp=np.array(gpmdprs_HS['precipRateNearSurface'])

gpmka = h5py.File(pfad_radar_Ka, 'r')
gpmka_HS=gpmka['HS']['SLV']
ka_lat=np.array(gpmka['HS']['Latitude'])			#(7934, 24)
ka_lon=np.array(gpmka['HS']['Longitude'])			#(7934, 24)
ka_pp=np.array(gpmka_HS['precipRateNearSurface'])

gpmku = h5py.File(pfad_radar_Ku, 'r')
gpmku_HS=gpmku['NS']['SLV']
ku_lat=np.array(gpmku['NS']['Latitude'])			#(7934, 24)
ku_lon=np.array(gpmku['NS']['Longitude'])			#(7934, 24)
ku_pp=np.array(gpmku_HS['precipRateNearSurface'])


gpmdprgmi = h5py.File(pfad_radar_DPRGMI, 'r')
gpmdprgmi_HS=gpmdprgmi['NS']
dprgmi_lat=np.array(gpmdprgmi['NS']['Latitude'])			#(7934, 24)
dprgmi_lon=np.array(gpmdprgmi['NS']['Longitude'])			#(7934, 24)
dpr_ppgmi=np.array(gpmdprgmi_HS['surfPrecipTotRate'])


# Lon Lat Bestimmung
# ------------------
radars = [dprs_pp, ka_pp, ku_pp, dpr_ppgmi]
rad_lat = [dprs_lat, ka_lat, ku_lat,dprgmi_lat]
rad_lon = [dprs_lon, ka_lon, ku_lon, dprgmi_lon]
rad_name = ['DPR','Ka','Ku','DPRGMI']

for ii in range(len(radars)):

    dpr_pp = radars[ii]
    dpr_lat = rad_lat[ii]
    dpr_lon = rad_lon[ii]
    radarname = rad_name[ii]

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
    gridded = wradlib.comp.togrid(xy, grid_xy, ranges[-1], np.array([x.mean(), y.mean()]), R.ravel(), ipoli[0],nnearest=30,p=2)
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

    mask = ~np.isnan(B) & ~np.isnan(A)
    slope, intercept, r_value, p_value, std_err = stats.linregress(B[mask], A[mask])

    line = slope*B+intercept
    plt.scatter(B,A, color='blue', label='RR [mm/h]')
    plt.plot(B,line,'r-')
    maxAB = np.nanmax([np.nanmax(A),np.nanmax(B)])
    plt.xlim(0,maxAB + 1)
    plt.ylim(0,maxAB + 1)
    legend = plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=2, fancybox=True, shadow=True,
                        fontsize='small', title="________"+str(radarname)+"_vs_BoxPol________" + "\n Slope: " + str(round(slope,3))
                                                + ', Intercept: '+  str(round(intercept,3)) + "\n Correlation: " +
                                                str(round(r_value,3)) + ', Std_err: '+  str(round(std_err,3)))
    plt.xlabel(str(radarname) + "RR [mm/h]")
    plt.ylabel("BoxPol RR [mm/h]")
    plt.title(" .")
    plt.grid(True)

    plt.subplot(222)  # ==== RainRate boxpol ==== #
    ax1, pm2 = wradlib.vis.plot_ppi(R,r,az,vmin=0.01,vmax=maxv, cmap=my_cmap)
    cbar = plt.colorbar(pm2, shrink=0.75)
    cbar.set_label("RainRate Boxpol [mm/h]")
    plt.xlim((-101000,101000))
    plt.ylim((-101000,101000))
    plt.xticks(())
    plt.yticks(())
    plt.xlabel("X Range")
    plt.ylabel("Y Range")
    plt.title(ppi_datapath[-28:-8])

    plt.subplot(224)
    pm2 = plt.pcolormesh(dpr_lon[latstart:latend], dpr_lat[latstart:latend], gridded,vmin=0.01,vmax=maxv, cmap=my_cmap)
    plt.xlim((bonn_lon1,bonn_lon2))
    plt.ylim((bonn_lat1,bonn_lat2))
    plt.title(ppi_datapath[-28:-8])
    cbar = plt.colorbar(pm2, shrink=0.75)
    cbar.set_label("Boxpol  interpolated [mm/h]")
    plt.xlabel("Easting")
    plt.ylabel("Northing")

    plt.subplot(223)
    pm2 = plt.pcolormesh(dpr_lon[latstart:latend], dpr_lat[latstart:latend],np.ma.masked_invalid(dpr_pp[latstart:latend])
                   ,vmin=0.01,vmax=maxv, cmap=my_cmap)
    plt.xlim((bonn_lon1,bonn_lon2))
    plt.ylim((bonn_lat1,bonn_lat2))
    plt.title(pfad_boxpol_rhi01[-28:-6])
    cbar = plt.colorbar(pm2, shrink=.75)
    cbar.set_label(str(radarname) + "[mm/h]")
    plt.xlabel("Easting")
    plt.ylabel("Northing")

    plt.tight_layout()
    plt.show()
    #plt.savefig('/user/velibor/SHKGPM/data/plot/DPR_'+str(radarname)+'_boxpol_'+ ppi_datapath[-28:-8] + 'Vergleich.png')
    #plt.close()



