
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
from pcc import get_miub_cmap as my_cmap
from pcc import plot_radar
from pcc import boxpol_pos
from pcc import plot_borders

import wradlib as wrl
from osgeo import osr

#Pos = boxpol_pos()
#blon0, blat0 = Pos['lon_ppi'], Pos['lat_ppi']
#bbx, bby = Pos['gkx_ppi'], Pos['gky_ppi']

# Pfad mit String
# ---------------

# Hohe von DPR
TH = 15 #Threshold um Nullen fuer Niederschlag raus zu filtern

ipoli = [wradlib.ipol.Idw, wradlib.ipol.Linear, wradlib.ipol.Nearest, wradlib.ipol.OrdinaryKriging]
offset = 2


ZP = '20141007023744' ; pfadnr=0# 0.47
#ZP = '20140826220500'; pfadnr=1 # 0.82
#ZP = '20141008094000'; pfadnr=1 # 0.82   #!!!!!!!!!!!!!!NICE
#ZP = '20141008094500'; pfadnr=1 # 0.679  #!!!!!!!!!!!!!!NICE
#ZP = '20150128171500'; pfadnr=0 #0.28
#ZP = '20150128172208'; pfadnr=0#0.321
#ZP = '20160209103500'; pfadnr=1 # 0.23
#ZP = '20151216024501'; pfadnr=0#0.589
#ZP = '20151216023500' ; pfadnr=0# 0.651
#ZP = '20160209103000'; pfadnr=1 ###PFAD=1


year = ZP[0:4]
m = ZP[4:6]
d = ZP[6:8]
ht = ZP[8:10]
mt = ZP[10:12]
st = ZP[12:14]


pfad_radar = glob.glob('/automount/ags/velibor/gpmdata/dpr/2A.GPM.DPR.V6-20160118.' + year + m + d + '*.HDF5')

print pfad_radar
pfad_radar = pfad_radar[pfadnr]
#pfad_radar_Ku = pfad_radar[0]

pfad_radar = glob.glob('/automount/ags/velibor/gpmdata/dpr/2A.GPM.DPR.V6-20160118.' + year + m + d + '*.HDF5')

print pfad_radar
pfad_radar = pfad_radar[pfadnr]
#pfad_radar_Ku = pfad_radar[0]

deg_scan =  ["/ppi_1p5deg/","/ppi_2p4deg/","/ppi_3p4deg/",
             "/n_ppi_010deg/","/n_ppi_045deg/",
             "/n_ppi_082deg/","/n_ppi_110deg/","/n_ppi_140deg/",
             "/n_ppi_180deg/","/n_ppi_280deg/","/n_vertical_scan/"][0]

#pfad_radar = glob.glob('/automount/radar-archiv/scans_juelich/2014/2014-10/2014-10-08/rainscanner.wuestebach/2014100809390000dBZ.azi')
#pfad_radar = glob.glob('/automount/radar-archiv/scans_juelich/2014/2014-10/2014-10-08/DWD_Vol_2/2014-10-08--09:40:00,00.mvol')
try:
    ppi_datapath=glob.glob('/automount/radar-archiv/scans/' + year+ "/" +
                           year +"-"+ m + "/" + year+ "-" + m +"-"+ d +
                           deg_scan+ year + "-" + m +"-"+ d + "--" +ht +
                           ":"+mt+":"+st+",*.mvol")
    print ppi_datapath
    ppi_datapath = ppi_datapath[0]

except:
    ppi_datapath=glob.glob('/automount/radar/scans/' + year+ "/" +
                           year +"-"+ m + "/" + year+ "-" + m +"-"+
                           d + deg_scan+ year + "-" + m +"-"+ d +
                           "--" +ht +":"+mt+":"+st+",*.mvol")
    print ppi_datapath
    ppi_datapath = ppi_datapath[0]


#pfad_radar = glob.glob('/automount/radar-archiv/scans_juelich/2014/2014-10/2014-10-08/rainscanner.wuestebach/2014100809390000dBZ.azi')
#pfad_radar = glob.glob('/automount/radar-archiv/scans_juelich/2014/2014-10/2014-10-08/DWD_Vol_2/2014-10-08--09:40:00,00.mvol')

ppi_datapath2 = '/automount/radar-archiv/scans_juelich/2014/2014-10/2014-10-07/DWD_Vol_2/2014-10-07--02:35:00,00.mvol'
# PPI BoxPol Daten einlesen
#---------------------------

ppi=h5py.File(ppi_datapath,'r')
data, attrs = wradlib.io.read_GAMIC_hdf5(ppi_datapath)

ZH0 = data['SCAN0']['ZH']['data']
PHIDP = data['SCAN0']['PHIDP']['data']
r = attrs['SCAN0']['r']
az = attrs['SCAN0']['az']
lon_ppi = attrs['VOL']['Longitude']
lat_ppi = attrs['VOL']['Latitude']
alt_ppi = attrs['VOL']['Height']
rho = data['SCAN0']['RHOHV']['data']
zdr = data['SCAN0']['ZDR']['data']



ppi2=h5py.File(ppi_datapath2,'r')
data2, attrs2 = wradlib.io.read_GAMIC_hdf5(ppi_datapath2)
scan2 = 'SCAN8'
ZH2 = data2[scan2]['ZH']['data']
PHIDP2 = data2[scan2]['PHIDP']['data']
r2 = attrs2[scan2]['r']
az2 = attrs2[scan2]['az']
lon_ppi2 = attrs2['VOL']['Longitude']
lat_ppi2 = attrs2['VOL']['Latitude']
alt_ppi2 = attrs2['VOL']['Height']
rho2 = data2[scan2]['RHOHV']['data']
zdr2 = data2[scan2]['ZDR']['data']

#R[151:165]=np.nan


### Threshold for DPR sensitivity
#R[R<TH]=np.nan
R, R2 = ZH0, ZH2
# DPR Einlesen
# ------------
gpmku = h5py.File(pfad_radar, 'r')
gpmku_HS = gpmku['NS']['SLV']
dpr_lat = np.array(gpmku['NS']['Latitude'])
dpr_lon = np.array(gpmku['NS']['Longitude'])
dpr_pp = np.array(gpmku_HS['zFactorCorrectedNearSurface'])
dpr_pp[dpr_pp < 0] = np.nan

# Cut the Swath
from pcc import cut_the_swath
dpr_lon, dpr_lat, dpr_pp = cut_the_swath(dpr_lon,dpr_lat,dpr_pp, eu=2)

# Koordinaten Projektion
# ------------------

proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)



dpr_lon, dpr_lat = wradlib.georef.reproject(dpr_lon, dpr_lat, projection_target=proj_stereo , projection_source=proj_wgs)
blon, blat = wradlib.georef.reproject(lon_ppi, lat_ppi, projection_target=proj_stereo , projection_source=proj_wgs)
blon2, blat2 = wradlib.georef.reproject(lon_ppi2, lat_ppi2, projection_target=proj_stereo , projection_source=proj_wgs)


# Dpr zuschneiden
#-----------------

lon0, lat0, radius = blon, blat, 100
rr = np.sqrt((dpr_lat - lat0)**2 + (dpr_lon - lon0)**2)
position = rr < radius

pp = dpr_pp.copy()

pp[np.where(rr > radius)] = np.nan
##################################################################


radar_location = (lon_ppi, lat_ppi, alt_ppi)
elevation = attrs['SCAN0']['el']

azimuths = az
ranges = r
polargrid = np.meshgrid(ranges, azimuths)
lon, lat, alt = wradlib.georef.polar2lonlatalt_n(polargrid[0], polargrid[1], elevation, radar_location)
lon, lat = wradlib.georef.reproject(lon, lat, projection_target=proj_stereo , projection_source=proj_wgs)

grid_xy = np.vstack((dpr_lon.ravel(), dpr_lat.ravel())).transpose()

xy=np.concatenate([lon.ravel()[:,None],lat.ravel()[:,None]], axis=1)

gridded = wradlib.comp.togrid(xy, grid_xy, ranges[-1], np.array([lon.mean(), lat.mean()]), R.ravel(), ipoli[0],nnearest=40,p=2)
gridded = np.ma.masked_invalid(gridded).reshape(dpr_lon.shape)
gridded[np.where(rr > radius)]=np.nan





#######juxpol
radar_location2 = (lon_ppi2, lat_ppi2, alt_ppi2)
elevation2 = attrs2[scan2]['el']

azimuths2 = az2
ranges2 = r2
polargrid2 = np.meshgrid(ranges2, azimuths2)
lon2, lat2, alt2 = wradlib.georef.polar2lonlatalt_n(polargrid2[0], polargrid2[1], elevation2, radar_location2)
lon2, lat2 = wradlib.georef.reproject(lon2, lat2, projection_target=proj_stereo , projection_source=proj_wgs)

grid_xy2 = np.vstack((dpr_lon.ravel(), dpr_lat.ravel())).transpose()

xy2=np.concatenate([lon2.ravel()[:,None],lat2.ravel()[:,None]], axis=1)

gridded2 = wradlib.comp.togrid(xy2, grid_xy2, ranges2[-1], np.array([lon2.mean(), lat2.mean()]), R2.ravel(), ipoli[0],nnearest=40,p=2)
gridded2 = np.ma.masked_invalid(gridded2).reshape(dpr_lon.shape)
gridded2[np.where(rr > radius)]=np.nan







fig = plt.figure(figsize=(18,12))
fig.suptitle('BoXPol vs DPR '+ZP+' Rho_th: ')#+str(rho_th))

###################
ax1 = fig.add_subplot(231, aspect='equal')
plt.pcolormesh(dpr_lon, dpr_lat,np.ma.masked_invalid(pp),vmin=0, vmax=40, cmap=my_cmap())

plt.colorbar()
plot_borders(ax1)
plot_radar(lon_ppi, lat_ppi, ax1, reproject=True, cband=False,col='black')
plt.plot(dpr_lon[:,0],dpr_lat[:,0], color='black',lw=1)
plt.plot(dpr_lon[:,-1],dpr_lat[:,-1], color='black',lw=1)
plt.plot(dpr_lon[:,dpr_lon.shape[1]/2],dpr_lat[:,dpr_lon.shape[1]/2], color='black',lw=1, ls='--')

plt.xlim(-350,-100)
plt.ylim(-4350, -4100)
plt.title('GPM - DPR')

plt.grid()

ax2 = fig.add_subplot(232, aspect='equal')
plt.pcolormesh(dpr_lon, dpr_lat,np.ma.masked_invalid(gridded),vmin=0, vmax=40, cmap=my_cmap())

plt.colorbar()
plot_borders(ax2)
plot_radar(lon_ppi, lat_ppi, ax2, reproject=True, cband=False,col='black')
plt.plot(dpr_lon[:,0],dpr_lat[:,0], color='black',lw=1)
plt.plot(dpr_lon[:,-1],dpr_lat[:,-1], color='black',lw=1)
plt.plot(dpr_lon[:,dpr_lon.shape[1]/2],dpr_lat[:,dpr_lon.shape[1]/2], color='black',lw=1, ls='--')

plt.xlim(-350,-100)
plt.ylim(-4350, -4100)
plt.title('BoXPol - onDPR\n' +attrs['SCAN0']['Time'])

plt.grid()

ax3 = fig.add_subplot(233, aspect='equal')
plt.pcolormesh(lon, lat,R,vmin=0, vmax=40, cmap=my_cmap())
plt.colorbar()
plot_borders(ax3)
plot_radar(lon_ppi, lat_ppi, ax3, reproject=True, cband=False,col='black')
plt.plot(dpr_lon[:,0],dpr_lat[:,0], color='black',lw=1)
plt.plot(dpr_lon[:,-1],dpr_lat[:,-1], color='black',lw=1)
plt.plot(dpr_lon[:,dpr_lon.shape[1]/2],dpr_lat[:,dpr_lon.shape[1]/2], color='black',lw=1, ls='--')

plt.xlim(-350,-100)
plt.ylim(-4350, -4100)
plt.title('BoXPol\n' +attrs['SCAN0']['Time'])

plt.grid()

ax4 = fig.add_subplot(236, aspect='equal')

plt.pcolormesh(lon2, lat2,R2,vmin=0, vmax=40, cmap=my_cmap())
plt.colorbar()
plot_borders(ax4)
plot_radar(lon_ppi2, lat_ppi2, ax4, reproject=True, cband=False,col='black')
plt.plot(dpr_lon[:,0],dpr_lat[:,0], color='black',lw=1)
plt.plot(dpr_lon[:,-1],dpr_lat[:,-1], color='black',lw=1)
plt.plot(dpr_lon[:,dpr_lon.shape[1]/2],dpr_lat[:,dpr_lon.shape[1]/2], color='black',lw=1, ls='--')

plt.xlim(-350,-100)
plt.ylim(-4350, -4100)
plt.title('JUXPOL \n' +attrs2[scan2]['Time'])
plt.grid()

ax5 = fig.add_subplot(235, aspect='equal')
plt.pcolormesh(dpr_lon, dpr_lat,np.ma.masked_invalid(gridded2),vmin=0, vmax=40, cmap=my_cmap())

plt.colorbar()
plot_borders(ax5)
plot_radar(lon_ppi2, lat_ppi2, ax5, reproject=True, cband=False,col='black')
plt.plot(dpr_lon[:,0],dpr_lat[:,0], color='black',lw=1)
plt.plot(dpr_lon[:,-1],dpr_lat[:,-1], color='black',lw=1)
plt.plot(dpr_lon[:,dpr_lon.shape[1]/2],dpr_lat[:,dpr_lon.shape[1]/2], color='black',lw=1, ls='--')

plt.xlim(-350,-100)
plt.ylim(-4350, -4100)
plt.title('JUXPOL - onDPR \n' +attrs2['SCAN8']['Time'])

plt.grid()

#plt.savefig('/automount/ags/velibor/plot/boxpol/boxpol_vs_DPR/boxdpr_'+ZP)


ax44 = fig.add_subplot(234, aspect='equal')

plt.pcolormesh(lon, lat,R,vmin=0, vmax=40, cmap=my_cmap(), zorder=2)

plt.pcolormesh(lon2, lat2,R2,vmin=0, vmax=40, cmap=my_cmap(),zorder=1)

plt.colorbar()
plot_borders(ax44)
plot_radar(lon_ppi, lat_ppi, ax44, reproject=True, cband=False,col='black')
plot_radar(lon_ppi2, lat_ppi2, ax44, reproject=True, cband=False,col='blue')

plt.plot(dpr_lon[:,0],dpr_lat[:,0], color='black',lw=1)
plt.plot(dpr_lon[:,-1],dpr_lat[:,-1], color='black',lw=1)
plt.plot(dpr_lon[:,dpr_lon.shape[1]/2],dpr_lat[:,dpr_lon.shape[1]/2], color='black',lw=1, ls='--')

plt.xlim(-350,-100)
plt.ylim(-4350, -4100)
plt.title('JUXPOL  BoXPol' )
plt.grid()

plt.show()


