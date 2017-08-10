
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

Pos = boxpol_pos()
blon0, blat0 = Pos['lon_ppi'], Pos['lat_ppi']
bbx, bby = Pos['gkx_ppi'], Pos['gky_ppi']

# Pfad mit String
# ---------------

# Hohe von DPR
TH = 18 #Threshold um Nullen fuer Niederschlag raus zu filtern

ipoli = [wradlib.ipol.Idw, wradlib.ipol.Linear, wradlib.ipol.Nearest, wradlib.ipol.OrdinaryKriging]
offset = 2


#ZP = '20141007023744'
#ZP = '20140826220500'
#ZP = '20140921071000'
#ZP = '20141008094000'
#ZP = '20141008094500'
#ZP = '20150128171500'
#ZP = '20150128172208'
#ZP = '20160209103500'
#ZP = '20151216024501'
#ZP = '20151216023500'
ZP = '20160209103000'


year = ZP[0:4]
m = ZP[4:6]
d = ZP[6:8]
ht = ZP[8:10]
mt = ZP[10:12]
st = ZP[12:14]


pfad_radar = glob.glob('/automount/ags/velibor/gpmdata/dpr/2A.GPM.DPR.V6-20160118.' + year + m + d + '*.HDF5')
print pfad_radar
pfad_radar = pfad_radar[1]
#pfad_radar_Ku = pfad_radar[0]

deg_scan =  ["/ppi_1p5deg/","/ppi_2p4deg/","/ppi_3p4deg/",
             "/n_ppi_010deg/","/n_ppi_045deg/",
             "/n_ppi_082deg/","/n_ppi_110deg/","/n_ppi_140deg/",
             "/n_ppi_180deg/","/n_ppi_280deg/","/n_vertical_scan/"][0]


try:
    ppi_datapath=glob.glob('/automount/radar-archiv/scans/' + year+ "/" + year +"-"+ m + "/" + year+ "-" + m +"-"+ d + deg_scan+ year + "-" + m +"-"+ d + "--" +ht +":"+mt+":"+st+",*.mvol")
    print ppi_datapath
    ppi_datapath = ppi_datapath[0]
except:
    ppi_datapath=glob.glob('/automount/radar/scans/' + year+ "/" + year +"-"+ m + "/" + year+ "-" + m +"-"+ d + deg_scan+ year + "-" + m +"-"+ d + "--" +ht +":"+mt+":"+st+",*.mvol")
    print ppi_datapath
    ppi_datapath = ppi_datapath[0]


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
R = ZH0
R[151:165]=np.nan

#R = (R0 + R1)/2.

# DPR Einlesen
# ------------
gpmku = h5py.File(pfad_radar, 'r')
gpmku_HS = gpmku['NS']['SLV']
dpr_lat = np.array(gpmku['NS']['Latitude'])			#(7934, 24)
dpr_lon = np.array(gpmku['NS']['Longitude'])			#(7934, 24)
dpr_pp = np.array(gpmku_HS['zFactorCorrectedNearSurface'])
dpr_pp[dpr_pp < 0] = np.nan


# Koordinaten Projektion
# ------------------

proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

#from pcc import boxpol_pos
#bonn_pos = boxpol_pos()

dpr_lon, dpr_lat = wradlib.georef.reproject(dpr_lon, dpr_lat, projection_target=proj_stereo , projection_source=proj_wgs)
blon, blat = wradlib.georef.reproject(blon0, blat0, projection_target=proj_stereo , projection_source=proj_wgs)


# Dpr zuschneiden
#-----------------

lon0, lat0, radius = blon, blat, 100
rr = np.sqrt((dpr_lat - lat0)**2 + (dpr_lon - lon0)**2)
position = rr < radius

pp = dpr_pp.copy()

pp[np.where(rr > radius)] = np.nan


radar_location = (lon_ppi, lat_ppi, alt_ppi)
elevation = 1.5
azimuths = az
ranges = r
polargrid = np.meshgrid(ranges, azimuths)
lon, lat, alt = wradlib.georef.polar2lonlatalt_n(polargrid[0], polargrid[1], elevation, radar_location)
lon, lat = wradlib.georef.reproject(lon, lat, projection_target=proj_stereo , projection_source=proj_wgs)

grid_xy = np.vstack((dpr_lon.ravel(), dpr_lat.ravel())).transpose()

xy=np.concatenate([lon.ravel()[:,None],lat.ravel()[:,None]], axis=1)

gridded = wradlib.comp.togrid(xy, grid_xy, ranges[-1], np.array([lon.mean(), lat.mean()]), R.ravel(), ipoli[0],nnearest=80,p=2)
gridded = np.ma.masked_invalid(gridded).reshape(dpr_lon.shape)
gridded[np.where(rr > radius)]=np.nan



fig = plt.figure(figsize=(14,12))
fig.suptitle('BoXPol vs DPR '+ZP)

###################
ax1 = fig.add_subplot(221, aspect='auto')
plt.pcolormesh(dpr_lon, dpr_lat,np.ma.masked_invalid(pp),vmin=0, vmax=40, cmap=my_cmap())

plt.colorbar()
plot_borders(ax1)
plot_radar(blon0, blat0, ax1, reproject=True, cband=False,col='black')
plt.plot(dpr_lon[:,0],dpr_lat[:,0], color='black',lw=1)
plt.plot(dpr_lon[:,-1],dpr_lat[:,-1], color='black',lw=1)
plt.plot(dpr_lon[:,dpr_lon.shape[1]/2],dpr_lat[:,dpr_lon.shape[1]/2], color='black',lw=1, ls='--')

plt.xlim(-350,-100)
plt.ylim(-4350, -4100)
plt.title('GPM - DPR')
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
plt.grid()

ax2 = fig.add_subplot(222, aspect='auto')
plt.pcolormesh(dpr_lon, dpr_lat,np.ma.masked_invalid(gridded),vmin=0, vmax=40, cmap=my_cmap())

plt.colorbar()
plot_borders(ax2)
plot_radar(blon0, blat0, ax2, reproject=True, cband=False,col='black')
plt.plot(dpr_lon[:,0],dpr_lat[:,0], color='black',lw=1)
plt.plot(dpr_lon[:,-1],dpr_lat[:,-1], color='black',lw=1)
plt.plot(dpr_lon[:,dpr_lon.shape[1]/2],dpr_lat[:,dpr_lon.shape[1]/2], color='black',lw=1, ls='--')

plt.xlim(-350,-100)
plt.ylim(-4350, -4100)
plt.title('BoXPol - onDPR')
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
plt.grid()

ax3 = fig.add_subplot(223, aspect='auto')
plt.pcolormesh(lon, lat,R,vmin=0, vmax=40, cmap=my_cmap())
plt.colorbar()
plot_borders(ax3)
plot_radar(blon0, blat0, ax3, reproject=True, cband=False,col='black')
plt.plot(dpr_lon[:,0],dpr_lat[:,0], color='black',lw=1)
plt.plot(dpr_lon[:,-1],dpr_lat[:,-1], color='black',lw=1)
plt.plot(dpr_lon[:,dpr_lon.shape[1]/2],dpr_lat[:,dpr_lon.shape[1]/2], color='black',lw=1, ls='--')

plt.xlim(-350,-100)
plt.ylim(-4350, -4100)
plt.title('BoXPol - onDPR')
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
plt.grid()

ax4 = fig.add_subplot(224, aspect='auto')
A, B = gridded.copy(), pp.copy()
A[A<=10]=np.nan
B[B<=10]=np.nan
maske = ~np.isnan(A) & ~np.isnan(B)

slope, intercept, r_value, p_value, std_err = stats.linregress(A[maske], B[maske])
line = slope * A +intercept

plt.scatter(A[maske],B[maske], color='black')
plt.hist2d(A[maske],B[maske], bins=10, cmap='PuBu')
plt.title('r:'+str(round(r_value,3))+'  +  '+str(round(std_err,3)))
plt.colorbar()
plt.grid()
plt.show()


# Plot
# ----




