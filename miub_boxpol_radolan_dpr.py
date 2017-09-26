
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
TH = 15 #Threshold um Nullen fuer Niederschlag raus zu filtern

ipoli = [wradlib.ipol.Idw, wradlib.ipol.Linear, wradlib.ipol.Nearest, wradlib.ipol.OrdinaryKriging]
offset = 2


#ZP = '20141007023744' ; pfadnr=0# 0.47
#ZP = '20140826220500'; pfadnr=1 # 0.82
ZP = '20141008094000'; pfadnr=1 # 0.82   #!!!!!!!!!!!!!!NICE
#ZP = '20141008094500'; pfadnr=1 # 0.679  #!!!!!!!!!!!!!!NICE
#ZP = '20150128171500'; pfadnr=0 #0.28
#ZP = '20150128172208'; pfadnr=0#0.321
#ZP = '20160209103500'; pfadnr=1 # 0.23
#ZP = '20151216024501'; pfadnr=0#0.589
#ZP = '20151216023500' ; pfadnr=0# 0.651
#ZP = '20160209103000'; pfadnr=1 ###PFAD=1


year = ZP[0:4]
ye = ZP[2:4]
m = ZP[4:6]
d = ZP[6:8]
ht = ZP[8:10]
mt = ZP[10:12]
st = ZP[12:14]


pfad_radar = glob.glob('/automount/ags/velibor/gpmdata/dpr/2A.GPM.DPR.V6-20160118.' + year + m + d + '*.HDF5')
print pfad_radar
pfad_radar = pfad_radar[pfadnr]
#pfad_radar_Ku = pfad_radar[0]

deg_scan =  ["/ppi_1p5deg/","/ppi_2p4deg/","/ppi_3p4deg/",
             "/n_ppi_010deg/","/n_ppi_045deg/",
             "/n_ppi_082deg/","/n_ppi_110deg/","/n_ppi_140deg/",
             "/n_ppi_180deg/","/n_ppi_280deg/","/n_vertical_scan/"][0]


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


#################################################### PPI BoxPol Daten einlesen
#------------------------------------------------------------------------------

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

R = ZH0
R[151:165]=np.nan


print ("________ATTCORR______")
pia_harrison = wrl.atten.correctAttenuationHB(
    R,
    coefficients = dict(a=4.57e-5, b=0.731, l=1.0),
    mode="warn",
    thrs=59.)
pia_harrison[pia_harrison > 4.8] = 4.8

print ("________ATTCORR2______")
R = R + pia_harrison

print ("________CLUTTER______")
rho_th  = 0.85
R[rho<= rho_th] = np.nan################WARUM GEHT DAS NICHT ?

print ("________ofset______")
#R = R + 2
#?
print ("________beambl.______")

print ("________DPR Threshold______")
### Threshold for DPR sensitivity
R[R<TH]=np.nan

################################################################# DPR Einlesen
# -----------------------------------------------------------------------------
gpmku = h5py.File(pfad_radar, 'r')
gpmku_HS = gpmku['NS']['SLV']
dpr_lat = np.array(gpmku['NS']['Latitude'])
dpr_lon = np.array(gpmku['NS']['Longitude'])
dpr_pp = np.array(gpmku_HS['zFactorCorrectedNearSurface'])
dpr_pp[dpr_pp < 0] = np.nan


############################################################## RADOLAN einlesen
# -----------------------------------------------------------------------------

mtt = mt
mtt = str(int(round(float(mtt)/5.0)*5.0))

if mtt == '0':
    mtt = '00'
elif mtt == '5':
    mtt = '05'
print mtt

r_pro = 'rx'

pfad = ('/automount/radar/dwd/'+ r_pro +'/'+str(year)+'/'+str(year)+'-'+
        str(m)+'/'+ str(year)+'-'+str(m)+'-'+str(d)+'/raa01-'+r_pro+'_10000-'+
        str(ye)+str(m)+ str(d)+str(ht)+str(mtt)+'-dwd---bin.gz')

pfad_radolan = pfad[:-3]

try:
    rw_filename = wradlib.util.get_wradlib_data_file(pfad)
except EnvironmentError:
    rw_filename = wradlib.util.get_wradlib_data_file(pfad_radolan)

rwdata, rwattrs = wradlib.io.read_RADOLAN_composite(rw_filename)

radolan_zeit = rwattrs['datetime'].strftime("%Y.%m.%d -- %H:%M:%S")
#Binaere Grid
rn = rwdata.copy()
rn[rn != -9999] = 1
rn[rn == -9999] = 0

radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
x = radolan_grid_xy[:,:,0]
y = radolan_grid_xy[:,:,1]
rwdata = np.ma.masked_equal(rwdata, -9999) / 2 - 32.5

### Threshold for DPR sensitivity
rwdata[rwdata<TH]=-9999



######################################################## Cut the Swath for Bonn
# -----------------------------------------------------------------------------

from pcc import cut_the_swath
dpr_lon, dpr_lat, dpr_pp = cut_the_swath(dpr_lon,dpr_lat,dpr_pp, eu=0)


######################################################## Koordinaten Projektion
# -----------------------------------------------------------------------------

proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)


dpr_lon, dpr_lat = wradlib.georef.reproject(dpr_lon, dpr_lat,
                                            projection_target=proj_stereo ,
                                            projection_source=proj_wgs)

blon, blat = wradlib.georef.reproject(blon0, blat0,
                                      projection_target=proj_stereo ,
                                      projection_source=proj_wgs)


############################################################### Dpr zuschneiden
#------------------------------------------------------------------------------

lon0, lat0, radius = blon, blat, 100
rr = np.sqrt((dpr_lat - lat0)**2 + (dpr_lon - lon0)**2)
position = rr < radius

pp = dpr_pp.copy()

pp[np.where(rr > radius)] = np.nan

########################################################### RADOLAN zuschneiden
#------------------------------------------------------------------------------

rr2 = np.sqrt((y - lat0)**2 + (x - lon0)**2)
position2 = rr2 < radius

pp2 = rwdata.copy()

pp2[np.where(rr2 > radius)] = np.nan


################################################# RADOLAN interpolieren auf DPR
#------------------------------------------------------------------------------
gk3 = wradlib.georef.epsg_to_osr(31467)

grid_gpm_xy = np.vstack((dpr_lon.ravel(), dpr_lat.ravel())).transpose()

xy = np.vstack((x.ravel(), y.ravel())).transpose()

result = wrl.ipol.interpolate(xy, grid_gpm_xy,
                              pp2.reshape(pp2.shape[0]*pp2.shape[1]),
                              wrl.ipol.Idw, nnearest=4)

result = np.ma.masked_invalid(result)

rrr = result.reshape(dpr_lon.shape)

rrr[np.where(rr > radius)] = np.nan


################################################## BoXPol interpolieren auf DPR
#------------------------------------------------------------------------------

from wradlib.trafo import idecibel
from wradlib.trafo import decibel
R = idecibel(R)

radar_location = (lon_ppi, lat_ppi, alt_ppi)
elevation = 1.5
azimuths = az
ranges = r
polargrid = np.meshgrid(ranges, azimuths)
lon, lat, alt = wradlib.georef.polar2lonlatalt_n(polargrid[0], polargrid[1],
                                                 elevation, radar_location)
lon, lat = wradlib.georef.reproject(lon, lat, projection_target=proj_stereo ,
                                    projection_source=proj_wgs)

grid_xy = np.vstack((dpr_lon.ravel(), dpr_lat.ravel())).transpose()

xy=np.concatenate([lon.ravel()[:,None],lat.ravel()[:,None]], axis=1)

gridded = wradlib.comp.togrid(xy, grid_xy, ranges[-1], np.array([lon.mean(),
                                                                 lat.mean()]),
                              R.ravel(), ipoli[0],nnearest=40,p=2)

gridded = np.ma.masked_invalid(gridded).reshape(dpr_lon.shape)
gridded[np.where(rr > radius)]=np.nan

R = decibel(R)
gridded = decibel(gridded)



# Hier entsteht beim Interpolieren ein Felhler der
#  sehr oft den Wert ~10.7 erstellt
# Unklar!
new_g = gridded.copy()


gridded[gridded<15]=np.nan
rrr[rrr<15]=np.nan
rwdata[rwdata<15]=np.nan


########################################################################## PLOT
# -----------------------------------------------------------------------------
ff = 15
cc = 0.5
fig = plt.figure(figsize=(12,12))
ax1 = fig.add_subplot(331, aspect='equal')#------------------------------------

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


ax2 = fig.add_subplot(335, aspect='equal')#------------------------------------

plt.pcolormesh(dpr_lon, dpr_lat,np.ma.masked_invalid(rrr),vmin=0, vmax=40, cmap=my_cmap())

plt.colorbar()
plot_borders(ax2)
plot_radar(blon0, blat0, ax2, reproject=True, cband=False,col='black')
plt.plot(dpr_lon[:,0],dpr_lat[:,0], color='black',lw=1)
plt.plot(dpr_lon[:,-1],dpr_lat[:,-1], color='black',lw=1)
plt.plot(dpr_lon[:,dpr_lon.shape[1]/2],dpr_lat[:,dpr_lon.shape[1]/2], color='black',lw=1, ls='--')

plt.xlim(-350,-100)
plt.ylim(-4350, -4100)
plt.title('RADOLAN')
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


ax3 = fig.add_subplot(339, aspect='equal')#------------------------------------

plt.pcolormesh(dpr_lon, dpr_lat,np.ma.masked_invalid(gridded),vmin=0, vmax=40, cmap=my_cmap())

plt.colorbar()
plot_borders(ax3)
plot_radar(blon0, blat0, ax3, reproject=True, cband=False,col='black')
plt.plot(dpr_lon[:,0],dpr_lat[:,0], color='black',lw=1)
plt.plot(dpr_lon[:,-1],dpr_lat[:,-1], color='black',lw=1)
plt.plot(dpr_lon[:,dpr_lon.shape[1]/2],dpr_lat[:,dpr_lon.shape[1]/2], color='black',lw=1, ls='--')

plt.xlim(-350,-100)
plt.ylim(-4350, -4100)
plt.title('BoXPol')
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

#Scatter
from satlib import corcor
ax4 = fig.add_subplot(332, aspect='equal')
plt.scatter(pp,rrr, label=corcor(pp,rrr))
plt.grid();plt.legend()

ax5 = fig.add_subplot(333, aspect='equal')#------------------------------------
plt.scatter(pp,gridded, label=corcor(pp,gridded))
plt.grid();plt.legend()
ax6 = fig.add_subplot(336, aspect='equal')#------------------------------------
plt.scatter(rrr,gridded, label=corcor(gridded,rrr))
plt.grid();plt.legend()


#Hist
bb, aa = 15, 1
ll=2
m1, m2, m3 = ~np.isnan(pp), ~np.isnan(rrr),~np.isnan(gridded)
ax4 = fig.add_subplot(334, aspect='auto')
plt.hist(pp[m1], bins=bb, alpha=aa, label='DPR',facecolor="None",
         edgecolor='blue', linewidth=ll, normed=1)
plt.hist(rrr[m2], bins=bb,alpha=aa, label='RADOLAN',facecolor="None",
         edgecolor='green', linewidth=ll, normed=1)
plt.legend()
plt.grid()
ax5 = fig.add_subplot(337, aspect='auto')#------------------------------------
plt.hist(pp[m1],  bins=bb,alpha=aa, label='DPR',facecolor="None",
         edgecolor='blue', linewidth=ll, normed=1)
plt.hist(gridded[m3],  bins=bb,alpha=aa, label='Boxpol',facecolor="None",
         edgecolor='red', linewidth=ll, normed=1)
plt.legend()
plt.grid()
ax6 = fig.add_subplot(338, aspect='auto')#------------------------------------
plt.hist(rrr[m2],  bins=bb,alpha=aa, label='RADOLAN',facecolor="None",
         edgecolor='green', linewidth=ll, normed=1)

mr = ~np.isnan(rwdata)
plt.hist(rwdata[mr],  bins=bb,alpha=aa, label='RADOLAN',facecolor="None",
         edgecolor='green', linewidth=ll, normed=1,linestyle='dashed')

plt.hist(gridded[m3],  bins=bb,alpha=aa, label='Boxpol',facecolor="None",
         edgecolor='red', linewidth=ll, normed=1)
mR = ~np.isnan(R)
plt.hist(R[mR],  bins=bb,alpha=aa, label='Boxpol',facecolor="None",
         edgecolor='red', linewidth=ll,linestyle='dashed', normed=1)
plt.legend()
plt.grid()
plt.show()

