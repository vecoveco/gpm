
'''

Darstellung der Overpasses ueber BoXPol

Plots fur die MIUB Website

'''

#!/usr/bin/env python

import h5py
import numpy as np
import matplotlib.pyplot as plt
import wradlib
import glob

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

pfad_radar2 = glob.glob('/automount/ags/velibor/gpmdata/dpr/2A.GPM.DPR.V6-20160118.' + year + m + d + '*.HDF5')

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


# PPI BoxPol Daten einlesen
#---------------------------

ppi=h5py.File(ppi_datapath,'r')
data, attrs = wradlib.io.read_GAMIC_hdf5(ppi_datapath)

ZH0 = data['SCAN0']['ZH']['data']
zdr = data['SCAN0']['ZDR']['data']

PHIDP = data['SCAN0']['PHIDP']['data']
r = attrs['SCAN0']['r']
az = attrs['SCAN0']['az']
lon_ppi = attrs['VOL']['Longitude']
lat_ppi = attrs['VOL']['Latitude']
alt_ppi = attrs['VOL']['Height']
rho = data['SCAN0']['RHOHV']['data']

boxpoltime = attrs['SCAN0']['Time']


R = ZH0
#R[151:165,:]=np.nan
#zdr[151:165,:]=np.nan
#rho[151:165,:]=np.nan


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
R[R<15]=np.nan

### Threshold for DPR sensitivity
R[R<TH]=np.nan

# DPR Einlesen
# ------------
gpmku = h5py.File(pfad_radar, 'r')
gpmku_NS = gpmku['NS']['SLV']
gpmka_MS = gpmku['MS']['SLV']
gpmka_HS = gpmku['HS']['SLV']


dpr_lat = np.array(gpmku['NS']['Latitude'])
dpr_lon = np.array(gpmku['NS']['Longitude'])
dpr_latms = np.array(gpmku['MS']['Latitude'])
dpr_lonms = np.array(gpmku['MS']['Longitude'])
dpr_laths = np.array(gpmku['HS']['Latitude'])
dpr_lonhs = np.array(gpmku['HS']['Longitude'])
dpr_pp = np.array(gpmku_NS['zFactorCorrectedNearSurface'])
dpr_pp_ms = np.array(gpmka_MS['zFactorCorrectedNearSurface'])
dpr_pp_hs = np.array(gpmka_HS['zFactorCorrectedNearSurface'])

dpr_pp[dpr_pp < 0] = np.nan
dpr_pp_ms[dpr_pp_ms < 0] = np.nan
dpr_pp_hs[dpr_pp_hs < 0] = np.nan

from satlib import get_time_of_gpm_over_boxpol
gpm_time = gpmku['NS']['ScanTime']
gt = get_time_of_gpm_over_boxpol(dpr_lon, dpr_lat, gpm_time)

gz = gt[0:4]+'-'+gt[4:6]+'-'+gt[6:8]+'T'+gt[8:10]+':'+gt[10:12]+':'+gt[12:14]+'Z UTC'

# Cut the Swath
from pcc import cut_the_swath
dpr_lon, dpr_lat, dpr_pp = cut_the_swath(dpr_lon,dpr_lat,dpr_pp, eu=0)
dpr_lonms, dpr_latms, dpr_pp_ms = cut_the_swath(dpr_lonms,dpr_latms,dpr_pp_ms, eu=0)
dpr_lonhs, dpr_laths, dpr_pp_hs = cut_the_swath(dpr_lonhs,dpr_laths,dpr_pp_hs, eu=0)

# Koordinaten Projektion
# ------------------

proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)



dpr_lon, dpr_lat = wradlib.georef.reproject(dpr_lon, dpr_lat, projection_target=proj_stereo , projection_source=proj_wgs)
dpr_lonms, dpr_latms = wradlib.georef.reproject(dpr_lonms, dpr_latms, projection_target=proj_stereo , projection_source=proj_wgs)
dpr_lonhs, dpr_laths = wradlib.georef.reproject(dpr_lonhs, dpr_laths, projection_target=proj_stereo , projection_source=proj_wgs)


blon, blat = wradlib.georef.reproject(blon0, blat0, projection_target=proj_stereo , projection_source=proj_wgs)

radar_location = (lon_ppi, lat_ppi, alt_ppi)
elevation = 1.5
azimuths = az
ranges = r
polargrid = np.meshgrid(ranges, azimuths)
lon, lat, alt = wradlib.georef.polar2lonlatalt_n(polargrid[0], polargrid[1], elevation, radar_location)
lon, lat = wradlib.georef.reproject(lon, lat, projection_target=proj_stereo , projection_source=proj_wgs)

# Dpr zuschneiden
#-----------------

lon0, lat0, radius = blon, blat, 100
rr = np.sqrt((dpr_lat - lat0)**2 + (dpr_lon - lon0)**2)
position = rr < radius

rrms = np.sqrt((dpr_latms - lat0)**2 + (dpr_lonms - lon0)**2)
positionms = rrms < radius

rrhs = np.sqrt((dpr_laths - lat0)**2 + (dpr_lonhs - lon0)**2)
positionhs = rrhs < radius


pp = dpr_pp.copy()
ppms = dpr_pp_ms.copy()
pphs = dpr_pp_hs.copy()

#pp[np.where(rr > radius)] = np.nan
#ppms[np.where(rrms > radius)] = np.nan
#pphs[np.where(rrhs > radius)] = np.nan



# Size Colorbar
cc = 0.6
# Reflectivity Ranges
zmin, zmax = -30, 65
step = 5
bounds = np.arange(zmin-5,zmax+5,step)
bounds2 = np.arange(-1.2, 3.4, .2)

# Fontsize Title
ff = 13
fig = plt.figure(figsize=(16,10))
#fig.suptitle('GPM Overpass BoXPol\n Time: '+ZP, fontsize=15)

###################
ax1 = fig.add_subplot(231, aspect='equal')
plt.pcolormesh(dpr_lon, dpr_lat,np.ma.masked_invalid(pp),vmin=zmin, vmax=zmax, cmap=my_cmap())


cb1 = plt.colorbar(shrink=cc,extend='both',ticks=bounds, boundaries=bounds)
cb1.set_label('Reflectivity (dBZ)')

plot_borders(ax1)
plot_radar(blon0, blat0, ax1, reproject=True, cband=False,col='black')
plt.plot(dpr_lon[:,0],dpr_lat[:,0], color='black',lw=1)
plt.plot(dpr_lon[:,-1],dpr_lat[:,-1], color='black',lw=1)
plt.plot(dpr_lon[:,dpr_lon.shape[1]/2],dpr_lat[:,dpr_lon.shape[1]/2], color='black',lw=1, ls='--')

plt.xlim(-350,-100)
plt.ylim(-4350, -4100)
plt.title(gz + ' DPR NS V06', loc='left', fontsize=ff)
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

ax2 = fig.add_subplot(232, aspect='equal')
plt.pcolormesh(dpr_lonms, dpr_latms,np.ma.masked_invalid(ppms),vmin=zmin, vmax=zmax, cmap=my_cmap())

cb2 = plt.colorbar(shrink=cc,extend='both',ticks=bounds, boundaries=bounds)
cb2.set_label('Reflectivity (dBZ)')

plot_borders(ax2)
plot_radar(blon0, blat0, ax2, reproject=True, cband=False,col='black')
plt.plot(dpr_lonms[:,0],dpr_latms[:,0], color='black',lw=1)
plt.plot(dpr_lonms[:,-1],dpr_latms[:,-1], color='black',lw=1)
plt.plot(dpr_lon[:,dpr_lon.shape[1]/2],dpr_lat[:,dpr_lon.shape[1]/2], color='black',lw=1, ls='--')

plt.xlim(-350,-100)
plt.ylim(-4350, -4100)
plt.title(gz + ' DPR MS V06', loc='left', fontsize=ff)
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

ax3 = fig.add_subplot(233, aspect='equal')
plt.pcolormesh(dpr_lonhs, dpr_laths,np.ma.masked_invalid(pphs),vmin=zmin, vmax=zmax, cmap=my_cmap())

cb3 =plt.colorbar(shrink=cc,extend='both',ticks=bounds, boundaries=bounds)
cb3.set_label('Reflectivity (dBZ)')

plot_borders(ax3)
plot_radar(blon0, blat0, ax3, reproject=True, cband=False,col='black')
plt.plot(dpr_lonhs[:,0],dpr_laths[:,0], color='black',lw=1)
plt.plot(dpr_lonhs[:,-1],dpr_laths[:,-1], color='black',lw=1)
plt.plot(dpr_lon[:,dpr_lon.shape[1]/2],dpr_lat[:,dpr_lon.shape[1]/2], color='black',lw=1, ls='--')

plt.xlim(-350,-100)
plt.ylim(-4350, -4100)
plt.title(gz + ' DPR HS V06', loc='left', fontsize=ff)
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


ax4 = fig.add_subplot(234, aspect='equal')
plt.pcolormesh(lon, lat,R,vmin=zmin, vmax=zmax, cmap=my_cmap())
cb4 = plt.colorbar(shrink=cc,extend='both',ticks=bounds, boundaries=bounds)
cb4.set_label('Reflectivity (dBZ)')
plot_borders(ax4)
plot_radar(blon0, blat0, ax4, reproject=True, cband=False,col='black')
plt.plot(dpr_lon[:,0],dpr_lat[:,0], color='black',lw=1)
plt.plot(dpr_lon[:,-1],dpr_lat[:,-1], color='black',lw=1)
plt.plot(dpr_lonhs[:,0],dpr_laths[:,0], color='black',lw=1)
plt.plot(dpr_lonhs[:,-1],dpr_laths[:,-1], color='black',lw=1)
plt.plot(dpr_lon[:,dpr_lon.shape[1]/2],dpr_lat[:,dpr_lon.shape[1]/2], color='black',lw=1, ls='--')

plt.xlim(-350,-100)
plt.ylim(-4350, -4100)
plt.title(boxpoltime + ' UTC - BoXPol - ZH', loc='left', fontsize=ff)
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

ax5 = fig.add_subplot(235, aspect='equal')
plt.pcolormesh(lon, lat,zdr, cmap=my_cmap(), vmin=-1.2, vmax=3.4)
cb5 = plt.colorbar(shrink=cc,extend='both',ticks=bounds2, boundaries=bounds2)
cb5.set_label('Differential Reflectivity (dB)')
plot_borders(ax5)
plot_radar(blon0, blat0, ax5, reproject=True, cband=False,col='black')
plt.plot(dpr_lon[:,0],dpr_lat[:,0], color='black',lw=1)
plt.plot(dpr_lon[:,-1],dpr_lat[:,-1], color='black',lw=1)
plt.plot(dpr_lonhs[:,0],dpr_laths[:,0], color='black',lw=1)
plt.plot(dpr_lonhs[:,-1],dpr_laths[:,-1], color='black',lw=1)
plt.plot(dpr_lon[:,dpr_lon.shape[1]/2],dpr_lat[:,dpr_lon.shape[1]/2], color='black',lw=1, ls='--')

plt.xlim(-350,-100)
plt.ylim(-4350, -4100)
plt.title(boxpoltime + ' UTC - BoXPol - ZDR', loc='left', fontsize=ff)
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

bounds3 = np.arange(0.6, 1, 0.02)
#bounds3 = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 0.92, 0.94, 0.96, 0.97, 0.975, 0.980, 0.985, 0.90, 0.95, 0.96, 0.97, 0.98, 0.99, 1])

ax6 = fig.add_subplot(236, aspect='equal')
plt.pcolormesh(lon, lat,rho, cmap=my_cmap(), vmin=0.6, vmax=.99)
cb6 = plt.colorbar(shrink=cc,extend='both',ticks=bounds3, boundaries=bounds3)
cb6.set_label('Crosscorrelation Coefficient ()')
plot_borders(ax6)
plot_radar(blon0, blat0, ax6, reproject=True, cband=False,col='black')
plt.plot(dpr_lon[:,0],dpr_lat[:,0], color='black',lw=1)
plt.plot(dpr_lon[:,-1],dpr_lat[:,-1], color='black',lw=1)
plt.plot(dpr_lonhs[:,0],dpr_laths[:,0], color='black',lw=1)
plt.plot(dpr_lonhs[:,-1],dpr_laths[:,-1], color='black',lw=1)
plt.plot(dpr_lon[:,dpr_lon.shape[1]/2],dpr_lat[:,dpr_lon.shape[1]/2], color='black',lw=1, ls='--')

plt.xlim(-350,-100)
plt.ylim(-4350, -4100)
plt.title(boxpoltime + ' UTC - BoXPol - RHOHV', loc='left', fontsize=ff)
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

plt.tight_layout()
#plt.savefig('/automount/ags/velibor/plot/boxpol/boxpol_vs_DPR/boxdpr_'+ZP)
plt.show()


