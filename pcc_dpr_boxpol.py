
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

Pos = boxpol_pos()
blon, blat = Pos['lon_ppi'], Pos['lat_ppi']
bbx, bby = Pos['gkx_ppi'], Pos['gky_ppi']

# Pfad mit String
# ---------------

# Hohe von DPR
TH = 18 #Threshold um Nullen fuer Niederschlag raus zu filtern

ipoli = [wradlib.ipol.Idw, wradlib.ipol.Linear, wradlib.ipol.Nearest, wradlib.ipol.OrdinaryKriging]
offset = 2


ZP = '20141007023744'
#ZP = '20140704134500'
#ZP = '20150225163500'
#ZP = '20150816070500'
#good#ZP = '20151208213500'
#good#ZP = '20151216024501'
#ZP = '20160503024500'

year = ZP[0:4]
m = ZP[4:6]
d = ZP[6:8]
ht = ZP[8:10]
mt = ZP[10:12]
st = ZP[12:14]



pfad_radar = glob.glob('/automount/ags/velibor/gpmdata/dpr/2A.GPM.DPR.V6-20160118.' + year + m + d + '*.HDF5')
print pfad_radar
#pfad_radar = sorted(glob.glob(pfad))
#print pfad_radar
pfad_radar_Ku = pfad_radar[0]

try:
    ppi_datapath=glob.glob('/automount/radar-archiv/scans/' + year+ "/" + year +"-"+ m + "/" + year+ "-" + m +"-"+ d + "/ppi_1p5deg/"+ year + "-" + m +"-"+ d + "--" +ht +":"+mt+":"+st+",*.mvol")[0]
except:
    ppi_datapath=glob.glob('/automount/radar/scans/' + year+ "/" + year +"-"+ m + "/" + year+ "-" + m +"-"+ d + "/ppi_1p5deg/"+ year + "-" + m +"-"+ d + "--" +ht +":"+mt+":"+st+",*.mvol")[0]

print ppi_datapath


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

R = ZH

#R[151:165]=np.nan



# DPR Einlesen
# ------------

gpmku = h5py.File(pfad_radar_Ku, 'r')
gpmku_HS=gpmku['NS']['SLV']
ku_lat=np.array(gpmku['NS']['Latitude'])			#(7934, 24)
ku_lon=np.array(gpmku['NS']['Longitude'])			#(7934, 24)
ku_pp=np.array(gpmku_HS['zFactorCorrectedNearSurface'])




# Lon Lat Bestimmung
# ------------------
radars = [ku_pp]
rad_lat = [ku_lat]
rad_lon = [ku_lon]
rad_name = ['DPR']

ii = 0

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
gridded = wradlib.comp.togrid(xy, grid_xy, ranges[-1], np.array([x.mean(), y.mean()]), R.ravel(), ipoli[0],nnearest=40,p=2)
gridded = np.ma.masked_invalid(gridded).reshape(xgrid.shape)




# Plot
# ----



################################################################Swap!
#rrr, ggg = ggg, rrr

ff = 15
cc = 0.5
fig = plt.figure(figsize=(12,12))
ax1 = fig.add_subplot(221, aspect='equal')#------------------------------------
ax1, pm1 = wradlib.vis.plot_ppi(R,r,az,vmin=0.01,vmax=50, cmap=my_cmap())

cb = plt.colorbar(pm1,shrink=cc)
cb.set_label("Reflectivity [dBZ]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.plot(xgrid[latstart:latend][:,0], ygrid[latstart:latend][:,0], color='black')
plt.plot(xgrid[latstart:latend][:,-1], ygrid[latstart:latend][:,-1], color='black')
#plot_borders(ax1)
plot_radar(bbx,bby, ax1, reproject=False, cband=False,col='black')
#plot_radar(blon, blat, ax1, reproject=True, cband=False,col='black')

plt.title('BoXPol Reflectivity:\n 2014-10-07--02:37:44',fontsize=ff)
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



ax2 = fig.add_subplot(222, aspect='equal')#------------------------------------

pm2 = plt.pcolormesh(dpr_lon[latstart:latend], dpr_lat[latstart:latend],np.ma.masked_invalid(dpr_pp[latstart:latend]),
                     cmap=my_cmap(), vmin=0.01, vmax=50, zorder=2)
plt.plot(dpr_lon[latstart:latend][:,0], dpr_lat[latstart:latend][:,0], color='black')
plt.plot(dpr_lon[latstart:latend][:,-1], dpr_lat[latstart:latend][:,-1], color='black')
plot_radar(blon, blat, ax2, reproject=False, cband=False,col='black')


cb = plt.colorbar(pm2,shrink=cc)
cb.set_label("Reflectivity [dBZ]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.title('GPM DPR Reflectivity \n 2014-10-07--02:35:27 ',fontsize=ff)
plt.grid(color='r')
plt.tick_params(
        axis='both',
        which='both',
        bottom='off',
        top='off',
        labelbottom='off',
        right='off',
        left='off',
        labelleft='off')
plt.xlim(5.5,9)
plt.ylim(49,52)


ax2 = fig.add_subplot(223, aspect='equal')#------------------------------------

pm3 = plt.pcolormesh(dpr_lon[latstart:latend], dpr_lat[latstart:latend], gridded,
                     cmap=my_cmap(), vmin=0.01, vmax=50,zorder=2)
plt.plot(dpr_lon[latstart:latend][:,0], dpr_lat[latstart:latend][:,0], color='black')
plt.plot(dpr_lon[latstart:latend][:,-1], dpr_lat[latstart:latend][:,-1], color='black')

cb = plt.colorbar(pm3, shrink=cc)
cb.set_label("Reflectivity [dBZ]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)

plt.title('BoXPol Reflectivity Interpoliert\n 2014-10-07--02:37:44',fontsize=ff) #RW Product Polar Stereo
plt.tick_params(
        axis='both',
        which='both',
        bottom='off',
        top='off',
        labelbottom='off',
        right='off',
        left='off',
        labelleft='off')
#lot_borders(ax2)
plot_radar(blon, blat, ax2, reproject=False, cband=False,col='black')
plt.grid(color='r')
plt.xlim(5.5,9)
plt.ylim(49,52)

ax4 = fig.add_subplot(224, aspect='equal')#------------------------------------

rrr = gridded.copy()
ggg = np.ma.masked_invalid(dpr_pp)[latstart:latend].copy()

rrr[rrr<TH]=np.nan
ggg[ggg<TH]=np.nan

maske = ~np.isnan(ggg) & ~np.isnan(rrr)
slope, intercept, r_value, p_value, std_err = stats.linregress(ggg[maske], rrr[maske])
line = slope * ggg +intercept

diffi = ggg[maske]-rrr[maske]
bias = np.nansum(diffi)/len(diffi)
rmse = np.sqrt(np.nansum(((diffi)**2.0)/len(diffi)))

ax4.scatter(ggg, rrr, label='Reflectivity [dBZ]', color='grey', alpha=0.6)

r_value_s, p_value_s = stats.spearmanr(ggg[maske],rrr[maske])

text = ('f(x) = ' + str(round(slope,3)) + 'x + ' + str(round(intercept,3)) +
           '\nCorr: ' + str(round(r_value,3)) + r'$\pm$: '+  str(round(std_err,3))+
        '\nbias:' +  str(round(bias,3))+
        '\nrmse:' +  str(round(rmse,3))+
        '\nCorrS:' +  str(round(r_value_s,3)))

ax4.annotate(text, xy=(0.01, 0.99), xycoords='axes fraction', fontsize=10,
                horizontalalignment='left', verticalalignment='top')
from scipy import stats, linspace

t1 = linspace(0,50,50)
plt.plot(t1,t1,'k-')
plt.plot(t1, t1*slope + intercept, 'r-', lw=3 ,label='Regression')

plt.legend(loc='lower right', fontsize=10, scatterpoints= 1, numpoints=1, shadow=True)

plt.xlim(0,50)
plt.ylim(0,50)


plt.xlabel('GPM DPR Reflectivity [dBZ]',fontsize=ff)
plt.ylabel('BoXPol Reflectivity [dBZ]',fontsize=ff)
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)
plt.grid(color='r')



plt.tight_layout()
#plt.savefig('/home/velibor/shkgpm/plot/gpm_dpr_radolan_v2_'+ZP + '.png' )
#plt.close()
plt.show()




