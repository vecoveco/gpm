
"""

Einlesen und darstellen von GPM und Radolan Dateien

Radolanpfad:

"""


import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import pandas as pd
import wradlib
import glob
import math
import pandas as pd
from scipy import stats
import matplotlib as mpl
import wradlib as wrl
from osgeo import osr
import os
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
#import mpl_toolkits.basemap.pyproj as pyproj


ipoli = [wradlib.ipol.Idw, wradlib.ipol.Linear, wradlib.ipol.Nearest, wradlib.ipol.OrdinaryKriging]
TH_ka, TH_ku = 0.2, 0.5

#DAS BESP'20141007023500'

# Zeitstempel nach YYYYMMDDhhmmss

ZP = '20141007023500'#'20160805054500'#'20141007023500''20161024232500'#'20150427223500' #'20141007023500'#'20161024232500'#'20140609132500'#'20160917102000'#'20160917102000'#
year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
ye = ZP[2:4]

## Read RADOLAN GK Koordinaten
## ----------------------------

pfad = ('/automount/radar/dwd/rx/'+str(year)+'/'+str(year)+'-'+str(m)+'/'+str(year)+'-'+str(m)+'-'+str(d)+'/raa01-rx_10000-'+str(ye)+str(m)+str(d)+str(ht)+str(mt)+'-dwd---bin.gz')

pfad_radolan = pfad[:-3]

####### pfad

rw_filename = wradlib.util.get_wradlib_data_file(pfad_radolan)
rwdata, rwattrs = wradlib.io.read_RADOLAN_composite(rw_filename)

rwdata = np.ma.masked_equal(rwdata, -9999) / 2 - 32.5

radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
x = radolan_grid_xy[:,:,0]
y = radolan_grid_xy[:,:,1]
#Z = wradlib.trafo.idecibel(rwdata)
#rwdata = wradlib.zr.z2r(Z, a=200., b=1.6)




## Read GPROF
## ------------
pfad2 = ('/home/velibor/shkgpm/data/'+str(year)+str(m)+str(d)+'/dpr/*.HDF5')
pfad_gprof = glob.glob(pfad2)
print pfad_gprof
pfad_gprof_g = pfad_gprof[0]


gpmdprs = h5py.File(pfad_gprof_g, 'r')
gprof_lat=np.array(gpmdprs['NS']['Latitude'])			#(7934, 24)
gprof_lon=np.array(gpmdprs['NS']['Longitude'])

#(7934, 24)
#gprof_pp=np.array(gpmdprs['NS']['SLV']['precipRateNearSurface'])

#gprof_pp=np.array(gpmdprs['NS']['SLV']['precipRate'])
#gprof_pp=np.array(gpmdprs['NS']['DSD']['phase'],dtype=float)
#print type(gprof_pp[1,2,0])
#gprof_pp = np.float(gprof_pp)
#gprof_pp= gprof_pp[:,:,:,0]
#gprof_pp=np.array(gpmdprs['NS']['DSD']['phase'])
##############################################Bei Regen
RR, ZZ = 'precipRate', 'zFactorCorrected' #precipRateNearSurface
gprof_pp=np.array(gpmdprs['NS']['SLV']['precipRateNearSurface'])

print gprof_pp.shape

gprof_pp[gprof_pp==-9999.9]= np.NaN


parameter3 = gpmdprs['NS']['DSD']['binNode']
Node = np.array(parameter3, dtype=float)
Node[Node<-1]= np.nan



##############################Bei Phase
#parameter2 = gpmdprs['NS']['cloudLiqWaterCont']
#parameter2 = gpmdprs['NS']['cloudIceWaterCont']

#parameter2 = gpmdprs['NS']['precipTotPSDparamHigh']
#parameter2 = gpmdprs['NS']['precipTotWaterCont']
#parameter2 = gpmdprs['NS']['correctedReflectFactor']
para_name = 'precipRate'
parameter2 = gpmdprs['NS']['SLV'][para_name]

dpr = np.array(parameter2, dtype=float)
dpr[dpr<-9998]=np.nan
#dpr[dpr<100]=dpr[dpr<100]-100
#dpr[dpr>=200]=dpr[dpr>=200]-200
#dpr[dpr==125]=0
#dpr[dpr==175]=0
#dpr[dpr==100]=0
#dpr[dpr==150]=0

print 'CloudIcemaxmin:', np.nanmin(dpr), np.nanmax(dpr)

######################################## Bei Dropsize
#parameter2 = gpmdprs['NS']['SLV']['paramDSD']
#dpr = np.array(parameter2, dtype=float)
#dpr = dpr[:,:,:,1]
#dpr[dpr<-9998]=np.nan




#####################################Parameter bestimmen
ip = 0
PV_vmin = [0.1,-15]
PV_vmax = [10,40]
PV_name = ['Rainrate (mm/h)','Z (dBZ)']


# Swath ueber Deutschland
from pcc import cut_the_swath
blon, blat, gprof_pp_b = cut_the_swath(gprof_lon,gprof_lat,gprof_pp)
ablon, ablat, dpr3 = cut_the_swath(gprof_lon,gprof_lat,dpr)

nblon, nblat, node = cut_the_swath(gprof_lon,gprof_lat, Node)
node = node[:,:,1:4]

dpr4 = np.copy(dpr3)


print('Shape: ', dpr3.shape)
#dpr3 = gprof_pp_b
#gprof_pp_b = gprof_pp_b[:,:,80]

#gprof_pp_b[gprof_pp_b==-9999.9]=np.nan

print 'gprof min max:' + str(np.nanmin(gprof_pp_b)), str(np.nanmax(gprof_pp_b)), gprof_pp_b.shape
## GPM lon/lat in GK

proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

gpm_x, gpm_y = wradlib.georef.reproject(blon, blat, projection_target=proj_stereo , projection_source=proj_wgs)
grid_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()



## Landgrenzenfunktion
## -------------------
from pcc import boxpol_pos
bonn_pos = boxpol_pos()
bx, by = bonn_pos['gkx_ppi'], bonn_pos['gky_ppi']
boxlat, boxlon = bonn_pos['lat_ppi'], bonn_pos['lon_ppi']

from pcc import plot_borders, plot_ocean

dataset1, inLayer1 = wradlib.io.open_shape('/automount/db01/python/data/ADM/germany/vg250_0101.gk3.shape.ebenen/vg250_ebenen/vg250_l.shp')

import matplotlib.cm as cm
my_cmap = cm.get_cmap('jet',40)
my_cmap.set_under('lightgrey')
my_cmap.set_over('darkred')


##################################################################INTERLOLATION
gk3 = wradlib.georef.epsg_to_osr(31467)

grid_gpm_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose() # GPM Grid erschaffen

xy = np.vstack((x.ravel(), y.ravel())).transpose()

mask = ~np.isnan(rwdata)

result = wrl.ipol.interpolate(xy, grid_gpm_xy, rwdata[mask].reshape(900*900,1), wrl.ipol.Idw, nnearest=4)  #Idw

result = np.ma.masked_invalid(result)

rrr = result.reshape(gpm_x.shape)



Z = wradlib.trafo.idecibel(rwdata)
rwdata = wradlib.zr.z2r(Z, a=200., b=1.6)


Zr = wradlib.trafo.idecibel(rrr)
rrr = wradlib.zr.z2r(Zr, a=200., b=1.6)
#rrr[rrr<=TH_ka]=np.NaN
#rrr[rrr==-9999.0]=np.nan



from pcc import plot_radar

######################################################################### PLOT
###########################################################################----


cut = 20
node[:,cut]
nn = (176-node[:,cut]) * 0.125


fig = plt.figure(figsize=(12,12))
ff = 13.1
fft = 10.0
ax1 = fig.add_subplot(223, aspect='equal')
plt.pcolormesh(x, y, rwdata, cmap=my_cmap,vmin=PV_vmin[ip],vmax=PV_vmax[ip], zorder=2)
#plt.scatter(x, y, rwdata, cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
cb = plt.colorbar(shrink=0.8)
cb.set_label("Rainrate (mm/h)",fontsize=ff)
cb.ax.tick_params(labelsize=ff)

plot_borders(ax1)
plot_radar(boxlon, boxlat, ax1, reproject=True)

plt.title('RADOLAN Rainrate: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff) #RW Product Polar Stereo
plt.xlabel("x [km] ",fontsize=ff)
plt.ylabel("y [km]  ",fontsize=ff)
plt.grid(color='r')
plt.xlim(-420,390)
plt.ylim(-4700, -3700)
plt.xticks(fontsize=fft)
plt.yticks(fontsize=fft)



ax2 = fig.add_subplot(222, aspect='equal')
pm2 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(gprof_pp_b),
                     cmap=my_cmap,vmin=PV_vmin[ip],vmax=PV_vmax[ip], zorder=2)

#pm2 = plt.pcolormesh(gprof_lon[latstart:latend], gprof_lat[latstart:latend],np.ma.masked_invalid(gprof_pp[latstart:latend]),
                     #cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
plt.plot(gpm_x[:,cut],gpm_y[:,cut], color='red',lw=1)
cb = plt.colorbar(shrink=0.8)
cb.set_label(PV_name[ip],fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.xlabel("x [km] ",fontsize=ff)
plt.ylabel("y [km]  ",fontsize=ff)
plt.title('GPM DPR Rainrate: \n'+ '2014-10-07 T: 023600 UTC',fontsize=ff)
plot_borders(ax2)
plot_radar(boxlon, boxlat, ax2, reproject=True)
plt.grid(color='r')
plt.tight_layout()
plt.xlim(-420,390)
plt.ylim(-4700, -3700)
plt.xticks(fontsize=fft)
plt.yticks(fontsize=fft)




ax3 = fig.add_subplot(221, aspect='equal')
pm3 = plt.pcolormesh(gpm_x, gpm_y,rrr,
                     cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)

cb = plt.colorbar(shrink=0.8)
cb.set_label("Rainrate (mm/h)",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.xlabel("x [km] ",fontsize=ff)
plt.ylabel("y [km]  ",fontsize=ff)
plt.title('RADOLAN Rainrate Interpoliert: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff) #RW Product Polar Stereo
plot_borders(ax3)
plot_radar(boxlon, boxlat, ax3, reproject=True)

plt.xlim(-420,390)
plt.ylim(-4700, -3700)
plt.grid(color='r')
plt.tight_layout()
plt.xticks(fontsize=fft)
plt.yticks(fontsize=fft)


ax4 = fig.add_subplot(224, aspect='auto')
h = np.arange(176,0,-1)*0.125 # Bei 88 500m und bei 176 ist es 250m
#level1 = np.arange(np.nanmin(dpr4[:,cut,:]),np.nanmax(dpr4[:,cut,:]),0.1)

level1 = np.arange(0.1,11.,0.1)
t_level = np.arange(1,11,1)

ax5 = plt.contourf(gpm_x[:,cut],h,dpr4[:,cut,:].transpose(),
             vmin=0.101,
             vmax=10,
             cmap=my_cmap,
             levels=level1)
plt.plot(gpm_x[:,cut], nn, '-k')

cb = plt.colorbar(shrink=0.8, ticks=t_level)
cb.set_label('Rainrate (mm/h)',fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.xlabel("x [km] ",fontsize=ff)
plt.ylabel("z [km]  ",fontsize=ff)
plt.xticks(fontsize=fft)
plt.yticks(fontsize=fft)
plt.title('GPM DPR Rainrate: \n'+ '2014-10-07 T: 023600 UTC',fontsize=ff)
plt.xlim(-420,390)

plt.grid(True)
plt.tight_layout()
plt.show()

