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


##############################Bei Phase
parameter2 = gpmdprs['NS']['DSD']['phase']
dpr = np.array(parameter2, dtype=float)
dpr[dpr==255]=np.nan
#dpr[dpr<100]=dpr[dpr<100]-100
#dpr[dpr>=200]=dpr[dpr>=200]-200
#dpr[dpr==125]=0
#dpr[dpr==175]=0
#dpr[dpr==100]=0
#dpr[dpr==150]=0

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

dpr4 = np.copy(dpr3)
dpr4[dpr4<100]=dpr4[dpr4<100]-100
dpr4[dpr4>=200]=dpr4[dpr4>=200]-200
dpr4[dpr4==125]=0
dpr4[dpr4==175]=0
dpr4[dpr4==100]=0
dpr4[dpr4==150]=0

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
blat, blon = bonn_pos['lat_ppi'], bonn_pos['lon_ppi']

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
'''
ff = 15
fig = plt.figure(figsize=(10,10))

ax1 = fig.add_subplot(223, aspect='equal')
plt.pcolormesh(x, y, rwdata, cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
#plt.scatter(x, y, rwdata, cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
cb = plt.colorbar(shrink=0.8)
cb.set_label("Rainrate (mm/h)",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plot_borders(ax1)

plot_radar(blon, blat, ax1, reproject=True)

plt.title('RADOLAN Rainrate: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff) #RW Product Polar Stereo
plt.xlabel("x [km] ",fontsize=ff)
plt.ylabel("y [km]  ",fontsize=ff)
#plt.xticks(fontsize=0)
#plt.yticks(fontsize=0)
plt.grid(color='r')
plt.xlim(-420,390)
plt.ylim(-4700, -3700)


ax2 = fig.add_subplot(222, aspect='equal')
pm2 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(gprof_pp_b),
                     cmap=my_cmap,vmin=PV_vmin[ip],vmax=PV_vmax[ip], zorder=2)

#pm2 = plt.pcolormesh(gprof_lon[latstart:latend], gprof_lat[latstart:latend],np.ma.masked_invalid(gprof_pp[latstart:latend]),
                     #cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
cb = plt.colorbar(shrink=0.8)
cb.set_label(PV_name[ip],fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.xlabel("x [km] ",fontsize=ff)
plt.ylabel("y [km]  ",fontsize=ff)
plt.title('GPM DPR Rainrate: \n'+ '2014-10-07 T: 023500 UTC',fontsize=ff)
plot_borders(ax2)
plot_radar(blon, blat, ax2, reproject=True)

#plt.xticks(fontsize=ff)
#plt.yticks(fontsize=ff)
#plt.xlim((bonn_lon1-1,bonn_lon2+1))
#plt.ylim((bonn_lat1-1,bonn_lat2+1))
plt.grid(color='r')
plt.tight_layout()
#Limits Setzen
#ax2.set_xlim(ax1.get_xlim())
#ax2.set_ylim(ax1.get_ylim())
plt.xlim(-420,390)
plt.ylim(-4700, -3700)
#plt.xticks(fontsize=0)
#plt.yticks(fontsize=0)




ax2 = fig.add_subplot(221, aspect='equal')
pm2 = plt.pcolormesh(gpm_x, gpm_y,rrr,
                     cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)

cb = plt.colorbar(shrink=0.8)
cb.set_label("Rainrate (mm/h)",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.xlabel("x [km] ",fontsize=ff)
plt.ylabel("y [km]  ",fontsize=ff)
plt.title('RADOLAN Rainrate Interpoliert: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff) #RW Product Polar Stereo
plot_borders(ax2)
plot_radar(blon, blat, ax2, reproject=True)

plt.xlim(-420,390)
plt.ylim(-4700, -3700)
plt.grid(color='r')
plt.tight_layout()



ax2 = fig.add_subplot(224, aspect='equal')

A = rrr
B = np.ma.masked_invalid(gprof_pp_b)
A[A<TH_rain] = np.nan
#B[B<TH_rain] = np.nan

ref = rrr
est = np.ma.masked_invalid(gprof_pp_b)

mask = ~np.isnan(B) & ~np.isnan(A)
slope, intercept, r_value, p_value, std_err = stats.linregress(B[mask], A[mask])
line = slope*B+intercept


xx = B[mask]
yy = A[mask]
xedges, yedges = np.linspace(-4, 4, 42), np.linspace(-25, 25, 42)
hist, xedges, yedges = np.histogram2d(xx, yy, (xedges, yedges))
xidx = np.clip(np.digitize(xx, xedges), 0, hist.shape[0]-1)
yidx = np.clip(np.digitize(yy, yedges), 0, hist.shape[1]-1)
c = hist[xidx, yidx]
plt.scatter(xx, yy, c=c, label='RR [mm/h]')
cb = plt.colorbar()
cb.set_label("Counts (#)",fontsize=ff)
plt.plot(B,line,'r-')
maxAB = np.nanmax([np.nanmax(xx),np.nanmax(yy)])
plt.xlim(0,maxAB + 1)
plt.ylim(0,maxAB + 1)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=2, fancybox=True, shadow=True,
                    fontsize='small', title= "Slope: " + str(round(slope,3))
                                            + ', Intercept: '+  str(round(intercept,3)) + "\n Correlation: " +
                                            str(round(r_value,3)) + ', Std_err: '+  str(round(std_err,3)))
plt.xlabel("DPR "+str(PV_name[ip]))
plt.ylabel("RADOLAN RR [mm/h]")
plt.title(" .")

plt.grid(True)
plt.tight_layout()
plt.show()
###############################################################################


cut = 15
#dpr3[dpr3 < 0]=np.nan
print ('----------dpr3-------')
print np.nanmin(dpr3[:,cut,:]),np.nanmax(dpr3[:,cut,:])

levels = np.arange((np.nanmin(dpr3[:,cut,:])),(np.nanmax(dpr3[:,cut,:])),0.1)
#levels = np.arange(PV_vmin[ip], PV_vmax[ip],0.1)

#levels = np.arange(0,10,0.01)

print np.nanmax(dpr3[:,cut,:])
print gprof_pp_b.shape, gpm_x.shape


## PLOT
## ----
ff = 15
fig = plt.figure(figsize=(10,10))
ax11 = fig.add_subplot(212, aspect='auto')
h = np.arange(88,0,-1)*0.125 # Bei 88 500m und bei 176 ist es 250m
#plt.contour(gpm_x[:,cut],h,dpr3[:,cut,:].transpose(),vmin=0.1,vmax=10, levels=levels2, colors='black')
plt.contourf(gpm_x[:,cut],h,dpr3[:,cut,:].transpose(),vmin=(np.nanmin(dpr3[:,cut,:])),vmax=(np.nanmax(dpr3[:,cut,:])), levels=levels,cmap=my_cmap)
#plt.pcolormesh(dpr3[:,cut,:].transpose())#,vmin=(np.nanmin(dpr3[:,cut,:])),vmax=(np.nanmax(dpr3[:,cut,:])), levels=levels,cmap=my_cmap)

cb = plt.colorbar(shrink=0.4)
cb.set_label(PV_name[ip],fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.xlabel("x [km] ",fontsize=ff)
plt.ylabel("z [km]  ",fontsize=ff)
plt.grid()
#plt.xlim((-400,0))
#plt.ylim(0,12)

ax12 = fig.add_subplot(211, aspect='equal')

pm12 = plt.pcolormesh(gpm_x[:,:], gpm_y[:,:],np.ma.masked_invalid(dpr3[:,:,80]),
                     cmap=my_cmap,vmin=PV_vmin[ip],vmax=PV_vmax[ip], zorder=2)
cb = plt.colorbar(shrink=0.8)

plt.plot(gpm_x[:,cut],gpm_y[:,cut], color='black',lw=1)
#plt.scatter(bx,by, color='red', marker='v', s=200, zorder=2)

plot_radar(blon, blat, ax12, reproject=True)

#plt.scatter(216,-4236, lw=3, color='magenta', marker='o')# BONN markieren
cb.set_label(PV_name[ip],fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.xlabel("x [km] ",fontsize=ff)
plt.ylabel("y [km]  ",fontsize=ff)
plt.title(str(PV_name[ip])+' GPM DPR: \n'+ '2014-10-07 T: 023500 UTC',fontsize=ff)
plot_borders(ax12)

#plt.xlim((-400,0))
plt.grid()

plt.show()







cut = 35
#dpr3[dpr3 < 0]=np.nan
print ('----------dpr3-------')
print np.nanmin(dpr3[cut,:,:]),np.nanmax(dpr3[cut,:,:])

#levels = np.arange(0,10,0.01)

#levels2 = np.arange((np.nanmin(dpr3[cut,:,:])),(np.nanmax(dpr3[cut,:,:])),0.1)
levels2 = np.arange(PV_vmin[ip],PV_vmax[ip],0.1)
print np.nanmax(dpr3[cut,:,:])
print gprof_pp_b.shape, gpm_x.shape


## PLOT
## ----
ff = 15
fig = plt.figure(figsize=(10,10))
ax21 = fig.add_subplot(212, aspect='auto')
h = np.arange(88,0,-1)*0.125 # Bei 88 250m und bei 176 ist es 125m

#plt.contour(gpm_x[:,cut],h,dpr3[:,cut,:].transpose(),vmin=0.1,vmax=10, levels=levels2, colors='black')
plt.contourf(gpm_x[cut,:],h,dpr3[cut,:,:].transpose(),vmin=PV_vmin[ip],vmax=PV_vmax[ip], levels=levels2,cmap=my_cmap)
#plt.pcolormesh(dpr3[:,cut,:].transpose())#,vmin=(np.nanmin(dpr3[:,cut,:])),vmax=(np.nanmax(dpr3[:,cut,:])), levels=levels,cmap=my_cmap)

cb = plt.colorbar(shrink=0.4)
cb.set_label(PV_name[ip],fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.xlabel("x [km] ",fontsize=ff)
plt.ylabel("z [km]  ",fontsize=ff)
plt.grid()
#plt.xlim((-400,0))
#plt.ylim(0,12)

ax22 = fig.add_subplot(211, aspect='equal')

pm22 = plt.pcolormesh(gpm_x[:,:], gpm_y[:,:],np.ma.masked_invalid(dpr3[:,:,80]),
                     cmap=my_cmap,vmin=PV_vmin[ip],vmax=PV_vmax[ip], zorder=2)
cb = plt.colorbar(shrink=0.8)

plt.plot(gpm_x[cut,:],gpm_y[cut,:], color='black',lw=1)
#plt.scatter(bx,by, color='red', marker='v', s=200, zorder=2)

plot_radar(blon, blat, ax22, reproject=True)

#plt.scatter(216,-4236, lw=3, color='magenta', marker='o')# BONN markieren
cb.set_label(PV_name[ip],fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.xlabel("x [km] ",fontsize=ff)
plt.ylabel("y [km]  ",fontsize=ff)
plt.title(str(PV_name[ip]) +' GPM DPR : \n'+ '2014-10-07 T: 023500 UTC',fontsize=ff)
plot_borders(ax22)

#plt.xlim((-400,0))
plt.grid()
plt.tight_layout()
plt.show()
'''


ff = 20
cut = 15

fig = plt.figure(figsize=(10,10))

ax1 = fig.add_subplot(223, aspect='equal')
plt.pcolormesh(x, y, rwdata, cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
#plt.scatter(x, y, rwdata, cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
cb = plt.colorbar(shrink=0.8)
cb.set_label("Rainrate (mm/h)",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plot_borders(ax1)

plot_radar(blon, blat, ax1, reproject=True)

plt.title('RADOLAN Rainrate: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff) #RW Product Polar Stereo
plt.xlabel("x [km] ",fontsize=ff)
plt.ylabel("y [km]  ",fontsize=ff)
#plt.xticks(fontsize=0)
#plt.yticks(fontsize=0)
plt.grid(color='r')
plt.xlim(-420,390)
plt.ylim(-4700, -3700)


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
plt.title('GPM DPR Rainrate: \n'+ '2014-10-07 T: 023500 UTC',fontsize=ff)
plot_borders(ax2)
plot_radar(blon, blat, ax2, reproject=True)

#plt.xticks(fontsize=ff)
#plt.yticks(fontsize=ff)
#plt.xlim((bonn_lon1-1,bonn_lon2+1))
#plt.ylim((bonn_lat1-1,bonn_lat2+1))
plt.grid(color='r')
plt.tight_layout()
#Limits Setzen
#ax2.set_xlim(ax1.get_xlim())
#ax2.set_ylim(ax1.get_ylim())
plt.xlim(-420,390)
plt.ylim(-4700, -3700)
#plt.xticks(fontsize=0)
#plt.yticks(fontsize=0)




ax2 = fig.add_subplot(221, aspect='equal')
pm2 = plt.pcolormesh(gpm_x, gpm_y,rrr,
                     cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)

cb = plt.colorbar(shrink=0.8)
cb.set_label("Rainrate (mm/h)",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.xlabel("x [km] ",fontsize=ff)
plt.ylabel("y [km]  ",fontsize=ff)
plt.title('RADOLAN Rainrate Interpoliert: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff) #RW Product Polar Stereo
plot_borders(ax2)
plot_radar(blon, blat, ax2, reproject=True)

plt.xlim(-420,390)
plt.ylim(-4700, -3700)
plt.grid(color='r')
plt.tight_layout()



ax2 = fig.add_subplot(224, aspect='auto')
h = np.arange(176,0,-1)*0.125 # Bei 88 500m und bei 176 ist es 250m
level1 = np.arange(np.nanmin(dpr4[:,cut,:]),np.nanmax(dpr4[:,cut,:]),1)
levels= [100,150,200]

print level1


#plt.contour(gpm_x[:,cut],h,dpr3[:,cut,:].transpose(), level = levels2)#,vmin=(np.nanmin(dpr3[:,cut,:])),vmax=(np.nanmax(dpr3[:,cut,:])),cmap=my_cmap)#, levels=levels
plt.contour(gpm_x[:,cut],h,dpr3[:,cut,:].transpose(), levels = levels, colors='black')#,vmin=(np.nanmin(dpr3[:,cut,:])),vmax=(np.nanmax(dpr3[:,cut,:])),cmap=my_cmap)#, levels=levels
plt.contourf(gpm_x[:,cut],h,dpr4[:,cut,:].transpose(), levels = level1, vmax=np.nanmax(dpr4[:,cut,:]))


cb = plt.colorbar(shrink=0.4)
cb.set_label('Temp [C]',fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.xlabel("x [km] ",fontsize=ff)
plt.ylabel("z [km]  ",fontsize=ff)

plt.grid(True)
plt.tight_layout()
plt.show()