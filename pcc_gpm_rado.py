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
TH_rain= 0.2

# Zeitstempel nach YYYYMMDDhhmmss
#ZP = '20141007023500'#'20141007023500'#'20161024232500'#'20150427223500' #'20141007023500'#'20161024232500'#'20140609132500'#'20160917102000'#'20160917102000'#'20160805054500'
#ZP = '20170203005500'
ZP = '20141007023500'
year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
ye = ZP[2:4]

## Read RADOLAN GK Koordinaten
## ----------------------------

pfad = ('/automount/radar/dwd/rx/'+str(year)+'/'+str(year)+'-'+str(m)+'/'+
        str(year)+'-'+str(m)+'-'+str(d)+'/raa01-rx_10000-'+str(ye)+str(m)+
        str(d)+str(ht)+str(mt)+'-dwd---bin.gz')

pfad_radolan = pfad[:-3]

print pfad
print pfad_radolan
####### pfad

rw_filename = wradlib.util.get_wradlib_data_file(pfad_radolan)
rwdata, rwattrs = wradlib.io.read_RADOLAN_composite(rw_filename)

rwdata = np.ma.masked_equal(rwdata, -9999) / 2 - 32.5

#sec = rwattrs['secondary']
#rwdata.flat[sec] = -9999
#rwdata = np.ma.masked_equal(rwdata, -9999)
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
gprof_lon=np.array(gpmdprs['NS']['Longitude'])			#(7934, 24)
#gprof_pp=np.array(gpmdprs['NS']['SLV']['precipRateNearSurface'])

#gprof_pp=np.array(gpmdprs['NS']['SLV']['precipRate'])
#gprof_pp=np.array(gpmdprs['NS']['DSD']['phase'],dtype=float)
gprof_pp=np.array(gpmdprs['NS']['SLV']['zFactorCorrectedNearSurface'])
print gprof_pp.shape
#print type(gprof_pp[1,2,0])
#gprof_pp = np.float(gprof_pp)
gprof_pp[gprof_pp==-9999.9]= 0



from pcc import cut_the_swath2
blon, blat, gprof_pp_b = cut_the_swath2(gprof_lon,gprof_lat,gprof_pp)

dpr3 = gprof_pp_b
#gprof_pp_b = gprof_pp_b[:,:,80]

gprof_pp_b[gprof_pp_b==-9999.9]=np.nan

print 'gprof min max:' + str(np.nanmin(gprof_pp_b)), str(np.nanmax(gprof_pp_b)), gprof_pp_b.shape


proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

#, ,np.ma.masked_invalid(gprof_pp[latstart:latend]
#gpm_x, gpm_y = wrl.georef.reproject(gprof_lon[latstart:latend], gprof_lat[latstart:latend], projection_source=proj_ll,projection_target=proj_gk)
gpm_x, gpm_y = wradlib.georef.reproject(blon, blat, projection_target=proj_stereo , projection_source=proj_wgs)
grid_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()





## Landgrenzenfunktion
## -------------------
from pcc import boxpol_pos
bonn_pos = boxpol_pos()
bx, by = bonn_pos['gkx_ppi'], bonn_pos['gky_ppi']
bonnlat, bonnlon = bonn_pos['lat_ppi'], bonn_pos['lon_ppi']
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



#Z = wradlib.trafo.idecibel(rwdata)
#rwdata = wradlib.zr.z2r(Z, a=200., b=1.6)


#Zr = wradlib.trafo.idecibel(rrr)
#rrr = wradlib.zr.z2r(Zr, a=200., b=1.6)
print np.nanmin(rrr)
rrr[rrr<=0]=0




from pcc import plot_radar

########################################################################## PLOT
###########################################################################----

ff = 15
fig = plt.figure(figsize=(10,10))

ax1 = fig.add_subplot(223, aspect='equal')
plt.pcolormesh(x, y, rwdata, cmap=my_cmap, vmin=0.01, vmax=50, zorder=2)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
#plt.scatter(x, y, rwdata, cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
cb = plt.colorbar(shrink=0.8)
cb.set_label("Ref [dbz]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plot_borders(ax1)

plot_radar(bonnlon, bonnlat, ax1, reproject=True)

plt.title('RADOLAN Rainrate: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff) #RW Product Polar Stereo
#plt.xlabel("x [km] ",fontsize=ff)
#plt.ylabel("y [km]  ",fontsize=ff)
#plt.xticks(fontsize=0)
#plt.yticks(fontsize=0)
plt.grid(color='r')
#plt.xlim(-1000, 850)
#plt.ylim(-5500, -3000)
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')

ax2 = fig.add_subplot(222, aspect='equal')
pm2 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(gprof_pp_b),
                     cmap=my_cmap, vmin=0.01, vmax=50, zorder=2)

cb = plt.colorbar(shrink=0.8)
cb.set_label("Ref [dbz]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
#plt.xlabel("x [km] ",fontsize=ff)
#plt.ylabel("y [km]  ",fontsize=ff)
plt.title('GPM DPR Ref: \n'+ '2014-10-07 T: 023500 UTC',fontsize=ff)
plot_borders(ax2)
plot_radar(bonnlon, bonnlat, ax2, reproject=True)

plt.grid(color='r')
plt.tight_layout()

#plt.xlim(-1000, 850)
#plt.ylim(-5500, -3000)
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')


ax2 = fig.add_subplot(221, aspect='equal')
pm3 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(rrr),
                     cmap=my_cmap, vmin=0.01, vmax=50,zorder=2)

cb = plt.colorbar(shrink=0.8)
cb.set_label("Ref [dbz]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
#plt.xlabel("x [km] ",fontsize=ff)
#plt.ylabel("y [km]  ",fontsize=ff)
plt.title('RADOLAN Ref Interpoliert: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff) #RW Product Polar Stereo
plot_borders(ax2)
plot_radar(bonnlon, bonnlat, ax2, reproject=True)

#plt.xlim(-420,390)
#plt.ylim(-4700, -3700)
plt.grid(color='r')
plt.tight_layout()
#plt.xlim(-1000, 850)
#plt.ylim(-5500, -3000)
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')

ax2 = fig.add_subplot(224, aspect='equal')


A = np.copy(rrr)
B = np.ma.masked_invalid(gprof_pp_b)
A[A<0.001] = np.nan
B[B<0.001] = np.nan

ref = np.copy(rrr)
est = np.ma.masked_invalid(gprof_pp_b)

mask = ~np.isnan(B) & ~np.isnan(A)
slope, intercept, r_value, p_value, std_err = stats.linregress(B[mask], A[mask])
line = slope*B+intercept

from pcc import skill_score
RR = skill_score(est,ref,0.)
from pcc import plot_score
plot_score(est,ref,RR)

#plt.scatter(B[mask],A[mask], label='RR [mm/h]')
plt.plot(B,line,'r-')
#maxAB = np.nanmax([np.nanmax(A[mask]),np.nanmax(B[mask])])
#plt.xlim(0,maxAB + 1)
#plt.ylim(0,maxAB + 1)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=2, fancybox=True, shadow=True,
                    fontsize='small', title= "Slope: " + str(round(slope,3))
                                            + ', Intercept: '+  str(round(intercept,3)) + "\n Correlation: " +
                                            str(round(r_value,3)) + ', Std_err: '+  str(round(std_err,3)))
plt.xlabel("DPR Ref [dbz]")
plt.ylabel("RADOLAN Ref [dbz]")
plt.title(" .")

plt.grid(True)
plt.tight_layout()
plt.show()



