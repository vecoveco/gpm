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

ZP = '20150110220500'

# Zeitstempel nach YYYYMMDDhhmmss
year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
ye = ZP[2:4]



## Read RADOLAN GK Koordinaten
## ----------------------------
iii = 0
r_pro = 'rx'
pfad = ('/automount/radar/dwd/'+r_pro+'/'+str(year)+'/'+str(year)+'-'+str(m)+
        '/'+str(year)+'-'+str(m)+'-'+str(d)+'/raa01-'+r_pro+'_10000-'+str(ye)+
        str(m)+str(d)+str(ht)+str(mt)+'-dwd---bin.gz')

pfad_radolan = pfad[:-3]

####### pfad

try:
    rw_filename = wradlib.util.get_wradlib_data_file(pfad)
except EnvironmentError:
    rw_filename = wradlib.util.get_wradlib_data_file(pfad_radolan)

rwdata, rwattrs = wradlib.io.read_RADOLAN_composite(rw_filename)

rwdata = np.ma.masked_equal(rwdata, -9999) / 2 - 32.5

radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
x = radolan_grid_xy[:,:,0]
y = radolan_grid_xy[:,:,1]





## Read  t
## ------------
pfadT = ('/home/velibor/shkgpm/data/20150110/T/*.HDF5')
pfad_T = glob.glob(pfadT)
print pfad_T
pfad_T_g = pfad_T[0]
T = h5py.File(pfad_T_g, 'r')
T_lat = np.array(T['S1']['Latitude'])
T_lon = np.array(T['S1']['Longitude'])
T_pp = np.array(T['S1']['Tc'])
T_pp[T_pp== -9999.9] = np.nan

T2_lat = np.array(T['S2']['Latitude'])
T2_lon = np.array(T['S2']['Longitude'])
T2_pp = np.array(T['S2']['Tc'])
T2_pp[T2_pp== -9999.9] = np.nan

pfad3 = ('/home/velibor/shkgpm/data/20150110/gprof/*.HDF5')
pfad_gprof = glob.glob(pfad3)
print pfad_gprof
pfad_gprof_g = pfad_gprof[0]
gpmgmi = h5py.File(pfad_gprof_g, 'r')
gpmgmi_S1=gpmgmi['S1']
gprof_lat=np.array(gpmgmi_S1['Latitude'])
gprof_lon=np.array(gpmgmi_S1['Longitude'])

gprof_pp=np.array(gpmgmi_S1['surfacePrecipitation'])
gprof_pp[gprof_pp<=0] = np.nan

gprof_snow=np.array(gpmgmi_S1['snowCoverIndex'], dtype=float)
gprof_snow[gprof_snow==-99] = np.nan


gprof_l=np.array(gpmgmi_S1['liquidPrecipFraction'], dtype=float)
gprof_l[gprof_l==-9999.9] = np.nan

pfad4 = ('/home/velibor/shkgpm/data/20150110/dpr/*.HDF5')
pfad_gprof4 = glob.glob(pfad4)
print pfad_gprof4
pfad_gprof_g = pfad_gprof4[0]
gpmgmi = h5py.File(pfad_gprof_g, 'r')
dpr_lat=np.array(gpmgmi['NS']['Latitude'])
dpr_lon=np.array(gpmgmi['NS']['Longitude'])
dpr_pp=np.array(gpmgmi['NS']['SLV']['precipRateNearSurface'])
dpr_pp[dpr_pp<=0] = np.nan


bonn_lat1 = 47.9400
bonn_lat2 = 55.3500
bonn_lon1 = 6.40000
bonn_lon2 = 14.10000

ilat= np.where((gprof_lat>bonn_lat1) & (gprof_lat<bonn_lat2))
ilon= np.where((gprof_lon>bonn_lon1) & (gprof_lon<bonn_lon2))
lonstart = ilon[0][0]
lonend = ilon[0][-1]
latstart = ilat[0][0]
latend = ilat[0][-1]


ilatd= np.where((dpr_lat>bonn_lat1) & (dpr_lat<bonn_lat2))
ilond= np.where((dpr_lon>bonn_lon1) & (dpr_lon<bonn_lon2))
lonstartd = ilond[0][0]
lonendd = ilond[0][-1]
latstartd = ilatd[0][0]
latendd = ilatd[0][-1]

proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)


gprof_x, gprof_y = wradlib.georef.reproject(gprof_lon[latstart:latend], gprof_lat[latstart:latend], projection_target=proj_stereo , projection_source=proj_wgs)
dpr_x, dpr_y = wradlib.georef.reproject(dpr_lon[latstartd:latendd], dpr_lat[latstartd:latendd], projection_target=proj_stereo , projection_source=proj_wgs)


itlat= np.where((T_lat>bonn_lat1) & (T_lat<bonn_lat2))
itlon= np.where((T_lon>bonn_lon1) & (T_lon<bonn_lon2))
tlonstart = itlon[0][0]
tlonend = itlon[0][-1]
tlatstart = itlat[0][0]
tlatend = itlat[0][-1]

T_x, T_y = wradlib.georef.reproject(T_lon[tlatstart:tlatend], T_lat[tlatstart:tlatend], projection_target=proj_stereo , projection_source=proj_wgs)


itlat2= np.where((T2_lat>bonn_lat1) & (T2_lat<bonn_lat2))
itlon2= np.where((T2_lon>bonn_lon1) & (T2_lon<bonn_lon2))
tlonstart2 = itlon2[0][0]
tlonend2 = itlon2[0][-1]
tlatstart2 = itlat2[0][0]
tlatend2 = itlat2[0][-1]

T2_x, T2_y = wradlib.georef.reproject(T2_lon[tlatstart:tlatend], T2_lat[tlatstart:tlatend], projection_target=proj_stereo , projection_source=proj_wgs)


## Landgrenzenfunktion
## -------------------
from pcc import plot_borders



dataset1, inLayer1 = wradlib.io.open_shape('/automount/db01/python/data/ADM/germany/vg250_0101.gk3.shape.ebenen/vg250_ebenen/vg250_l.shp')

import matplotlib.cm as cm
my_cmap = cm.get_cmap('jet',40)
my_cmap.set_under('lightgrey')
my_cmap.set_over('darkred')




#Z = wradlib.trafo.idecibel(rwdata)
#rwdata = wradlib.zr.z2r(Z, a=200., b=1.6)

#rwdata = rwdata*8

######################################################################################## PLOT
from pcc import get_miub_cmap
#ff = 15
#fig = plt.figure(figsize=(12,12))
def plot_all():
    ax1 = fig.add_subplot(231, aspect='equal')
    plt.pcolormesh(x, y, rwdata, cmap=get_miub_cmap(),
                   vmin=0.1,vmax=50, zorder=2)
    cb = plt.colorbar(shrink=0.5)
    cb.set_label("Ref (dbz)",fontsize=ff)
    #cb.set_label("Rainrate (mm/h)",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    plot_borders(ax1)
    plt.title('RADOLAN Ref: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
           ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff) #RW Product Polar Stereo

    plt.xlabel("x [km] ",fontsize=ff)
    plt.ylabel("y [km]  ",fontsize=ff)
    plt.grid(color='r')
    plt.xlim(-420,390)
    plt.ylim(-4700, -3700)


    ax2 = fig.add_subplot(232, aspect='equal')
    pm2 = plt.pcolormesh(gprof_x, gprof_y,np.ma.masked_invalid(gprof_pp[latstart:latend]),
                         cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)

    cb = plt.colorbar(shrink=0.5)
    cb.set_label("Rainrate (mm/h)",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    plt.xlabel("x [km] ",fontsize=ff)
    plt.ylabel("y [km]  ",fontsize=ff)
    plt.title('GPM GPROF Rainrate: \n' + str(pfad_gprof_g[66:70]) + '-' +str(pfad_gprof_g[70:72])+ '-' +
              str(pfad_gprof_g[72:74]) + ' T: ' +str(pfad_gprof_g[76:82]) + '-' + str(pfad_gprof_g[84:90]) + ' UTC',fontsize=ff)
    plot_borders(ax2)
    plt.xticks(fontsize=ff)
    plt.yticks(fontsize=ff)
    plt.grid(color='r')
    plt.tight_layout()
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_ylim(ax1.get_ylim())



    ax37 = fig.add_subplot(233, aspect='equal')
    plt.pcolormesh(dpr_x, dpr_y,np.ma.masked_invalid(dpr_pp[latstartd:latendd]),
                         cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)

    cb = plt.colorbar(shrink=0.5)
    cb.set_label("Rainrate (mm/h)",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    plt.xlabel("x [km] ",fontsize=ff)
    plt.ylabel("y [km]  ",fontsize=ff)
    plt.title('GPM DPR Rainrate: \n' + str(pfad_gprof_g[66:70]) + '-' +str(pfad_gprof_g[70:72])+ '-' +
              str(pfad_gprof_g[72:74]) + ' T: ' +str(pfad_gprof_g[76:82]) + '-' + str(pfad_gprof_g[84:90]) + ' UTC',fontsize=ff)
    plot_borders(ax37)
    plt.xticks(fontsize=ff)
    plt.yticks(fontsize=ff)
    plt.grid(color='r')
    ax37.set_xlim(ax1.get_xlim())
    ax37.set_ylim(ax1.get_ylim())

    ##################################
    ax29 = fig.add_subplot(234, aspect='equal')
    #pm2 = plt.pcolormesh(gprof_x, gprof_y,np.ma.masked_invalid(gprof_pp[latstart:latend]),
    #                     cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
    pm2 = plt.pcolormesh(gprof_x, gprof_y,np.ma.masked_invalid(gprof_snow[latstart:latend]),
                         cmap=get_miub_cmap(),vmin=0,vmax=5, zorder=1)


    cb = plt.colorbar(shrink=0.5)
    cb.set_label("SnowCoverIndex ",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    plt.xlabel("x [km] ",fontsize=ff)
    plt.ylabel("y [km]  ",fontsize=ff)
    plt.title('GPM GPROF SnowCoverIndex: \n' + str(pfad_gprof_g[66:70]) + '-' +str(pfad_gprof_g[70:72])+ '-' +
              str(pfad_gprof_g[72:74]) + ' T: ' +str(pfad_gprof_g[76:82]) + '-' + str(pfad_gprof_g[84:90]) + ' UTC',fontsize=ff)
    plot_borders(ax29)
    plt.xticks(fontsize=ff)
    plt.yticks(fontsize=ff)
    plt.grid(color='r')
    #plt.tight_layout()
    ax29.set_xlim(ax1.get_xlim())
    ax29.set_ylim(ax1.get_ylim())

    from pcc import get_2_cmap
    ax28 = fig.add_subplot(235, aspect='equal')
    pm2 = plt.pcolormesh(gprof_x, gprof_y,np.ma.masked_invalid(gprof_l[latstart:latend]),
                         cmap=get_2_cmap(),vmin=0,vmax=1, zorder=1)

    cb = plt.colorbar(shrink=0.5)
    cb.set_label("LiquidWaterFraction",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    plt.xlabel("x [km] ",fontsize=ff)
    plt.ylabel("y [km]  ",fontsize=ff)
    plt.title('GPM GPROF LiquidWaterFraction: \n' + str(pfad_gprof_g[66:70]) + '-' +str(pfad_gprof_g[70:72])+ '-' +
              str(pfad_gprof_g[72:74]) + ' T: ' +str(pfad_gprof_g[76:82]) + '-' + str(pfad_gprof_g[84:90]) + ' UTC',fontsize=ff)
    plot_borders(ax28)
    plt.xticks(fontsize=ff)
    plt.yticks(fontsize=ff)
    plt.grid(color='r')
    #plt.tight_layout()
    ax28.set_xlim(ax1.get_xlim())
    ax28.set_ylim(ax1.get_ylim())





#'''

S1 = ['10.65 GHz V-Pol', '10.65 GHz H-Pol','18.7 GHz V-Pol' , '18.7 GHz H-Pol',
      '23.8 GHz V-Pol','36.64 GHz V-Pol',  '36.64 GHz H-Pol','89.0 GHz V-Pol',
      '89.0 GHz H-Pol']

S2 = ['166.0 GHz V-Pol', '166.0 GHz H-Pol','183.31 +/-3 GHz V-Pol', '183.31 +/-7 GHz V-Pol']

qmin,qmax  = 10,90
maxi, mini = 300, 100

s2qmin = []
s2qmax = []



for jj in range(9):
    ff = 15
    fig = plt.figure(figsize=(20,20))
    #ax3 = plt.subplot(3,3,jj+1)
    plot_all()

    ax35 = plt.subplot(236, aspect='equal')
    TT = T_pp[:,:,jj][latstart:latend]
    plt.pcolormesh(T_x, T_y, np.ma.masked_invalid(TT), cmap='jet',#,vmin = mini, vmax = maxi)
                   vmin=np.nanpercentile(TT,qmin), vmax=np.nanpercentile(TT,qmax))
    cb = plt.colorbar(shrink=0.5)
    cb.set_label("Tb (K)",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    plot_borders(ax35)
    plt.xlabel("x [km] ",fontsize=ff)
    plt.ylabel("y [km]  ",fontsize=ff)
    plt.xlim(-420,390)
    plt.ylim(-4700, -3700)
    plt.title(str(S1[jj]))
    plt.tight_layout()
    plt.savefig("/home/velibor/shkgpm/plot/T/Tb_"+ str(jj)[0:5] +".png")
    plt.close()
#plt.show()




for jj in range(4):
    ff = 15
    fig = plt.figure(figsize=(20,20))
    plot_all()
    ax34 = plt.subplot(236, aspect='equal')
    #ax3 = plt.subplot(2,2,jj+1)
    TT2 = T2_pp[:,:,jj][latstart:latend]
    plt.pcolormesh(T2_x, T2_y, np.ma.masked_invalid(TT2), cmap='jet',
                   vmin=np.nanpercentile(TT2,qmin), vmax=np.nanpercentile(TT2,qmax))
    cb = plt.colorbar(shrink=0.5)
    cb.set_label("Tb (K)",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    plot_borders(ax34)
    plt.xlabel("x [km] ",fontsize=ff)
    plt.ylabel("y [km]  ",fontsize=ff)
    plt.xlim(-420,390)
    plt.ylim(-4700, -3700)
    plt.title(str(S2[jj]))

    plt.tight_layout()
    plt.savefig("/home/velibor/shkgpm/plot/T/Tb2_"+ str(jj) +".png")
    plt.close()

#plt.show()
#'''
'''
for ii in range(9):
    A = T_pp[:,:,ii][latstart:latend]
    B = np.ma.masked_invalid(gprof_pp[latstart:latend])
    #A[A<TH_rain] = np.nan
    B[B<2] = np.nan
    #A[A==10] = np.nan

    plt.subplot(3,3,ii+1)
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
    #plt.colorbar()
    plt.plot(B,line,'r-')
    maxA = np.nanmax(yy)
    maxB = np.nanmax(xx)
    minA = np.nanmin(yy)
    minB = np.nanmin(xx)

    plt.xlim(minB,maxB + 1)
    plt.ylim(minA,maxA + 10)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=2, fancybox=True, shadow=True,
                        fontsize='small', title= str(S1[ii]) + "\n Correlation: " +
                                                str(round(r_value,3)) + ', Std_err: '+  str(round(std_err,3)))
    plt.xlabel("GPROF RR [mm/h]")
    plt.ylabel("Tb [K]")
    plt.title(" .")
    plt.tight_layout()
    plt.grid(True)

plt.show()
'''
'''
for ii in range(4):
    A2 = T2_pp[:,:,ii][latstart:latend]
    B2 = np.ma.masked_invalid(gprof_pp[latstart:latend])
    #A[A<TH_rain] = np.nan
    B2[B2<2] = np.nan
    #A[A==10] = np.nan

    plt.subplot(2,2,ii+1)
    mask = ~np.isnan(B2) & ~np.isnan(A2)
    slope, intercept, r_value, p_value, std_err = stats.linregress(B2[mask], A2[mask])
    line = slope*B2+intercept

    xx2 = B2[mask]
    yy2 = A2[mask]
    xedges2, yedges2 = np.linspace(-4, 4, 42), np.linspace(-25, 25, 42)
    hist2, xedges2, yedges2 = np.histogram2d(xx2, yy2, (xedges2, yedges2))
    xidx2 = np.clip(np.digitize(xx2, xedges2), 0, hist.shape[0]-1)
    yidx2 = np.clip(np.digitize(yy2, yedges2), 0, hist.shape[1]-1)
    c = hist[xidx2, yidx2]
    plt.scatter(xx2, yy2, c=c, label='RR [mm/h]')
    #plt.colorbar()
    plt.plot(B,line,'r-')
    maxA2 = np.nanmax(yy2)
    maxB2 = np.nanmax(xx2)
    minA2 = np.nanmin(yy2)
    minB2 = np.nanmin(xx2)

    plt.xlim(minB2,maxB2 + 1)
    plt.ylim(minA2,maxA2 + 10)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=2, fancybox=True, shadow=True,
                        fontsize='small', title= str(S2[ii]) + "\n Correlation: " +
                                                str(round(r_value,3)) + ', Std_err: '+  str(round(std_err,3)))
    plt.xlabel("GPROF RR [mm/h]")
    plt.ylabel("Tb [K]")
    plt.title(" .")
    plt.tight_layout()
    plt.grid(True)

plt.show()
'''
'''
# Nach TRMM TMI
#SI = (451.9 - 0.044 * T_pp[:,:,2]-1.775 * T_pp[:,:,4] + 0.00575 * (T_pp[:,:,4]**2) - T_pp[:,:,7])
#RR = 0.00513 * (SI **1.9468)

# Nach GPM GMI 2010
RR_conv = -0.000011769 * (T_pp[:,:,7]**3)  + 0.0080267 * (T_pp[:,:,7]**2) +1.9461* T_pp[:,:,7]+182.677
RR_strat = -0.0708* T_pp[:,:,7]+19.7034

CSP = 0.5
RR = RR_conv * CSP + RR_strat * (1-CSP)

ax3 = plt.subplot(1,1,1)
TT = RR[latstart:latend]
plt.pcolormesh(T_x, T_y, np.ma.masked_invalid(TT), cmap=my_cmap, vmin=np.nanmin(TT), vmax=np.nanmax(TT))
cb = plt.colorbar(shrink=0.5)
cb.set_label("RR (mm/h)",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plot_borders(ax3)
plt.xlabel("x [km] ",fontsize=ff)
plt.ylabel("y [km]  ",fontsize=ff)
plt.xlim(-420,390)
plt.ylim(-4700, -3700)
plt.tight_layout()

plt.show()
'''