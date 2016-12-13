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

ZP = '20141007023500'

# Zeitstempel nach YYYYMMDDhhmmss
year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
ye = ZP[2:4]



## Read RADOLAN GK Koordinaten
## ----------------------------
iii = 0
pfad = ('/automount/radar/dwd/rx/'+str(year)+'/'+str(year)+'-'+str(m)+'/'+str(year)+'-'+str(m)+'-'+str(d)+'/raa01-rx_10000-'+str(ye)+str(m)+str(d)+str(ht)+str(mt)+'-dwd---bin.gz')

pfad_radolan = pfad[:-3]

####### pfad

rw_filename = wradlib.util.get_wradlib_data_file(pfad_radolan)
rwdata, rwattrs = wradlib.io.read_RADOLAN_composite(rw_filename)

rwdata = np.ma.masked_equal(rwdata, -9999) / 2 - 32.5

radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
x = radolan_grid_xy[:,:,0]
y = radolan_grid_xy[:,:,1]





## Read Corra
## ------------
pfadT = ('/home/velibor/shkgpm/data/example/T/*.HDF5')
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

pfad3 = ('/home/velibor/shkgpm/data/example/gprof/*.HDF5')
pfad_gprof = glob.glob(pfad3)
print pfad_gprof
pfad_gprof_g = pfad_gprof[0]
gpmgmi = h5py.File(pfad_gprof_g, 'r')
gpmgmi_S1=gpmgmi['S1']
gprof_lat=np.array(gpmgmi_S1['Latitude'])
gprof_lon=np.array(gpmgmi_S1['Longitude'])
gprof_pp=np.array(gpmgmi_S1['surfacePrecipitation'])
gprof_pp[gprof_pp<=0] = np.nan




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



proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)


gprof_x, gprof_y = wradlib.georef.reproject(gprof_lon[latstart:latend], gprof_lat[latstart:latend], projection_target=proj_stereo , projection_source=proj_wgs)


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

def plot_borders(ax):

    from osgeo import osr
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(4326)
    india = osr.SpatialReference()
    # asia south albers equal area conic
    india.ImportFromEPSG(102028)

    proj_gk = osr.SpatialReference()
    proj_gk.ImportFromEPSG(31466)
    proj_ll = osr.SpatialReference()
    proj_ll.ImportFromEPSG(4326)
    gk3 = wradlib.georef.epsg_to_osr(31467)
    proj_stereo = wrl.georef.create_osr("dwd-radolan")
    proj_wgs = osr.SpatialReference()
    proj_wgs.ImportFromEPSG(4326)

    # country list
    countries = ['Germany']#,'France','Denmark', 'Netherlands', 'Poland']
    # open the input data source and get the layer
    filename = wradlib.util.get_wradlib_data_file('/automount/db01/python/data/NED/10m/cultural/10m_cultural/10m_cultural/ne_10m_admin_0_countries.shp')
    dataset, inLayer = wradlib.io.open_shape(filename)
    # iterate over countries, filter accordingly, get coordinates and plot
    for item in countries:
        #print item
        # SQL-like selection syntax
        fattr = "(name='"+item+"')"
        inLayer.SetAttributeFilter(fattr)
        # get borders and names
        borders, keys = wradlib.georef.get_shape_coordinates(inLayer, key='name')

        for j in range(borders.shape[0]):
            bu = np.array(borders[j].shape)
            a = np.array(bu.shape)

            if a==1:
                for i in range(0,borders[j].shape[0],1):
                    bordx, bordy = wrl.georef.reproject(borders[j][i][:,0], borders[j][i][:,1], projection_source=proj_wgs, projection_target=proj_stereo)
                    bord_xy = np.vstack((bordx.ravel(), bordy.ravel())).transpose()

                    wradlib.vis.add_lines(ax, bord_xy, color='black', lw=2, zorder=3)
            if a==2:    #richtig
                bordx, bordy = wrl.georef.reproject(borders[j][:,0], borders[j][:,1], projection_source=proj_wgs, projection_target=proj_stereo)
                bord_xy = np.vstack((bordx.ravel(), bordy.ravel())).transpose()

                wradlib.vis.add_lines(ax, bord_xy, color='black', lw=2, zorder=3)

            bord_xy = np.vstack((bordx.ravel(), bordy.ravel())).transpose()

            wradlib.vis.add_lines(ax, bord_xy, color='black', lw=2, zorder=3)
    ax.autoscale()


dataset1, inLayer1 = wradlib.io.open_shape('/automount/db01/python/data/ADM/germany/vg250_0101.gk3.shape.ebenen/vg250_ebenen/vg250_l.shp')

import matplotlib.cm as cm
my_cmap = cm.get_cmap('jet',40)
my_cmap.set_under('lightgrey')
my_cmap.set_over('darkred')




Z = wradlib.trafo.idecibel(rwdata)
rwdata = wradlib.zr.z2r(Z, a=200., b=1.6)


######################################################################################## PLOT

ff = 15
fig = plt.figure(figsize=(10,10))
ax1 = fig.add_subplot(121, aspect='equal')
plt.pcolormesh(x, y, rwdata, cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
cb = plt.colorbar(shrink=0.5)
cb.set_label("Rainrate (mm/h)",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plot_borders(ax1)
plt.title('RADOLAN Rainrate: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff) #RW Product Polar Stereo

plt.xlabel("x [km] ",fontsize=ff)
plt.ylabel("y [km]  ",fontsize=ff)
plt.grid(color='r')
plt.xlim(-420,390)
plt.ylim(-4700, -3700)


ax2 = fig.add_subplot(122, aspect='equal')
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

plt.show()



S1 = ['10.65 GHz V-Pol', '10.65 GHz H-Pol','18.7 GHz V-Pol' , '18.7 GHz H-Pol',
      '23.8 GHz V-Pol','36.64 GHz V-Pol',  '36.64 GHz H-Pol','89.0 GHz V-Pol',
      '89.0 GHz H-Pol']

S2 = ['166.0 GHz V-Pol', '166.0 GHz H-Pol','183.31 +/-3 GHz V-Pol', '183.31 +/-7 GHz V-Pol']


ff = 15
fig = plt.figure(figsize=(20,15))

for jj in range(9):
    ax3 = plt.subplot(3,3,jj+1)
    TT = T_pp[:,:,jj][latstart:latend]
    plt.pcolormesh(T_x, T_y, np.ma.masked_invalid(TT), cmap=my_cmap, vmin=np.nanmin(TT), vmax=np.nanmax(TT))
    cb = plt.colorbar(shrink=0.5)
    cb.set_label("Tb (K)",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    plot_borders(ax3)
    plt.xlabel("x [km] ",fontsize=ff)
    plt.ylabel("y [km]  ",fontsize=ff)
    plt.xlim(-420,390)
    plt.ylim(-4700, -3700)
    plt.title(str(S1[jj]))
    plt.tight_layout()

plt.show()


ff = 15
fig = plt.figure(figsize=(20,15))

for jj in range(4):
    ax3 = plt.subplot(2,2,jj+1)
    TT = T2_pp[:,:,jj][latstart:latend]
    plt.pcolormesh(T2_x, T2_y, np.ma.masked_invalid(TT), cmap=my_cmap, vmin=np.nanmin(TT), vmax=np.nanmax(TT))
    cb = plt.colorbar(shrink=0.5)
    cb.set_label("Tb (K)",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    plot_borders(ax3)
    plt.xlabel("x [km] ",fontsize=ff)
    plt.ylabel("y [km]  ",fontsize=ff)
    plt.xlim(-420,390)
    plt.ylim(-4700, -3700)
    plt.title(str(S2[jj]))
    plt.tight_layout()

plt.show()

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

    ##############################################################################
    #plt.hist2d(B[mask],A[mask], bins=100, norm=LogNorm())
    #plt.hexbin(B[mask],A[mask])
    #plt.colorbar()
    #plt.show()

    xx = B[mask]
    yy = A[mask]
    xedges, yedges = np.linspace(-4, 4, 42), np.linspace(-25, 25, 42)
    hist, xedges, yedges = np.histogram2d(xx, yy, (xedges, yedges))
    xidx = np.clip(np.digitize(xx, xedges), 0, hist.shape[0]-1)
    yidx = np.clip(np.digitize(yy, yedges), 0, hist.shape[1]-1)
    c = hist[xidx, yidx]
    plt.scatter(xx, yy, c=c, label='RR [mm/h]')
    plt.colorbar()
    plt.plot(B,line,'r-')
    maxA = np.nanmax(yy)
    maxB = np.nanmax(xx)
    minA = np.nanmin(yy)
    minB = np.nanmin(xx)

    plt.xlim(minB,maxB + 1)
    plt.ylim(minA,maxA + 1)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=2, fancybox=True, shadow=True,
                        fontsize='small', title= "Slope: " + str(round(slope,3))
                                                + ', Intercept: '+  str(round(intercept,3)) + "\n Correlation: " +
                                                str(round(r_value,3)) + ', Std_err: '+  str(round(std_err,3)))
    plt.xlabel("GPROF RR [mm/h]")
    plt.ylabel("Tb [K]")
    plt.title(" .")

    plt.grid(True)

plt.show()