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

## Read RADOLAN GK Koordinaten
## ----------------------------
iii = 0
#pfad = ('/user/velibor/SHKGPM/data/radolan/*bin')
pfad = ('/automount/radar/dwd/rx/2014/2014-10/2014-10-16/raa01-rx_10000-1410162310-dwd---bin')
pfad_radolan= sorted(glob.glob(pfad))
pfad_radolan = pfad_radolan[iii]


rw_filename = wradlib.util.get_wradlib_data_file(pfad_radolan)
rwdata, rwattrs = wradlib.io.read_RADOLAN_composite(rw_filename)

for key, value in rwattrs.items():
    print(key + ':', value)

rwdata = np.ma.masked_equal(rwdata, -9999) / 2 - 32.5

#sec = rwattrs['secondary']
#rwdata.flat[sec] = -9999
#rwdata = np.ma.masked_equal(rwdata, -9999)
radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
x = radolan_grid_xy[:,:,0]
y = radolan_grid_xy[:,:,1]
Z = wradlib.trafo.idecibel(rwdata)
rwdata = wradlib.zr.z2r(Z, a=200., b=1.6)

## Radolan lat lon Koordinaten
# ------------------------------
'''
## mask data
sec = rwattrs['secondary']
rwdata.flat[sec] = -9999
rwdata = np.ma.masked_equal(rwdata, -9999)

## create radolan projection object
proj_stereo = wradlib.georef.create_osr("dwd-radolan")

## create wgs84 projection object
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

## get radolan grid
radolan_grid_xy = wradlib.georef.get_radolan_grid(900, 900)
x1 = radolan_grid_xy[:, :, 0]
y1 = radolan_grid_xy[:, :, 1]

## convert to lonlat
radolan_grid_ll = wradlib.georef.reproject(radolan_grid_xy,
                                       projection_source=proj_stereo,
                                       projection_target=proj_wgs)
lon1 = radolan_grid_ll[:, :, 0]
lat1 = radolan_grid_ll[:, :, 1]
'''


## Read GPROF
## ------------
#pfad2 = ('/home/velibor/shkgpm/data/20140921/gprof/*.HDF5')
pfad2 = ('/home/velibor/shkgpm/data/20141016/gprof/*.HDF5')

pfad_gprof = glob.glob(pfad2)
pfad_gprof_g = pfad_gprof[0]

print pfad_gprof_g

gpmgmi = h5py.File(pfad_gprof_g, 'r')

gpmgmi.keys()
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

## GPM lon/lat in GK
## ------------
#proj_gk = osr.SpatialReference()
#proj_gk.ImportFromEPSG(31466)
#proj_ll = osr.SpatialReference()
#proj_ll.ImportFromEPSG(4326)
#gk3 = wradlib.georef.epsg_to_osr(31467)

proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

#, ,np.ma.masked_invalid(gprof_pp[latstart:latend]
#gpm_x, gpm_y = wrl.georef.reproject(gprof_lon[latstart:latend], gprof_lat[latstart:latend], projection_source=proj_ll,projection_target=proj_gk)
gpm_x, gpm_y = wradlib.georef.reproject(gprof_lon[latstart:latend], gprof_lat[latstart:latend], projection_target=proj_stereo , projection_source=proj_wgs)
grid_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()





## Landgrenzenfunktion
## -------------------
'''
def plot_ocean(ax):
    # open the input data source and get the layer
    filename = os.path.join('/automount/db01/python/data/NED/10m/physical/10m_physical/ne_10m_ocean.shp')
    dataset, inLayer = wradlib.io.open_shape(filename)
    ocean, keys = wradlib.georef.get_shape_coordinates(inLayer)
    wradlib.vis.add_lines(ax, ocean, color='black', lw=2 , zorder=4)
'''
def plot_ocean(ax):

    filename = os.path.join('/automount/db01/python/data/NED/10m/physical/10m_physical/ne_10m_ocean.shp')
    dataset, inLayer = wradlib.io.open_shape(filename)
    inLayer.SetSpatialFilterRect(1, 45, 19, 56.5)
    borders, keys = wradlib.georef.get_shape_coordinates(inLayer)
    proj_gk = osr.SpatialReference()
    proj_gk.ImportFromEPSG(31466)
    proj_ll = osr.SpatialReference()
    proj_ll.ImportFromEPSG(4326)
    gk3 = wradlib.georef.epsg_to_osr(31467)
    proj_stereo = wrl.georef.create_osr("dwd-radolan")
    proj_wgs = osr.SpatialReference()
    proj_wgs.ImportFromEPSG(4326)
    print borders.shape

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
'''
def plot_dem(ax):
    filename = wradlib.util.get_wradlib_data_file('/home/velibor/cosmo_de_4326.tif')
    # pixel_spacing is in output units (lonlat)
    rastercoords, rastervalues = wradlib.io.read_raster_data(filename,
                                                         spacing=0.005)
    # specify kwargs for plotting, using terrain colormap and LogNorm
    dem = ax.pcolormesh(rastercoords[..., 0], rastercoords[..., 1],
                        rastervalues, cmap=plt.cm.gist_earth, norm=LogNorm(),
                        vmin=1, vmax=3000, zorder=1)
    # make some space on the right for colorbar axis
    div1 = make_axes_locatable(ax)
    #cax1 = div1.append_axes("right", size="5%", pad=0.1)
    # add colorbar and title
    # we use LogLocator for colorbar
    #cb = plt.gcf().colorbar(dem, cax=cax1,
                           #ticks=ticker.LogLocator(subs=range(10)))
    #cb.set_label('terrain height [m]')

def get_miub_cmap ():
    startcolor = 'white' # a dark olive
    color1 = '#8ec7ff' #'cyan' # a bright yellow
    color2 = 'dodgerblue'
    color3 = 'lime'
    color4 = 'yellow'
    color5 = 'darkorange'
    color6 = 'red'
    color7 = 'purple'
    #color6 = 'grey'
    endcolor = 'darkmagenta' # medium dark red
    colors = [startcolor, color1, color2, color3, color4, color5, color6, endcolor]
    return col.LinearSegmentedColormap.from_list('miub1' ,colors)
'''

dataset1, inLayer1 = wradlib.io.open_shape('/automount/db01/python/data/ADM/germany/vg250_0101.gk3.shape.ebenen/vg250_ebenen/vg250_l.shp')

import matplotlib.cm as cm
my_cmap = cm.get_cmap('jet',40)
my_cmap.set_under('lightgrey')
my_cmap.set_over('darkred')



#INTERLOLATION
gk3 = wradlib.georef.epsg_to_osr(31467)

grid_gpm_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose() # GPM Grid erschaffen

xy = np.vstack((x.ravel(), y.ravel())).transpose()

result = wrl.ipol.interpolate(xy, grid_gpm_xy, rwdata.reshape(900*900,1), wrl.ipol.Idw, nnearest=40)

result = np.ma.masked_invalid(result)
#Todo: rwdata muss vonn 900x900 reshaped werden auf die gleiche Form wie xy (810000, 2) FEHLER IWO
#SCHEMA http://wradlib.org/wradlib-docs/latest/notebooks/interpolation/wradlib_ipol_example.html



## PLOT
## ----
ff = 15
fig = plt.figure()
ax1 = fig.add_subplot(221, aspect='equal')
plt.pcolormesh(x, y, rwdata, cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
cb = plt.colorbar(shrink=0.8)
cb.set_label("Rainrate (mm/h)",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plot_borders(ax1)
plt.title('RADOLAN Rainrate: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff) #RW Product Polar Stereo

#plot_ocean(ax1)
plt.xlabel("Longitude ",fontsize=ff)
plt.ylabel("Latitude  ",fontsize=ff)
#plt.xticks(fontsize=0)
#plt.yticks(fontsize=0)
plt.grid(color='r')
plt.xlim(-420,390)
plt.ylim(-4700, -3700)


ax2 = fig.add_subplot(222, aspect='equal')
pm2 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(gprof_pp[latstart:latend]),
                     cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
#pm2 = plt.pcolormesh(gprof_lon[latstart:latend], gprof_lat[latstart:latend],np.ma.masked_invalid(gprof_pp[latstart:latend]),
                     #cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
cb = plt.colorbar(shrink=0.8)
cb.set_label("Rainrate (mm/h)",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.xlabel("Longitude ",fontsize=ff)
plt.ylabel("Latitude  ",fontsize=ff)
plt.title('GPM GPROF Rainrate: \n' + str(pfad_gprof_g[66:70]) + '-' +str(pfad_gprof_g[70:72])+ '-' +
          str(pfad_gprof_g[72:74]) + ' T: ' +str(pfad_gprof_g[76:82]) + '-' + str(pfad_gprof_g[84:90]) + ' UTC',fontsize=ff)
#plot_ocean(ax2)
plot_borders(ax2)
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)
#plt.xlim((bonn_lon1-1,bonn_lon2+1))
#plt.ylim((bonn_lat1-1,bonn_lat2+1))
plt.grid(color='r')
plt.tight_layout()
#Limits Setzen
ax2.set_xlim(ax1.get_xlim())
ax2.set_ylim(ax1.get_ylim())
#plt.xticks(fontsize=0)
#plt.yticks(fontsize=0)


rrr = result.reshape(gpm_x.shape)
rrr[rrr==10] = np.nan

ax2 = fig.add_subplot(223, aspect='equal')
pm2 = plt.pcolormesh(gpm_x, gpm_y,rrr,
                     cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)

cb = plt.colorbar(shrink=0.8)
cb.set_label("Rainrate (mm/h)",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.xlabel("Longitude ",fontsize=ff)
plt.ylabel("Latitude  ",fontsize=ff)
plt.title('RADOLAN Rainrate Interpoliert: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff) #RW Product Polar Stereo
#plot_ocean(ax2)
plot_borders(ax2)
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)
#plt.xlim((bonn_lon1-1,bonn_lon2+1))
#plt.ylim((bonn_lat1-1,bonn_lat2+1))
plt.grid(color='r')
plt.tight_layout()
#Limits Setzen
ax2.set_xlim(ax1.get_xlim())
ax2.set_ylim(ax1.get_ylim())
#plt.xticks(fontsize=0)
#plt.yticks(fontsize=0)


ax2 = fig.add_subplot(224, aspect='equal')

A = rrr
B = np.ma.masked_invalid(gprof_pp[latstart:latend])
A[A<0.1] = np.nan
B[B<0.1] = np.nan
#A[A==10] = np.nan

mask = ~np.isnan(B) & ~np.isnan(A)
slope, intercept, r_value, p_value, std_err = stats.linregress(B[mask], A[mask])
line = slope*B+intercept
plt.scatter(B,A, color='blue', label='RR [mm/h]')
plt.plot(B,line,'r-')
maxAB = np.nanmax([np.nanmax(A),np.nanmax(B)])
plt.xlim(0,maxAB + 1)
plt.ylim(0,maxAB + 1)
legend = plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=2, fancybox=True, shadow=True,
                    fontsize='small', title="________"+"_vs_BoxPol________" + "\n Slope: " + str(round(slope,3))
                                            + ', Intercept: '+  str(round(intercept,3)) + "\n Correlation: " +
                                            str(round(r_value,3)) + ', Std_err: '+  str(round(std_err,3)))
plt.xlabel("GPROF RR [mm/h]")
plt.ylabel("RADOLAN RR [mm/h]")
plt.title(" .")

plt.grid(True)
plt.show()


