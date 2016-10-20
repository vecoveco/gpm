"""Einlesen und darstellen von GPM und Radolan Dateien"""


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
iii = 1
pfad = ('/user/velibor/SHKGPM/data/radolan/*.gz')
pfad_radolan= sorted(glob.glob(pfad))
pfad_radolan = pfad_radolan[iii]


rw_filename = wradlib.util.get_wradlib_data_file(pfad_radolan)
rwdata, rwattrs = wradlib.io.read_RADOLAN_composite(rw_filename)

sec = rwattrs['secondary']
rwdata.flat[sec] = -9999
rwdata = np.ma.masked_equal(rwdata, -9999)
radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
x = radolan_grid_xy[:,:,0]
y = radolan_grid_xy[:,:,1]


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
pfad2 = ('/home/velibor/shkgpm/data/20140921/gprof/*.HDF5')
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
proj_gk = osr.SpatialReference()
proj_gk.ImportFromEPSG(31466)
proj_ll = osr.SpatialReference()
proj_ll.ImportFromEPSG(4326)
gk3 = wradlib.georef.epsg_to_osr(31467)
proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

#, ,np.ma.masked_invalid(gprof_pp[latstart:latend]
#gpm_x, gpm_y = wrl.georef.reproject(gprof_lon[latstart:latend], gprof_lat[latstart:latend], projection_source=proj_ll,projection_target=proj_gk)
gpm_x, gpm_y = wradlib.georef.reproject(gprof_lon[latstart:latend], gprof_lat[latstart:latend], projection_target=proj_stereo , projection_source=proj_wgs)
grid_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()

print gpm_x.shape, gpm_y.shape



## Landgrenzenfunktion
## -------------------

def plot_ocean(ax):
    # open the input data source and get the layer
    filename = os.path.join('/automount/db01/python/data/NED/10m/physical/10m_physical/ne_10m_ocean.shp')
    dataset, inLayer = wradlib.io.open_shape(filename)
    ocean, keys = wradlib.georef.get_shape_coordinates(inLayer)
    wradlib.vis.add_lines(ax, ocean, color='black', lw=2 , zorder=4)

def plot_borders(ax):
    # plot country borders from esri vector shape, filter by attribute
    # create wgs84 and india osr objects (spatial reference system)
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
    countries = ['Germany']
    # open the input data source and get the layer
    filename = wradlib.util.get_wradlib_data_file('geo/ne_10m_admin_0_boundary_'
                                              'lines_land.shp')
    dataset, inLayer = wradlib.io.open_shape(filename)
    # iterate over countries, filter accordingly, get coordinates and plot
    for item in countries:
        print item
        # SQL-like selection syntax
        fattr = "(adm0_left = '" + item + "' or adm0_right = '" + item + "')"
        inLayer.SetAttributeFilter(fattr)
        # get borders and names
        borders, keys = wradlib.georef.get_shape_coordinates(inLayer, key='name')
        #print borders[0].shape
        #print borders[12][:,1]
        #print borders[12][:,0]
        for j in range(14):
            bu = np.array(borders[j].shape)
            #print bu.shape
            a = np.array(bu.shape)
            print borders[j].shape
            #print borders[0]
            #print a
            print borders[j].shape[0]
            if a==1:
                for i in range(0,borders[j].shape[0],1):
                    bordx, bordy = wrl.georef.reproject(borders[j][i][:,0], borders[j][i][:,1], projection_source=proj_wgs, projection_target=proj_stereo)
                    bord_xy = np.vstack((bordx.ravel(), bordy.ravel())).transpose()
                    #print bord_xy
                    wradlib.vis.add_lines(ax, bord_xy, color='black', lw=2, zorder=3)
            if a==2:    #richtig
                bordx, bordy = wrl.georef.reproject(borders[j][:,0], borders[j][:,1], projection_source=proj_wgs, projection_target=proj_stereo)
                bord_xy = np.vstack((bordx.ravel(), bordy.ravel())).transpose()
            #print bord_xy
                wradlib.vis.add_lines(ax, bord_xy, color='black', lw=2, zorder=3)
            #print j
            #print bordx
            bord_xy = np.vstack((bordx.ravel(), bordy.ravel())).transpose()
            #print bord_xy
            wradlib.vis.add_lines(ax, bord_xy, color='black', lw=2, zorder=3)
        #print bordx.shape, bordy.shape, bord_xy.shape
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



#veli = np.ma.masked_where(rwdata<=0, rwdata)
print pfad_radolan
## PLOT
## ----
ff = 15
fig = plt.figure()
ax1 = fig.add_subplot(121, aspect='equal')
plt.pcolormesh(x, y, rwdata, cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
cb = plt.colorbar(shrink=0.3)
cb.set_label("Rainrate (mm/h)",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plot_borders(ax1)
plt.title('RADOLAN Rainrate: \n'+'20' + str(pfad_radolan[-23:-21])+'-'+str(pfad_radolan[-21:-19])+'-'+str(pfad_radolan[-19:-17])+' T: '+str(pfad_radolan[-17:-13]) + '00 UTC',fontsize=ff) #RW Product Polar Stereo
#plot_dem(ax1)
#plot_ocean(ax1)
plt.xlabel("Longitude",fontsize=ff)
plt.ylabel("Latitude ",fontsize=ff)
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)

plt.grid(color='r')
plt.xlim(-420,390)
plt.ylim(-4700, -3700)


ax2 = fig.add_subplot(122, aspect='equal')
pm2 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(gprof_pp[latstart:latend]),
                     cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
#pm2 = plt.pcolormesh(gprof_lon[latstart:latend], gprof_lat[latstart:latend],np.ma.masked_invalid(gprof_pp[latstart:latend]),
                     #cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
cb = plt.colorbar(shrink=0.3)
cb.set_label("Rainrate (mm/h)",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.xlabel("Longitude",fontsize=ff)
plt.ylabel("Latitude ",fontsize=ff)
plt.title('GPM GPROF Rainrate: \n' + str(pfad_gprof_g[66:70]) + '-' +str(pfad_gprof_g[70:72])+ '-' +
          str(pfad_gprof_g[72:74]) + ' T: ' +str(pfad_gprof_g[76:82]) + '-' + str(pfad_gprof_g[84:90]) + ' UTC',fontsize=ff)
#plot_ocean(ax2)
#plot_dem(ax2)
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
plt.show()





#INTERLOLATION
'''
gk3 = wradlib.georef.epsg_to_osr(31467)
x, y = wradlib.georef.reproject(lon1, lat1, projection_target=gk3)
xgrid, ygrid = wradlib.georef.reproject(gprof_lon[latstart:latend], gprof_lat[latstart:latend],
                                        projection_target=gk3)

grid_xy = np.vstack((xgrid.ravel(), ygrid.ravel())).transpose()

xy=np.concatenate([x.ravel()[:,None],y.ravel()[:,None]], axis=1)

print (x.shape, y.shape, xgrid.shape, ygrid.shape, grid_xy.shape, xy.shape )

#gridded = wradlib.comp.togrid(xy, grid_xy, ranges[-1], np.array([x.mean(), y.mean()]), R.ravel(),ipoli[0],nnearest=5,p=2)
#gridded = np.ma.masked_invalid(gridded).reshape(xgrid.shape)
grid = wradlib.ipol.Idw(xy, grid_xy, nnearest=4, p=2.)
print grid
# Interpolation objects
#idw = ipol.Idw(src, trg)
'''