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
ZP = '20160917102000'#'20161024232500'#'20140609132500'#'20160917102000'#'20160917102000'#'20160805054500'#'20141007023500'
# '20160904135000'
year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
ye = ZP[2:4]


#ToDo: Veralgemeinern bzgl aller Parameter:
## Read GPROF
## ------------
#pfad2 = ('/home/velibor/shkgpm/data/'+str(year)+str(m)+str(d)+'/corra/*.HDF5')

pfad2 = ('/home/velibor/shkgpm/data/example/dpr/*.HDF5')

pfad_gprof = glob.glob(pfad2)
pfad_gprof_g = pfad_gprof[0]


gpmdprs = h5py.File(pfad_gprof_g, 'r')
gprof_lat=np.array(gpmdprs['MS']['Latitude'])			#(7934, 24)
gprof_lon=np.array(gpmdprs['MS']['Longitude'])			#(7934, 24)
#gprof_pp=np.array(gpmdprs['NS']['surfPrecipTotRate'])
#gprof_pp=np.array(gpmdprs['MS']['tenMeterWindSpeed'])

gprof_pp=np.array(gpmdprs['MS']['DSD']['binNode'])
#gprof_pp=np.array(gpmdprs['MS']['pia'])

#gprof_pp = gprof_pp[:,:,:,0]

#gprof_pp[gprof_pp==-9999.9] = np.nan

print gprof_pp.shape, gprof_lat.shape



bonn_lat1 = 47.9400
bonn_lat2 = 55.3500
bonn_lon1 = 6.40000
bonn_lon2 = 14.10000

#Europa
#bonn_lat1 = 30.9400
#bonn_lat2 = 70.3500
#bonn_lon1 = 0.40000
#bonn_lon2 = 50.10000


ilat= np.where((gprof_lat>bonn_lat1) & (gprof_lat<bonn_lat2))
ilon= np.where((gprof_lon>bonn_lon1) & (gprof_lon<bonn_lon2))
lonstart = ilon[0][0]
lonend = ilon[0][-1]
latstart = ilat[0][0]
latend = ilat[0][-1]

alon = gprof_lon[latstart:latend]
alat = gprof_lat[latstart:latend]
gprof_pp_a = gprof_pp[latstart:latend]


ailat= np.where((alat>bonn_lat1) & (alat<bonn_lat2))
ailon= np.where((alon>bonn_lon1) & (alon<bonn_lon2))
alonstart = ailon[0][0]
alonend = ailon[0][-1]
alatstart = ailat[0][0]
alatend = ailat[0][-1]

blon = alon[alonstart:alonend]
blat = alat[alonstart:alonend]
gprof_pp_b = gprof_pp_a[alonstart:alonend]


## Nur wenn genzer Swath genutzt wird
#gprof_pp_b = gprof_pp
#blat, blon = gprof_lat, gprof_lon

###########PROJECTION

proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

gpm_x, gpm_y = wradlib.georef.reproject(blon, blat, projection_target=proj_stereo , projection_source=proj_wgs)
grid_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()


print np.nanmin(gprof_pp_b)
print np.nanmax(gprof_pp_b)



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

cut = 13

print np.nanmax(gprof_pp_b[:,cut,:])
print gprof_pp_b.shape, gpm_x.shape
## PLOT
## ----
ff = 15
fig = plt.figure(figsize=(10,10))
ax2 = fig.add_subplot(121, aspect='auto')
pm1 = plt.imshow(np.transpose(gprof_pp_b[:,cut,:]),vmin=0, vmax=176, aspect='auto')

cb = plt.colorbar(shrink=0.4)
cb.set_label("Ref (dBZ)",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.xlabel("x [km] ",fontsize=ff)
plt.ylabel("z [km]  ",fontsize=ff)
plt.grid()

ax2 = fig.add_subplot(122, aspect='equal')
pm2 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(gprof_pp_b[:,:,0]),
                     cmap=my_cmap,vmin=0,vmax=176, zorder=2)

cb = plt.colorbar(shrink=0.4)
cb.set_label("Rainrate (mm/h)",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.xlabel("x [km] ",fontsize=ff)
plt.ylabel("y [km]  ",fontsize=ff)
plot_borders(ax2)
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)
plt

plt.tight_layout()

plt.show()




#plt.pcolormesh(l,h ,np.transpose(np.ma.masked_invalid(gprof_pp_b[:,cut,:])),
 #                    cmap=my_cmap,vmin=0,vmax=1, zorder=2)
#plt.show()

