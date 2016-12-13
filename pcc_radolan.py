"""

Das Program dient der Veranschaulichung der 5 Minutigen RX Radolan Daten!

"""


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



# Zeitstempel nach YYYYMMDDhhmmss
ZP = '20160707155500'#'20161024232500'#'20140609132500'#'20160917102000'#'20160917102000'#'20160805054500'#'20141007023500'
year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
ye = ZP[2:4]



## Read RADOLAN GK Koordinaten
## ----------------------------
iii = 0
pfad = ('/automount/radar/dwd/rx/'+str(year)+'/'+str(year)+'-'+str(m)+'/'+
        str(year)+'-'+str(m)+'-'+str(d)+'/raa01-rx_10000-'+str(ye)+str(m)+
        str(d)+str(ht)+str(mt)+'-dwd---bin.gz')

pfad_radolan = pfad[:-3]

rw_filename = wradlib.util.get_wradlib_data_file(pfad)
rwdata, rwattrs = wradlib.io.read_RADOLAN_composite(rw_filename)

rwdata = np.ma.masked_equal(rwdata, -9999) / 2 - 32.5

radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
x = radolan_grid_xy[:,:,0]
y = radolan_grid_xy[:,:,1]


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

dataset1, inLayer1 = wradlib.io.open_shape('/automount/db01/python/data/ADM/'
                                           'germany/vg250_0101.gk3.shape.ebenen'
                                           '/vg250_ebenen/vg250_l.shp')

import matplotlib.cm as cm
my_cmap = cm.get_cmap('jet',40)
my_cmap.set_under('lightgrey')
my_cmap.set_over('darkred')



Z = wradlib.trafo.idecibel(rwdata)
rwdata = wradlib.zr.z2r(Z, a=200., b=1.6)


########################################################################## PLOT

ff = 15
fig = plt.figure(figsize=(10,10))
ax1 = fig.add_subplot(111, aspect='equal')
plt.pcolormesh(x, y, rwdata, cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
#plt.scatter(x, y, rwdata, cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
cb = plt.colorbar(shrink=0.8)
cb.set_label("Rainrate (mm/h)",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plot_borders(ax1)
plt.title('RADOLAN Rainrate: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+
          str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff)

#plot_ocean(ax1)
plt.xlabel("x [km] ",fontsize=ff)
plt.ylabel("y [km]  ",fontsize=ff)
#plt.xticks(fontsize=0)
#plt.yticks(fontsize=0)
plt.grid(color='r')
plt.xlim(-420,390)
plt.ylim(-4700, -3700)



plt.show()
