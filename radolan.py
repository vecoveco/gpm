"""Einlesen und darstellen von GPM und Radolan Dateien"""


import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import wradlib
import glob
import math
import pandas as pd
from scipy import stats
import matplotlib as mpl

ipoli = [wradlib.ipol.Idw, wradlib.ipol.Linear, wradlib.ipol.Nearest, wradlib.ipol.OrdinaryKriging]

## Read RADOLAN
## ------------

pfad = ('/user/velibor/SHKGPM/data/radolan/*.gz')
pfad_radolan= sorted(glob.glob(pfad))
pfad_radolan = pfad_radolan[0]

rw_filename = wradlib.util.get_wradlib_data_file(pfad_radolan)
rwdata, rwattrs = wradlib.io.read_RADOLAN_composite(rw_filename)

sec = rwattrs['secondary']
rwdata.flat[sec] = -9999
rwdata = np.ma.masked_equal(rwdata, -9999)

radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
x = radolan_grid_xy[:,:,0]
y = radolan_grid_xy[:,:,1]

## Read GPROF
## ------------

pfad2 = ('/home/velibor/shkgpm/data/20140921/gprof/*.HDF5')
pfad_gprof = glob.glob(pfad2)
pfad_gprof_g = pfad_gprof[0]

gpmgmi = h5py.File(pfad_gprof_g, 'r')

gpmgmi.keys()
gpmgmi_S1=gpmgmi['S1']
gprof_lat=np.array(gpmgmi_S1['Latitude'])
gprof_lon=np.array(gpmgmi_S1['Longitude'])
gprof_pp=np.array(gpmgmi_S1['surfacePrecipitation'])
gprof_pp[gprof_pp==-9999] = np.nan


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

def plot_borders(ax):
    # plot country borders from esri vector shape, filter by attribute
    # create wgs84 and india osr objects (spatial reference system)
    from osgeo import osr
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(4326)
    india = osr.SpatialReference()
    # asia south albers equal area conic
    india.ImportFromEPSG(102028)

    # country list
    countries = ['Germany']
    # open the input data source and get the layer
    filename = wradlib.util.get_wradlib_data_file('geo/ne_10m_admin_0_boundary_'
                                              'lines_land.shp')
    dataset, inLayer = wradlib.io.open_shape(filename)
    # iterate over countries, filter accordingly, get coordinates and plot
    for item in countries:
        # SQL-like selection syntax
        fattr = "(adm0_left = '" + item + "' or adm0_right = '" + item + "')"
        inLayer.SetAttributeFilter(fattr)
        # get borders and names
        borders, keys = wradlib.georef.get_shape_coordinates(inLayer, key='name')
        wradlib.vis.add_lines(ax, borders, color='black', lw=2, zorder=4)
    ax.autoscale()


#radar_location = (lon_ppi, lat_ppi, alt_ppi)
#elevation = 1.5
#azimuths = az
#ranges = r
#polargrid = np.meshgrid(ranges, azimuths)
#lon, lat, alt = wradlib.georef.polar2lonlatalt_n(polargrid[0], polargrid[1], elevation, radar_location)

#gk3 = wradlib.georef.epsg_to_osr(31467)
#x, y = wradlib.georef.reproject(x, y, projection_target=gk3)
#xgrid, ygrid = wradlib.georef.reproject(gprof_lon[latstart:latend], gprof_lat[latstart:latend], projection_target=gk3)

#grid_xy = np.vstack((xgrid.ravel(), ygrid.ravel())).transpose()

#xy=np.concatenate([x.ravel()[:,None],y.ravel()[:,None]], axis=1)
#gridded = wradlib.comp.togrid(xy, grid_xy, ranges[-1], np.array([x.mean(), y.mean()]), rwdata.ravel(), ipoli[0],nnearest=30,p=2)
#gridded = np.ma.masked_invalid(gridded).reshape(xgrid.shape)


## PLOT
## ----
fig = plt.figure()
ax = fig.add_subplot(211, aspect='equal')
plt.pcolormesh(x, y, rwdata, cmap="spectral",vmin=0,vmax=10)
cb = plt.colorbar(shrink=0.75)
cb.set_label("mm/h")
plt.title('RADOLAN RW Product Polar Stereo \n' + rwattrs['datetime'].isoformat())
plt.grid(color='r')

ax = fig.add_subplot(212, aspect='equal')
pm2 = plt.pcolormesh(gprof_lon[latstart:latend], gprof_lat[latstart:latend],np.ma.masked_invalid(gprof_pp[latstart:latend]),
                     cmap="spectral",vmin=0,vmax=10)
cb = plt.colorbar(shrink=0.75)
cb.set_label("mm/h")
plt.xlim((bonn_lon1,bonn_lon2))
plt.ylim((bonn_lat1,bonn_lat2))
plt.title('GPM GPROF: \n' + str(pfad_gprof_g[66:74]))
#plot_borders(ax)
plt.grid(color='r')
plt.tight_layout()


plt.show()
